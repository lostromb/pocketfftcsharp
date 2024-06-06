using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    internal class PlanBluesteinFloat64 : IReal1DFFTPlanFloat64, IComplex1DFFTPlanFloat64
    {
        internal int n;
        internal int n2;
        internal IComplex1DFFTPlanFloat64 plan;
        internal cmplx[] mem;
        internal int bk_idx; // indexes into mem
        internal int bkf_idx; // indexes into mem

        public int Length => plan.Length;

        public PlanBluesteinFloat64(int length, Func<int, IComplex1DFFTPlanFloat64> innerTransformFactory)
        {
            n = length;
            n2 = Intrinsics.good_size(n * 2 - 1);
            mem = ArrayPool<cmplx>.Shared.Rent(n + n2);
            bk_idx = 0;
            bkf_idx = bk_idx + n;

            /* initialize b_k */
            Span<double> tmp = new double[4 * n];
            Intrinsics.sincos_2pibyn(2 * n, tmp);
            Span<cmplx> bk = mem.AsSpan(bk_idx);
            Span<cmplx> bkf = mem.AsSpan(bkf_idx);
            bk[0].r = 1;
            bk[0].i = 0;

            int coeff = 0;
            for (int m = 1; m < n; ++m)
            {
                coeff += 2 * m - 1;
                if (coeff >= 2 * n)
                {
                    coeff -= 2 * n;
                }

                bk[m].r = tmp[2 * coeff];
                bk[m].i = tmp[2 * coeff + 1];
            }

            /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
            double xn2 = 1.0 / n2;
            bkf[0].r = bk[0].r * xn2;
            bkf[0].i = bk[0].i * xn2;
            for (int m = 1; m < n; m++)
            {
                bkf[m].r = bkf[n2 - m].r = bk[m].r * xn2;
                bkf[m].i = bkf[n2 - m].i = bk[m].i * xn2;
            }

            // OPT Span.Clear(....)
            for (int m = n; m <= n2 - n; ++m)
            {
                bkf[m].r = 0.0;
                bkf[m].i = 0.0;
            }

            plan = innerTransformFactory(n2);
            plan.Forward(bkf, 1.0);
        }

        public void Dispose()
        {
            if (mem != null)
            {
                ArrayPool<cmplx>.Shared.Return(mem);
            }
        }

        public void Forward(Span<double> c, double fct)
        {
            cmplx[] tmp = ArrayPool<cmplx>.Shared.Rent(n);
            for (int m = 0; m < n; ++m)
            {
                tmp[m].r = c[m];
                tmp[m].i = 0.0;
            }

            fftblue_fft(tmp, -1, fct);

            c[0] = tmp[0].r;

            // Jank method
            MemoryMarshal.Cast<cmplx, double>(tmp.AsSpan(1, n - 1)).CopyTo(c.Slice(1));
            
            // Safe method
            //for (int i = 1; i < n; i++)
            //{
            //    c[2 * i - 1] = tmp[i].r;
            //    c[2 * i] = tmp[i].i;
            //}

            ArrayPool<cmplx>.Shared.Return(tmp);
        }

        public void Forward(Span<cmplx> c, double fct)
        {
            fftblue_fft(c, -1, fct);
        }

        public void Backward(Span<double> c, double fct)
        {
            cmplx[] tmp = ArrayPool<cmplx>.Shared.Rent(n);
            tmp[0].r = c[0];
            tmp[1].i = 0.0;
            
            // Jank method
            c.Slice(1, (n - 1)).CopyTo(MemoryMarshal.Cast<cmplx, double>(tmp.AsSpan(1)));

            // Safe method
            //for (int i = 1; i < n; i++)
            //{
            //    tmp[i].r = c[2 * i - 1];
            //    tmp[i].i = c[2 * i];
            //}

            if ((n & 1) == 0)
            {
                tmp[n - 1].i = 0.0;
            }

            for (int m = 1; m < n; m++)
            {
                Intrinsics.BLUESTEINSTEP0(ref tmp[n - m], ref tmp[m]);
            }

            fftblue_fft(tmp, 1, fct);
            for (int m = 0; m < n; ++m)
            {
                c[m] = tmp[m].r;
            }

            ArrayPool<cmplx>.Shared.Return(tmp);
        }

        public void Backward(Span<cmplx> c, double fct)
        {
            fftblue_fft(c, 1, fct);
        }

        private void fftblue_fft(Span<cmplx> c, int isign, double fct)
        {
            Span<cmplx> bk = mem.AsSpan(bk_idx);
            Span<cmplx> bkf = mem.AsSpan(bkf_idx);
            cmplx[] akf = ArrayPool<cmplx>.Shared.Rent(n2);

            /* initialize a_k and FFT it */
            if (isign > 0)
            {
                for (int m = 0; m < n; m++)
                {
                    Intrinsics.BLUESTEINSTEP1A(ref akf[m], ref bk[m], ref c[m]);
                }
            }
            else
            {
                for (int m = 0; m < n; m++)
                {
                    Intrinsics.BLUESTEINSTEP1B(ref akf[m], ref c[m], ref bk[m]);
                }
            }

            MemoryMarshal.Cast<cmplx, double>(akf.AsSpan(n, n2 - n)).Clear();
            //for (int m = n; m < n2; ++m)
            //{
            //    akf[m].r = 0;
            //    akf[m].i = 0;
            //}

            plan.Forward(akf, fct);

            /* do the convolution */
            if (isign > 0)
            {
                for (int m = 0; m < n2; m++)
                {
                    Intrinsics.BLUESTEINSTEP2A(ref akf[m], ref bkf[m]);
                }
            }
            else
            {
                for (int m = 0; m < n2; m++)
                {
                    Intrinsics.BLUESTEINSTEP2B(ref akf[m], ref bkf[m]);
                }
            }

            /* inverse FFT */
            plan.Backward(akf, 1.0);

            /* multiply by b_k */
            if (isign > 0)
            {
                for (int m = 0; m < n; m++)
                {
                    Intrinsics.BLUESTEINSTEP3A(ref c[m], ref bk[m], ref akf[m]);
                }
            }
            else
            {
                for (int m = 0; m < n; m++)
                {
                    Intrinsics.BLUESTEINSTEP3B(ref c[m], ref bk[m], ref akf[m]);
                }
            }

            ArrayPool<cmplx>.Shared.Return(akf);
        }
    }
}
