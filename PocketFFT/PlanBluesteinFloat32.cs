using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    internal class PlanBluesteinFloat32 : IReal1DFFTPlanFloat32, IComplex1DFFTPlanFloat32
    {
        internal int n;
        internal int n2;
        internal IComplex1DFFTPlanFloat32 plan;
        internal cmplxF[] mem;
        internal int bk_idx; // indexes into mem
        internal int bkf_idx; // indexes into mem

        public int Length => plan.Length;

        public PlanBluesteinFloat32(int length, Func<int, IComplex1DFFTPlanFloat32> innerTransformFactory)
        {
            n = length;
            n2 = Intrinsics.good_size(n * 2 - 1);
            mem = ArrayPool<cmplxF>.Shared.Rent(n + n2);
            bk_idx = 0;
            bkf_idx = bk_idx + n;

            /* initialize b_k */
            Span<float> tmp = new float[4 * n];
            Intrinsics.sincos_2pibyn(2 * n, tmp);
            Span<cmplxF> bk = mem.AsSpan(bk_idx);
            Span<cmplxF> bkf = mem.AsSpan(bkf_idx);
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
            float xn2 = 1.0f / n2;
            bkf[0].r = bk[0].r * xn2;
            bkf[0].i = bk[0].i * xn2;
            for (int m = 1; m < n; m++)
            {
                bkf[m].r = bkf[n2 - m].r = bk[m].r * xn2;
                bkf[m].i = bkf[n2 - m].i = bk[m].i * xn2;
            }

            for (int m = n; m <= n2 - n; ++m)
            {
                bkf[m].r = 0.0f;
                bkf[m].i = 0.0f;
            }

            plan = innerTransformFactory(n2);
            plan.Forward(bkf, 1.0f);
        }

        public void Dispose()
        {
            if (mem != null)
            {
                ArrayPool<cmplxF>.Shared.Return(mem);
            }
        }

        public void Forward(Span<float> c, float fct)
        {
            cmplxF[] tmp = ArrayPool<cmplxF>.Shared.Rent(n);
            for (int m = 0; m < n; ++m)
            {
                tmp[m].r = c[m];
                tmp[m].i = 0.0f;
            }

            fftblue_fft(tmp, -1, fct);

            c[0] = tmp[0].r;

            // Jank method
            MemoryMarshal.Cast<cmplxF, float>(tmp.AsSpan(1, n - 1)).CopyTo(c.Slice(1));

            // Safe method
            //for (int i = 1; i < n; i++)
            //{
            //    c[2 * i - 1] = tmp[i].r;
            //    c[2 * i] = tmp[i].i;
            //}

            ArrayPool<cmplxF>.Shared.Return(tmp);
        }

        public void Forward(Span<cmplxF> c, float fct)
        {
            fftblue_fft(c, -1, fct);
        }

        public void Backward(Span<float> c, float fct)
        {
            cmplxF[] tmp = ArrayPool<cmplxF>.Shared.Rent(n);
            tmp[0].r = c[0];
            tmp[1].i = 0.0f;

            // Jank method
            c.Slice(1, (n - 1)).CopyTo(MemoryMarshal.Cast<cmplxF, float>(tmp.AsSpan(1)));

            // Safe method
            //for (int i = 1; i < n; i++)
            //{
            //    tmp[i].r = c[2 * i - 1];
            //    tmp[i].i = c[2 * i];
            //}

            if ((n & 1) == 0)
            {
                tmp[n - 1].i = 0.0f;
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

            ArrayPool<cmplxF>.Shared.Return(tmp);
        }

        public void Backward(Span<cmplxF> c, float fct)
        {
            fftblue_fft(c, 1, fct);
        }

        private void fftblue_fft(Span<cmplxF> c, int isign, float fct)
        {
            Span<cmplxF> bk = mem.AsSpan(bk_idx);
            Span<cmplxF> bkf = mem.AsSpan(bkf_idx);
            cmplxF[] akf = ArrayPool<cmplxF>.Shared.Rent(n2);

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

            MemoryMarshal.Cast<cmplxF, float>(akf.AsSpan(n, n2 - n)).Clear();
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
            plan.Backward(akf, 1.0f);

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

            ArrayPool<cmplxF>.Shared.Return(akf);
        }
    }
}
