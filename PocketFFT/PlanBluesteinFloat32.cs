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
        internal PlanComplexPackedFloat32 plan;
        internal cmplxF[] mem;
        internal int bk; // indexes into mem
        internal int bkf; // indexes into mem

        public int Length => plan.length;

        public PlanBluesteinFloat32(int length)
        {
            this.n = length;
            this.n2 = Intrinsics.good_size(this.n * 2 - 1);
            this.mem = ArrayPool<cmplxF>.Shared.Rent(this.n + this.n2);
            this.bk = 0;
            this.bkf = this.bk + this.n;

            /* initialize b_k */
            Span<float> tmp = new float[4 * this.n];
            Intrinsics.sincos_2pibyn(2 * this.n, tmp);
            Span<cmplxF> bk = this.mem.AsSpan(this.bk);
            Span<cmplxF> bkf = this.mem.AsSpan(this.bkf);
            bk[0].r = 1;
            bk[0].i = 0;

            int coeff = 0;
            for (int m = 1; m < this.n; ++m)
            {
                coeff += 2 * m - 1;
                if (coeff >= 2 * this.n)
                {
                    coeff -= 2 * this.n;
                }

                bk[m].r = tmp[2 * coeff];
                bk[m].i = tmp[2 * coeff + 1];
            }

            /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
            float xn2 = 1.0f / this.n2;
            bkf[0].r = bk[0].r * xn2;
            bkf[0].i = bk[0].i * xn2;
            for (int m = 1; m < this.n; m++)
            {
                bkf[m].r = bkf[this.n2 - m].r = bk[m].r * xn2;
                bkf[m].i = bkf[this.n2 - m].i = bk[m].i * xn2;
            }

            for (int m = this.n; m <= this.n2 - this.n; ++m)
            {
                bkf[m].r = 0.0f;
                bkf[m].i = 0.0f;
            }

            this.plan = new PlanComplexPackedFloat32(this.n2);
            this.plan.Forward(bkf, 1.0f);
        }

        public void Dispose()
        {
            if (this.mem != null)
            {
                ArrayPool<cmplxF>.Shared.Return(this.mem);
            }
        }

        public void Forward(Span<float> c, float fct)
        {
            int n = this.n;
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
            int n = this.n;
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
                tmp[n - m].r = tmp[m].r;
                tmp[n - m].i = -tmp[m].i;
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
            int n = this.n;
            int n2 = this.n2;
            Span<cmplxF> bk = this.mem.AsSpan(this.bk);
            Span<cmplxF> bkf = this.mem.AsSpan(this.bkf);
            cmplxF[] akf = ArrayPool<cmplxF>.Shared.Rent(n2);

            /* initialize a_k and FFT it */
            if (isign > 0)
            {
                for (int m = 0; m < n; m++)
                {
                    akf[m].r = c[m].r * bk[m].r - c[m].i * bk[m].i;
                    akf[m].i = c[m].r * bk[m].i + c[m].i * bk[m].r;
                }
            }
            else
            {
                for (int m = 0; m < n; m++)
                {
                    akf[m].r = c[m].r * bk[m].r + c[m].i * bk[m].i;
                    akf[m].i = -c[m].r * bk[m].i + c[m].i * bk[m].r;
                }
            }

            for (int m = n; m < n2; ++m)
            {
                akf[m].r = 0;
                akf[m].i = 0;
            }

            this.plan.Forward(akf, fct);

            /* do the convolution */
            if (isign > 0)
            {
                for (int m = 0; m < n2; m++)
                {
                    float im = -akf[m].r * bkf[m].i + akf[m].i * bkf[m].r;
                    akf[m].r = akf[m].r * bkf[m].r + akf[m].i * bkf[m].i;
                    akf[m].i = im;
                }
            }
            else
            {
                for (int m = 0; m < n2; m++)
                {
                    float im = akf[m].r * bkf[m].i + akf[m].i * bkf[m].r;
                    akf[m].r = akf[m].r * bkf[m].r - akf[m].i * bkf[m].i;
                    akf[m].i = im;
                }
            }

            /* inverse FFT */
            this.plan.Backward(akf, 1.0f);

            /* multiply by b_k */
            if (isign > 0)
            {
                for (int m = 0; m < n; m++)
                {
                    c[m].r = bk[m].r * akf[m].r - bk[m].i * akf[m].i;
                    c[m].i = bk[m].i * akf[m].r + bk[m].r * akf[m].i;
                }
            }
            else
            {
                for (int m = 0; m < n; m++)
                {
                    c[m].r = bk[m].r * akf[m].r + bk[m].i * akf[m].i;
                    c[m].i = -bk[m].i * akf[m].r + bk[m].r * akf[m].i;
                }
            }

            ArrayPool<cmplxF>.Shared.Return(akf);
        }
    }
}
