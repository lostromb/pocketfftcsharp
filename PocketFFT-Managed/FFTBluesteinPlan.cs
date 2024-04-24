using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    internal class FFTBluesteinPlan : IRealFFTPlan, IComplexFFTPlan
    {
        internal int n;
        internal int n2;
        internal ComplexFFTPackedPlan plan;
        internal cmplx[] mem;
        internal int bk; // indexes into mem
        internal int bkf; // indexes into mem

        public int Length => plan.length;

        public FFTBluesteinPlan(int length)
        {
            this.n = length;
            this.n2 = Intrinsics.good_size(this.n * 2 - 1);
            this.mem = new cmplx[this.n + this.n2];
            this.bk = 0;
            this.bkf = this.bk + this.n;

            /* initialize b_k */
            double[] tmp = new double[4 * this.n];
            Intrinsics.sincos_2pibyn(2 * this.n, tmp);
            Span<cmplx> bk = this.mem.AsSpan(this.bk);
            Span<cmplx> bkf = this.mem.AsSpan(this.bkf);
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
            double xn2 = 1.0 / this.n2;
            bkf[0].r = bk[0].r * xn2;
            bkf[0].i = bk[0].i * xn2;
            for (int m = 1; m < this.n; m++)
            {
                bkf[m].r = bkf[this.n2 - m].r = bk[m].r * xn2;
                bkf[m].i = bkf[this.n2 - m].i = bk[m].i * xn2;
            }

            for (int m = this.n; m <= this.n2 - this.n; ++m)
            {
                bkf[m].r = 0.0;
                bkf[m].i = 0.0;
            }

            this.plan = new ComplexFFTPackedPlan(this.n2);
            this.plan.Forward(bkf, 1.0);
        }

        public void Dispose() { }

        public void Forward(Span<double> c, double fct)
        {
            int n = this.n;
            double[] tmp = new double[2 * n];
            for (int m = 0; m < n; ++m)
            {
                tmp[2 * m] = c[m];
                tmp[2 * m + 1] = 0.0;
            }

            // It's possible to remove this span cast, but we would just move it to a different cast lower down in the function so I guess it stays
            fftblue_fft(MemoryMarshal.Cast<double, cmplx>(tmp), -1, fct);

            c[0] = tmp[0];
            tmp.AsSpan(2, (n - 1)).CopyTo(c.Slice(1));
            //memcpy(c + 1, tmp + 2, (n - 1) * sizeof(double));

            //DEALLOC(tmp);
        }

        public void Forward(Span<cmplx> c, double fct)
        {
            fftblue_fft(c, -1, fct);
        }

        public void Backward(Span<double> c, double fct)
        {
            int n = this.n;
            double[] tmp = new double[2 * n];
            tmp[0] = c[0];
            tmp[1] = 0.0;
            c.Slice(1, (n - 1)).CopyTo(tmp.AsSpan(2));
            //memcpy(tmp + 2, c + 1, (n - 1) * sizeof(double));

            if ((n & 1) == 0) tmp[n + 1] = 0.0;
            for (int m = 2; m < n; m += 2)
            {
                tmp[2 * n - m] = tmp[m];
                tmp[2 * n - m + 1] = -tmp[m + 1];
            }

            fftblue_fft(MemoryMarshal.Cast<double, cmplx>(tmp), 1, fct);
            for (int m = 0; m < n; ++m)
            {
                c[m] = tmp[2 * m];
            }

            //DEALLOC(tmp);
        }

        public void Backward(Span<cmplx> c, double fct)
        {
            fftblue_fft(c, 1, fct);
        }

        private void fftblue_fft(Span<cmplx> c, int isign, double fct)
        {
            int n = this.n;
            int n2 = this.n2;
            Span<cmplx> bk = this.mem.AsSpan(this.bk);
            Span<cmplx> bkf = this.mem.AsSpan(this.bkf);
            Span<cmplx> akf = new cmplx[n2];

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
                    double im = -akf[m].r * bkf[m].i + akf[m].i * bkf[m].r;
                    akf[m].r = akf[m].r * bkf[m].r + akf[m].i * bkf[m].i;
                    akf[m].i = im;
                }
            }
            else
            {
                for (int m = 0; m < n2; m++)
                {
                    double im = akf[m].r * bkf[m].i + akf[m].i * bkf[m].r;
                    akf[m].r = akf[m].r * bkf[m].r - akf[m].i * bkf[m].i;
                    akf[m].i = im;
                }
            }

            /* inverse FFT */
            this.plan.Backward(akf, 1.0);

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

            //DEALLOC(akf);
        }
    }
}
