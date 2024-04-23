using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    public class fftblue_plan : IRealFFTPlan, IComplexFFTPlan
    {
        internal int n;
        internal int n2;
        internal cfftp_plan plan;
        internal double[] mem;
        internal int bk; // indexes into mem
        internal int bkf; // indexes into mem

        public int Length => plan.length;

        public fftblue_plan(int length)
        {
            this.n = length;
            this.n2 = Intrinsics.good_size(this.n * 2 - 1);
            this.mem = new double[2 * this.n + 2 * this.n2];
            this.bk = 0;
            this.bkf = this.bk + 2 * this.n;

            /* initialize b_k */
            double[] tmp = new double[4 * this.n];
            Intrinsics.sincos_2pibyn(2 * this.n, tmp);
            Span<double> bk = this.mem.AsSpan(this.bk);
            Span<double> bkf = this.mem.AsSpan(this.bkf);
            bk[0] = 1;
            bk[1] = 0;

            int coeff = 0;
            for (int m = 1; m < this.n; ++m)
            {
                coeff += 2 * m - 1;
                if (coeff >= 2 * this.n) coeff -= 2 * this.n;
                bk[2 * m] = tmp[2 * coeff];
                bk[2 * m + 1] = tmp[2 * coeff + 1];
            }

            /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
            double xn2 = 1.0 / this.n2;
            bkf[0] = bk[0] * xn2;
            bkf[1] = bk[1] * xn2;
            for (int m = 2; m < 2 * this.n; m += 2)
            {
                bkf[m] = bkf[2 * this.n2 - m] = bk[m] * xn2;
                bkf[m + 1] = bkf[2 * this.n2 - m + 1] = bk[m + 1] * xn2;
            }

            for (int m = 2 * this.n; m <= (2 * this.n2 - 2 * this.n + 1); ++m)
            {
                bkf[m] = 0.0;
            }

            this.plan = new cfftp_plan(this.n2);
            this.plan.Forward(MemoryMarshal.Cast<double, cmplx>(bkf), 1.0);
        }

        public void Dispose() { }

        public void Forward(Span<double> c, double fct)
        {
            fftblue_fft(c, -1, fct);
        }

        public void Forward(Span<cmplx> c, double fct)
        {
            fftblue_fft(MemoryMarshal.Cast<cmplx, double>(c), -1, fct);
        }

        public void Backward(Span<double> c, double fct)
        {
            fftblue_fft(c, 1, fct);
        }

        public void Backward(Span<cmplx> c, double fct)
        {
            fftblue_fft(MemoryMarshal.Cast<cmplx, double>(c), 1, fct);
        }

        private void fftblue_fft(Span<double> c, int isign, double fct)
        {
            int n = this.n;
            int n2 = this.n2;
            Span<double> bk = this.mem.AsSpan(this.bk);
            Span<double> bkf = this.mem.AsSpan(this.bkf);
            Span<double> akf = new double[2 * n2];

            /* initialize a_k and FFT it */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
                    akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
                }
            }
            else
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
                    akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
                }
            }
            for (int m = 2 * n; m < 2 * n2; ++m)
            {
                akf[m] = 0;
            }

            this.plan.Forward(MemoryMarshal.Cast<double, cmplx>(akf), fct);

            /* do the convolution */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n2; m += 2)
                {
                    double im = -akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
                    akf[m] = akf[m] * bkf[m] + akf[m + 1] * bkf[m + 1];
                    akf[m + 1] = im;
                }
            }
            else
            {
                for (int m = 0; m < 2 * n2; m += 2)
                {
                    double im = akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
                    akf[m] = akf[m] * bkf[m] - akf[m + 1] * bkf[m + 1];
                    akf[m + 1] = im;
                }
            }

            /* inverse FFT */
            this.plan.Backward(MemoryMarshal.Cast<double, cmplx>(akf), 1.0);

            /* multiply by b_k */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
                    c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
                }
            }
            else
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
                    c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
                }
            }

            //DEALLOC(akf);
        }

        private void rfftblue_forward(Span<double> c, double fct)
        {
            int n = this.n;
            double[] tmp = new double[2 * n];
            for (int m = 0; m < n; ++m)
            {
                tmp[2 * m] = c[m];
                tmp[2 * m + 1] = 0.0;
            }

            fftblue_fft(tmp, -1, fct);

            c[0] = tmp[0];
            tmp.AsSpan(2, (n - 1)).CopyTo(c.Slice(1));
            //memcpy(c + 1, tmp + 2, (n - 1) * sizeof(double));

            //DEALLOC(tmp);
        }

        private void rfftblue_backward(fftblue_plan plan, Span<double> c, double fct)
        {
            int n = plan.n;
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

            fftblue_fft(tmp, 1, fct);
            for (int m = 0; m < n; ++m)
            {
                c[m] = tmp[2 * m];
            }

            //DEALLOC(tmp);
        }
    }
}
