﻿using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;

namespace PocketFFT
{
    internal struct cfftp_fctdata
    {
        internal int fct;
        internal int tw; // originally these were cmplx* pointers, but they just indexed into plan.mem
        internal int tws; // so they have been replaced by basic integer indexes

        internal cfftp_fctdata(int fct, int tw, int tws)
        {
            this.fct = fct;
            this.tw = tw;
            this.tws = tws;
        }
    }

#if NET8_0_OR_GREATER
    [System.Runtime.CompilerServices.InlineArray(Constants.NFCT)]
    internal struct cfftp_fctdata_array
    {
        internal cfftp_fctdata data;
    }
#endif

    internal class ComplexFFTPackedPlan : IComplexFFTPlan
    {
        internal int length;
        internal int nfct;
        internal cmplx[] mem;
#if NET8_0_OR_GREATER
        internal cfftp_fctdata_array fct;
#else
        internal cfftp_fctdata[] fct; // [Constants.NFCT]
#endif

        public int Length => length;

        public ComplexFFTPackedPlan(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

            this.length = length;
            this.nfct = 0;
#if !NET8_0_OR_GREATER
            if (this.fct == null) // if it's not an inline struct array...
            {
                this.fct = new cfftp_fctdata[Constants.NFCT];
            }
#endif

            for (int i = 0; i < Constants.NFCT; ++i)
            {
                this.fct[i] = new cfftp_fctdata(0, 0, 0);
            }

            this.mem = null;
            if (length == 1)
            {
                return;
            }

            if (!cfftp_factorize())
            {
                throw new ArithmeticException($"Could not factorize FFT of length {length}");
            }

            int tws = cfftp_twsize();
            this.mem = new cmplx[tws];

            cfftp_comp_twiddle();
        }

        public void Dispose() { }

        public void Forward(Span<cmplx> c, double fct)
        {
            pass_all(c, fct, -1);
        }

        public void Backward(Span<cmplx> c, double fct)
        {
            pass_all(c, fct, 1);
        }

        private bool cfftp_factorize()
        {
            int length = this.length;
            int nfct = 0;

            while ((length % 4) == 0)
            {
                if (nfct >= Constants.NFCT)
                {
                    return false;
                }

                this.fct[nfct++].fct = 4;
                length >>= 2;
            }
            if ((length % 2) == 0)
            {
                length >>= 1;
                // factor 2 should be at the front of the factor list
                if (nfct >= Constants.NFCT)
                {
                    return false;
                }

                this.fct[nfct++].fct = 2;
                Intrinsics.Swap(ref this.fct[0].fct, ref this.fct[nfct - 1].fct);
            }

            int maxl = (int)(Math.Sqrt((double)length)) + 1;
            for (int divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
            {
                if ((length % divisor) == 0)
                {
                    while ((length % divisor) == 0)
                    {
                        if (nfct >= Constants.NFCT)
                        {
                            return false;
                        }

                        this.fct[nfct++].fct = divisor;
                        length /= divisor;
                    }

                    maxl = (int)(Math.Sqrt((double)length)) + 1;
                }
            }

            if (length > 1)
            {
                this.fct[nfct++].fct = length;
            }

            this.nfct = nfct;
            return true;
        }

        private int cfftp_twsize()
        {
            int twsize = 0, l1 = 1;
            for (int k = 0; k < this.nfct; ++k)
            {
                int ip = this.fct[k].fct, ido = this.length / (l1 * ip);
                twsize += (ip - 1) * (ido - 1);
                if (ip > 11)
                {
                    twsize += ip;
                }

                l1 *= ip;
            }

            return twsize;
        }

        private void cfftp_comp_twiddle()
        {
            int length = this.length;
            double[] twid = new double[2 * length];
            Intrinsics.sincos_2pibyn(length, twid);
            int l1 = 1;
            int memofs = 0;
            for (int k = 0; k < this.nfct; ++k)
            {
                int ip = this.fct[k].fct, ido = length / (l1 * ip);
                this.fct[k].tw = memofs;
                Span<cmplx> tw = this.mem.AsSpan(memofs);
                memofs += (ip - 1) * (ido - 1);
                for (int j = 1; j < ip; ++j)
                {
                    for (int i = 1; i < ido; ++i)
                    {
                        tw[(j - 1) * (ido - 1) + i - 1].r = twid[2 * j * l1 * i];
                        tw[(j - 1) * (ido - 1) + i - 1].i = twid[2 * j * l1 * i + 1];
                    }
                }
                if (ip > 11)
                {
                    this.fct[k].tws = memofs;
                    Span<cmplx> tws = this.mem.AsSpan(memofs);
                    memofs += ip;
                    for (int j = 0; j < ip; ++j)
                    {
                        tws[j].r = twid[2 * j * l1 * ido];
                        tws[j].i = twid[2 * j * l1 * ido + 1];
                    }
                }
                l1 *= ip;
            }
        }

        private void pass_all(Span<cmplx> c, double fct, int sign)
        {
            if (this.length == 1)
            {
                return;
            }

            int len = this.length;
            int l1 = 1, nf = nfct;
            Span<cmplx> ch = new cmplx[len];
            Span<cmplx> p1 = c, p2 = ch;

            for (int k1 = 0; k1 < nf; k1++)
            {
                int ip = this.fct[k1].fct;
                int l2 = ip * l1;
                int ido = len / l2;
                if (ip == 4)
                {
                    if (sign > 0)
                    {
                        pass4b(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                    else
                    {
                        pass4f(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                }
                else if (ip == 2)
                {
                    if (sign > 0)
                    {
                        pass2b(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                    else
                    {
                        pass2f(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                }
                else if (ip == 3)
                {
                    if (sign > 0)
                    {
                        pass3b(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                    else
                    {
                        pass3f(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                }
                else if (ip == 5)
                {
                    if (sign > 0)
                    {
                        pass5b(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                    else
                    {
                        pass5f(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw));
                    }
                }
                else if (ip == 7)
                {
                    pass7(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw), sign);
                }
                else if (ip == 11)
                {
                    pass11(ido, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw), sign);
                }
                else
                {
                    passg(ido, ip, l1, p1, p2, this.mem.AsSpan(this.fct[k1].tw), this.mem.AsSpan(this.fct[k1].tws), sign);
                    {
                        Span<cmplx> tmp = p1;
                        p1 = p2;
                        p2 = tmp;
                    }
                }

                {
                    Span<cmplx> tmp = p1;
                    p1 = p2;
                    p2 = tmp;
                }
                l1 = l2;
            }

            if (p1 != c) // FIXME span address comparison?
            {
                if (fct != 1.0)
                    for (int i = 0; i < len; ++i)
                    {
                        c[i].r = ch[i].r * fct;
                        c[i].i = ch[i].i * fct;
                    }
                else
                {
                    p1.Slice(0, len).CopyTo(c);
                    //memcpy(c, p1, len * sizeof(cmplx));
                }
            }
            else
            {
                if (fct != 1.0)
                {
                    for (int i = 0; i < len; ++i)
                    {
                        c[i].r *= fct;
                        c[i].i *= fct;
                    }
                }
            }

            //DEALLOC(ch);
        }

        private static void pass2b_opt(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;

            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    Intrinsics.PMC(ref ch[(0) + ido * ((k) + l1 * (0))],
                        ref ch[(0) + ido * ((k) + l1 * (1))],
                        ref cc[(0) + ido * ((0) + cdim * (k))],
                        ref cc[(0) + ido * ((1) + cdim * (k))]);
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    Intrinsics.PMC(ref ch[(0) + ido * ((k) + l1 * (0))],
                        ref ch[(0) + ido * ((k) + l1 * (1))],
                        ref cc[(0) + ido * ((0) + cdim * (k))],
                        ref cc[(0) + ido * ((1) + cdim * (k))]);

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t = default;
                        Intrinsics.PMC(ref ch[(i) + ido * ((k) + l1 * (0))],
                            ref t,
                            ref cc[(i) + ido * ((0) + cdim * (k))],
                            ref cc[(i) + ido * ((1) + cdim * (k))]);
                        Intrinsics.A_EQ_B_MUL_C(
                            ref ch[(i) + ido * ((k) + l1 * (1))],
                            ref wa[(i) - 1 + (0) * (ido - 1)],
                            ref t);
                    }
                }
            }
        }

        private static void pass2b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;
            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t;
                        ch[(i) + ido * ((k) + l1 * (0))].r = cc[(i) + ido * ((0) + cdim * (k))].r + cc[(i) + ido * ((1) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = cc[(i) + ido * ((0) + cdim * (k))].i + cc[(i) + ido * ((1) + cdim * (k))].i;
                        t.r = cc[(i) + ido * ((0) + cdim * (k))].r - cc[(i) + ido * ((1) + cdim * (k))].r;
                        t.i = cc[(i) + ido * ((0) + cdim * (k))].i - cc[(i) + ido * ((1) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (0) * (ido - 1)].r * t.r - wa[(i) - 1 + (0) * (ido - 1)].i * t.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (0) * (ido - 1)].r * t.i + wa[(i) - 1 + (0) * (ido - 1)].i * t.r;
                    }
                }
            }
        }

        private static void pass2f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;
            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t;
                        ch[(i) + ido * ((k) + l1 * (0))].r = cc[(i) + ido * ((0) + cdim * (k))].r + cc[(i) + ido * ((1) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = cc[(i) + ido * ((0) + cdim * (k))].i + cc[(i) + ido * ((1) + cdim * (k))].i;
                        t.r = cc[(i) + ido * ((0) + cdim * (k))].r - cc[(i) + ido * ((1) + cdim * (k))].r;
                        t.i = cc[(i) + ido * ((0) + cdim * (k))].i - cc[(i) + ido * ((1) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (0) * (ido - 1)].r * t.r + wa[(i) - 1 + (0) * (ido - 1)].i * t.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (0) * (ido - 1)].r * t.i - wa[(i) - 1 + (0) * (ido - 1)].i * t.r;
                    }
                }
            }
        }

        private static void pass3b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 3;
            const double tw1r = -0.5, tw1i = 0.86602540378443864676;
            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r; ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                        t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                        t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                        ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                        cmplx ca, cb;
                        ca.r = t0.r + tw1r * t1.r;
                        ca.i = t0.i + tw1r * t1.i;
                        cb.i = tw1i * t2.r;
                        cb.r = -(tw1i * t2.i);
                        ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                        ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                        ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                        ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                    }

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2;
                        t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((2) + cdim * (k))].r;
                        t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((2) + cdim * (k))].i;
                        t2.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((2) + cdim * (k))].r;
                        t2.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((2) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                        cmplx ca, cb, da, db;
                        ca.r = t0.r + tw1r * t1.r;
                        ca.i = t0.i + tw1r * t1.i;
                        cb.i = tw1i * t2.r;
                        cb.r = -(tw1i * t2.i);
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }

        private static void pass3f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 3;
            const double tw1r = -0.5, tw1i = -0.86602540378443864676;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t3 = cc[(i) + ido * ((0) + cdim * (k))], t4, t5;
                        t4.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((2) + cdim * (k))].r;
                        t4.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((2) + cdim * (k))].i;
                        t5.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((2) + cdim * (k))].r;
                        t5.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((2) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (0))].r = t3.r + t4.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t3.i + t4.i;
                        cmplx ca2, cb2, da, db;
                        ca2.r = t3.r + tw1r * t4.r;
                        ca2.i = t3.i + tw1r * t4.i;
                        cb2.i = tw1i * t5.r;
                        cb2.r = -(tw1i * t5.i);
                        da.r = ca2.r + cb2.r;
                        da.i = ca2.i + cb2.i;
                        db.r = ca2.r - cb2.r;
                        db.i = ca2.i - cb2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }

        private static void pass4b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 4;
            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1, t2, t3, t4;
                    t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
                    t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
                    t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
                    t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    Intrinsics.ROT90(ref t4);

                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1, t2, t3, t4;
                    t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
                    t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
                    t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
                    t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    Intrinsics.ROT90(ref t4);

                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx c2, c3, c4, t5, t6, t7, t8;
                        cmplx cc0 = cc[(i) + ido * ((0) + cdim * (k))],
                            cc1 = cc[(i) + ido * ((1) + cdim * (k))],
                            cc2 = cc[(i) + ido * ((2) + cdim * (k))],
                            cc3 = cc[(i) + ido * ((3) + cdim * (k))];
                        t6.r = cc0.r + cc2.r;
                        t6.i = cc0.i + cc2.i;
                        t5.r = cc0.r - cc2.r;
                        t5.i = cc0.i - cc2.i;
                        t7.r = cc1.r + cc3.r;
                        t7.i = cc1.i + cc3.i;
                        t8.r = cc1.r - cc3.r;
                        t8.i = cc1.i - cc3.i;
                        Intrinsics.ROT90(ref t8);

                        cmplx wa0 = wa[(i) - 1 + (0) * (ido - 1)],
                        wa1 = wa[(i) - 1 + (1) * (ido - 1)],
                        wa2 = wa[(i) - 1 + (2) * (ido - 1)];
                        ch[(i) + ido * ((k) + l1 * (0))].r = t6.r + t7.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t6.i + t7.i;
                        c3.r = t6.r - t7.r;
                        c3.i = t6.i - t7.i;
                        c2.r = t5.r + t8.r;
                        c2.i = t5.i + t8.i;
                        c4.r = t5.r - t8.r;
                        c4.i = t5.i - t8.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa0.r * c2.r - wa0.i * c2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa0.r * c2.i + wa0.i * c2.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa1.r * c3.r - wa1.i * c3.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa1.r * c3.i + wa1.i * c3.r;
                        ch[(i) + ido * ((k) + l1 * (3))].r = wa2.r * c4.r - wa2.i * c4.i;
                        ch[(i) + ido * ((k) + l1 * (3))].i = wa2.r * c4.i + wa2.i * c4.r;
                    }
                }
            }
        }

        private static void pass4f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 4;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1, t2, t3, t4;
                    t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
                    t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
                    t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
                    t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    Intrinsics.ROTM90(ref t4);

                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1, t2, t3, t4;
                    t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
                    t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
                    t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
                    t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    Intrinsics.ROTM90(ref t4);

                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx c2, c3, c4, t5, t6, t7, t8;
                        cmplx cc0 = cc[(i) + ido * ((0) + cdim * (k))],
                            cc1 = cc[(i) + ido * ((1) + cdim * (k))],
                            cc2 = cc[(i) + ido * ((2) + cdim * (k))],
                            cc3 = cc[(i) + ido * ((3) + cdim * (k))];
                        t6.r = cc0.r + cc2.r;
                        t6.i = cc0.i + cc2.i;
                        t5.r = cc0.r - cc2.r;
                        t5.i = cc0.i - cc2.i;
                        t7.r = cc1.r + cc3.r;
                        t7.i = cc1.i + cc3.i;
                        t8.r = cc1.r - cc3.r;
                        t8.i = cc1.i - cc3.i;
                        Intrinsics.ROTM90(ref t8);

                        cmplx wa0 = wa[(i) - 1 + (0) * (ido - 1)],
                        wa1 = wa[(i) - 1 + (1) * (ido - 1)],
                        wa2 = wa[(i) - 1 + (2) * (ido - 1)];
                        ch[(i) + ido * ((k) + l1 * (0))].r = t6.r + t7.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t6.i + t7.i;
                        c3.r = t6.r - t7.r;
                        c3.i = t6.i - t7.i;
                        c2.r = t5.r + t8.r;
                        c2.i = t5.i + t8.i;
                        c4.r = t5.r - t8.r;
                        c4.i = t5.i - t8.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa0.r * c2.r + wa0.i * c2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa0.r * c2.i - wa0.i * c2.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa1.r * c3.r + wa1.i * c3.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa1.r * c3.i - wa1.i * c3.r;
                        ch[(i) + ido * ((k) + l1 * (3))].r = wa2.r * c4.r + wa2.i * c4.i;
                        ch[(i) + ido * ((k) + l1 * (3))].i = wa2.r * c4.i - wa2.i * c4.r;
                    }
                }
            }
        }

        static void pass5b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {

            const int cdim = 5;
            const double tw1r = 0.3090169943749474241,
                tw1i = 0.95105651629515357212,
                tw2r = -0.8090169943749474241,
                tw2i = 0.58778525229247312917;
            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));

            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                    {
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    {
                        t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                    {
                        cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                        }
                    }
                }

            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t4.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((3) + cdim * (k))].r; t2.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((3) + cdim * (k))].i; t3.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((3) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass5f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 5;
            const double tw1r = 0.3090169943749474241,
                tw1i = -0.95105651629515357212,
                tw2r = -0.8090169943749474241,
                tw2i = -0.58778525229247312917;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                    {
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    {
                        t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                    {
                        cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t4.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((3) + cdim * (k))].r; t2.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((3) + cdim * (k))].i; t3.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((3) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass7(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, int sign)
        {
            const int cdim = 7;
            double tw1r = 0.623489801858733530525,
                tw1i = sign * 0.7818314824680298087084,
                tw2r = -0.222520933956314404289,
                tw2i = sign * 0.9749279121818236070181,
                tw3r = -0.9009688679024191262361,
                tw3i = sign * 0.4338837391175581204758;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                    {
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                    }
                    {
                        t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
                    }
                    {
                        t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                    {
                        cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                        {
                            t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                        {
                            cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t1 = cc[(i) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                        {
                            t2.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((6) + cdim * (k))].r; t2.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((6) + cdim * (k))].i; t7.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((6) + cdim * (k))].r; t7.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((6) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((5) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((5) + cdim * (k))].i; t6.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((5) + cdim * (k))].r; t6.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((5) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(i) + ido * ((3) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((3) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t5.r = cc[(i) + ido * ((3) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t5.i = cc[(i) + ido * ((3) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(i) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (6))].r = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (6))].i = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (5))].r = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (5))].i = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass11(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, int sign)
        {
            const int cdim = 11;
            double tw1r = 0.8412535328311811688618,
                tw1i = sign * 0.5406408174555975821076,
                tw2r = 0.4154150130018864255293,
                tw2i = sign * 0.9096319953545183714117,
                tw3r = -0.1423148382732851404438,
                tw3i = sign * 0.9898214418809327323761,
                tw4r = -0.6548607339452850640569,
                tw4i = sign * 0.755749574354258283774,
                tw5r = -0.9594929736144973898904,
                tw5i = sign * 0.2817325568414296977114;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                    {
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
                    }
                    {
                        t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
                    }
                    {
                        t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
                    }
                    {
                        t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
                    }
                    {
                        t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                    {
                        cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (10))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (9))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (8))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (7))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                        {
                            t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
                        }
                        {
                            t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
                        }
                        {
                            t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                        {
                            cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (10))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (9))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (8))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (7))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t1 = cc[(i) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                        {
                            t2.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((10) + cdim * (k))].r; t2.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((10) + cdim * (k))].i; t11.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((10) + cdim * (k))].r; t11.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((10) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((9) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((9) + cdim * (k))].i; t10.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((9) + cdim * (k))].r; t10.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((9) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(i) + ido * ((3) + cdim * (k))].r + cc[(i) + ido * ((8) + cdim * (k))].r; t4.i = cc[(i) + ido * ((3) + cdim * (k))].i + cc[(i) + ido * ((8) + cdim * (k))].i; t9.r = cc[(i) + ido * ((3) + cdim * (k))].r - cc[(i) + ido * ((8) + cdim * (k))].r; t9.i = cc[(i) + ido * ((3) + cdim * (k))].i - cc[(i) + ido * ((8) + cdim * (k))].i;
                        }
                        {
                            t5.r = cc[(i) + ido * ((4) + cdim * (k))].r + cc[(i) + ido * ((7) + cdim * (k))].r; t5.i = cc[(i) + ido * ((4) + cdim * (k))].i + cc[(i) + ido * ((7) + cdim * (k))].i; t8.r = cc[(i) + ido * ((4) + cdim * (k))].r - cc[(i) + ido * ((7) + cdim * (k))].r; t8.i = cc[(i) + ido * ((4) + cdim * (k))].i - cc[(i) + ido * ((7) + cdim * (k))].i;
                        }
                        {
                            t6.r = cc[(i) + ido * ((5) + cdim * (k))].r + cc[(i) + ido * ((6) + cdim * (k))].r; t6.i = cc[(i) + ido * ((5) + cdim * (k))].i + cc[(i) + ido * ((6) + cdim * (k))].i; t7.r = cc[(i) + ido * ((5) + cdim * (k))].r - cc[(i) + ido * ((6) + cdim * (k))].r; t7.i = cc[(i) + ido * ((5) + cdim * (k))].i - cc[(i) + ido * ((6) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(i) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (10))].r = wa[(i) - 1 + (10 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (10 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (10))].i = wa[(i) - 1 + (10 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (10 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (9))].r = wa[(i) - 1 + (9 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (9 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (9))].i = wa[(i) - 1 + (9 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (9 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (8))].r = wa[(i) - 1 + (8 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (8 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (8))].i = wa[(i) - 1 + (8 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (8 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (7))].r = wa[(i) - 1 + (7 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (7 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (7))].i = wa[(i) - 1 + (7 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (7 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (5))].r = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (5))].i = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (6))].r = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (6))].i = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void passg(int ido, int ip, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, Span<cmplx> csarr, int sign)
        {
            int cdim = ip;
            int ipph = (ip + 1) / 2;
            int idl1 = ido * l1;

            if (Constants.NAIL_TEST) NailTest.PrintComplexArray(cc.Slice(0, cdim * l1));
            Span<cmplx> wal = new cmplx[ip];
            wal[0] = new cmplx(1.0, 0.0);

            for (int i = 1; i < ip; ++i)
            {
                wal[i] = new cmplx(csarr[i].r, sign * csarr[i].i);
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 0; i < ido; ++i)
                {
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))];
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
            {
                for (int k = 0; k < l1; ++k)
                {
                    for (int i = 0; i < ido; ++i)
                    {
                        ch[(i) + ido * ((k) + l1 * (j))].r = cc[(i) + ido * ((j) + cdim * (k))].r + cc[(i) + ido * ((jc) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (j))].i = cc[(i) + ido * ((j) + cdim * (k))].i + cc[(i) + ido * ((jc) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (jc))].r = cc[(i) + ido * ((j) + cdim * (k))].r - cc[(i) + ido * ((jc) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (jc))].i = cc[(i) + ido * ((j) + cdim * (k))].i - cc[(i) + ido * ((jc) + cdim * (k))].i;
                    }
                }
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 0; i < ido; ++i)
                {
                    cmplx tmp = ch[(i) + ido * ((k) + l1 * (0))];
                    for (int j = 1; j < ipph; ++j)
                    {
                        tmp.r = tmp.r + ch[(i) + ido * ((k) + l1 * (j))].r;
                        tmp.i = tmp.i + ch[(i) + ido * ((k) + l1 * (j))].i;
                    }

                    cc[(i) + ido * ((k) + l1 * (0))] = tmp;
                }
            }

            for (int l = 1, lc = ip - 1; l < ipph; ++l, --lc)
            {
                // j=0
                for (int ik = 0; ik < idl1; ++ik)
                {
                    cc[(ik) + idl1 * (l)].r = ch[(ik) + idl1 * (0)].r + wal[l].r * ch[(ik) + idl1 * (1)].r + wal[2 * l].r * ch[(ik) + idl1 * (2)].r;
                    cc[(ik) + idl1 * (l)].i = ch[(ik) + idl1 * (0)].i + wal[l].r * ch[(ik) + idl1 * (1)].i + wal[2 * l].r * ch[(ik) + idl1 * (2)].i;
                    cc[(ik) + idl1 * (lc)].r = -wal[l].i * ch[(ik) + idl1 * (ip - 1)].i - wal[2 * l].i * ch[(ik) + idl1 * (ip - 2)].i;
                    cc[(ik) + idl1 * (lc)].i = wal[l].i * ch[(ik) + idl1 * (ip - 1)].r + wal[2 * l].i * ch[(ik) + idl1 * (ip - 2)].r;
                }

                int iwal = 2 * l;
                int j = 3, jc = ip - 3;
                for (; j < ipph - 1; j += 2, jc -= 2)
                {
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal = wal[iwal];
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal2 = wal[iwal];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)].r += ch[(ik) + idl1 * (j)].r * xwal.r + ch[(ik) + idl1 * (j + 1)].r * xwal2.r;
                        cc[(ik) + idl1 * (l)].i += ch[(ik) + idl1 * (j)].i * xwal.r + ch[(ik) + idl1 * (j + 1)].i * xwal2.r;
                        cc[(ik) + idl1 * (lc)].r -= ch[(ik) + idl1 * (jc)].i * xwal.i + ch[(ik) + idl1 * (jc - 1)].i * xwal2.i;
                        cc[(ik) + idl1 * (lc)].i += ch[(ik) + idl1 * (jc)].r * xwal.i + ch[(ik) + idl1 * (jc - 1)].r * xwal2.i;
                    }
                }

                for (; j < ipph; ++j, --jc)
                {
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal = wal[iwal];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)].r += ch[(ik) + idl1 * (j)].r * xwal.r;
                        cc[(ik) + idl1 * (l)].i += ch[(ik) + idl1 * (j)].i * xwal.r;
                        cc[(ik) + idl1 * (lc)].r -= ch[(ik) + idl1 * (jc)].i * xwal.i;
                        cc[(ik) + idl1 * (lc)].i += ch[(ik) + idl1 * (jc)].r * xwal.i;
                    }
                }
            }
            //DEALLOC(wal);

            // shuffling and twiddling
            if (ido == 1)
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                {
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cmplx t1 = cc[(ik) + idl1 * (j)], t2 = cc[(ik) + idl1 * (jc)];
                        cc[(ik) + idl1 * (j)].r = t1.r + t2.r; cc[(ik) + idl1 * (j)].i = t1.i + t2.i;
                        cc[(ik) + idl1 * (jc)].r = t1.r - t2.r; cc[(ik) + idl1 * (jc)].i = t1.i - t2.i;
                    }
                }
            }
            else
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        cmplx t1 = cc[(0) + ido * ((k) + l1 * (j))], t2 = cc[(0) + ido * ((k) + l1 * (jc))];
                        cc[(0) + ido * ((k) + l1 * (j))].r = t1.r + t2.r;
                        cc[(0) + ido * ((k) + l1 * (j))].i = t1.i + t2.i;
                        cc[(0) + ido * ((k) + l1 * (jc))].r = t1.r - t2.r;
                        cc[(0) + ido * ((k) + l1 * (jc))].i = t1.i - t2.i;

                        for (int i = 1; i < ido; ++i)
                        {
                            cmplx x1, x2;
                            x1.r = cc[(i) + ido * ((k) + l1 * (j))].r + cc[(i) + ido * ((k) + l1 * (jc))].r;
                            x1.i = cc[(i) + ido * ((k) + l1 * (j))].i + cc[(i) + ido * ((k) + l1 * (jc))].i;
                            x2.r = cc[(i) + ido * ((k) + l1 * (j))].r - cc[(i) + ido * ((k) + l1 * (jc))].r;
                            x2.i = cc[(i) + ido * ((k) + l1 * (j))].i - cc[(i) + ido * ((k) + l1 * (jc))].i;

                            int idij = (j - 1) * (ido - 1) + i - 1;
                            cc[(i) + ido * ((k) + l1 * (j))].r = wa[idij].r * x1.r - sign * wa[idij].i * x1.i;
                            cc[(i) + ido * ((k) + l1 * (j))].i = wa[idij].r * x1.i + sign * wa[idij].i * x1.r;
                            idij = (jc - 1) * (ido - 1) + i - 1;
                            cc[(i) + ido * ((k) + l1 * (jc))].r = wa[idij].r * x2.r - sign * wa[idij].i * x2.i;
                            cc[(i) + ido * ((k) + l1 * (jc))].i = wa[idij].r * x2.i + sign * wa[idij].i * x2.r;
                        }
                    }
                }
            }
        }
    }
}