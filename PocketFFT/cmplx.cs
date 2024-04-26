using System;
using System.Collections.Generic;
using System.Text;

namespace PocketFFT
{
    public struct cmplx
    {
        public double r;
        public double i;

        public cmplx(double r, double i)
        {
            this.r = r;
            this.i = i;
        }
    }

    public struct cmplxF
    {
        public float r;
        public float i;

        public cmplxF(float r, float i)
        {
            this.r = r;
            this.i = i;
        }
    }
}
