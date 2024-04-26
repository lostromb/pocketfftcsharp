using System.Runtime.InteropServices.ComTypes;

namespace PocketFFT
{
    using System;
    using System.Buffers.Binary;
    using System.Collections;
    using System.Collections.Generic;
    using System.Runtime.InteropServices;

    // https://gitlab.mpcdf.mpg.de/mtr/pocketfft

    // compare also https://www.shoup.net/PGFFT/
    // https://bitbucket.org/jpommier/pffft/src/master/
    public static class FFTPlanFactory
    {
        public static IComplexFFTPlan Create1DComplexFFTPlan(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new NativeComplexFFTPlan(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new ComplexFFTPackedPlan(length);
            }

            double comp1 = Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */

            if (comp2 < comp1) // use Bluestein
            {
                return new FFTBluesteinPlan(length);
            }
            else
            {
                return new ComplexFFTPackedPlan(length);
            }
        }

        public static IRealFFTPlan Create1DRealFFTPlan(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new NativeRealFFTPlan(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new RealFFTPackedPlan(length);
            }

            double comp1 = 0.5 * Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
            if (comp2 < comp1) // use Bluestein
            {
                return new FFTBluesteinPlan(length);
            }
            else
            {
                return new RealFFTPackedPlan(length);
            }
        }
    }
}