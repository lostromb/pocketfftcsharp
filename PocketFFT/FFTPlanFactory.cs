namespace PocketFFT
{
    using System;

    // https://gitlab.mpcdf.mpg.de/mtr/pocketfft

    // compare also https://www.shoup.net/PGFFT/
    // https://bitbucket.org/jpommier/pffft/src/master/
    public static class FFTPlanFactory
    {
        public static IComplex1DFFTPlanFloat32 Create1DComplexFFTPlanFloat32(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new PlanNativeComplexFloat32(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new PlanComplexPackedFloat32(length);
            }

            double comp1 = Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */

            if (comp2 < comp1) // use Bluestein
            {
                return new PlanBluesteinFloat32(length);
            }
            else
            {
                return new PlanComplexPackedFloat32(length);
            }
        }

        public static IComplex1DFFTPlanFloat64 Create1DComplexFFTPlanFloat64(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new PlanNativeComplexFloat64(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new PlanComplexPackedFloat64(length);
            }

            double comp1 = Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */

            if (comp2 < comp1) // use Bluestein
            {
                return new PlanBluesteinFloat64(length);
            }
            else
            {
                return new PlanComplexPackedFloat64(length);
            }
        }

        public static IReal1DFFTPlanFloat32 Create1DRealFFTPlanFloat32(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new PlanNativeRealFloat32(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new PlanRealPackedFloat32(length);
            }

            double comp1 = 0.5 * Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
            if (comp2 < comp1) // use Bluestein
            {
                return new PlanBluesteinFloat32(length);
            }
            else
            {
                return new PlanRealPackedFloat32(length);
            }
        }

        public static IReal1DFFTPlanFloat64 Create1DRealFFTPlanFloat64(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

#if NETCOREAPP
            return new PlanNativeRealFloat64(length);
#endif

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new PlanRealPackedFloat64(length);
            }

            double comp1 = 0.5 * Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
            if (comp2 < comp1) // use Bluestein
            {
                return new PlanBluesteinFloat64(length);
            }
            else
            {
                return new PlanRealPackedFloat64(length);
            }
        }
    }
}