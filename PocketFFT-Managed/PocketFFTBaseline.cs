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
    public static class PocketFFTBaseline
    {
        
        public static IComplexFFTPlan make_cfft_plan(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new cfftp_plan(length);
            }

            double comp1 = Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */

            if (comp2 < comp1) // use Bluestein
            {
                return new fftblue_plan(length);
            }
            else
            {
                return new cfftp_plan(length);
            }
        }

        public static IRealFFTPlan make_rfft_plan(int length)
        {
            if (length == 0)
            {
                throw new ArgumentOutOfRangeException("FFT length must be greater than zero");
            }

            if ((length < 50) || (Intrinsics.largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                return new rfftp_plan(length);
            }

            double comp1 = 0.5 * Intrinsics.cost_guess(length);
            double comp2 = 2 * Intrinsics.cost_guess(Intrinsics.good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
            if (comp2 < comp1) // use Bluestein
            {
                return new fftblue_plan(length);
            }
            else
            {
                return new rfftp_plan(length);
            }
        }
    }
}