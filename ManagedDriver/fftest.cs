using PocketFFT;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ManagedDriver
{
    public static class fftest
    {
        private const int maxlen = 8196;

        private static void fill_random(Span<double> data, int length)
        {
            Random rand = new Random();
            for (int m = 0; m < length; ++m)
            {
                data[m] = rand.NextDouble() - 0.5;
            }
        }

        private static void fill_random(Span<cmplx> data, int length)
        {
            Random rand = new Random();
            for (int m = 0; m < length; ++m)
            {
                data[m].r = rand.NextDouble() - 0.5;
                data[m].i = rand.NextDouble() - 0.5;
            }
        }

        private static double errcalc(Span<double> data, Span<double> odata, int length)
        {
            double sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m] - odata[m]) * (data[m] - odata[m]);
                sum += odata[m] * odata[m];
            }

            return Math.Sqrt(errsum / sum);
        }

        private static double errcalc(Span<cmplx> data, Span<cmplx> odata, int length)
        {
            double sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m].r - odata[m].r) * (data[m].r - odata[m].r);
                errsum += (data[m].i - odata[m].i) * (data[m].i - odata[m].i);
                sum += odata[m].r * odata[m].r;
                sum += odata[m].i * odata[m].i;
            }

            return Math.Sqrt(errsum / sum);
        }

        public static int test_real()
        {
            double[] data = new double[2 * maxlen];
            double[] odata = new double[2 * maxlen];
            const double epsilon = 2e-15;
            int ret = 0;
            fill_random(odata, maxlen);
            double errsum = 0;
            for (int length = 1; length <= maxlen; ++length)
            {
                //Console.WriteLine("testing real " + length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IRealFFTPlan plan = FFTPlanFactory.Create1DRealFFTPlan(length))
                {
                    plan.Forward(data, 1.0);
                    plan.Backward(data, 1.0 / length);
                }

                double err = errcalc(data, odata, length);
                if (err > epsilon)
                {
                    Console.WriteLine($"problem at real length {length}: {err}");
                    ret = 1;
                }
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
            return ret;
        }

        public static int test_complex()
        {
            cmplx[] data = new cmplx[maxlen];
            cmplx[] odata = new cmplx[maxlen];
            fill_random(odata, maxlen);
            const double epsilon = 2e-15;
            int ret = 0;
            double errsum = 0;

            for (int length = 1; length <= maxlen; ++length)
            {
                //Console.WriteLine("testing cmplx " + length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IComplexFFTPlan plan = FFTPlanFactory.Create1DComplexFFTPlan(length))
                {
                    plan.Forward(data, 1.0);
                    plan.Backward(data, 1.0 / length);
                }
                double err = errcalc(data, odata, length);
                if (err > epsilon)
                {
                    Console.WriteLine($"problem at complex length {length}: {err}");
                    ret = 1;
                }
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
            return ret;
        }

        public static void TestReal(int length)
        {
            double[] data = new double[length];
            double[] odata = new double[length];
            fill_random(odata, length);
            odata.AsSpan(0, length).CopyTo(data);
            using (IRealFFTPlan plan = FFTPlanFactory.Create1DRealFFTPlan(length))
            {
                plan.Forward(data, 1.0);
                plan.Backward(data, 1.0 / length);
            }

            double err = errcalc(data, odata, length);
            Console.WriteLine($"err: {err}");
        }
    }
}
