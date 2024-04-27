using BenchmarkDotNet.Running;
using PocketFFT;
using System.Numerics;
using System.Text;

namespace Driver
{
    internal static class Program
    {
        public static void Main(string[] args)
        {
            //using (FileStream tsvOut = new FileStream("out.csv", FileMode.Create, FileAccess.Write))
            //using (StreamWriter writer = new StreamWriter(tsvOut))
            //{
            //    double[] data = new double[1000];
            //    writer.WriteLine(string.Join(",", Enumerable.Range(0, 1000)));
            //    //AddNoise(0.01, data.AsSpan());
            //    //AddSineWave(0.2, 0.0, 500, data.AsSpan());
            //    //AddSineWave(0.2, 0.0, 100, data.AsSpan());
            //    //AddSineWave(0.2, 0.0, 20, data.AsSpan());
            //    //AddSineWave(1.0, 0.25, 2, data.AsSpan());
            //    writer.WriteLine(string.Join(",", data));
            //    var realPlan = FFTPlanFactory.Create1DRealFFTPlanFloat64(1000);
            //    realPlan.Forward(data, 1.0);
            //    writer.WriteLine(string.Join(",", data));
            //}

            BenchmarkRunner.Run<Benchmarks>();
        }

        private static void AddSineWave(double amplitude, double phase, double period, Span<double> target)
        {
            for (int c = 0; c < target.Length; c++)
            {
                double thisPhase = (Math.PI * 2) * (((double)c / period) + phase);
                while (thisPhase > (Math.PI * 2))
                {
                    thisPhase -= (Math.PI * 2);
                }

                target[c] += amplitude * Math.Sin(thisPhase);
            }
        }

        private static void AddNoise(double amplitude, Span<double> target)
        {
            for (int c = 0; c < target.Length; c++)
            {
                target[c] += amplitude * (Random.Shared.NextDouble() + -0.5) * 2;
            }
        }
    }
}
