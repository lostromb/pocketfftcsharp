using BenchmarkDotNet.Running;

namespace Driver
{
    internal static class Program
    {
        public static void Main(string[] args)
        {
            BenchmarkRunner.Run<Benchmarks>();
        }
    }
}
