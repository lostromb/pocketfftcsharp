using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using static PocketFFT.PocketFFTBaseline;

namespace PocketFFT
{
    internal static class NailTest
    {
        public static void PrintComplexArray(ReadOnlySpan<cmplx> val)
        {
            for (int c = 0; c < val.Length; c++)
            {
                PrintComplex(val[c]);
            }
        }

        public static void PrintComplex(cmplx val)
        {
#if NET8_0_OR_GREATER
            Console.WriteLine(string.Format("0x{0:X16}r|0x{1:X16}i",
                new DoubleLayout(val.r).LVal,
                new DoubleLayout(val.i).LVal).ToLowerInvariant());
#endif
        }

        public static void PrintDoubleArray(ReadOnlySpan<double> val)
        {
#if NET8_0_OR_GREATER
            int col = 0;
            for (int c = 0; c < val.Length; c++)
            {
                PrintDouble(val[c]);
                if (c != (val.Length - 1))
                {
                    Console.Write(", ");
                }
                if (++col > 4)
                {
                    Console.WriteLine();
                    col = 0;
                }
            }

           Console.WriteLine();
#endif
        }

        public static void PrintDouble(double val)
        {
#if NET8_0_OR_GREATER
            Console.Write(string.Format("0x{0:X}", new DoubleLayout(val).LVal));
#endif
        }

        [StructLayout(LayoutKind.Explicit)]
        private struct DoubleLayout
        {
            [FieldOffset(0)]
            public double DVal;

            [FieldOffset(0)]
            public long LVal;

            public DoubleLayout(double d)
            {
                LVal = 0;
                DVal = d;
            }
        }
    }
}
