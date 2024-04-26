#if NETCOREAPP

namespace PocketFFT
{
    using Microsoft.Win32.SafeHandles;
    using System;
    using System.Buffers;
    using System.Collections.Generic;
    using System.Numerics;
    using System.Reflection.Metadata;
    using System.Runtime.InteropServices;
    using System.Text;

    public class PlanNativeRealFloat32 : IReal1DFFTPlanFloat32
    {
        private int _length;
        private NativePocketFFT.NativeRFFTPlanHandle _handle;

        public PlanNativeRealFloat32(int length)
        {
            _length = length;
            _handle = NativePocketFFT.make_rfft_plan(length);
        }

        public int Length => _length;

        public unsafe void Forward(Span<float> c, float fct)
        {
            // literally just faster to convert float32 -> float64 and run it through the native code
            double[] scratch = ArrayPool<double>.Shared.Rent(c.Length);
            Intrinsics.CastSingleToDouble(c, scratch.AsSpan(0, c.Length));

            fixed (double* ptr = scratch)
            {
                NativePocketFFT.rfft_forward(_handle, ptr, fct);
            }

            Intrinsics.CastDoubleToSingle(scratch.AsSpan(0, c.Length), c);
            ArrayPool<double>.Shared.Return(scratch);
        }

        public unsafe void Backward(Span<float> c, float fct)
        {
            double[] scratch = ArrayPool<double>.Shared.Rent(c.Length);
            Intrinsics.CastSingleToDouble(c, scratch.AsSpan(0, c.Length));

            fixed (double* ptr = scratch)
            {
                NativePocketFFT.rfft_backward(_handle, ptr, fct);
            }

            Intrinsics.CastDoubleToSingle(scratch.AsSpan(0, c.Length), c);
            ArrayPool<double>.Shared.Return(scratch);
        }

        public void Dispose()
        {
            _handle.Dispose();
        }
    }
}

#endif