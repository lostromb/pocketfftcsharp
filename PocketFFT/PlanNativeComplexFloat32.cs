#if NETCOREAPP

namespace PocketFFT
{
    using Microsoft.Win32.SafeHandles;
    using System;
    using System.Buffers;
    using System.Collections.Generic;
    using System.Runtime.InteropServices;
    using System.Text;

    public class PlanNativeComplexFloat32 : IComplex1DFFTPlanFloat32
    {
        private int _length;
        private NativePocketFFT.NativeCFFTPlanHandle _handle;

        public PlanNativeComplexFloat32(int length)
        {
            _length = length;
            _handle = NativePocketFFT.make_cfft_plan(length);
        }

        public int Length => _length;

        public unsafe void Forward(Span<cmplxF> c, float fct)
        {
            // literally just faster to convert float32 -> float64 and run it through the native code
            cmplx[] scratch = ArrayPool<cmplx>.Shared.Rent(c.Length);
            Span<float> singleSpan = MemoryMarshal.Cast<cmplxF, float>(c);
            Span<double> doubleSpan = MemoryMarshal.Cast<cmplx, double>(scratch.AsSpan(0, c.Length));
            Intrinsics.CastSingleToDouble(singleSpan, doubleSpan);

            fixed (cmplx* ptr = scratch)
            {
                NativePocketFFT.cfft_forward(_handle, ptr, fct);
            }

            Intrinsics.CastDoubleToSingle(doubleSpan, singleSpan);
            //for (int i = 0; i < c.Length; i++)
            //{
            //    c[i] = new cmplxF((float)scratch[i].r, (float)scratch[i].i);
            //}

            ArrayPool<cmplx>.Shared.Return(scratch);
        }

        public unsafe void Backward(Span<cmplxF> c, float fct)
        {
            cmplx[] scratch = ArrayPool<cmplx>.Shared.Rent(c.Length);
            Span<float> singleSpan = MemoryMarshal.Cast<cmplxF, float>(c);
            Span<double> doubleSpan = MemoryMarshal.Cast<cmplx, double>(scratch.AsSpan(0, c.Length));
            Intrinsics.CastSingleToDouble(singleSpan, doubleSpan);

            fixed (cmplx* ptr = scratch)
            {
                NativePocketFFT.cfft_backward(_handle, ptr, fct);
            }

            Intrinsics.CastDoubleToSingle(doubleSpan, singleSpan);
            //for (int i = 0; i < c.Length; i++)
            //{
            //    c[i] = new cmplxF((float)scratch[i].r, (float)scratch[i].i);
            //}

            ArrayPool<cmplx>.Shared.Return(scratch);
        }

        public void Dispose()
        {
            _handle.Dispose();
        }
    }
}

#endif