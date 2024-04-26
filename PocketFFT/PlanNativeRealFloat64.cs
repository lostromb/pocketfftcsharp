#if NETCOREAPP

namespace PocketFFT
{
    using Microsoft.Win32.SafeHandles;
    using System;
    using System.Collections.Generic;
    using System.Reflection.Metadata;
    using System.Runtime.InteropServices;
    using System.Text;

    public class PlanNativeRealFloat64 : IReal1DFFTPlanFloat64
    {
        private int _length;
        private NativePocketFFT.NativeRFFTPlanHandle _handle;

        public PlanNativeRealFloat64(int length)
        {
            _length = length;
            _handle = NativePocketFFT.make_rfft_plan(length);
        }

        public int Length => _length;

        public unsafe void Backward(Span<double> c, double fct)
        {
            fixed (double* ptr = c)
            {
                NativePocketFFT.rfft_backward(_handle, ptr, fct);
            }
        }

        public void Dispose()
        {
            _handle.Dispose();
        }

        public unsafe void Forward(Span<double> c, double fct)
        {
            fixed (double* ptr = c)
            {
                NativePocketFFT.rfft_forward(_handle, ptr, fct);
            }
        }
    }
}

#endif