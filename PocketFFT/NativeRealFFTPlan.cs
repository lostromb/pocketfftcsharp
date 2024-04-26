using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
#if NETCOREAPP
    public class NativeRealFFTPlan : IRealFFTPlan
    {
        private int _length;
        private IntPtr _handle;

        public NativeRealFFTPlan(int length)
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
            NativePocketFFT.destroy_rfft_plan(_handle);
        }

        public unsafe void Forward(Span<double> c, double fct)
        {
            fixed (double* ptr = c)
            {
                NativePocketFFT.rfft_forward(_handle, ptr, fct);
            }
        }
    }
#endif
}
