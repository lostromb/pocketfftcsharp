using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
#if NETCOREAPP
    public class NativeComplexFFTPlan : IComplexFFTPlan
    {
        private int _length;
        private IntPtr _handle;

        public NativeComplexFFTPlan(int length)
        {
            _length = length;
            _handle = NativePocketFFT.make_cfft_plan(length);
        }

        public int Length => _length;

        public unsafe void Backward(Span<cmplx> c, double fct)
        {
            fixed (cmplx* ptr = c)
            {
                NativePocketFFT.cfft_backward(_handle, ptr, fct);
            }
        }

        public void Dispose()
        {
            NativePocketFFT.destroy_cfft_plan(_handle);
        }

        public unsafe void Forward(Span<cmplx> c, double fct)
        {
            fixed (cmplx* ptr = c)
            {
                NativePocketFFT.cfft_forward(_handle, ptr, fct);
            }
        }
    }
#endif
}
