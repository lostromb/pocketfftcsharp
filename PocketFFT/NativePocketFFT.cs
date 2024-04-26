
#if NETCOREAPP

using Microsoft.Win32.SafeHandles;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    internal unsafe static class NativePocketFFT
    {
        private const string DLL_NAME = "libpocketfft";

        internal class NativeCFFTPlanHandle : SafeHandleZeroOrMinusOneIsInvalid
        {
            public NativeCFFTPlanHandle() : base(ownsHandle: true)
            {
            }

            protected override bool ReleaseHandle()
            {
                destroy_cfft_plan(handle);
                return true;
            }
        }

        [DllImport(DLL_NAME)]
        internal static extern NativeCFFTPlanHandle make_cfft_plan(int length);

        [DllImport(DLL_NAME)]
        internal static extern void destroy_cfft_plan(IntPtr plan);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_backward(NativeCFFTPlanHandle plan, cmplx* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_forward(NativeCFFTPlanHandle plan, cmplx* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_length(NativeCFFTPlanHandle plan);



        internal class NativeRFFTPlanHandle : SafeHandleZeroOrMinusOneIsInvalid
        {
            public NativeRFFTPlanHandle() : base(ownsHandle: true)
            {
            }

            protected override bool ReleaseHandle()
            {
                destroy_rfft_plan(handle);
                return true;
            }
        }

        [DllImport(DLL_NAME)]
        internal static extern NativeRFFTPlanHandle make_rfft_plan(int length);

        [DllImport(DLL_NAME)]
        internal static extern void destroy_rfft_plan(IntPtr plan);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_backward(NativeRFFTPlanHandle plan, double* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_forward(NativeRFFTPlanHandle plan, double* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_length(NativeRFFTPlanHandle plan);
    }
}

#endif