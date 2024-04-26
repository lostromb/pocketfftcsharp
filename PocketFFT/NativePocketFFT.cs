using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
#if NETCOREAPP
    internal unsafe static class NativePocketFFT
    {
        private const string DLL_NAME = "libpocketfft";

        [DllImport(DLL_NAME)]
        internal static extern IntPtr make_cfft_plan(int length);

        [DllImport(DLL_NAME)]
        internal static extern void destroy_cfft_plan(IntPtr plan);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_backward(IntPtr plan, cmplx* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_forward(IntPtr plan, cmplx* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int cfft_length(IntPtr plan);


        [DllImport(DLL_NAME)]
        internal static extern IntPtr make_rfft_plan(int length);

        [DllImport(DLL_NAME)]
        internal static extern void destroy_rfft_plan(IntPtr plan);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_backward(IntPtr plan, double* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_forward(IntPtr plan, double* c, double fct);

        [DllImport(DLL_NAME)]
        internal static extern int rfft_length(IntPtr plan);
    }
#endif
}
