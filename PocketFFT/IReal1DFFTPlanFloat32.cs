using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    public interface IReal1DFFTPlanFloat32 : I1DFFTPlan
    {
        /// <summary>
        /// Performs a forward Fourier transform on a real input array. I do not know what the intended output is.
        /// </summary>
        /// <param name="c">The array of real numbers to be transformed.</param>
        /// <param name="fct">A scaling parameter to scale the results of the transform.</param>
        void Forward(Span<float> c, float fct);

        /// <summary>
        /// Performs a forward Fourier transform on a real input array. I do not know what the intended output is.
        /// </summary>
        /// <param name="c">The array of real numbers to be transformed.</param>
        /// <param name="fct">A scaling parameter to scale the results of the transform.</param>
        void Backward(Span<float> c, float fct);
    }
}
