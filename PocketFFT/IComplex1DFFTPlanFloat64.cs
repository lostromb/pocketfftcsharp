﻿using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    public interface IComplex1DFFTPlanFloat64 : I1DFFTPlan
    {
        /// <summary>
        /// Performs a forward Fourier transform on a complex input array.
        /// </summary>
        /// <param name="c">The array of complex numbers to be transformed.</param>
        /// <param name="fct">A scaling parameter to scale the results of the transform.</param>
        void Forward(Span<cmplx> c, double fct);

        /// <summary>
        /// Performs a reverse Fourier transform on a complex input array.
        /// </summary>
        /// <param name="c">The array of complex numbers to be transformed.</param>
        /// <param name="fct">A scaling parameter to scale the results of the transform.</param>
        void Backward(Span<cmplx> c, double fct);
    }
}