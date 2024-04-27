using System;
using System.Collections.Generic;
using System.Text;

namespace PocketFFT
{
    /// <summary>
    /// Represents a precomputed state engine for running 1-dimensional Fourier transforms.
    /// </summary>
    public interface I1DFFTPlan : IDisposable
    {
        /// <summary>
        /// Gets the length of the transform which this plan operates on.
        /// </summary>
        int Length { get; }
    }
}
