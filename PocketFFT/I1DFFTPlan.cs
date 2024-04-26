using System;
using System.Collections.Generic;
using System.Text;

namespace PocketFFT
{
    public interface I1DFFTPlan : IDisposable
    {
        int Length { get; }
    }
}
