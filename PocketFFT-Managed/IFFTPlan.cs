using System;
using System.Collections.Generic;
using System.Text;

namespace PocketFFT
{
    public interface IFFTPlan : IDisposable
    {
        int Length { get; }
    }
}
