using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    public interface IRealFFTPlan : IFFTPlan
    {
        void Forward(Span<double> c, double fct);
        void Backward(Span<double> c, double fct);
    }
}
