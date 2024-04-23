using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PocketFFT
{
    public interface IComplexFFTPlan : IFFTPlan
    {
        void Forward(Span<cmplx> c, double fct);
        void Backward(Span<cmplx> c, double fct);
    }
}
