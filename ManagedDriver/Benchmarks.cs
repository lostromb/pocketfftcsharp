using BenchmarkDotNet.Attributes;
using PocketFFT;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using static Driver.Benchmarks;

namespace Driver
{
    [MemoryDiagnoser]
    public class Benchmarks
    {
        private IReal1DFFTPlanFloat64 _real64;
        private IReal1DFFTPlanFloat32 _real32;
        private IComplex1DFFTPlanFloat64 _complex64;
        private IComplex1DFFTPlanFloat32 _complex32;

        private double[] _realInput64;
        private float[] _realInput32;
        private cmplx[] _complexInput64;
        private cmplxF[] _complexInput32;

        [Params(64, 191, 1000, 2310, 5980)]
        public int TransformLength { get; set; }

        [GlobalSetup]
        public void Setup()
        {
            _realInput64 = new double[TransformLength * 2];
            _real64 = FFTPlanFactory.Create1DRealFFTPlanFloat64(TransformLength);
            _realInput32 = new float[TransformLength * 2];
            _real32 = FFTPlanFactory.Create1DRealFFTPlanFloat32(TransformLength);
            _complexInput64 = new cmplx[TransformLength];
            _complex64 = FFTPlanFactory.Create1DComplexFFTPlanFloat64(TransformLength);
            _complexInput32 = new cmplxF[TransformLength];
            _complex32 = FFTPlanFactory.Create1DComplexFFTPlanFloat32(TransformLength);

            for (int c = 0; c < _realInput64.Length; c++)
            {
                _realInput64[c] = Random.Shared.NextDouble() - 0.5;
                _realInput32[c] = (float)_realInput64[c];
            }

            for (int c = 0; c < _complexInput64.Length; c++)
            {
                _complexInput64[c].r = Random.Shared.NextDouble() - 0.5;
                _complexInput64[c].i = Random.Shared.NextDouble() - 0.5;
                _complexInput32[c].r = (float)_complexInput64[c].r;
                _complexInput32[c].i = (float)_complexInput64[c].i;
            }
        }

        [Benchmark]
        public void Real32()
        {
            _real32.Forward(_realInput32.AsSpan(), 3.0f);
            _real32.Backward(_realInput32.AsSpan(), 3.0f);
        }

        [Benchmark]
        public void Real64()
        {
            _real64.Forward(_realInput64.AsSpan(), 3.0);
            _real64.Backward(_realInput64.AsSpan(), 3.0);
        }

        [Benchmark]
        public void Complex32()
        {
            _complex32.Forward(_complexInput32.AsSpan(), 3.0f);
            _complex32.Backward(_complexInput32.AsSpan(), 3.0f);
        }

        [Benchmark]
        public void Complex64()
        {
            _complex64.Forward(_complexInput64.AsSpan(), 3.0);
            _complex64.Backward(_complexInput64.AsSpan(), 3.0);
        }
    }
}
