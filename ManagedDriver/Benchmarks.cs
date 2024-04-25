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
        private IRealFFTPlan _real;
        private IComplexFFTPlan _complex;

        private double[] _realInput;
        private cmplx[] _complexInput;

        [Params(64, 191, 1000, 2310, 5980)]
        public int TransformLength { get; set; }

        [GlobalSetup]
        public void Setup()
        {
            _realInput = new double[TransformLength * 2];
            _real = FFTPlanFactory.Create1DRealFFTPlan(TransformLength);
            _complexInput = new cmplx[TransformLength];
            _complex = FFTPlanFactory.Create1DComplexFFTPlan(TransformLength);

            for (int c = 0; c < _realInput.Length; c++)
            {
                _realInput[c] = Random.Shared.NextDouble() - 0.5;
            }

            for (int c = 0; c < _complexInput.Length; c++)
            {
                _complexInput[c].r = Random.Shared.NextDouble() - 0.5;
                _complexInput[c].i = Random.Shared.NextDouble() - 0.5;
            }
        }

        [Benchmark]
        public void Real()
        {
            _real.Forward(_realInput.AsSpan(), 3.0);
            _real.Backward(_realInput.AsSpan(), 3.0);
        }

        //[Benchmark]
        //public void Complex()
        //{
        //    _complex.Forward(_complexInput.AsSpan(), 3.0);
        //    _complex.Backward(_complexInput.AsSpan(), 3.0);
        //}
    }
}
