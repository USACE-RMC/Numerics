using BenchmarkDotNet.Attributes;
using Numerics.Distributions.Copulas;
using Numerics.Mathematics.RootFinding;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The Archimedean copula generator members against the 2.2.0 members (kept as references), each over 64 points of (0, 1)
    /// (or the matching generator and slope values) at a representative θ.
    /// </summary>
    /// <remarks>
    /// The 2.2.0 members are copied into the reference classes below, each reading θ from an instance field as the 2.2.0
    /// members read the copula's Theta property, so the JIT cannot fold θ into the reference formulas. Each inverse is
    /// evaluated at the image of the same 64 points under its own forward member: the 2.2.0 Clayton generator lacks
    /// Nelsen's factor 1/θ, so its values differ from the library generator's.
    /// </remarks>
    [ShortRunJob]
    public class ArchimedeanGeneratorBenchmarks
    {
        private readonly double[] _t = new double[64];
        private double[] _claytonPhi = null!, _claytonSlope = null!, _frankPhi = null!, _frankSlope = null!;
        private double[] _gumbelPhi = null!, _gumbelSlope = null!, _joePhi = null!, _joeSlope = null!, _amhPhi = null!, _amhSlope = null!;
        private double[] _claytonPhiPrevious = null!, _claytonSlopePrevious = null!, _frankPhiPrevious = null!, _frankSlopePrevious = null!;
        private double[] _gumbelPhiPrevious = null!, _gumbelSlopePrevious = null!, _joePhiPrevious = null!, _joeSlopePrevious = null!;
        private double[] _amhPhiPrevious = null!, _amhSlopePrevious = null!;
        private readonly ClaytonCopula _clayton = new ClaytonCopula(2d);
        private readonly FrankCopula _frank = new FrankCopula(5d);
        private readonly GumbelCopula _gumbel = new GumbelCopula(2d);
        private readonly JoeCopula _joe = new JoeCopula(2d);
        private readonly AMHCopula _amh = new AMHCopula(0.5);
        private readonly PreviousClayton _claytonPrevious = new PreviousClayton(2d);
        private readonly PreviousFrank _frankPrevious = new PreviousFrank(5d);
        private readonly PreviousGumbel _gumbelPrevious = new PreviousGumbel(2d);
        private readonly PreviousJoe _joePrevious = new PreviousJoe(2d);
        private readonly PreviousAmh _amhPrevious = new PreviousAmh(0.5);

        /// <summary>Fills the evaluation points and the generator and slope values at them.</summary>
        [GlobalSetup]
        public void Setup()
        {
            for (int i = 0; i < _t.Length; i++) _t[i] = (i + 0.5) / _t.Length;
            _claytonPhi = Map(_clayton.Generator); _claytonSlope = Map(_clayton.GeneratorPrime);
            _frankPhi = Map(_frank.Generator); _frankSlope = Map(_frank.GeneratorPrime);
            _gumbelPhi = Map(_gumbel.Generator); _gumbelSlope = Map(_gumbel.GeneratorPrime);
            _joePhi = Map(_joe.Generator); _joeSlope = Map(_joe.GeneratorPrime);
            _amhPhi = Map(_amh.Generator); _amhSlope = Map(_amh.GeneratorPrime);
            _claytonPhiPrevious = Map(_claytonPrevious.Generator); _claytonSlopePrevious = Map(_claytonPrevious.GeneratorPrime);
            _frankPhiPrevious = Map(_frankPrevious.Generator); _frankSlopePrevious = Map(_frankPrevious.GeneratorPrime);
            _gumbelPhiPrevious = Map(_gumbelPrevious.Generator); _gumbelSlopePrevious = Map(_gumbelPrevious.GeneratorPrime);
            _joePhiPrevious = Map(_joePrevious.Generator); _joeSlopePrevious = Map(_joePrevious.GeneratorPrime);
            _amhPhiPrevious = Map(_amhPrevious.Generator); _amhSlopePrevious = Map(_amhPrevious.GeneratorPrime);
        }

        private double[] Map(Func<double, double> f)
        {
            var values = new double[_t.Length];
            for (int i = 0; i < _t.Length; i++) values[i] = f(_t[i]);
            return values;
        }

        private static double Sum(double[] x, Func<double, double> f)
        {
            double s = 0d;
            for (int i = 0; i < x.Length; i++) s += f(x[i]);
            return s;
        }

        // Clayton, θ = 2.
        /// <summary>Previous Clayton generator.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonGeneratorPrevious() => Sum(_t, _claytonPrevious.Generator);
        /// <summary>Clayton generator.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonGenerator() => Sum(_t, _clayton.Generator);
        /// <summary>Previous Clayton inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonInversePrevious() => Sum(_claytonPhiPrevious, _claytonPrevious.GeneratorInverse);
        /// <summary>Clayton inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonInverse() => Sum(_claytonPhi, _clayton.GeneratorInverse);
        /// <summary>Previous Clayton slope.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrimePrevious() => Sum(_t, _claytonPrevious.GeneratorPrime);
        /// <summary>Clayton slope.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrime() => Sum(_t, _clayton.GeneratorPrime);
        /// <summary>Previous Clayton second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrime2Previous() => Sum(_t, _claytonPrevious.GeneratorPrime2);
        /// <summary>Clayton second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrime2() => Sum(_t, _clayton.GeneratorPrime2);
        /// <summary>Previous Clayton slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrimeInversePrevious() => Sum(_claytonSlopePrevious, _claytonPrevious.GeneratorPrimeInverse);
        /// <summary>Clayton slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double ClaytonPrimeInverse() => Sum(_claytonSlope, _clayton.GeneratorPrimeInverse);

        // Frank, θ = 5.
        /// <summary>Previous Frank generator.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankGeneratorPrevious() => Sum(_t, _frankPrevious.Generator);
        /// <summary>Frank generator.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankGenerator() => Sum(_t, _frank.Generator);
        /// <summary>Previous Frank inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankInversePrevious() => Sum(_frankPhiPrevious, _frankPrevious.GeneratorInverse);
        /// <summary>Frank inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankInverse() => Sum(_frankPhi, _frank.GeneratorInverse);
        /// <summary>Previous Frank slope.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrimePrevious() => Sum(_t, _frankPrevious.GeneratorPrime);
        /// <summary>Frank slope.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrime() => Sum(_t, _frank.GeneratorPrime);
        /// <summary>Previous Frank second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrime2Previous() => Sum(_t, _frankPrevious.GeneratorPrime2);
        /// <summary>Frank second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrime2() => Sum(_t, _frank.GeneratorPrime2);
        /// <summary>Previous Frank slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrimeInversePrevious() => Sum(_frankSlopePrevious, _frankPrevious.GeneratorPrimeInverse);
        /// <summary>Frank slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double FrankPrimeInverse() => Sum(_frankSlope, _frank.GeneratorPrimeInverse);

        // Gumbel, θ = 2.
        /// <summary>Previous Gumbel generator.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelGeneratorPrevious() => Sum(_t, _gumbelPrevious.Generator);
        /// <summary>Gumbel generator.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelGenerator() => Sum(_t, _gumbel.Generator);
        /// <summary>Previous Gumbel inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelInversePrevious() => Sum(_gumbelPhiPrevious, _gumbelPrevious.GeneratorInverse);
        /// <summary>Gumbel inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelInverse() => Sum(_gumbelPhi, _gumbel.GeneratorInverse);
        /// <summary>Previous Gumbel slope.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrimePrevious() => Sum(_t, _gumbelPrevious.GeneratorPrime);
        /// <summary>Gumbel slope.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrime() => Sum(_t, _gumbel.GeneratorPrime);
        /// <summary>Previous Gumbel second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrime2Previous() => Sum(_t, _gumbelPrevious.GeneratorPrime2);
        /// <summary>Gumbel second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrime2() => Sum(_t, _gumbel.GeneratorPrime2);
        /// <summary>Previous Gumbel slope inverse (Brent, tolerance 1E-8).</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrimeInversePrevious() => Sum(_gumbelSlopePrevious, _gumbelPrevious.GeneratorPrimeInverse);
        /// <summary>Gumbel slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double GumbelPrimeInverse() => Sum(_gumbelSlope, _gumbel.GeneratorPrimeInverse);

        // Joe, θ = 2.
        /// <summary>Previous Joe generator.</summary><returns>The sum.</returns>
        [Benchmark] public double JoeGeneratorPrevious() => Sum(_t, _joePrevious.Generator);
        /// <summary>Joe generator.</summary><returns>The sum.</returns>
        [Benchmark] public double JoeGenerator() => Sum(_t, _joe.Generator);
        /// <summary>Previous Joe inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double JoeInversePrevious() => Sum(_joePhiPrevious, _joePrevious.GeneratorInverse);
        /// <summary>Joe inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double JoeInverse() => Sum(_joePhi, _joe.GeneratorInverse);
        /// <summary>Previous Joe slope.</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrimePrevious() => Sum(_t, _joePrevious.GeneratorPrime);
        /// <summary>Joe slope.</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrime() => Sum(_t, _joe.GeneratorPrime);
        /// <summary>Previous Joe second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrime2Previous() => Sum(_t, _joePrevious.GeneratorPrime2);
        /// <summary>Joe second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrime2() => Sum(_t, _joe.GeneratorPrime2);
        /// <summary>Previous Joe slope inverse (Brent, tolerance 1E-8).</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrimeInversePrevious() => Sum(_joeSlopePrevious, _joePrevious.GeneratorPrimeInverse);
        /// <summary>Joe slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double JoePrimeInverse() => Sum(_joeSlope, _joe.GeneratorPrimeInverse);

        // AMH, θ = 0.5.
        /// <summary>Previous AMH generator.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhGeneratorPrevious() => Sum(_t, _amhPrevious.Generator);
        /// <summary>AMH generator.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhGenerator() => Sum(_t, _amh.Generator);
        /// <summary>Previous AMH inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhInversePrevious() => Sum(_amhPhiPrevious, _amhPrevious.GeneratorInverse);
        /// <summary>AMH inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhInverse() => Sum(_amhPhi, _amh.GeneratorInverse);
        /// <summary>Previous AMH slope.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrimePrevious() => Sum(_t, _amhPrevious.GeneratorPrime);
        /// <summary>AMH slope.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrime() => Sum(_t, _amh.GeneratorPrime);
        /// <summary>Previous AMH second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrime2Previous() => Sum(_t, _amhPrevious.GeneratorPrime2);
        /// <summary>AMH second derivative.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrime2() => Sum(_t, _amh.GeneratorPrime2);
        /// <summary>Previous AMH slope inverse (Brent, tolerance 1E-8).</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrimeInversePrevious() => Sum(_amhSlopePrevious, _amhPrevious.GeneratorPrimeInverse);
        /// <summary>AMH slope inverse.</summary><returns>The sum.</returns>
        [Benchmark] public double AmhPrimeInverse() => Sum(_amhSlope, _amh.GeneratorPrimeInverse);

        /// <summary>The 2.2.0 Clayton generator members.</summary>
        private sealed class PreviousClayton(double theta)
        {
            private readonly double _theta = theta;

            public double Generator(double t)
            {
                double a = t;
                return Math.Sign(a) * Math.Pow(Math.Abs(a), -_theta) - 1.0d;
            }

            public double GeneratorInverse(double t)
            {
                double a = 1.0d + t;
                return Math.Sign(a) * Math.Pow(Math.Abs(a), -1.0d / _theta);
            }

            public double GeneratorPrime(double t)
            {
                double a = t;
                return -_theta * Math.Sign(a) * Math.Pow(Math.Abs(a), -_theta - 1.0d);
            }

            public double GeneratorPrime2(double t)
            {
                double a = t;
                return -_theta * Math.Sign(a) * Math.Pow(Math.Abs(a), -_theta - 2.0d) * (-_theta - 1.0d);
            }

            public double GeneratorPrimeInverse(double t)
            {
                double a = t / -_theta;
                return Math.Sign(a) * Math.Pow(Math.Abs(a), -1.0d / (_theta + 1.0d));
            }
        }

        /// <summary>The 2.2.0 Frank generator members.</summary>
        private sealed class PreviousFrank(double theta)
        {
            private readonly double _theta = theta;

            public double Generator(double t)
            {
                if (Math.Abs(_theta) < 1e-10) return t;
                return -Math.Log((Math.Exp(-_theta * t) - 1d) / (Math.Exp(-_theta) - 1d));
            }

            public double GeneratorInverse(double t) => -Math.Log(Math.Exp(-_theta - t) - Math.Exp(-t) + 1d) / _theta;

            public double GeneratorPrime(double t) => _theta / (1d - Math.Exp(_theta * t));

            public double GeneratorPrime2(double t) => _theta * _theta * Math.Exp(_theta * t) / Math.Pow(1d - Math.Exp(_theta * t), 2d);

            public double GeneratorPrimeInverse(double t) => Math.Log((_theta - t) / -t) / _theta;
        }

        /// <summary>
        /// The 2.2.0 Gumbel generator members. The 2.2.0 slope takes the sign of log t rather than −log t, which makes it
        /// positive on (0, 1) and leaves its Brent inverse without a bracket; the slope here takes the sign of −log t,
        /// which costs the same.
        /// </summary>
        private sealed class PreviousGumbel(double theta)
        {
            private readonly double _theta = theta;

            public double Generator(double t)
            {
                double a = -Math.Log(t);
                return Math.Sign(a) * Math.Pow(Math.Abs(a), _theta);
            }

            public double GeneratorInverse(double t)
            {
                double a = -t;
                return Math.Exp(Math.Sign(a) * Math.Pow(Math.Abs(a), 1.0d / _theta));
            }

            public double GeneratorPrime(double t)
            {
                double a = -Math.Log(t);
                return -_theta * Math.Sign(a) * Math.Pow(Math.Abs(a), _theta - 1.0d) / t;
            }

            public double GeneratorPrime2(double t)
            {
                double a = -Math.Log(t);
                return _theta * Math.Sign(a) * Math.Pow(Math.Abs(a), _theta - 2.0d) * (-_theta + Math.Log(t) + 1.0d) / (t * t);
            }

            public double GeneratorPrimeInverse(double t) => Brent.Solve(x => GeneratorPrime(x) - t, 0d, 1d);
        }

        /// <summary>The 2.2.0 Joe generator members.</summary>
        private sealed class PreviousJoe(double theta)
        {
            private readonly double _theta = theta;

            public double Generator(double t)
            {
                double a = 1.0d - t;
                return -Math.Log(1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), _theta));
            }

            public double GeneratorInverse(double t)
            {
                double a = 1.0d - Math.Exp(-t);
                return 1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), 1.0d / _theta);
            }

            public double GeneratorPrime(double t)
            {
                double a = 1.0d - t;
                return -(_theta * Math.Sign(a) * Math.Pow(Math.Abs(a), _theta - 1.0d)) / (1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), _theta));
            }

            public double GeneratorPrime2(double t)
            {
                double a = 1.0d - t;
                double num = _theta * (_theta + Math.Sign(a) * Math.Pow(Math.Abs(a), _theta) - 1.0d) * Math.Sign(a) * Math.Pow(Math.Abs(a), _theta - 2.0d);
                double aa = 1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), _theta);
                double den = Math.Sign(aa) * Math.Pow(Math.Abs(aa), 2d);
                return num / den;
            }

            public double GeneratorPrimeInverse(double t) => Brent.Solve(x => GeneratorPrime(x) - t, 0d, 1d);
        }

        /// <summary>The 2.2.0 Ali-Mikhail-Haq generator members.</summary>
        private sealed class PreviousAmh(double theta)
        {
            private readonly double _theta = theta;

            public double Generator(double t) => Math.Log((1.0d - _theta * (1.0d - t)) / t);

            public double GeneratorInverse(double t) => (1.0d - _theta) / (Math.Exp(t) - _theta);

            public double GeneratorPrime(double t) => (_theta - 1.0d) / (t * (_theta * (t - 1.0d) + 1.0d));

            public double GeneratorPrime2(double t)
            {
                double num = (_theta - 1.0d) * (_theta * (2.0d * t - 1.0d) + 1.0d);
                double a = _theta * (t - 1.0d) * t + t;
                double den = Math.Sign(a) * Math.Pow(Math.Abs(a), 2d);
                return -num / den;
            }

            public double GeneratorPrimeInverse(double t) => Brent.Solve(x => GeneratorPrime(x) - t, 0d, 1d);
        }
    }
}
