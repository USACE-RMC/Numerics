using System;
using System.Globalization;
using System.IO;

namespace Numerics.Sampling
{
    /// <summary>
    /// A class for generating a Sobol sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// A Sobol sequence is a low-discrepancy sequence with the property that for all values of N,
    /// its subsequence (x1, ... xN) has a low discrepancy. It can be used to generate pseudo-random
    /// points in a space S, which are equi-distributed.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> The implementation already comes with support for up to 21201 dimensions with direction numbers
    /// calculated from <see href="http://web.maths.unsw.edu.au/~fkuo/sobol/" />  </item>
    /// <item> This code was converted from the Apache Math Commons.  
    /// <see href = "https://commons.apache.org/proper/commons-math/apidocs/src-html/org/apache/commons/math4/random/SobolSequenceGenerator.html" /> </item>
    /// <item> <see href = "http://en.wikipedia.org/wiki/Sobol_sequence" /> </item>
    /// <item> <see href = "http://web.maths.unsw.edu.au/~fkuo/sobol/" /> </item>
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017. </item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class SobolSequence
    {
        /// <summary>
        /// Constructs a new Sobol Sequence.
        /// </summary>
        /// <param name="dimension">Optional. The spatial dimension. Default = 1.</param>
        /// <exception cref="ArgumentException"></exception>
        public SobolSequence(int dimension = 1)
        {
            if (dimension < 1 || dimension > MAX_DIMENSION)
            {
                throw new ArgumentException("The dimension must be between 1 and " + MAX_DIMENSION);
            }
            Dimension = dimension;
            _direction = new long[dimension, BITS + 1];
            x = new long[dimension];
            _shift = new long[dimension];

            initialize();
        }

        /// <summary>
        /// Constructs a new randomized (scrambled) Sobol Sequence with a reproducible seed.
        /// </summary>
        /// <param name="dimension">The spatial dimension.</param>
        /// <param name="seed">The pseudorandom seed for the scrambling.</param>
        /// <exception cref="ArgumentException"></exception>
        /// <remarks>
        /// <para>
        /// The randomization applies Matousek's linear matrix scrambling followed by a random
        /// digital shift - the practical realization of Owen-style scrambling. For each dimension a
        /// random unit lower-triangular bit matrix is applied once to the direction numbers (the
        /// generator recursion is linear over GF(2), so pre-scrambling the direction numbers
        /// scrambles every generated point exactly), and a random 52-bit digital shift is applied to
        /// each emitted coordinate. Both preserve the sequence's dyadic equidistribution, so each
        /// scrambled dimension remains a base-2 (0,1)-sequence.
        /// </para>
        /// <para>
        /// The same seed always reproduces the same sequence. The seeded draw order is part of that
        /// contract: for each dimension in order, the sub-diagonal matrix bits row by row, then the
        /// shift bits. The parameterless-seed constructor generates the original unrandomized
        /// sequence and is unaffected.
        /// </para>
        /// <b> References: </b>
        /// <list type="bullet">
        /// <item>Matousek, J. (1998). On the L2-discrepancy for anchored boxes. Journal of Complexity, 14(4), 527-556.</item>
        /// <item>Owen, A. B. (2021). On dropping the first Sobol' point. arXiv:2008.08051.</item>
        /// </list>
        /// </remarks>
        public SobolSequence(int dimension, int seed) : this(dimension)
        {
            Seed = seed;
            var prng = new MersenneTwister(seed);
            ScrambleDirections(prng);
        }

        /// <summary>
        /// The number of bits to use. 
        /// </summary>
        private static int BITS = 52;

        /// <summary>
        /// The scaling factor.
        /// </summary>
        private static double SCALE = Math.Pow(2, BITS);

        /// <summary>
        /// The maximum supported space dimension. 
        /// </summary>
        private static int MAX_DIMENSION = 21201;

        /// <summary>
        /// The current index in the sequence.
        /// </summary>
        private int count;

        /// <summary>
        /// Space dimension.
        /// </summary>
        public int Dimension { get; private set; }

        /// <summary>
        /// The pseudorandom seed of the scrambling, or null for the original unrandomized sequence.
        /// </summary>
        public int? Seed { get; }

        /// <summary>
        /// The direction vector for each component.
        /// </summary>
        private long[,] _direction;

        /// <summary>
        /// The current state.
        /// </summary>
        private long[] x;

        /// <summary>
        /// The per-dimension digital shift applied at emission; all zero for the unrandomized sequence.
        /// </summary>
        private readonly long[] _shift;

        /// <summary>
        /// Initialize the Sobol Sequence.
        /// </summary>
        private void initialize()
        {

            // special case: dimension 1 -> use unit initialization
            for (int i = 1; i <= BITS; i++)
            {
                _direction[0, i] = 1L << (BITS - i);
            }

            using (StreamReader reader = new StreamReader(new MemoryStream(Properties.Resources.new_joe_kuo_6)))
            {
                // ignore the first line
                reader.ReadLine();

                int index = 1;
                string? line = null;
                while ((line = reader.ReadLine()) != null)
                {
                    var st = line.Split(' ');

                    int dim;
                    int.TryParse(st[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out dim);

                   if (dim >= 2 && dim <= Dimension)
                   {
                        // we have found the right dimension
                        int i, s = 0, a = 0;
                        for (i = 1; i < st.Length; i++)
                        {
                            if (st[i] != "")
                            {
                                int.TryParse(st[i], NumberStyles.Integer, CultureInfo.InvariantCulture, out s);
                                break;
                            }
                        }
                        i++;
                        for (; i < st.Length; i++)
                        {
                            if (st[i] != "")
                            {
                                int.TryParse(st[i], NumberStyles.Integer, CultureInfo.InvariantCulture, out a);
                                break;
                            }
                        }
                        i++;
                        int[] m = new int[s + 1];
                        for (; i < st.Length; i++)
                        {
                            if (st[i] != "")
                            {
                                for (int j = 1; j <= s; j++)
                                {
                                    int.TryParse(st[i + j - 1], NumberStyles.Integer, CultureInfo.InvariantCulture, out m[j]);
                                }
                                break;
                            }
                        }
                        initDirectionVector(index++, a, m);
                    }

                    if (dim > Dimension)
                    {
                       return;
                    }
                }

            }

        }


        private void initDirectionVector(int d,  int a,  int[] m)
        {
            int s = m.Length - 1;
            for (int i = 1; i <= s; i++)
            {
                _direction[d, i] = ((long)m[i]) << (BITS - i);
            }
            for (int i = s + 1; i <= BITS; i++)
            {
                _direction[d, i] = _direction[d, i - s] ^ (_direction[d, i - s] >> s);
                for (int k = 1; k <= s - 1; k++)
                {
                    _direction[d, i] ^= ((a >> (s - 1 - k)) & 1) * _direction[d, i - k];
                }
            }
        }

        /// <summary>
        /// Applies the seeded linear matrix scramble to the direction numbers and draws the digital
        /// shifts. For each dimension in order: the sub-diagonal bits of a random unit
        /// lower-triangular bit matrix, row by row, then the 52 shift bits.
        /// </summary>
        /// <param name="prng">The seeded pseudorandom number generator.</param>
        private void ScrambleDirections(Random prng)
        {
            var rowMasks = new long[BITS];
            for (int d = 0; d < Dimension; d++)
            {
                // Row i keeps digit i on the diagonal and mixes a random subset of earlier digits.
                for (int i = 1; i <= BITS; i++)
                {
                    long mask = 1L << (BITS - i);
                    for (int j = 1; j < i; j++)
                    {
                        if (prng.Next(2) == 1) mask |= 1L << (BITS - j);
                    }
                    rowMasks[i - 1] = mask;
                }
                // The generator XORs direction numbers, and the scramble is linear over GF(2), so
                // scrambling the direction numbers once scrambles every generated point exactly.
                for (int c = 1; c <= BITS; c++)
                {
                    _direction[d, c] = ApplyLinearScramble(rowMasks, _direction[d, c]);
                }
                long shift = 0;
                for (int i = 1; i <= BITS; i++)
                {
                    if (prng.Next(2) == 1) shift |= 1L << (BITS - i);
                }
                _shift[d] = shift;
            }
        }

        /// <summary>
        /// Multiplies a direction number's digit vector by the unit lower-triangular bit matrix over GF(2).
        /// </summary>
        /// <param name="rowMasks">The matrix rows as digit masks.</param>
        /// <param name="vector">The direction number.</param>
        /// <returns>The scrambled direction number.</returns>
        private static long ApplyLinearScramble(long[] rowMasks, long vector)
        {
            long result = 0;
            for (int i = 1; i <= BITS; i++)
            {
                if (Parity(rowMasks[i - 1] & vector)) result |= 1L << (BITS - i);
            }
            return result;
        }

        /// <summary>
        /// Returns the bit parity of a value (true when the number of set bits is odd).
        /// </summary>
        /// <param name="value">The value.</param>
        /// <returns>True for odd parity.</returns>
        private static bool Parity(long value)
        {
            value ^= value >> 32;
            value ^= value >> 16;
            value ^= value >> 8;
            value ^= value >> 4;
            value ^= value >> 2;
            value ^= value >> 1;
            return (value & 1L) != 0;
        }

        /// <summary>
        /// Returns a double-precision number that is greater than or equal to 0.0, and less than 1.0.
        /// </summary>
        /// <returns>A double-precision number that is greater than or equal to 0.0, and less than 1.0.</returns>
        public double[] NextDouble()
        {
            double[] v = new double[Dimension];
            if (count == 0)
            {
                count++;
                //return v;
            }

            // find the index c of the rightmost 0
            int c = 1;
            int value = count - 1;
            while ((value & 1) == 1)
            {
                value >>= 1;
                c++;
            }

            for (int i = 0; i < Dimension; i++)
            {
                x[i] ^= _direction[i, c];
                // The digital shift is all zero for the unrandomized sequence, so the emission is
                // arithmetically identical there.
                v[i] = (double)(x[i] ^ _shift[i]) / SCALE;
            }
            count++;
            return v;
        }

        /// <summary>
        /// Skip to a specific index in the sequence.
        /// </summary>
        /// <param name="index">The index in the sequence.</param>
        /// <returns>A double-precision number that is greater than or equal to 0.0, and less than 1.0.</returns>
        public double[] SkipTo(int index)
        {
            if (index == 0)
            {
                // reset x vector
                x.Fill(0);
            }
            else
            {
                int i = index - 1;
                long grayCode = i ^ (i >> 1); // compute the gray code of i = i XOR floor(i / 2)
                for (int j = 0; j < Dimension; j++)
                {
                    long result = 0;
                    for (int k = 1; k <= BITS; k++)
                    {
                        long shift = grayCode >> (k - 1);
                        if (shift == 0)
                        {
                            // stop, as all remaining bits will be zero
                            break;
                        }
                        // the k-th bit of i
                         long ik = shift & 1;
                        result ^= ik * _direction[j, k];
                    }
                    x[j] = result;
                }
            }
            count = index;
            return NextDouble();
        }

    }
}
