/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System;
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

            initialize();
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
        /// The direction vector for each component.
        /// </summary>
        private long[,] _direction;

        /// <summary>
        /// The current state.
        /// </summary>
        private long[] x;

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
                string line = null;
                while ((line = reader.ReadLine()) != null)
                {
                    var st = line.Split(' ');

                    try
                    {
                        int dim;
                        int.TryParse(st[0], out dim);
                        
                       if (dim >= 2 && dim <= Dimension)
                       {
                            // we have found the right dimension
                            int i, s = 0, a = 0;
                            for (i = 1; i < st.Length; i++)
                            {
                                if (st[i] != "")
                                {
                                    int.TryParse(st[i], out s);
                                    break;
                                }
                            }
                            i++;
                            for (; i < st.Length; i++)
                            {
                                if (st[i] != "")
                                {
                                    int.TryParse(st[i], out a);
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
                                        int.TryParse(st[i + j - 1], out m[j]);
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
                    catch (Exception ex)
                    {
                       throw ex;
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
                v[i] = (double)x[i] / SCALE;
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
