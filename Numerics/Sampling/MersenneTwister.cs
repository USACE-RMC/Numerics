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

//
//   A C-program for MT19937, with initialization improved 2002/1/26.
//  Coded by Takuji Nishimura and Makoto Matsumoto.
//
//  Before using, initialize the state by using init_genrand(seed)
//   or init_by_array(init_key, key_length).
//
//  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//  are met:
//
//    1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote
//      products derived from this software without specific prior written
//      permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//   Any feedback is very welcome.
//  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
//

using System;

namespace Numerics.Sampling
{

    /// <summary>
    /// The Mersenne Twister pseudo-random number generator.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// The source code can be found at:
    /// <see href="http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class MersenneTwister : Random
    {

        /// <summary>
        /// Construct a Mersenne Twister PRNG using the clock to create a random seed.
        /// </summary>
        public MersenneTwister()
        {
            Initialize((uint)DateTime.UtcNow.Ticks);
        }

        /// <summary>
        /// Construct a Mersenne Twister PRNG given a seed. 
        /// </summary>
        /// <param name="seed">The PRNG seed.</param>
        public MersenneTwister(int seed)
        {
            Initialize((uint)seed);
        }

        /// <summary>
        /// Construct a Mersenne Twister PRNG given an array of seeds.
        /// </summary>
        /// <param name="seeds">The array of PRNG seeds.</param>
        public MersenneTwister(int[] seeds)
        {
            var useeds = new uint[seeds.Length];
            for (int i = 0; i < seeds.Length; i++)
            {
                useeds[i] = (uint)seeds[i];
            }
            Initialize(useeds, useeds.Length);
        }


        private const int N = 624;
        private const int M = 397;
        private const uint matrixA = 0x9908b0dfU;     // constant vector a 
        private const uint upperMask = 0x80000000U;   // most significant w-r bits 
        private const uint lowerMask = 0x7fffffffU;   // least significant r bits 
        private uint[] mt = new uint[N]; // the array for the state vector
        int mti = N + 1; // mti==N+1 means mt[N] is not initialized
        private static readonly uint[] mag01 = { 0x0U, matrixA }; // mag01[x] = x * MATRIX_A for x=0,1


        // The following constants, are used within the ADDITIONAL functions genrand_real2b() and
        // genrand_real3b(), equivalent to genrand_real() and genrand_real3(), but that return
        // evenly distributed values in the ranges [0, 1-kMT_Gap] and [0+kMT_Gap, 1-kMT_Gap],
        // respectively. A similar statement is valid also for genrand_real2c(), genrand_real4b()
        // and genrand_real5b(). See the section "Functions and procedures implemented" above,
        // for more details.
        // 
        // If you want to change the value of kMT_Gap, it is suggested to do it so that:
        // 5e-15 <= kMT_Gap <= 5e-2

        const double k2_31 = 2147483648.0;     // 2^31   ==  2147483648 == 80000000
        const double k2_32 = 2.0 * k2_31;      // 2^32   ==  4294967296 == 0
        const double k2_32b = k2_32 - 1.0;     // 2^32-1 ==  4294967295 == FFFFFFFF == -1

        const double kMT_Gap = 0.0000000000005;       // 5.0E-13
        const double kMT_Gap2 = 2.0 * kMT_Gap;         // 1.0E-12
        const double kMT_GapInterval = 1.0 - kMT_Gap2; // 0.9999999999990
        const double kMT_2b = kMT_GapInterval / k2_32b;
        const double kMT_2c = kMT_2b;
        const double kMT_3b = kMT_2b;

        /// <summary>
        /// Initialize with seed.
        /// </summary>
        /// <param name="seed">The PRNG seed.</param>
        public void Initialize(uint seed)     
        {
            mt[0] = seed & 0xffffffffU;
            for (mti = 1; mti < N; mti++)
            {
                mt[mti] = (uint)(1812433253U * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
                /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
                /* In the previous versions, MSBs of the seed affect   */
                /* only MSBs of the array _mt[].                       */
                /* 2002/01/09 modified by Makoto Matsumoto             */
                mt[mti] &= 0xffffffffU;
                /* for >32 bit machines */
            }
        }

        /// <summary>
        /// Initialize by an array with array-length.
        /// </summary>
        /// <param name="init_key">The array of PRNG seeds.</param>
        /// <param name="key_length">The array length.</param>
        private void Initialize(uint[] init_key, int key_length)
        {
            Initialize(19650218U);
            int i = 1, j = 0;
            int k = (N > key_length) ? N : key_length;
            for (; k > 0; k--)
            {
                mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525U)) + init_key[j] + (uint)j;
                mt[i] &= 0xffffffffU; // For WORDSIZE > 32 machines
                i++; j++;
                if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
                if (j >= key_length) j = 0;
            }
            for (k = N - 1; k > 0; k--)
            {
                mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941U)) - (uint)i;
                mt[i] &= 0xffffffffU; // For WORDSIZE > 32 machines
                i++;
                if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
            }

            mt[0] = upperMask; // MSB is 1; assuring non-zero initial array
        }

        /// <summary>
        /// Generates a random number on [0,0xffffffff]-interval.
        /// </summary>
        public uint GenRandInt32()  
        {

            uint y;

            /* mag01[x] = x * MATRIX_A  for x=0,1 */

            if (mti >= N)
            {
                /* generate _n words at one time */
                int kk;

                if (mti == N + 1) /* if init_genrand() has not been called, */
                {
                    Initialize(5489); /* a default initial seed is used */
                }

                for (kk = 0; kk < N - M; kk++)
                {
                    y = (mt[kk] & upperMask) | (mt[kk + 1] & lowerMask);
                    mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
                }

                for (; kk < N - 1; kk++)
                {
                    y = (mt[kk] & upperMask) | (mt[kk + 1] & lowerMask);
                    mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
                }

                y = (mt[N - 1] & upperMask) | (mt[0] & lowerMask);
                mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

                mti = 0;
            }

            y = mt[mti++];

            /* Tempering */
            y ^= y >> 11;
            y ^= (y << 7) & 0x9d2c5680;
            y ^= (y << 15) & 0xefc60000;
            y ^= y >> 18;

            return y;
        }

        /// <summary>
        /// Generates a random number on [0,0x7fffffff]-interval.
        /// </summary>
        public int GenRandInt31()
        {
            return (int)(GenRandInt32() >> 1);
        }

        /// <summary>
        /// Generates a random number on [0,1]-real-interval.
        /// </summary>
        /// <returns></returns>
        public double GenRandReal1()
        {
            return GenRandInt32() * (1.0 / 4294967295.0); // Divided by 2^32-1
        }

        /// <summary>
        /// Generates a random number on [0,1)-real-interval.
        /// </summary>
        public double GenRandReal2()
        {
            return GenRandInt32() * (1.0 / 4294967296.0); // Divided by 2^32
        }

        /// <summary>
        /// Generates a random number on (0,1)-real-interval.
        /// </summary>
        /// <returns></returns>
        public double GenRandReal3()
        {
            return (GenRandInt32() + 0.5) * (1.0 / 4294967296.0); // Divided by 2^32
        }

        /// <summary>
        /// Generates a random number on [0,1) with 53-bit resolution.
        /// </summary>
        public double GenRandRes53()
        {
            uint a = GenRandInt32() >> 5, b = GenRandInt32() >> 6;
            return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
        }

        /// <inheritdoc/>
        public sealed override int Next()
        {
            int result;
            do
            {
                result = GenRandInt31();
            } while (result == int.MaxValue);

            return result;
        }

        /// <inheritdoc/>
        public sealed override int Next(int maxExclusive)
        {
            if (maxExclusive <= 0)
                throw new ArgumentOutOfRangeException(nameof(maxExclusive), "Must be positive.");

            uint threshold = (uint.MaxValue - (uint.MaxValue % (uint)maxExclusive));

            uint r;
            do
            {
                r = GenRandInt32();
            } while (r >= threshold);

            return (int)(r % (uint)maxExclusive);
        }

        /// <inheritdoc/>
        public sealed override int Next(int minInclusive, int maxExclusive)
        {
            if (minInclusive >= maxExclusive)
                throw new ArgumentOutOfRangeException(nameof(minInclusive), "minInclusive must be less than maxExclusive.");

            return Next(maxExclusive - minInclusive) + minInclusive;
        }

        /// <inheritdoc/>
        public sealed override double NextDouble()
        {
            return GenRandReal2();
        }

    }
}
