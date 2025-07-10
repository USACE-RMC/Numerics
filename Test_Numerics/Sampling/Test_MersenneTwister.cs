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

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Sampling;

namespace Sampling
{
    /// <summary>
    /// Unit tests for the Mersenne Twister class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// Mersenne Twister with improved initialization. <see href="http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html"/>
    /// </remarks>
    [TestClass]
    public class Test_MersenneTwister
    {
        /// <summary>
        /// Test against output file found here: <see href="http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html"/> 
        /// </summary>
        [TestMethod]
        public void Test_MT19937()
        {
            // Initialize from array
            int[] seeds = { 291, 564, 837, 1110 };
            var prng = new MersenneTwister(seeds);

            // Test Int32
            var trueInt32 = new uint[] { 1067595299, 955945823, 477289528, 4107218783, 4228976476, 3344332714, 3355579695, 227628506, 810200273, 2591290167 };
            for (int i = 0; i < 10; i++)
            {
               Assert.AreEqual(trueInt32[i], prng.GenRandInt32(), 1E-8);
            }
            // run to 1000
            for (int i = 10; i < 1000; i++)
            {
                prng.GenRandInt32();
            }
            // Test NextDouble()
            var trueDouble = new double[] { 0.76275443, 0.99000644, 0.98670464, 0.10143112, 0.27933125, 0.69867227, 0.94218740, 0.03427201, 0.78842173, 0.28180608 };
            for (int i = 0; i < 10; i++)
            {
                Assert.AreEqual(trueDouble[i], prng.NextDouble(), 1E-8);
            }

        }
    }
}
