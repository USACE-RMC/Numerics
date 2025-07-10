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
    /// Unit test for the Sobol Sequence class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// R Core Team (2024). _R: A Language and Environment for Statistical Computing_.R Foundation for Statistical Computing, Vienna,
    /// Austria. <see href="https://www.R-project.org/"/>
    /// </remarks>
    [TestClass]
    public class Test_SobolSequence
    {

        /// <summary>
        /// Tested against the 'sobol' method in the 'randtoolbox' R package.
        /// </summary>
        [TestMethod]
        public void Test_Sobol()
        {
            // the results from R
            var trueResult = new double[,]
            { 
                { 0.5000, 0.5000 },
                { 0.7500, 0.2500 },
                { 0.2500, 0.7500 },
                { 0.3750, 0.3750 },
                { 0.8750, 0.8750 },
                { 0.6250, 0.1250 },
                { 0.1250, 0.6250 },
                { 0.1875, 0.3125 },
                { 0.6875, 0.8125 },
                { 0.9375, 0.0625 }
            };

            var sobol = new SobolSequence(2);
            for (int i = 0; i < 10; i++)
            {
                var rnd = sobol.NextDouble();
                for (int j = 0; j < 2; j++)
                {
                    Assert.AreEqual(trueResult[i, j], rnd[j]);
                }
            }
        
        }
    }
}
