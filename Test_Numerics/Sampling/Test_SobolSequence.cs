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
