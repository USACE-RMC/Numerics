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
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for Box-Cox transformation. These methods were tested against R's "boxcox()" and "boxcoxTransform()" methods
    /// from the "EnvStats" package and the "bxcx()" method from the "sae" package.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// Millard SP (2013). EnvStats: An R Package for Environmental Statistics. Springer, New York. ISBN 978-1-4614-8455-4, <see href="https://www.springer.com"/>
    /// </item>
    /// <item>
    /// Molina I, Marhuenda Y (2015). “sae: An R Package for Small Area Estimation.” The R Journal, 7(1), 81–98. 
    /// <see href="https://journal.r-project.org/archive/2015/RJ-2015-007/RJ-2015-007.pdf"/>
    /// </item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_BoxCox
    {
        /// <summary>
        /// Test the fitting method
        /// </summary>
        [TestMethod()]
        public void Test_Fit()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            double l1 = 0d;
            BoxCox.FitLambda(sample, out l1);
            Assert.AreEqual(1.670035d, l1,  1E-4);
        }

        /// <summary>
        /// Test the fitting method against R's "boxcox()" method
        /// </summary>
        [TestMethod]
        public void Test_Fit_R()
        {
            var sample = new double[] { 1, 3, 59, 1, 6, 4, 9, 13, 5, 84, 35, 8, 31, 34, 9, 1, 35, 66, 1, 65, 4, 68, 46 };
            double l1 = 0;
            BoxCox.FitLambda(sample, out l1);
            Assert.AreEqual(0.1314482, l1, 1E-6);
        }

        /// <summary>
        /// Test the transform and reverse transform methods
        /// </summary>
        [TestMethod()]
        public void Test_Transform()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var trueVals = new[] { 2359.592d, 2331.397d, 2334.155d, 2320.102d, 2333.603d, 2313.5d, 2341.056d, 2356.821d, 2382.357d, 2376.241d, 2385.139d, 2367.633d, 2382.357d, 2359.869d, 2379.576d, 2404.372d, 2368.466d, 2357.93d, 2362.364d, 2370.131d, 2412.474d, 2428.711d, 2427.87d, 2435.442d, 2439.653d, 2430.112d, 2437.968d, 2406.606d, 2401.86d, 2430.673d, 2455.118d, 2448.646d, 2452.022d, 2453.429d, 2465.826d, 2457.089d, 2459.624d, 2404.931d, 2308.279d, 2279.513d, 2298.949d, 2241.111d, 2264.23d, 2286.349d, 2312.4d, 2333.603d, 2309.103d, 2302.514d, 2254.699d, 2269.41d, 2278.693d, 2277.327d, 2299.771d, 2327.813d, 2393.772d, 2388.201d, 2394.33d, 2383.748d, 2356.545d, 2367.633d, 2359.315d, 2357.099d, 2379.576d, 2399.628d, 2399.07d, 2426.468d, 2417.508d, 2431.234d, 2401.023d, 2433.478d, 2459.624d, 2491.827d, 2481.357d, 2479.943d, 2527.33d, 2531.603d, 2527.614d, 2541.873d, 2563.607d, 2552.158d, 2567.905d, 2533.029d, 2546.728d, 2577.372d, 2599.803d, 2603.549d, 2594.907d, 2596.635d, 2599.226d, 2571.059d, 2602.972d, 2598.362d, 2601.243d, 2617.977d, 2636.202d, 2652.735d, 2656.221d, 2663.489d, 2651.283d, 2621.444d, 2634.174d, 2620.289d, 2684.466d, 2690.597d, 2680.673d, 2686.509d, 2660.871d, 2623.179d, 2566.185d, 2604.99d, 2597.498d, 2592.605d, 2637.65d, 2669.018d, 2676.009d, 2651.864d, 2681.84d, 2609.316d, 2622.023d, 2585.127d, 2584.265d, 2524.767d, 2589.153d, 2604.125d, 2603.26d, 2641.128d, 2647.509d, 2646.929d, 2669.6d, 2644.318d, 2599.803d, 2646.639d, 2710.78d, 2723.095d, 2721.334d, 2700.827d, 2729.26d, 2720.16d, 2699.364d, 2668.436d, 2636.492d, 2583.69d, 2527.045d, 2462.442d, 2530.464d, 2441.338d, 2485.882d, 2500.613d, 2414.152d, 2462.724d, 2532.744d, 2504.868d, 2419.187d, 2478.53d, 2455.962d, 2372.63d, 2291.823d, 2450.334d, 2441.057d, 2425.068d, 2464.697d, 2502.315d, 2477.965d, 2516.8d, 2468.929d, 2418.628d, 2448.927d, 2511.115d, 2505.436d, 2526.76d, 2519.644d, 2477.965d, 2479.095d, 2471.187d, 2500.896d, 2536.736d, 2519.929d, 2529.894d, 2544.728d, 2676.883d, 2674.843d, 2656.512d, 2648.09d, 2601.531d, 2643.158d, 2662.326d, 2660.581d, 2653.316d, 2711.952d, 2688.552d, 2697.318d, 2724.269d, 2734.844d, 2741.021d, 2753.98d, 2780.268d, 2736.903d, 2760.764d, 2706.094d, 2722.801d, 2678.049d, 2667.562d, 2545.585d, 2620.289d, 2613.068d, 2625.781d, 2666.689d };
            var vals = BoxCox.Transform(sample, 1.670035d);
            for (int i = 0; i < vals.Count; i++)
                Assert.AreEqual(vals[i], trueVals[i], 0.001d);

            var revVals = BoxCox.InverseTransform(vals, 1.670035d);
            for (int i = 0; i < revVals.Count; i++)
                Assert.AreEqual(revVals[i], sample[i], 0.001d);
        }

        /// <summary>
        /// Test the transform and reverse transform methods against R's "boxcoxTransform()" and "bxcx()" methods
        /// </summary>
        [TestMethod]
        public void Test_Transform_R()
        {
            double x = BoxCox.Transform(5d, 1d);
            double y = BoxCox.InverseTransform(x, 1d);
            double trueX = 4;
            double trueY = 5;
            Assert.AreEqual(trueX, x);
            Assert.AreEqual(trueY, y);

            var sample = new double[] { 1, 3, 59, 1, 6, 4, 9, 13, 5, 84, 35, 8, 31, 34, 9, 1, 35, 66, 1, 65, 4, 68, 46 };
            var trueVals = new double[] { 0.000000, 1.181898, 5.394755, 0.000000, 2.020349, 1.520639, 2.547415, 3.050330, 1.792351, 6.012795, 4.532207, 2.391402, 4.340082,
                4.486038, 2.547415, 0.000000, 4.532207, 5.587797, 0.000000, 5.561342, 1.520639, 5.639679, 4.976243 };
            var vals = BoxCox.Transform(sample, 0.1314482);
            for (int i = 0; i < vals.Count; i++)
                Assert.AreEqual(vals[i], trueVals[i], 1E-6);

            var revVals = BoxCox.InverseTransform(vals, 0.1314482);
            for (int i = 0; i < revVals.Count; i++)
                Assert.AreEqual(revVals[i], sample[i], 0.001d);
        }
    }
}
