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
using Numerics.Data;

namespace Data.Interpolation
{
    /// <summary>
    /// Unit tests for the Bilinear class. The output values of the class were tested against
    /// R's "bilinear( )" function of the "akima" package.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> References: </b>
    /// Akima H, Gebhardt A (2022). akima: Interpolation of Irregularly and Regularly Spaced Data. 
    /// R package version 0.6-3.4, https://CRAN.R-project.org/package=akima.
    /// </remarks>
    [TestClass]
    public class Test_Bilinear
    {
        /// <summary>
        /// Test the most basic implementation of the bilinear class and its interpolation function
        /// </summary>
        [TestMethod]
        public void Test_BiLinear()
        {
            var x1Array = new double[] { 100d, 200d, 300d, 400d, 500d };
            var x2Array = new double[] { 50d, 100d, 150d, 200d, 250d };
            var yArray = new double[5, 5];
            yArray[0, 0] = 850.36d;
            yArray[1, 0] = 865.39d;
            yArray[2, 0] = 867.27d;
            yArray[3, 0] = 869.71d;
            yArray[4, 0] = 871.84d;
            yArray[0, 1] = 868.45d;
            yArray[1, 1] = 878.66d;
            yArray[2, 1] = 880.47d;
            yArray[3, 1] = 881.91d;
            yArray[4, 1] = 883.19d;
            yArray[0, 2] = 896.71d;
            yArray[1, 2] = 901.2d;
            yArray[2, 2] = 901.92d;
            yArray[3, 2] = 902.4d;
            yArray[4, 2] = 903.22d;
            yArray[0, 3] = 914.78d;
            yArray[1, 3] = 918.54d;
            yArray[2, 3] = 919.14d;
            yArray[3, 3] = 919.9d;
            yArray[4, 3] = 920.18d;
            yArray[0, 4] = 928.87d;
            yArray[1, 4] = 929.46d;
            yArray[2, 4] = 929.54d;
            yArray[3, 4] = 929.62d;
            yArray[4, 4] = 929.68d;

            var bilinear = new Bilinear(x1Array, x2Array, yArray);
            double x1 = 350d;
            double x2 = 75d;
            double y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 874.84d, 1E-6);
        }

        /// <summary>
        /// Test all of the possible log transforms within the bilinear class and its interpolation function
        /// </summary>
        [TestMethod]
        public void Test_Log()
        {
            var x1Array = new double[] { 100d, 200d, 300d, 400d, 500d };
            var x2Array = new double[] { 50d, 100d, 150d, 200d, 250d };
            var yArray = new double[5, 5];
            yArray[0, 0] = 850.36d;
            yArray[1, 0] = 865.39d;
            yArray[2, 0] = 867.27d;
            yArray[3, 0] = 869.71d;
            yArray[4, 0] = 871.84d;
            yArray[0, 1] = 868.45d;
            yArray[1, 1] = 878.66d;
            yArray[2, 1] = 880.47d;
            yArray[3, 1] = 881.91d;
            yArray[4, 1] = 883.19d;
            yArray[0, 2] = 896.71d;
            yArray[1, 2] = 901.2d;
            yArray[2, 2] = 901.92d;
            yArray[3, 2] = 902.4d;
            yArray[4, 2] = 903.22d;
            yArray[0, 3] = 914.78d;
            yArray[1, 3] = 918.54d;
            yArray[2, 3] = 919.14d;
            yArray[3, 3] = 919.9d;
            yArray[4, 3] = 920.18d;
            yArray[0, 4] = 928.87d;
            yArray[1, 4] = 929.46d;
            yArray[2, 4] = 929.54d;
            yArray[3, 4] = 929.62d;
            yArray[4, 4] = 929.68d;

            double x1 = 350d;
            double x2 = 75d;

            var LogLinLin = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.Logarithmic };
            double y1 = LogLinLin.Interpolate(x1, x2);
            Assert.AreEqual(y1, 874.909523653025d, 1E-6);

            var LinLogLin = new Bilinear(x1Array, x2Array, yArray) { X2Transform = Transform.Logarithmic };
            double y2 = LinLogLin.Interpolate(x1, x2);
            Assert.AreEqual(y2, 875.919023759159d, 1E-6);

            var LinLinLog = new Bilinear(x1Array, x2Array, yArray) { YTransform = Transform.Logarithmic };
            double y3 = LinLinLog.Interpolate(x1, x2);
            Assert.AreEqual(y3, 874.8164, 1E-4);

            var LinLogLog = new Bilinear(x1Array, x2Array, yArray) { X2Transform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            double y4 = LinLogLog.Interpolate(x1, x2);
            Assert.AreEqual(y4, 875.896104342695d, 1E-6);

            var LogLogLin = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.Logarithmic, X2Transform = Transform.Logarithmic };
            double y5 = LogLogLin.Interpolate(x1, x2);
            Assert.AreEqual(y5, 875.9855, 1E-4);

            var LogLinLog = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            double y6 = LogLinLog.Interpolate(x1, x2);
            Assert.AreEqual(y6, 874.886, 1E-4);

            var LogLogLog = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.Logarithmic, X2Transform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            double y7 = LogLogLog.Interpolate(x1, x2);
            Assert.AreEqual(y7, 875.962713889793d, 1E-6);
        }


        /// <summary>
        /// Test all of the possible z variate transforms of x1 and x2 of the bilinear function
        /// </summary>
        [TestMethod]
        public void Test_Z_X()
        {
            var x1Array = new double[] { 0.7, 0.8, 0.88, 0.9, 0.99 };
            var x2Array = new double[] { 0.65, 0.78, 0.8, 0.86, 0.9 };
            var yArray = new double[5, 5];
            yArray[0, 0] = 850.36d;
            yArray[1, 0] = 865.39d;
            yArray[2, 0] = 867.27d;
            yArray[3, 0] = 869.71d;
            yArray[4, 0] = 871.84d;
            yArray[0, 1] = 868.45d;
            yArray[1, 1] = 878.66d;
            yArray[2, 1] = 880.47d;
            yArray[3, 1] = 881.91d;
            yArray[4, 1] = 883.19d;
            yArray[0, 2] = 896.71d;
            yArray[1, 2] = 901.2d;
            yArray[2, 2] = 901.92d;
            yArray[3, 2] = 902.4d;
            yArray[4, 2] = 903.22d;
            yArray[0, 3] = 914.78d;
            yArray[1, 3] = 918.54d;
            yArray[2, 3] = 919.14d;
            yArray[3, 3] = 919.9d;
            yArray[4, 3] = 920.18d;
            yArray[0, 4] = 928.87d;
            yArray[1, 4] = 929.46d;
            yArray[2, 4] = 929.54d;
            yArray[3, 4] = 929.62d;
            yArray[4, 4] = 929.68d;

            double x1 = 0.86d;
            double x2 = 0.79d;

            var ZLinLin = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.NormalZ };
            double y1 = ZLinLin.Interpolate(x1, x2);
            Assert.AreEqual(y1, 890.8358, 1E-4);

            var LinZLin = new Bilinear(x1Array, x2Array, yArray) { X2Transform = Transform.NormalZ };
            double y2 = LinZLin.Interpolate(x1, x2);
            Assert.AreEqual(y2, 890.7267, 1E-4);

            var ZZLin = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.NormalZ, X2Transform = Transform.NormalZ };
            double y3 = ZZLin.Interpolate(x1, x2);
            Assert.AreEqual(y3, 890.6835, 1E-4);
        }

        /// <summary>
        /// Test all of the possible z variate transforms of y within the bilinear function
        /// </summary>
        [TestMethod]
        public void Test_Z_Y()
        {
            var x1Array = new double[] { 0.7, 0.8, 0.88, 0.9, 0.99 };
            var x2Array = new double[] { 0.65, 0.78, 0.8, 0.86, 0.9 };
            var yArray = new double[5, 5];
            yArray[0, 0] = 0.91468032d;
            yArray[1, 0] = 0.930847173d;
            yArray[2, 0] = 0.932869374d;
            yArray[3, 0] = 0.935493933d;
            yArray[4, 0] = 0.937785044d;
            yArray[0, 1] = 0.934138628d;
            yArray[1, 1] = 0.945120902d;
            yArray[2, 1] = 0.947067808d;
            yArray[3, 1] = 0.948616728d;
            yArray[4, 1] = 0.949993546d;
            yArray[0, 2] = 0.964536184d;
            yArray[1, 2] = 0.969365803d;
            yArray[2, 2] = 0.970140263d;
            yArray[3, 2] = 0.97065657;
            yArray[4, 2] = 0.971538594d;
            yArray[0, 3] = 0.98397298d;
            yArray[1, 3] = 0.988017382d;
            yArray[2, 3] = 0.988662766d;
            yArray[3, 3] = 0.989480251d;
            yArray[4, 3] = 0.98978143d;
            yArray[0, 4] = 0.999128732d;
            yArray[1, 4] = 0.999763359d;
            yArray[2, 4] = 0.999849411d;
            yArray[3, 4] = 0.999935462d;
            yArray[4, 4] = 0.999999999d;

            double x1 = 0.86d;
            double x2 = 0.79d;

            var LinLinZ = new Bilinear(x1Array, x2Array, yArray) { YTransform = Transform.NormalZ };
            double y1 = LinLinZ.Interpolate(x1, x2);
            Assert.AreEqual(y1, 0.9596228, 1E-6);

            var LinZZ = new Bilinear(x1Array, x2Array, yArray) { X2Transform = Transform.NormalZ, YTransform = Transform.NormalZ };
            double y2 = LinZZ.Interpolate(x1, x2);
            Assert.AreEqual(y2, 0.95946, 1E-6);

            var ZLinZ = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.NormalZ, YTransform = Transform.NormalZ };
            double y3 = ZLinZ.Interpolate(x1, x2);
            Assert.AreEqual(y3, 0.9595799, 1E-6);

            var ZZZ = new Bilinear(x1Array, x2Array, yArray) { X1Transform = Transform.NormalZ, X2Transform = Transform.NormalZ, YTransform = Transform.NormalZ };
            double y4 = ZZZ.Interpolate(x1, x2);
            Assert.AreEqual(y4, 0.9594168, 1E-6);
        }

        /// <summary>
        /// Test the edge case of the bilinear function
        /// </summary>
        [TestMethod]
        public void Test_BilinearEdgeCases()
        {
            var x1Array = new double[] { 100d, 200d, 300d, 400d, 500d };
            var x2Array = new double[] { 50d, 100d, 150d, 200d, 250d };
            var yArray = new double[5, 5];
            yArray[0, 0] = 850.36d;
            yArray[1, 0] = 865.39d;
            yArray[2, 0] = 867.27d;
            yArray[3, 0] = 869.71d;
            yArray[4, 0] = 871.84d;
            yArray[0, 1] = 868.45d;
            yArray[1, 1] = 878.66d;
            yArray[2, 1] = 880.47d;
            yArray[3, 1] = 881.91d;
            yArray[4, 1] = 883.19d;
            yArray[0, 2] = 896.71d;
            yArray[1, 2] = 901.2d;
            yArray[2, 2] = 901.92d;
            yArray[3, 2] = 902.4d;
            yArray[4, 2] = 903.22d;
            yArray[0, 3] = 914.78d;
            yArray[1, 3] = 918.54d;
            yArray[2, 3] = 919.14d;
            yArray[3, 3] = 919.9d;
            yArray[4, 3] = 920.18d;
            yArray[0, 4] = 928.87d;
            yArray[1, 4] = 929.46d;
            yArray[2, 4] = 929.54d;
            yArray[3, 4] = 929.62d;
            yArray[4, 4] = 929.68d;
            var bilinear = new Bilinear(x1Array, x2Array, yArray);

            // Ascending - Both out of range
            // Top left
            double x1 = 50;
            double x2 = 28;
            double y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 850.36, 1E-6);
            // Top Right
            x1 = 50;
            x2 = 300;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 928.87, 1E-6);
            // Bottom Left
            x1 = 600;
            x2 = 25;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 871.84, 1E-6);
            // Bottom Right
            x1 = 600;
            x2 = 300;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 929.68, 1E-6);

            // Ascending - x1 out
            // Top
            x1 = 50;
            x2 = 75;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 859.405, 1E-6);
            // Bottom
            x1 = 600;
            x2 = 225;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 924.93, 1E-6);

            // Ascending - x2 out
            // Top
            x1 = 125;
            x2 = 25;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 854.11750, 1E-6);
            // Bottom
            x1 = 450;
            x2 = 300;
            y = bilinear.Interpolate(x1, x2);
            Assert.AreEqual(y, 929.65000, 1E-6);

        }
    }
}
