using System;
using System.Collections.Generic;

namespace Numerics.Sampling
{
    /// <summary>
    /// A class to perform Latin hypercube sampling (LHS).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para> 
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Stratified sampling is used to generate uniform random numbers by dividing the interval [0,1) into n bins 
    /// of equal probability, where n is the total number of samples required.
    /// <para>
    /// For each column:
    /// </para>
    /// <list type="number">
    /// <item>
    /// One value is placed in every bin, at a uniform random offset within the bin (<see cref="Random(int, int, int)"/>)
    /// or at the bin midpoint (<see cref="Median(int, int, int)"/>).
    /// </item>
    /// <item>
    /// The n values are put in random order with a Fisher-Yates shuffle.
    /// </item>
    /// </list>
    /// <para>
    /// Each bin is therefore used exactly once per column. Every column draws from its own
    /// <see cref="MersenneTwister"/>, seeded from one value of the master generator.
    /// </para>
    /// </para>
    /// </remarks>
    public class LatinHypercube
    {

        /// <summary>
        /// Generate Latin Hypercube samples with uniform stratified sampling.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="seed">Random seed (optional).</param>
        /// <remarks>
        /// A seed of zero or below, including the default -1, seeds the generator from the clock, so the
        /// sample differs on every call. Pass a positive seed, or a seeded generator through
        /// <see cref="Random(int, int, System.Random)"/>, for a reproducible sample.
        /// </remarks>
        public static double[,] Random(int sampleSize, int dimension = 1, int seed = -1)
        {
            return Random(sampleSize, dimension, seed > 0 ? new MersenneTwister(seed) : new MersenneTwister());
        }

        /// <summary>
        /// Generate Latin Hypercube samples with uniform stratified sampling, drawing from a supplied generator.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="random">The master generator. One value is drawn from it per column, to seed that column's own <see cref="MersenneTwister"/>.</param>
        /// <returns>A <paramref name="sampleSize"/> by <paramref name="dimension"/> array of values in [0, 1).</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="random"/> is null.</exception>
        /// <remarks>
        /// <para>
        /// The sample is a deterministic function of the generator's state, with no clock fallback, so a
        /// seeded generator always gives the same sample. Passing <c>new MersenneTwister(seed)</c> with a
        /// positive seed reproduces <see cref="Random(int, int, int)"/> with that seed exactly.
        /// </para>
        /// <para>
        /// The generator advances by <paramref name="dimension"/> calls to <see cref="System.Random.Next()"/>,
        /// so one generator can be shared by consecutive calls to draw independent samples.
        /// </para>
        /// </remarks>
        public static double[,] Random(int sampleSize, int dimension, Random random)
        {
            if (random == null) throw new ArgumentNullException(nameof(random));

            double[,] lhs = new double[sampleSize, dimension];
            var rndMaster = random;

            for (int col = 0; col < dimension; col++)
            {
                // Create bins
                double[] bins = new double[sampleSize];
                var rnd = new MersenneTwister(rndMaster.Next());

                for (int i = 0; i < sampleSize; i++)
                {
                    double delta = rnd.NextDouble();
                    bins[i] = (i + delta) / sampleSize;
                }

                // Shuffle bins
                Shuffle(bins, rnd);

                // Assign to LHS matrix
                for (int row = 0; row < sampleSize; row++)
                    lhs[row, col] = bins[row];
            }

            return lhs;
        }

        /// <summary>
        /// Generate Latin Hypercube samples using median bin locations.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="seed">Random seed (optional).</param>
        /// <remarks>
        /// A seed of zero or below, including the default -1, seeds the generator from the clock, so the
        /// sample differs on every call. Pass a positive seed, or a seeded generator through
        /// <see cref="Median(int, int, System.Random)"/>, for a reproducible sample.
        /// </remarks>
        public static double[,] Median(int sampleSize, int dimension = 1, int seed = -1)
        {
            return Median(sampleSize, dimension, seed > 0 ? new MersenneTwister(seed) : new MersenneTwister());
        }

        /// <summary>
        /// Generate Latin Hypercube samples using median bin locations, shuffling with a supplied generator.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="random">The master generator. One value is drawn from it per column, to seed that column's shuffle.</param>
        /// <returns>A <paramref name="sampleSize"/> by <paramref name="dimension"/> array of bin midpoints in (0, 1).</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="random"/> is null.</exception>
        /// <remarks>
        /// The sample is a deterministic function of the generator's state, with no clock fallback. Passing
        /// <c>new MersenneTwister(seed)</c> with a positive seed reproduces <see cref="Median(int, int, int)"/>
        /// with that seed exactly. The generator advances by <paramref name="dimension"/> calls to
        /// <see cref="System.Random.Next()"/>.
        /// </remarks>
        public static double[,] Median(int sampleSize, int dimension, Random random)
        {
            if (random == null) throw new ArgumentNullException(nameof(random));

            double[,] lhs = new double[sampleSize, dimension];
            var rndMaster = random;

            for (int col = 0; col < dimension; col++)
            {
                // Create median-centered bins
                double[] bins = new double[sampleSize];
                for (int i = 0; i < sampleSize; i++)
                    bins[i] = (i + 0.5) / sampleSize;

                // Shuffle bins
                Shuffle(bins, new MersenneTwister(rndMaster.Next()));

                // Assign to LHS matrix
                for (int row = 0; row < sampleSize; row++)
                    lhs[row, col] = bins[row];
            }

            return lhs;
        }

        /// <summary>
        /// Performs in-place Fisher–Yates shuffle.
        /// </summary>
        private static void Shuffle(double[] array, Random rnd)
        {
            for (int i = array.Length - 1; i > 0; i--)
            {
                int j = rnd.Next(i + 1);
                (array[i], array[j]) = (array[j], array[i]);
            }
        }
    }
}
