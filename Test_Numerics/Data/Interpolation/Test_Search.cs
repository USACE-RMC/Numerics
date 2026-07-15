using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;

namespace Data.Interpolation
{
    /// <summary>
    /// Regression tests for interpolation search helpers.
    /// </summary>
    [TestClass]
    public class Test_Search
    {
        /// <summary>
        /// Gets the interior, endpoint, and out-of-range values exercised by every search overload.
        /// </summary>
        private static readonly double[] Queries = { -5d, 0d, 5d, 25d, 35d, 40d, 45d };

        /// <summary>
        /// Verifies that ascending searches retain their existing behavior for every overload.
        /// </summary>
        [TestMethod]
        public void AscendingSearchesAgreeWithSequential()
        {
            VerifyAllOverloads(new[] { 0d, 10d, 20d, 30d, 40d }, SortOrder.Ascending);
        }

        /// <summary>
        /// Verifies descending interiors, endpoints, and out-of-range values for every overload.
        /// </summary>
        [TestMethod]
        public void DescendingSearchesAgreeWithSequential()
        {
            VerifyAllOverloads(new[] { 40d, 30d, 20d, 10d, 0d }, SortOrder.Descending);
        }

        /// <summary>
        /// Verifies nonzero bisection starts and hunt guesses in both directions.
        /// </summary>
        [TestMethod]
        public void NonzeroStartsReturnTheSequentialBracket()
        {
            VerifyStarts(new[] { 0d, 10d, 20d, 30d, 40d }, SortOrder.Ascending, 35d, 15d);
            VerifyStarts(new[] { 40d, 30d, 20d, 10d, 0d }, SortOrder.Descending, 15d, 25d);
        }

        /// <summary>
        /// Verifies that bisection and hunt agree with sequential search for every supported data container.
        /// </summary>
        /// <param name="values">The ordered values to search.</param>
        /// <param name="order">The sort direction of <paramref name="values"/>.</param>
        private static void VerifyAllOverloads(double[] values, SortOrder order)
        {
            var paired = new OrderedPairedData(values, values, true, order, true, order);
            IList<Ordinate> ordinates = values.Select(value => new Ordinate(value, value)).ToList();

            foreach (double query in Queries)
            {
                int expected = Search.Sequential(query, values, 0, order);
                Assert.AreEqual(expected, Search.Bisection(query, values, 0, order));
                Assert.AreEqual(expected, Search.Hunt(query, values, 0, order));

                expected = Search.Sequential(query, paired);
                Assert.AreEqual(expected, Search.Bisection(query, paired));
                Assert.AreEqual(expected, Search.Hunt(query, paired));

                expected = Search.Sequential(query, ordinates, 0, order);
                Assert.AreEqual(expected, Search.Bisection(query, ordinates, 0, order));
                Assert.AreEqual(expected, Search.Hunt(query, ordinates, 0, order));
            }
        }

        /// <summary>
        /// Verifies bisection and hunt from nonzero starting indices for every supported data container.
        /// </summary>
        /// <param name="values">The ordered values to search.</param>
        /// <param name="order">The sort direction of <paramref name="values"/>.</param>
        /// <param name="bisectionQuery">The value used to verify bisection.</param>
        /// <param name="huntQuery">The value used to verify hunt.</param>
        private static void VerifyStarts(
            double[] values,
            SortOrder order,
            double bisectionQuery,
            double huntQuery)
        {
            var paired = new OrderedPairedData(values, values, true, order, true, order);
            IList<Ordinate> ordinates = values.Select(value => new Ordinate(value, value)).ToList();

            int expected = Search.Sequential(bisectionQuery, values, 0, order);
            Assert.AreEqual(expected, Search.Bisection(bisectionQuery, values, 2, order));
            Assert.AreEqual(expected, Search.Bisection(bisectionQuery, paired, 2));
            Assert.AreEqual(expected, Search.Bisection(bisectionQuery, ordinates, 2, order));

            expected = Search.Sequential(huntQuery, values, 0, order);
            Assert.AreEqual(expected, Search.Hunt(huntQuery, values, 3, order));
            Assert.AreEqual(expected, Search.Hunt(huntQuery, paired, 3));
            Assert.AreEqual(expected, Search.Hunt(huntQuery, ordinates, 3, order));
        }
    }
}
