using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests the compiled-topology network wrapper over the Dijkstra solvers.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_Network
    {
        /// <summary>
        /// Builds the two-way five-node test network used throughout: a short chain 0-1-2 with a
        /// slow bypass 0-3-2 and a spur 2-4.
        /// </summary>
        /// <returns>The edges and the network destinations (node 4).</returns>
        private static (Edge[] Edges, int[] Destinations) BuildFixture()
        {
            var edges = new[]
            {
                new Edge(0, 1, 1, 0),
                new Edge(1, 0, 1, 0),
                new Edge(1, 2, 1, 1),
                new Edge(2, 1, 1, 1),
                new Edge(0, 3, 5, 2),
                new Edge(3, 0, 5, 2),
                new Edge(3, 2, 5, 3),
                new Edge(2, 3, 5, 3),
                new Edge(2, 4, 1, 4),
                new Edge(4, 2, 1, 4),
            };
            return (edges, new[] { 4 });
        }

        /// <summary>
        /// Construction succeeds for a non-empty edge set, groups the adjacency correctly with
        /// empty lists for edge-free directions, and reports the node count.
        /// </summary>
        [TestMethod]
        public void CtorBuildsAdjacencyForNonEmptyEdges()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            Assert.AreEqual(5, network.NodeCount);
            Assert.HasCount(5, network.IncomingEdges);
            Assert.HasCount(5, network.OutgoingEdges);
            Assert.HasCount(2, network.OutgoingEdges[0]);   // 0->1 and 0->3
            Assert.HasCount(1, network.IncomingEdges[4]);
            Assert.AreEqual(2, network.IncomingEdges[4][0].FromIndex);
            CollectionAssert.AreEqual(new[] { 4 }, network.DestinationIndices);

            // An isolated node (index 5 present only via the node count) is impossible here;
            // instead verify a node with no incoming edges gets an empty list, not null.
            var oneWay = new Network(new[] { new Edge(0, 1, 1, 0) }, new[] { 1 });
            Assert.IsNotNull(oneWay.IncomingEdges[0]);
            Assert.IsEmpty(oneWay.IncomingEdges[0]);
            Assert.IsNotNull(oneWay.OutgoingEdges[1]);
            Assert.IsEmpty(oneWay.OutgoingEdges[1]);
        }

        /// <summary>
        /// The constructor copies both input arrays: mutating them afterward changes nothing.
        /// </summary>
        [TestMethod]
        public void CtorCopiesInputs()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var before = network.Solve(4);

            edges[0] = new Edge(0, 1, 100, 0);
            destinations[0] = 0;

            var after = network.Solve(4);
            Assert.AreEqual(before[0, 2], after[0, 2], 0f);
            CollectionAssert.AreEqual(new[] { 4 }, network.DestinationIndices);
        }

        /// <summary>
        /// Malformed construction inputs fail with clear argument exceptions.
        /// </summary>
        [TestMethod]
        public void CtorValidationThrows()
        {
            var edges = new[] { new Edge(0, 1, 1, 0) };

            Assert.Throws<ArgumentNullException>(() => new Network(null!, [0]));
            Assert.Throws<ArgumentNullException>(() => new Network(edges, null!));
            Assert.Throws<ArgumentException>(() => new Network(new Edge[0], [0]));
            Assert.Throws<ArgumentException>(() => new Network(edges, new int[0]));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Network(edges, [5]));
            Assert.Throws<ArgumentException>(() => new Network(new[] { new Edge(-1, 1, 1, 0) }, [0]));
            Assert.Throws<ArgumentException>(() => new Network(new[] { new Edge(0, 1, 1, -1) }, [0]));
        }

        /// <summary>
        /// The network solves equal the static Dijkstra solves on the same inputs, cell for
        /// cell, for single and multiple destinations.
        /// </summary>
        [TestMethod]
        public void SolveMatchesDijkstra()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var networkSingle = network.Solve(4);
            var directSingle = Dijkstra.Solve(edges, 4, network.NodeCount);
            var networkMulti = network.Solve(new[] { 0, 4 });
            var directMulti = Dijkstra.Solve(edges, new[] { 0, 4 }, network.NodeCount);

            for (int i = 0; i < network.NodeCount; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(directSingle[i, c], networkSingle[i, c], 0f, $"Single cell [{i},{c}].");
                    Assert.AreEqual(directMulti[i, c], networkMulti[i, c], 0f, $"Multi cell [{i},{c}].");
                }
            }

            Assert.Throws<ArgumentOutOfRangeException>(() => network.Solve(9));
            Assert.Throws<ArgumentNullException>(() => network.Solve((int[])null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => network.Solve(new[] { 9 }));
        }

        /// <summary>
        /// Custom edge weights are actually applied: raising the short chain's weights reroutes
        /// traffic over the bypass.
        /// </summary>
        [TestMethod]
        public void SolveCustomWeightsAreApplied()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            // Base network: node 0 routes 0-1-2-4 at cost 3.
            var baseTable = network.Solve(4);
            Assert.AreEqual(3f, baseTable[0, 2], 0f);
            Assert.AreEqual(1f, baseTable[0, 0], 0f);

            // Flood the 0-1 and 1-2 links: node 0 must reroute 0-3-2-4 at cost 11.
            var weights = new float[] { 100, 100, 100, 100, 5, 5, 5, 5, 1, 1 };
            var flooded = network.Solve(weights);
            Assert.AreEqual(11f, flooded[0, 2], 0f);
            Assert.AreEqual(3f, flooded[0, 0], 0f);

            Assert.Throws<ArgumentNullException>(() => network.Solve((float[])null!));
            Assert.Throws<ArgumentException>(() => network.Solve(new float[3]));
        }

        /// <summary>
        /// An infinite custom weight severs its edge: with every route to the destination
        /// flooded to infinity, the start node becomes unreachable.
        /// </summary>
        [TestMethod]
        public void SolveCustomWeightsInfinityBlocksEdge()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            float inf = float.PositiveInfinity;
            // Sever 1->2 and 3->2 (positions 2 and 6): node 0 cannot reach node 4.
            var weights = new float[] { 1, 1, inf, 1, 5, 5, inf, 5, 1, 1 };
            var severed = network.Solve(weights);

            Assert.IsTrue(float.IsPositiveInfinity(severed[0, 2]));
            Assert.IsTrue(float.IsPositiveInfinity(severed[1, 2]));
            Assert.AreEqual(1f, severed[2, 2], 0f);
        }

        /// <summary>
        /// The table-reuse overloads produce cell-identical results to the allocating solves,
        /// repeated reuse into the same buffer stays identical, and wrong dimensions throw.
        /// </summary>
        [TestMethod]
        public void SolveReuseOverloadIsBitIdentical()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var weights = new float[] { 2, 2, 3, 3, 4, 4, 6, 6, 1, 1 };

            var allocated = network.Solve(weights);
            var reused = new float[network.NodeCount, 3];
            network.Solve(weights, reused);
            network.Solve(weights, reused); // steady-state second pass into the same buffer

            var nearestAllocated = network.SolveNearest();
            var nearestReused = new float[network.NodeCount, 3];
            var baseWeights = new float[] { 1, 1, 1, 1, 5, 5, 5, 5, 1, 1 };
            network.SolveNearest(baseWeights, nearestReused);

            for (int i = 0; i < network.NodeCount; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(allocated[i, c], reused[i, c], 0f, $"Reuse cell [{i},{c}].");
                    Assert.AreEqual(nearestAllocated[i, c], nearestReused[i, c], 0f, $"Nearest reuse cell [{i},{c}].");
                }
            }

            Assert.Throws<ArgumentNullException>(() => network.Solve(weights, null!));
            Assert.Throws<ArgumentException>(() => network.Solve(weights, new float[2, 3]));
            Assert.Throws<ArgumentException>(() => network.SolveNearest(weights, new float[network.NodeCount, 2]));
        }

        /// <summary>
        /// The single-pass nearest solve matches the merged solve costs over the network's
        /// destinations.
        /// </summary>
        [TestMethod]
        public void SolveNearestMatchesMergedNetworkSolve()
        {
            var (edges, _) = BuildFixture();
            var network = new Network(edges, new[] { 0, 4 });

            var merged = network.Solve(new[] { 0, 4 });
            var nearest = network.SolveNearest();

            for (int i = 0; i < network.NodeCount; i++)
            {
                Assert.AreEqual(merged[i, 2], nearest[i, 2], 0f, $"Cost mismatch at node {i}.");
            }
        }

        /// <summary>
        /// With nothing removed, GetPath returns the exact unblocked edge sequence to the
        /// nearest destination through both overloads.
        /// </summary>
        [TestMethod]
        public void GetPathHappyPath()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var direct = network.GetPath(new int[0], 0);
            Assert.IsNotNull(direct);
            CollectionAssert.AreEqual(new List<int> { 0, 1, 4 }, direct);

            var table = network.Solve(destinations);
            var viaTable = network.GetPath(new int[0], 0, table);
            Assert.IsNotNull(viaTable);
            CollectionAssert.AreEqual(new List<int> { 0, 1, 4 }, viaTable);
        }

        /// <summary>
        /// Removing an edge index on the recorded route forces the detour over the bypass, and
        /// every edge bearing the removed index is excluded.
        /// </summary>
        [TestMethod]
        public void GetPathDetoursAroundRemovedEdges()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var detour = network.GetPath(new[] { 1 }, 0);
            Assert.IsNotNull(detour);
            CollectionAssert.AreEqual(new List<int> { 2, 3, 4 }, detour);

            var table = network.Solve(destinations);
            var detourViaTable = network.GetPath(new[] { 1 }, 0, table);
            Assert.IsNotNull(detourViaTable);
            CollectionAssert.AreEqual(new List<int> { 2, 3, 4 }, detourViaTable);
        }

        /// <summary>
        /// When every route is severed, the first overload returns null and the second returns
        /// an empty list, per their respective contracts.
        /// </summary>
        [TestMethod]
        public void GetPathReturnsNullOrEmptyWhenSevered()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            Assert.IsNull(network.GetPath(new[] { 1, 3 }, 0));

            var table = network.Solve(destinations);
            var severed = network.GetPath(new[] { 1, 3 }, 0, table);
            Assert.IsNotNull(severed);
            Assert.IsEmpty(severed!);
        }

        /// <summary>
        /// Starting at a destination returns an empty path through both overloads.
        /// </summary>
        [TestMethod]
        public void GetPathStartIsDestinationReturnsEmpty()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var direct = network.GetPath(new[] { 1 }, 4);
            Assert.IsNotNull(direct);
            Assert.IsEmpty(direct!);

            var table = network.Solve(destinations);
            var viaTable = network.GetPath(new[] { 1 }, 4, table);
            Assert.IsNotNull(viaTable);
            Assert.IsEmpty(viaTable!);
        }

        /// <summary>
        /// The caller's removal array is never sorted or otherwise mutated.
        /// </summary>
        [TestMethod]
        public void GetPathDoesNotMutateEdgesToRemove()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var removals = new[] { 4, 1, 3 };
            network.GetPath(removals, 0);
            CollectionAssert.AreEqual(new[] { 4, 1, 3 }, removals);
        }

        /// <summary>
        /// A valid route passing through node zero is returned — zero is a real node index,
        /// not an unreachable sentinel.
        /// </summary>
        [TestMethod]
        public void GetPathThroughNodeZero()
        {
            var edges = new[]
            {
                new Edge(1, 0, 1, 0),
                new Edge(0, 2, 1, 1),
            };
            var network = new Network(edges, new[] { 2 });

            var path = network.GetPath(new int[0], 1);
            Assert.IsNotNull(path);
            CollectionAssert.AreEqual(new List<int> { 0, 1 }, path);
        }

        /// <summary>
        /// With several destinations, GetPath routes to the nearest one, and an exact tie
        /// resolves deterministically across repeated calls.
        /// </summary>
        [TestMethod]
        public void GetPathPicksNearestDestination()
        {
            // Chain 0-1-2-3 with destinations 0 and 3: node 1 is nearer to 0.
            var edges = new[]
            {
                new Edge(1, 0, 1, 0),
                new Edge(0, 1, 1, 0),
                new Edge(1, 2, 1, 1),
                new Edge(2, 1, 1, 1),
                new Edge(2, 3, 1, 2),
                new Edge(3, 2, 1, 2),
            };
            var network = new Network(edges, new[] { 0, 3 });

            var fromOne = network.GetPath(new int[0], 1);
            Assert.IsNotNull(fromOne);
            CollectionAssert.AreEqual(new List<int> { 0 }, fromOne);

            // The middle of an even chain ties exactly; the choice must repeat.
            var tieNetwork = new Network(new[]
            {
                new Edge(1, 0, 1, 0),
                new Edge(1, 2, 1, 1),
            }, new[] { 0, 2 });
            var first = tieNetwork.GetPath(new int[0], 1);
            var second = tieNetwork.GetPath(new int[0], 1);
            Assert.IsNotNull(first);
            Assert.HasCount(1, first!);
            CollectionAssert.AreEqual(first, second);
        }

        /// <summary>
        /// Removing the only edge index into the destination severs the route even though two
        /// directed edges bear that index.
        /// </summary>
        [TestMethod]
        public void GetPathHandlesDuplicateEdgeIndices()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            Assert.IsNull(network.GetPath(new[] { 4 }, 0));
        }

        /// <summary>
        /// The table fast path returns the identical route to the full solve when the removals
        /// miss the recorded route.
        /// </summary>
        [TestMethod]
        public void GetPathFastPathMatchesFullSolve()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var table = network.Solve(destinations);

            // Removing the bypass leaves the recorded route 0-1-2-4 untouched.
            var viaTable = network.GetPath(new[] { 3 }, 0, table);
            var direct = network.GetPath(new[] { 3 }, 0);
            Assert.IsNotNull(viaTable);
            Assert.IsNotNull(direct);
            CollectionAssert.AreEqual(direct, viaTable);
        }

        /// <summary>
        /// A detour from a custom-weight result table is solved with the same positional weights;
        /// using the constructor weights would choose the other available route.
        /// </summary>
        [TestMethod]
        public void GetPathCustomWeightTableUsesCustomWeightsForDetour()
        {
            var edges = new[]
            {
                new Edge(0, 1, 1, 0),
                new Edge(1, 4, 1, 1),
                new Edge(0, 2, 1, 2),
                new Edge(2, 4, 1, 3),
                new Edge(0, 3, 10, 4),
                new Edge(3, 4, 10, 5),
            };
            var customWeights = new float[] { 1, 1, 100, 100, 2, 2 };
            var network = new Network(edges, new[] { 4 });
            var customTable = network.Solve(customWeights);

            CollectionAssert.AreEqual(
                new List<int> { 4, 5 },
                network.GetPath(new[] { 1 }, 0, customTable, customWeights));
        }

        /// <summary>
        /// GetPath rejects null inputs, out-of-range start nodes, and wrong table dimensions
        /// with clear argument exceptions.
        /// </summary>
        [TestMethod]
        public void GetPathValidationThrows()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var table = network.Solve(destinations);

            Assert.Throws<ArgumentNullException>(() => network.GetPath(null!, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => network.GetPath(new int[0], 9));
            Assert.Throws<ArgumentNullException>(() => network.GetPath(null!, 0, table));
            Assert.Throws<ArgumentNullException>(() => network.GetPath(new int[0], 0, null!));
            Assert.Throws<ArgumentException>(() => network.GetPath(new int[0], 0, new float[2, 3]));
            Assert.Throws<ArgumentOutOfRangeException>(() => network.GetPath(new int[0], 9, table));
            Assert.Throws<ArgumentNullException>(() => network.GetPath(new int[0], 0, table, null!));
            Assert.Throws<ArgumentException>(() => network.GetPath(new int[0], 0, table, new float[1]));
        }

        /// <summary>
        /// Performance smoke: one thousand custom-weight re-solves through the table-reuse
        /// overload on a 100 by 100 grid.
        /// </summary>
        [TestMethod]
        [TestCategory("Performance")]
        public void P3_RepeatedWeightResolves()
        {
            const int side = 100;
            var edges = new List<Edge>();
            int index = 0;
            for (int r = 0; r < side; r++)
            {
                for (int c = 0; c < side; c++)
                {
                    int node = r * side + c;
                    if (c + 1 < side)
                    {
                        edges.Add(new Edge(node, node + 1, 1, index++));
                        edges.Add(new Edge(node + 1, node, 1, index++));
                    }
                    if (r + 1 < side)
                    {
                        edges.Add(new Edge(node, node + side, 1, index++));
                        edges.Add(new Edge(node + side, node, 1, index++));
                    }
                }
            }
            var network = new Network(edges.ToArray(), new[] { 0, side * side - 1 });
            var weights = new float[edges.Count];
            for (int j = 0; j < weights.Length; j++) weights[j] = 1f;
            var table = new float[network.NodeCount, 3];

            var stopwatch = Stopwatch.StartNew();
            for (int pass = 0; pass < 1000; pass++)
            {
                weights[pass % weights.Length] = 1f + (pass % 5);
                network.SolveNearest(weights, table);
            }
            stopwatch.Stop();
            Console.WriteLine($"P3 1000 re-solves on {side}x{side}: {stopwatch.Elapsed.TotalSeconds:F3} s ({stopwatch.Elapsed.TotalMilliseconds / 1000.0:F3} ms per solve)");

            Assert.AreEqual(0f, table[0, 2], 0f);
        }
    }
}
