using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests the internal indexed min-heap that backs the shortest-path solvers.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_IndexedMinHeap
    {
        /// <summary>
        /// Fixed-seed fuzz of Add, RemoveMin, and DecreaseKey against a linear-scan reference
        /// model, across three Clear-and-reuse rounds with a full ordered drain per round.
        /// </summary>
        [TestMethod]
        public void IndexedHeapFuzzMatchesReferenceModel()
        {
            var randy = new MersenneTwister(24680);
            var heap = new IndexedMinHeap(48);

            for (int round = 0; round < 3; round++)
            {
                heap.Clear();
                var model = new Dictionary<int, float>();

                for (int op = 0; op < 1200; op++)
                {
                    int node = randy.Next(0, 48);
                    double action = randy.NextDouble();
                    float weight = (float)randy.NextDouble();

                    if (action < 0.45)
                    {
                        if (!model.ContainsKey(node))
                        {
                            heap.Add(node, weight);
                            model[node] = weight;
                        }
                    }
                    else if (action < 0.65)
                    {
                        if (model.Count == 0) continue;
                        heap.RemoveMin(out int popped, out double poppedWeight);
                        float modelMin = float.MaxValue;
                        foreach (var entry in model)
                        {
                            if (entry.Value < modelMin) modelMin = entry.Value;
                        }
                        Assert.AreEqual(modelMin, poppedWeight, 0f);
                        Assert.AreEqual(modelMin, model[popped], 0f);
                        model.Remove(popped);
                    }
                    else
                    {
                        heap.DecreaseKey(node, weight);
                        if (!model.ContainsKey(node) || weight < model[node]) model[node] = weight;
                    }

                    Assert.AreEqual(model.Count, heap.Count);
                }

                while (model.Count > 0)
                {
                    heap.RemoveMin(out int popped, out double poppedWeight);
                    float modelMin = float.MaxValue;
                    foreach (var entry in model)
                    {
                        if (entry.Value < modelMin) modelMin = entry.Value;
                    }
                    Assert.AreEqual(modelMin, poppedWeight, 0f);
                    Assert.AreEqual(modelMin, model[popped], 0f);
                    model.Remove(popped);
                }
                Assert.AreEqual(0, heap.Count);
            }
        }

        /// <summary>
        /// The heap fills to exactly its capacity and drains fully sorted; the count never
        /// exceeds the capacity under solver-shaped operation sequences.
        /// </summary>
        [TestMethod]
        public void IndexedHeapCapacityInvariant()
        {
            const int capacity = 64;
            var randy = new MersenneTwister(13579);
            var heap = new IndexedMinHeap(capacity);

            // Exact fill: every node enters once.
            var weights = new float[capacity];
            for (int node = 0; node < capacity; node++)
            {
                weights[node] = (float)randy.NextDouble();
                heap.Add(node, weights[node]);
                Assert.IsLessThanOrEqualTo(capacity, heap.Count);
            }
            Assert.AreEqual(capacity, heap.Count);
            Assert.Throws<InvalidOperationException>(() => heap.Add(0, 0f));

            // Decrease keys while full never grows the heap.
            for (int node = 0; node < capacity; node += 3)
            {
                heap.DecreaseKey(node, weights[node] / 2f);
                Assert.AreEqual(capacity, heap.Count);
            }

            // Full drain returns ascending weights.
            double previous = double.NegativeInfinity;
            for (int i = 0; i < capacity; i++)
            {
                heap.RemoveMin(out _, out double weight);
                Assert.IsGreaterThanOrEqualTo(previous, weight, "The drain must be non-decreasing.");
                previous = weight;
            }
            Assert.AreEqual(0, heap.Count);
            Assert.Throws<InvalidOperationException>(() => heap.RemoveMin(out _, out _));
        }

        /// <summary>
        /// Fixed-seed fuzz with heavy weight ties: the indexed heap pops the same node sequence as the
        /// swap-based public <see cref="BinaryHeap{T}"/> under the same operation sequence, so tie order is
        /// part of the contract the shortest-path solvers rely on.
        /// </summary>
        [TestMethod]
        public void IndexedHeapMatchesBinaryHeapPopOrderOnTies()
        {
            const int capacity = 40;
            var randy = new MersenneTwister(97531);
            var heap = new IndexedMinHeap(capacity);
            var reference = new BinaryHeap<int>(capacity);
            var present = new bool[capacity];

            for (int op = 0; op < 4000; op++)
            {
                int node = randy.Next(0, capacity);
                float weight = randy.Next(0, 4);
                double action = randy.NextDouble();
                if (action < 0.4)
                {
                    if (present[node]) continue;
                    heap.Add(node, weight);
                    reference.Add(new BinaryHeap<int>.Node(weight, node, node));
                    present[node] = true;
                }
                else if (action < 0.7)
                {
                    if (heap.Count == 0) continue;
                    Assert.AreEqual(reference.RemoveMin().Index, heap.MinNode);
                    heap.RemoveMin(out int popped, out _);
                    present[popped] = false;
                }
                else
                {
                    heap.DecreaseKey(node, weight);
                    reference.DecreaseKey(new BinaryHeap<int>.Node(weight, node, node));
                    present[node] = true;
                }
                Assert.AreEqual(reference.Count, heap.Count);
            }
            while (heap.Count > 0)
            {
                heap.RemoveMin(out int popped, out _);
                Assert.AreEqual(reference.RemoveMin().Index, popped);
            }
        }

        /// <summary>
        /// The per-node record: <see cref="IndexedMinHeap.WeightOf"/> is positive infinity before a node is added,
        /// follows decreases, keeps the removal weight afterwards, and resets on Clear.
        /// </summary>
        [TestMethod]
        public void IndexedHeapWeightOfFollowsTheNode()
        {
            var heap = new IndexedMinHeap(4);
            Assert.AreEqual(double.PositiveInfinity, heap.WeightOf(2));

            heap.Add(2, 5d);
            heap.Add(1, 3d);
            heap.DecreaseKey(2, 1d);
            heap.DecreaseKey(2, 4d); // not smaller: ignored
            Assert.AreEqual(1d, heap.WeightOf(2));
            Assert.AreEqual(2, heap.MinNode);

            heap.RemoveMin(out int node, out double weight);
            Assert.AreEqual(2, node);
            Assert.AreEqual(1d, weight);
            Assert.AreEqual(1d, heap.WeightOf(2), "The removal weight is kept.");

            heap.DecreaseKey(2, 0.5d); // absent: added again
            Assert.AreEqual(0.5d, heap.WeightOf(2));
            Assert.AreEqual(2, heap.Count);

            heap.Clear();
            Assert.AreEqual(0, heap.Count);
            Assert.AreEqual(double.PositiveInfinity, heap.WeightOf(2));
            Assert.AreEqual(double.PositiveInfinity, heap.WeightOf(1));
        }
    }
}
