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
                        heap.RemoveMin(out int popped, out float poppedWeight);
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
                    heap.RemoveMin(out int popped, out float poppedWeight);
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
            float previous = float.NegativeInfinity;
            for (int i = 0; i < capacity; i++)
            {
                heap.RemoveMin(out _, out float weight);
                Assert.IsGreaterThanOrEqualTo(previous, weight, "The drain must be non-decreasing.");
                previous = weight;
            }
            Assert.AreEqual(0, heap.Count);
            Assert.Throws<InvalidOperationException>(() => heap.RemoveMin(out _, out _));
        }
    }
}
