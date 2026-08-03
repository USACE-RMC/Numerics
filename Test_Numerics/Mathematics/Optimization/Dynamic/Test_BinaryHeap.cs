using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests binary heap behavior used by dynamic optimization routines.
    /// </summary>
    [TestClass]
    public class Test_BinaryHeap
    {
        /// <summary>
        /// Checking weights on heap.
        /// </summary>
        [TestMethod]
        public void HeapTest1()
        {
            //Node weights
            float[] weights = new float[] { .3f, .5f, 32f, 15f, 12f, .01f, -4f };

            BinaryHeap<double> heap = new BinaryHeap<double>(30);
            for (int i = 0; i < weights.Length; i++)
            {
                heap.Add(new BinaryHeap<double>.Node(weights[i], i, i * .5));
            }

            Assert.AreEqual(heap.RemoveMin().Weight, weights[6]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[5]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[0]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[1]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[4]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[3]);
            Assert.AreEqual(heap.RemoveMin().Weight, weights[2]);
        }

        /// <summary>
        /// Random heap and weight allotment.
        /// </summary>
        [TestMethod]
        public void HeapTest2()
        {
            //Random Node weights
            float[] weights = new float[1000]; //= new float[] { .3f, .5f, 32f, 15f, 12f, .01f, -4f };

            Random randy = new MersenneTwister(12345);
            for (int i = 0; i < weights.Length; i++)
            {
                weights[i] = (float)randy.NextDouble();
            }
            //Add to the heap
            BinaryHeap<double> heap = new BinaryHeap<double>(weights.Length);
            for (int i = 0; i < weights.Length; i++)
            {
                heap.Add(new BinaryHeap<double>.Node(weights[i], i, i * .5));
            }

            Array.Sort(weights);
            //Compare
            for (int i = 0; i < weights.Length; i++)
            {
                Assert.AreEqual(weights[i],heap.RemoveMin().Weight);
            }
        }

        /// <summary>
        /// Making sure heap is ordering correctly.
        /// </summary>
        [TestMethod]
        public void HeapTest3()
        {
            //Node weights
            float[] weights = new float[] { .3f, .5f, 32f, 15f, 12f, .01f, -4f };

            BinaryHeap<double> heap = new BinaryHeap<double>(30);
            for (int i = 0; i < weights.Length; i++)
            {
                heap.Add(new BinaryHeap<double>.Node(weights[i], i, i * .5));
            }

            for (int i = 0; i < weights.Length; i++)
            {
                heap.Replace(new BinaryHeap<double>.Node(weights[i], i, i));
            }

            Array.Sort(weights);
            //Compare
            for (int i = 0; i < weights.Length; i++)
            {
                Assert.AreEqual(weights[i],heap.RemoveMin().Weight );
            }
        }

        /// <summary>
        /// Testing edge case.
        /// </summary>
        [TestMethod]
        public void HeapTest4()
        {
            //Node weights
            float[] weights = new float[] { .3f, .5f, 32f, 15f, 12f, .01f, -4f };

            BinaryHeap<double> heap = new BinaryHeap<double>(30);
            for (int i = 0; i < weights.Length; i++)
            {
                heap.Add(new BinaryHeap<double>.Node(weights[i], i, i * .5));
            }

            for (int i = 0; i < weights.Length; i++)
            {
                heap.Replace(new BinaryHeap<double>.Node(weights[i], i, i * 5));
            }

            //Compare
            for (int i = 0; i < weights.Length; i++)
            {
                Assert.AreNotEqual(heap.RemoveMin().Value , weights[i]);
            }
        }

        /// <summary>
        /// Decrease key for ordering
        /// </summary>
        [TestMethod]
        public void DecreaseKeyTest()
        {
            var heap = new BinaryHeap<string>(10);
            heap.Add(new BinaryHeap<string>.Node(10f, 1, "A"));
            heap.Add(new BinaryHeap<string>.Node(20f, 2, "B"));
            heap.Add(new BinaryHeap<string>.Node(30f, 3, "C"));

            heap.DecreaseKey(new BinaryHeap<string>.Node(5f, 3, "C")); // now should be top

            var node = heap.RemoveMin();
            Assert.AreEqual(3, node.Index);
            Assert.AreEqual(5f, node.Weight);
        }

        /// <summary>
        /// Edge case for heap capacity.
        /// </summary>
        [TestMethod]
        public void HeapCapacityExceededTest()
        {
            var ex = Assert.Throws<Exception>(() =>
            {
                var heap = new BinaryHeap<int>(3);
                heap.Add(new BinaryHeap<int>.Node(1f, 1, 1));
                heap.Add(new BinaryHeap<int>.Node(2f, 2, 2));
                heap.Add(new BinaryHeap<int>.Node(3f, 3, 3));
                heap.Add(new BinaryHeap<int>.Node(4f, 4, 4)); // exceeds capacity
            });
        }

        /// <summary>
        /// Testing RemoveMin() is getting called correctly on heap.
        /// </summary>
        [TestMethod]
        public void RemoveMinFromEmptyHeapTest()
        {
            var ex = Assert.Throws<Exception>(() =>
            {
                var heap = new BinaryHeap<int>(10);
                heap.RemoveMin(); // should throw
            });
        }

        /// <summary>
        /// Using heap and Decrease key and name of the test explains it.
        /// </summary>
        [TestMethod]
        public void ReplaceWithHigherWeightShouldDoNothingTest()
        {
            var heap = new BinaryHeap<int>(10);
            heap.Add(new BinaryHeap<int>.Node(5f, 1, 1));
            heap.DecreaseKey(new BinaryHeap<int>.Node(10f, 1, 1)); // Should not update weight

            Assert.AreEqual(5f, heap.RemoveMin().Weight); // weight should still be 5
        }

        /// <summary>
        /// Edge case to check ordering is correct with negative weights.
        /// </summary>
        [TestMethod]
        public void OrderingWithNegativeWeightsTest()
        {
            var heap = new BinaryHeap<string>(5);
            heap.Add(new BinaryHeap<string>.Node(-10f, 1, "A"));
            heap.Add(new BinaryHeap<string>.Node(-20f, 2, "B"));
            heap.Add(new BinaryHeap<string>.Node(0f, 3, "C"));

            Assert.AreEqual("B", heap.RemoveMin().Value);
            Assert.AreEqual("A", heap.RemoveMin().Value);
            Assert.AreEqual("C", heap.RemoveMin().Value);
        }

        /// <summary>
        /// Replacing the minimum node with a heavier weight must sift it down, not leave it at
        /// the root: the drain returns every node in weight order.
        /// </summary>
        [TestMethod]
        public void ReplaceWithHigherWeightReordersHeap()
        {
            var heap = new BinaryHeap<string>(10);
            heap.Add(new BinaryHeap<string>.Node(1f, 1, "A"));
            heap.Add(new BinaryHeap<string>.Node(2f, 2, "B"));
            heap.Add(new BinaryHeap<string>.Node(3f, 3, "C"));
            heap.Add(new BinaryHeap<string>.Node(4f, 4, "D"));

            heap.Replace(new BinaryHeap<string>.Node(10f, 1, "A")); // the minimum becomes the maximum

            Assert.AreEqual("B", heap.RemoveMin().Value);
            Assert.AreEqual("C", heap.RemoveMin().Value);
            Assert.AreEqual("D", heap.RemoveMin().Value);
            var last = heap.RemoveMin();
            Assert.AreEqual("A", last.Value);
            Assert.AreEqual(10f, last.Weight);
        }

        /// <summary>
        /// Replacing an index that is not in the heap is a silent no-op.
        /// </summary>
        [TestMethod]
        public void ReplaceUnknownIndexIsNoOp()
        {
            var heap = new BinaryHeap<int>(5);
            heap.Add(new BinaryHeap<int>.Node(5f, 1, 1));

            heap.Replace(new BinaryHeap<int>.Node(1f, 99, 99));

            Assert.AreEqual(1, heap.Count);
            var node = heap.RemoveMin();
            Assert.AreEqual(1, node.Index);
            Assert.AreEqual(5f, node.Weight);
        }

        /// <summary>
        /// Decreasing the key of an index that is not in the heap adds the node.
        /// </summary>
        [TestMethod]
        public void DecreaseKeyUnknownIndexAdds()
        {
            var heap = new BinaryHeap<int>(5);
            heap.Add(new BinaryHeap<int>.Node(5f, 1, 1));

            heap.DecreaseKey(new BinaryHeap<int>.Node(2f, 7, 7));

            Assert.AreEqual(2, heap.Count);
            Assert.AreEqual(7, heap.RemoveMin().Index);
        }

        /// <summary>
        /// Fixed-seed fuzz of Add, RemoveMin, DecreaseKey, and Replace (including weight
        /// increases and unknown indices) against a linear-scan reference model, with a full
        /// ordered drain at the end.
        /// </summary>
        [TestMethod]
        public void HeapFuzzMatchesReferenceModel()
        {
            var randy = new MersenneTwister(12345);
            var heap = new BinaryHeap<double>(64);
            var model = new Dictionary<int, float>();

            for (int op = 0; op < 2000; op++)
            {
                int index = randy.Next(0, 64);
                double action = randy.NextDouble();
                float weight = (float)randy.NextDouble();

                if (action < 0.4)
                {
                    if (!model.ContainsKey(index))
                    {
                        heap.Add(new BinaryHeap<double>.Node(weight, index, index));
                        model[index] = weight;
                    }
                }
                else if (action < 0.6)
                {
                    if (model.Count == 0) continue;
                    var popped = heap.RemoveMin();
                    float modelMin = float.MaxValue;
                    foreach (var entry in model)
                    {
                        if (entry.Value < modelMin) modelMin = entry.Value;
                    }
                    Assert.AreEqual(modelMin, popped.Weight, 0f);
                    Assert.AreEqual(modelMin, model[popped.Index], 0f);
                    model.Remove(popped.Index);
                }
                else if (action < 0.8)
                {
                    heap.DecreaseKey(new BinaryHeap<double>.Node(weight, index, index));
                    if (!model.ContainsKey(index) || weight < model[index]) model[index] = weight;
                }
                else
                {
                    heap.Replace(new BinaryHeap<double>.Node(weight, index, index));
                    if (model.ContainsKey(index)) model[index] = weight;
                }

                Assert.AreEqual(model.Count, heap.Count);
            }

            while (model.Count > 0)
            {
                var popped = heap.RemoveMin();
                float modelMin = float.MaxValue;
                foreach (var entry in model)
                {
                    if (entry.Value < modelMin) modelMin = entry.Value;
                }
                Assert.AreEqual(modelMin, popped.Weight, 0f);
                Assert.AreEqual(modelMin, model[popped.Index], 0f);
                model.Remove(popped.Index);
            }
            Assert.AreEqual(0, heap.Count);
        }

    }
}
