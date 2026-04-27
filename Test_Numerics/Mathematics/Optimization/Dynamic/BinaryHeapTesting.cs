using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    [TestClass]
    public class BinaryHeapTesting
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


    }
}