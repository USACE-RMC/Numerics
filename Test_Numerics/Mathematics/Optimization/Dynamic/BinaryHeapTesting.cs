/**
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
* **/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    [TestClass]
    public class BinaryHeapTesting
    {
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

            Assert.AreEqual(heap.RemoveMin().Weight == weights[6], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[5], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[0], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[1], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[4], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[3], true);
            Assert.AreEqual(heap.RemoveMin().Weight == weights[2], true);
        }

        [TestMethod]
        public void HeapTest2()
        {
            //Random Node weights
            float[] weights = new float[1000]; //= new float[] { .3f, .5f, 32f, 15f, 12f, .01f, -4f };

            Random randy = new Random();
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
                Assert.AreEqual(heap.RemoveMin().Weight == weights[i], true);
            }
        }

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
                Assert.AreEqual(heap.RemoveMin().Weight == weights[i], true);
            }
        }

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
                Assert.AreEqual(heap.RemoveMin().Value == weights[i], false);
            }
        }
    }
}