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
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// This is an implementation of the binary heap data structure. The binary heap is especially convenient for shortest path algorithms 
    /// such as Djikstra's shortest path.
    /// source of inspiration: http://opendatastructures.org/versions/edition-0.1e/ods-java/10_1_BinaryHeap_Implicit_Bi.html
    /// </summary>
    /// <typeparam name="T">Generic variable to store with each node. Typically used to store important data associated with the network that isn't required for the binary heap.</typeparam>
    public class BinaryHeap<T>
    {
        public struct Node
        {
            public float Weight;
            public int Index;
            public T Value;

            public Node(float nodeWeight, int nodeIndex, T nodeValue)
            {
                Weight = nodeWeight;
                Index = nodeIndex;
                Value = nodeValue;
            }
        }

        private readonly Node[] _heap;
        private readonly Dictionary<int, int> _positionMap = new();

        private int _n = 0; // Number of nodes.
        //private int _p = 0; // Parent Index

        public int Count => _n;

        public BinaryHeap(int heapSize)
        {
            _heap = new Node[heapSize];
        }

        /// <summary>
        /// Putting the new inem in the first vacant cell in the array. 
        /// Then move it up in the heap based on its value compared to its parent.
        /// </summary>
        /// <param name="i">Index of the node.</param>
        private void BubbleUp(int i)
        {
            while (i > 0)
            {
                int parent = (i - 1) / 2;
                if (_heap[i].Weight >= _heap[parent].Weight) break;

                //Swap
                (_heap[i], _heap[parent]) = (_heap[parent], _heap[i]);
                

                _positionMap[_heap[i].Index] = i;
                _positionMap[_heap[parent].Index] = parent;
                i = parent;
            }
        }

        /// <summary>
        /// Used in heap deletion. Compares the parent nodes with child nodes in subtree.
        /// </summary>
        /// <param name="i"></param>
        private void BubbleDown(int i)
        {
            while (true)
            {
                int left = 2 * i + 1;
                int right = 2 * i + 2;
                int smallest = i; 

                if (left <_n && _heap[left].Weight < _heap[smallest].Weight)
                    smallest = left;
                if (right < _n && _heap[right].Weight < _heap[smallest].Weight)
                    smallest = right;
                if (smallest == i) break; 

                (_heap[i], _heap[smallest]) = (_heap[smallest], _heap[i]);
                _positionMap[_heap[i].Index] = i;
                _positionMap[_heap[smallest].Index] = smallest;

                i = smallest;
            }
        }

        /// <summary>
        /// Updates the distance (priority) of a node if a shorter path is found.
        /// </summary>
        /// <param name="newNode"></param>
        public void DecreaseKey(Node newNode)
        {
            if(!_positionMap.TryGetValue(newNode.Index, out int position))
            {
                Add(newNode);
                return;
            }
            if (newNode.Weight >= _heap[position].Weight) return; // No need to decrease key if the new weight is not smaller.

            _heap[position] = newNode;
            BubbleUp(position);
        }

        /// <summary>
        /// Add a node to the heap.
        /// </summary>
        public void Add(Node node)
        {
            if (_n >= _heap.Length) 
                throw new InvalidOperationException("Heap is full.");

            _heap[_n] = node;
            _positionMap[node.Index] = _n; // Map the index to the position in the heap array (for Replace method)
            BubbleUp(_n);
            _n++;
        }

        /// <summary>
        /// Remove the minimum (top) node from the heap.
        /// </summary>
        /// <returns></returns>
        public Node RemoveMin()
        {
            if (_n == 0) 
                throw new InvalidOperationException("Heap is empty.");

            Node min = _heap[0];
            _positionMap.Remove(min.Index);

            _n--;

            if (_n > 0)
            {
                _heap[0] = _heap[_n];
                _positionMap[_heap[0].Index] = 0; // Update position map
                BubbleDown(0);
            }
            return min;
        }

        /// <summary>
        /// Replace a node that has the same index value as the new node.
        /// </summary>
        /// <param name="newNode"></param>
        public void Replace(Node newNode)
        {
            for (int i = 0; i < _n; i++)
            {
                if (_heap[i].Index == newNode.Index)
                {
                    _heap[i] = newNode;
                    BubbleUp(i);
                    break;
                }
            }
        }


    }
}