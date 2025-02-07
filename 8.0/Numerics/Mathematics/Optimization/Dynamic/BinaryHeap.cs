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
        private int _n = 0; // Number of nodes.
        private int _p = 0; // Parent Index

        public int Count => _n;

        public BinaryHeap(int heapSize)
        {
            _heap = new Node[heapSize];
        }

        /// <summary>
        /// Add a node to the heap.
        /// </summary>
        /// <param name="value">Node to add.</param>
        public void Add(Node value)
        {
            _heap[_n] = value;
            _n += 1;

            //Push the node up the tree.
            if (_n <= 1) return;
            int i = _n - 1;
            _p = (i - 1) / 2; // parent index (use integer division).
            Node parentNode = _heap[_p];
            Node currentNode = value;

            while (i > 0 && currentNode.Weight < parentNode.Weight)
            {
                _heap[i] = parentNode;
                _heap[_p] = currentNode;
                i = _p;
                _p = (i - 1) / 2; // parent index (use integer division).
                parentNode = _heap[_p];
                currentNode = _heap[i];
            }
        }

        /// <summary>
        /// Remove the minimum (top) node from the heap.
        /// </summary>
        /// <returns></returns>
        public Node RemoveMin()
        {
            Node result = _heap[0];

            //Push the empty node down the tree then replace with last node in array
            int r;
            int l;
            int i = 0;
            int pIndex = 0;
            Node lNode;
            Node rNode;
            do
            {
                r = 2 * i + 2; //right child index
                l = r - 1; //left child index
                if (r < _n)
                {
                    lNode = _heap[l];
                    rNode = _heap[r];
                    if (lNode.Weight < rNode.Weight)
                    {
                        _heap[i] = lNode;
                        pIndex = i;
                        i = l;
                    }
                    else
                    {
                        _heap[i] = rNode;
                        pIndex = i;
                        i = r;
                    }
                }
                else
                {
                    if (l < _n)
                    {
                        _heap[i] = _heap[l];
                        pIndex = i;
                        i = l;
                    }
                    else
                    {
                        break;
                    }
                }
            } while (i >= 0);
            //
            _n -= 1;
            //Replace with last node in array and push up.
            if (i != _n)
            {
                _heap[i] = _heap[_n];
                if (_heap[pIndex].Weight > _heap[i].Weight)
                {
                    //Push the node up.
                    if (i <= 0) return result;
                    _p = pIndex;
                    Node parentNode = _heap[_p];
                    Node currentNode = _heap[i];
                    while (i > 0 && currentNode.Weight < parentNode.Weight)
                    {
                        _heap[i] = parentNode;
                        _heap[_p] = currentNode;
                        i = _p;
                        _p = (i - 1) / 2; // parent index (use integer division).
                        parentNode = _heap[_p];
                        currentNode = _heap[i];
                    }
                }
            }
            //
            return result;
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
                    //Push the node up the tree
                    if (i <= 0) return;
                    _p = (i - 1) / 2; // parent index (use integer division).
                    Node parentNode = _heap[_p];
                    Node currentNode = newNode;
                    while (i > 0 && currentNode.Weight < parentNode.Weight)
                    {
                        _heap[i] = parentNode;
                        _heap[_p] = currentNode;
                        i = _p;
                        _p = (i - 1) / 2;
                        parentNode = _heap[_p];
                        currentNode = _heap[i];
                    }
                    break;
                }
            }
        }


    }
}