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
        /// <summary>
        /// Represents a node in the binary heap with a weight, index, and value.
        /// </summary>
        public struct Node
        {
            /// <summary>
            /// The weight (priority) of the node.
            /// </summary>
            public float Weight;
            /// <summary>
            /// The index identifier of the node.
            /// </summary>
            public int Index;
            /// <summary>
            /// The value stored in the node.
            /// </summary>
            public T Value;

            /// <summary>
            /// Creates a new node with the specified weight, index, and value.
            /// </summary>
            /// <param name="nodeWeight">The weight (priority) of the node.</param>
            /// <param name="nodeIndex">The index identifier of the node.</param>
            /// <param name="nodeValue">The value to store in the node.</param>
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

        /// <summary>
        /// The number of nodes in the heap.
        /// </summary>
        public int Count => _n;

        /// <summary>
        /// Creates a new binary heap with the specified maximum size.
        /// </summary>
        /// <param name="heapSize">The maximum number of nodes the heap can hold.</param>
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