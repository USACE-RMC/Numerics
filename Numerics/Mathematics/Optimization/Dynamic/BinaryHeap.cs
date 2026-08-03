using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// An array-backed binary min-heap keyed by node weight, with an index-to-position map that
    /// supports keyed decrease-key and replace operations. The heap is especially convenient for
    /// shortest path algorithms such as Dijkstra's method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The heap has a fixed capacity set at construction; adding past the capacity throws. Node
    /// weights may be any finite ordering value, including negatives — the heap orders by weight
    /// and imposes no algorithmic precondition of its own. The keyed operations
    /// (<see cref="DecreaseKey"/> and <see cref="Replace"/>) address nodes by their
    /// <see cref="Node.Index"/> and assume at most one live node per index; behavior with
    /// duplicate indices is unsupported for those operations.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href="http://opendatastructures.org/versions/edition-0.1e/ods-java/10_1_BinaryHeap_Implicit_Bi.html"/>
    /// </description></item>
    /// </list>
    /// </remarks>
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

        /// <summary>
        /// The heap slots in implicit binary-tree order.
        /// </summary>
        private readonly Node[] _heap;

        /// <summary>
        /// Maps a node index to its current slot in <see cref="_heap"/>.
        /// </summary>
        private readonly Dictionary<int, int> _positionMap = new();

        /// <summary>
        /// The number of nodes currently in the heap.
        /// </summary>
        private int _n = 0;

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
        /// Moves the node at the given slot up the heap until its parent is no heavier.
        /// </summary>
        /// <param name="i">The slot of the node to sift up.</param>
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
        /// Moves the node at the given slot down the heap until no child is lighter.
        /// </summary>
        /// <param name="i">The slot of the node to sift down.</param>
        private void BubbleDown(int i)
        {
            while (true)
            {
                int left = 2 * i + 1;
                int right = 2 * i + 2;
                int smallest = i;

                if (left < _n && _heap[left].Weight < _heap[smallest].Weight)
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
        /// Updates the weight (priority) of the node with the same index if the new weight is
        /// smaller; adds the node when its index is not in the heap. Larger weights are ignored.
        /// </summary>
        /// <param name="newNode">The node carrying the index to address and the candidate weight.</param>
        public void DecreaseKey(Node newNode)
        {
            if (!_positionMap.TryGetValue(newNode.Index, out int position))
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
        /// <param name="node">The node to add.</param>
        /// <exception cref="InvalidOperationException">Thrown when the heap is at capacity.</exception>
        public void Add(Node node)
        {
            if (_n >= _heap.Length)
                throw new InvalidOperationException("Heap is full.");

            _heap[_n] = node;
            _positionMap[node.Index] = _n; // Map the index to the position in the heap array for the keyed operations.
            BubbleUp(_n);
            _n++;
        }

        /// <summary>
        /// Remove the minimum (top) node from the heap.
        /// </summary>
        /// <returns>The node with the smallest weight.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the heap is empty.</exception>
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
        /// Replaces the node that has the same index value as the new node, restoring heap order
        /// in either direction. Does nothing when the index is not in the heap.
        /// </summary>
        /// <param name="newNode">The node carrying the index to address and the replacement weight and value.</param>
        public void Replace(Node newNode)
        {
            if (!_positionMap.TryGetValue(newNode.Index, out int position)) return;

            float previousWeight = _heap[position].Weight;
            _heap[position] = newNode;
            if (newNode.Weight < previousWeight)
            {
                BubbleUp(position);
            }
            else
            {
                BubbleDown(position);
            }
        }

    }
}
