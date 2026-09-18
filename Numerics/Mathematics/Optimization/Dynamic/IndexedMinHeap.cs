using System;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// An array-backed indexed binary min-heap specialized for the shortest-path solvers: node
    /// identifiers are dense integers in [0, capacity), so the index-to-slot map is a flat
    /// array rather than a dictionary, and the heap slots store the weight and node identifier
    /// in parallel arrays.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The sift operations use hole shifting (one write per level instead of a swap) and stop on
    /// equal weights, so the slot layout — and therefore the pop order among equal weights —
    /// evolves exactly as it does for the swap-based public <see cref="BinaryHeap{T}"/> given
    /// the same operation sequence. Capacity is fixed at construction; the shortest-path solvers
    /// size it at the node count, which the solver state machine makes an exact bound (at most
    /// one live entry per node). <see cref="Clear"/> resets only the occupied slots, so a heap
    /// can be reused across solves without reallocation.
    /// </para>
    /// </remarks>
    internal sealed class IndexedMinHeap
    {
        /// <summary>The heap-ordered weights.</summary>
        private readonly float[] _weights;

        /// <summary>The heap-ordered node identifiers.</summary>
        private readonly int[] _nodes;

        /// <summary>Maps a node identifier to its heap slot; -1 when the node is not in the heap.</summary>
        private readonly int[] _position;

        /// <summary>The number of nodes currently in the heap.</summary>
        private int _count;

        /// <summary>
        /// The number of nodes in the heap.
        /// </summary>
        internal int Count => _count;

        /// <summary>
        /// Creates a new indexed min-heap for node identifiers in [0, capacity).
        /// </summary>
        /// <param name="capacity">The maximum number of nodes the heap can hold, and the exclusive upper bound on node identifiers.</param>
        internal IndexedMinHeap(int capacity)
        {
            _weights = new float[capacity];
            _nodes = new int[capacity];
            _position = new int[capacity];
            for (int i = 0; i < capacity; i++) _position[i] = -1;
        }

        /// <summary>
        /// Empties the heap, resetting only the occupied position entries.
        /// </summary>
        internal void Clear()
        {
            for (int i = 0; i < _count; i++) _position[_nodes[i]] = -1;
            _count = 0;
        }

        /// <summary>
        /// Adds a node that is not currently in the heap.
        /// </summary>
        /// <param name="node">The node identifier.</param>
        /// <param name="weight">The node's weight (priority).</param>
        /// <exception cref="InvalidOperationException">Thrown when the heap is at capacity.</exception>
        internal void Add(int node, float weight)
        {
            if (_count >= _nodes.Length)
                throw new InvalidOperationException("Heap is full.");

            int slot = _count;
            _count++;
            SiftUp(slot, node, weight);
        }

        /// <summary>
        /// Removes the minimum-weight node from the heap.
        /// </summary>
        /// <param name="node">The removed node identifier.</param>
        /// <param name="weight">The removed node's weight.</param>
        /// <exception cref="InvalidOperationException">Thrown when the heap is empty.</exception>
        internal void RemoveMin(out int node, out float weight)
        {
            if (_count == 0)
                throw new InvalidOperationException("Heap is empty.");

            node = _nodes[0];
            weight = _weights[0];
            _position[node] = -1;

            _count--;
            if (_count > 0)
            {
                SiftDown(0, _nodes[_count], _weights[_count]);
            }
        }

        /// <summary>
        /// Lowers the weight of a node already in the heap; adds the node when it is absent.
        /// Weights that are not smaller than the current weight are ignored.
        /// </summary>
        /// <param name="node">The node identifier.</param>
        /// <param name="weight">The candidate weight.</param>
        internal void DecreaseKey(int node, float weight)
        {
            int slot = _position[node];
            if (slot < 0)
            {
                Add(node, weight);
                return;
            }
            if (weight >= _weights[slot]) return;

            SiftUp(slot, node, weight);
        }

        /// <summary>
        /// Moves the given entry up from the given slot until its parent is no heavier, using
        /// hole shifting, and writes the entry at its final slot.
        /// </summary>
        /// <param name="slot">The starting slot (treated as a hole).</param>
        /// <param name="node">The entry's node identifier.</param>
        /// <param name="weight">The entry's weight.</param>
        private void SiftUp(int slot, int node, float weight)
        {
            while (slot > 0)
            {
                int parent = (slot - 1) >> 1;
                if (weight >= _weights[parent]) break;

                _weights[slot] = _weights[parent];
                _nodes[slot] = _nodes[parent];
                _position[_nodes[slot]] = slot;
                slot = parent;
            }
            _weights[slot] = weight;
            _nodes[slot] = node;
            _position[node] = slot;
        }

        /// <summary>
        /// Moves the given entry down from the given slot until no child is lighter, using hole
        /// shifting, and writes the entry at its final slot. Child ties prefer the left child,
        /// matching the swap-based sift.
        /// </summary>
        /// <param name="slot">The starting slot (treated as a hole).</param>
        /// <param name="node">The entry's node identifier.</param>
        /// <param name="weight">The entry's weight.</param>
        private void SiftDown(int slot, int node, float weight)
        {
            while (true)
            {
                int left = 2 * slot + 1;
                if (left >= _count) break;
                int right = left + 1;
                int smallest = (right < _count && _weights[right] < _weights[left]) ? right : left;
                if (_weights[smallest] >= weight) break;

                _weights[slot] = _weights[smallest];
                _nodes[slot] = _nodes[smallest];
                _position[_nodes[slot]] = slot;
                slot = smallest;
            }
            _weights[slot] = weight;
            _nodes[slot] = node;
            _position[node] = slot;
        }
    }
}
