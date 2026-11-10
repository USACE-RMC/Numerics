using System;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// An array-backed indexed binary min-heap specialized for the shortest-path solvers: node
    /// identifiers are dense integers in [0, capacity), so the index-to-slot map is a flat
    /// array rather than a dictionary, and each heap slot stores its weight and node identifier
    /// together.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The sift operations use hole shifting (one write per level instead of a swap) and stop on
    /// equal weights, so the slot layout, and therefore the pop order among equal weights,
    /// evolves exactly as it does for the swap-based public <see cref="BinaryHeap{T}"/> given
    /// the same operation sequence. Capacity is fixed at construction; the shortest-path solvers
    /// size it at the node count, which the solver state machine makes an exact bound (at most
    /// one live entry per node).
    /// </para>
    /// <para>
    /// Each node also keeps the last weight it held in the heap next to its slot, so a solver
    /// can use the heap as its distance store: one cache line then answers both whether a path
    /// is shorter and where the node sits in the heap. <see cref="Clear"/> resets every node,
    /// so a heap can be reused across solves without reallocation.
    /// </para>
    /// </remarks>
    internal sealed class IndexedMinHeap
    {
        /// <summary>A heap slot: a node identifier and its weight, stored together so a sift step reads one cache line.</summary>
        private struct Entry
        {
            /// <summary>The node's weight (priority).</summary>
            public double Weight;

            /// <summary>The node identifier.</summary>
            public int Node;
        }

        /// <summary>The per-node record: the last weight the node held in the heap, and its current slot.</summary>
        private struct NodeRecord
        {
            /// <summary>The last weight the node held in the heap; positive infinity when it has not been added.</summary>
            public double Weight;

            /// <summary>The node's heap slot; -1 when the node is not in the heap.</summary>
            public int Slot;
        }

        /// <summary>The heap-ordered entries.</summary>
        private readonly Entry[] _heap;

        /// <summary>The per-node records, indexed by node identifier.</summary>
        private readonly NodeRecord[] _records;

        /// <summary>The number of nodes currently in the heap.</summary>
        private int _count;

        /// <summary>
        /// The number of nodes in the heap.
        /// </summary>
        internal int Count => _count;

        /// <summary>
        /// The node identifier at the top of the heap, the next one <see cref="RemoveMin"/> returns.
        /// The heap must not be empty.
        /// </summary>
        internal int MinNode => _heap[0].Node;

        /// <summary>
        /// Creates a new indexed min-heap for node identifiers in [0, capacity).
        /// </summary>
        /// <param name="capacity">The maximum number of nodes the heap can hold, and the exclusive upper bound on node identifiers.</param>
        internal IndexedMinHeap(int capacity)
        {
            _heap = new Entry[capacity];
            _records = new NodeRecord[capacity];
            ResetRecords();
        }

        /// <summary>
        /// Empties the heap and forgets every node's last weight.
        /// </summary>
        internal void Clear()
        {
            _count = 0;
            ResetRecords();
        }

        /// <summary>
        /// Returns the last weight the node held in the heap: its current weight while it is in
        /// the heap, the weight it was removed with afterwards, and positive infinity when it has
        /// not been added since construction or the last <see cref="Clear"/>.
        /// </summary>
        /// <param name="node">The node identifier.</param>
        /// <returns>The node's last weight.</returns>
        internal double WeightOf(int node) => _records[node].Weight;

        /// <summary>
        /// Adds a node that is not currently in the heap.
        /// </summary>
        /// <param name="node">The node identifier.</param>
        /// <param name="weight">The node's weight (priority).</param>
        /// <exception cref="InvalidOperationException">Thrown when the heap is at capacity.</exception>
        internal void Add(int node, double weight)
        {
            if (_count >= _heap.Length)
                throw new InvalidOperationException("Heap is full.");

            int slot = _count;
            _count++;
            _records[node].Weight = weight;
            SiftUp(slot, node, weight);
        }

        /// <summary>
        /// Removes the minimum-weight node from the heap.
        /// </summary>
        /// <param name="node">The removed node identifier.</param>
        /// <param name="weight">The removed node's weight.</param>
        /// <exception cref="InvalidOperationException">Thrown when the heap is empty.</exception>
        internal void RemoveMin(out int node, out double weight)
        {
            if (_count == 0)
                throw new InvalidOperationException("Heap is empty.");

            Entry top = _heap[0];
            node = top.Node;
            weight = top.Weight;
            _records[node].Slot = -1;

            _count--;
            if (_count > 0)
            {
                Entry last = _heap[_count];
                SiftDown(0, last.Node, last.Weight);
            }
        }

        /// <summary>
        /// Lowers the weight of a node already in the heap; adds the node when it is absent.
        /// Weights that are not smaller than the current weight of a node in the heap are ignored.
        /// </summary>
        /// <param name="node">The node identifier.</param>
        /// <param name="weight">The candidate weight.</param>
        internal void DecreaseKey(int node, double weight)
        {
            ref NodeRecord record = ref _records[node];
            int slot = record.Slot;
            if (slot < 0)
            {
                Add(node, weight);
                return;
            }
            if (weight >= record.Weight) return;

            record.Weight = weight;
            SiftUp(slot, node, weight);
        }

        /// <summary>
        /// Marks every node absent, with positive infinite weight.
        /// </summary>
        private void ResetRecords()
        {
            NodeRecord[] records = _records;
            for (int i = 0; i < records.Length; i++)
            {
                records[i].Weight = double.PositiveInfinity;
                records[i].Slot = -1;
            }
        }

        /// <summary>
        /// Moves the given entry up from the given slot until its parent is no heavier, using
        /// hole shifting, and writes the entry at its final slot.
        /// </summary>
        /// <param name="slot">The starting slot (treated as a hole).</param>
        /// <param name="node">The entry's node identifier.</param>
        /// <param name="weight">The entry's weight.</param>
        private void SiftUp(int slot, int node, double weight)
        {
            Entry[] heap = _heap;
            NodeRecord[] records = _records;
            while (slot > 0)
            {
                int parent = (slot - 1) >> 1;
                Entry above = heap[parent];
                if (weight >= above.Weight) break;

                heap[slot] = above;
                records[above.Node].Slot = slot;
                slot = parent;
            }
            heap[slot].Weight = weight;
            heap[slot].Node = node;
            records[node].Slot = slot;
        }

        /// <summary>
        /// Moves the given entry down from the given slot until no child is lighter, using hole
        /// shifting, and writes the entry at its final slot. Child ties prefer the left child,
        /// matching the swap-based sift.
        /// </summary>
        /// <param name="slot">The starting slot (treated as a hole).</param>
        /// <param name="node">The entry's node identifier.</param>
        /// <param name="weight">The entry's weight.</param>
        /// <remarks>
        /// The sift is bottom-up (Floyd's form): the hole first descends along the lighter-child
        /// path to a leaf without comparing against the entry, then climbs back while the
        /// entry is no heavier than the element above it. Weights along the lighter-child path never
        /// decrease, so the entry lands on the same slot as the top-down sift, which stops at the
        /// first lighter child that is no lighter than the entry, and every slot ends up
        /// holding the same element. The descent has no data-dependent branch, and the climb is
        /// usually short because the entry, taken from the last slot, is usually heavy.
        /// </remarks>
        private void SiftDown(int slot, int node, double weight)
        {
            Entry[] heap = _heap;
            NodeRecord[] records = _records;
            int count = _count;
            int hole = slot;
            while (true)
            {
                int left = 2 * hole + 1;
                if (left >= count) break;
                int child = left;
                int right = left + 1;
                if (right < count) child += heap[right].Weight < heap[left].Weight ? 1 : 0;
                Entry below = heap[child];
                heap[hole] = below;
                records[below.Node].Slot = hole;
                hole = child;
            }
            while (hole > slot)
            {
                int parent = (hole - 1) >> 1;
                Entry above = heap[parent];
                if (above.Weight < weight) break;

                heap[hole] = above;
                records[above.Node].Slot = hole;
                hole = parent;
            }
            heap[hole].Weight = weight;
            heap[hole].Node = node;
            records[node].Slot = hole;
        }
    }
}
