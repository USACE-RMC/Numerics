"""Shortest paths on directed networks from ``Numerics.Mathematics.Optimization``.

:class:`Edge` and :class:`Network` describe a directed network; the ``dijkstra_*`` functions and
:meth:`Network.solve` return a routing table with one row per node: column 0 is the next node
toward the destination, column 1 the index of the edge to take and column 2 the cumulative cost
(float32; the costs accumulate in float64 and are rounded once). The destination row is
``(itself, -1, 0)`` and an unreachable node's row is ``(-1, -1, inf)``. :func:`shortest_distances` gives
float64 distances from several sources on a CSR graph, as ``scipy.sparse.csgraph.dijkstra`` does.

:class:`BinaryHeap` is the keyed min-heap the shortest-path solvers are built on. The knapsack
and dynamic programming solvers of the same namespace are in :mod:`pynumerics.optimization`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pynumerics._util import _opt
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _net_rectangular,
    _numpy_rectangular,
    to_net_doubles,
    to_net_ints,
    to_net_singles,
    to_numpy_doubles,
    to_numpy_ints,
)


def _table_to_numpy(table: Any) -> np.ndarray:
    """Copy a .NET ``float[,]`` routing table into a ``(nodes, 3)`` float32 array, bit for bit, in one copy."""
    return _numpy_rectangular(table, np.float32)


def _table_to_net(table: Any) -> Any:
    """Copy a ``(nodes, 3)`` array into a new .NET ``float[,]`` (rounded to float32) in one copy."""
    array = np.asarray(table, dtype=np.float32)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2-D routing table; got shape {array.shape}.")
    return _net_rectangular(array)


def _singles(values: Any) -> Any:
    """Copy edge weights into a ``System.Single[]`` (rounded to float32)."""
    array = np.asarray(values, dtype=np.float32).reshape(1, 1, -1)
    return to_net_singles(array)[0]


def _path(path: Any) -> np.ndarray | None:
    """Copy a .NET ``List<int>`` of edge indices, or keep None."""
    return None if path is None else to_numpy_ints(path.ToArray())


@dataclass(frozen=True, slots=True)
class Edge:
    """A directed, weighted edge of a network (``Numerics.Mathematics.Optimization.Edge``).

    Parameters
    ----------
    from_index : int
        The start node.
    to_index : int
        The end node.
    weight : float
        The cost of the edge (stored as float32 by Numerics); infinity makes it impassable.
    index : int
        The edge index reported in paths and routing tables.
    """

    from_index: int
    to_index: int
    weight: float
    index: int

    @classmethod
    def from_net(cls, net: Any) -> Edge:
        """Copy a .NET ``Edge``.

        Parameters
        ----------
        net : Numerics.Mathematics.Optimization.Edge
            The .NET edge.

        Returns
        -------
        Edge
            The copy.
        """
        return cls(int(net.FromIndex), int(net.ToIndex), float(net.Weight), int(net.Index))

    def to_net(self) -> Any:
        """Return the .NET ``Edge``.

        Returns
        -------
        Numerics.Mathematics.Optimization.Edge
            A new .NET edge.
        """
        return _opt().Edge(int(self.from_index), int(self.to_index), float(self.weight), int(self.index))


def _edges(edges: Any) -> Any:
    """Return the edges as a .NET ``Edge[]``.

    Accepts :class:`Edge` values, ``(from, to, weight)`` or ``(from, to, weight, index)`` tuples (a
    missing index is the position), or a .NET ``Edge[]`` passed through.
    """
    opt = _opt()
    import System  # noqa: PLC0415 - the runtime is running

    if isinstance(edges, System.Array):
        return edges
    items = []
    for position, edge in enumerate(edges):
        if isinstance(edge, Edge):
            items.append(edge.to_net())
            continue
        fields = tuple(edge)
        index = int(fields[3]) if len(fields) > 3 else position
        items.append(opt.Edge(int(fields[0]), int(fields[1]), float(fields[2]), index))
    return System.Array[opt.Edge](items)


def _destinations(destinations: Any) -> Any:
    """Return one destination as an int or several as a ``System.Int32[]``."""
    if np.ndim(destinations) == 0:
        return int(destinations)
    return to_net_ints(destinations)


def dijkstra_solve(edges: Any, destinations: Any, node_count: int = -1) -> np.ndarray:
    """Return the shortest-path routing table toward one or more destinations (``Dijkstra.Solve``).

    Parameters
    ----------
    edges : sequence
        :class:`Edge` values or ``(from, to, weight[, index])`` tuples.
    destinations : int or array_like of int
        One destination node, or several; each node routes to its cheapest one, earlier destinations
        winning exact ties.
    node_count : int, optional
        Number of nodes; -1, the default, derives it from the edges.

    Returns
    -------
    numpy.ndarray
        The ``(nodes, 3)`` float32 routing table.
    """
    return _table_to_numpy(_opt().Dijkstra.Solve(_edges(edges), _destinations(destinations), int(node_count)))


def dijkstra_solve_nearest(edges: Any, destinations: Any, node_count: int = -1) -> np.ndarray:
    """Return the routing table to the nearest destination in one multi-source pass (``Dijkstra.SolveNearest``).

    Costs equal those of :func:`dijkstra_solve`; exact ties may route differently.

    Parameters
    ----------
    edges : sequence
        :class:`Edge` values or ``(from, to, weight[, index])`` tuples.
    destinations : array_like of int
        The destination nodes, at least one.
    node_count : int, optional
        Number of nodes; -1, the default, derives it from the edges.

    Returns
    -------
    numpy.ndarray
        The ``(nodes, 3)`` float32 routing table.
    """
    return _table_to_numpy(_opt().Dijkstra.SolveNearest(_edges(edges), to_net_ints(np.atleast_1d(destinations)),
                                                        int(node_count)))


def dijkstra_path_exists(table: Any, node: int) -> bool:
    """Return whether a node reaches the destination in a routing table (``Dijkstra.PathExists``).

    Parameters
    ----------
    table : array_like
        A ``(nodes, 3)`` routing table.
    node : int
        The node.

    Returns
    -------
    bool
        True when the node's cost is finite.
    """
    return bool(_opt().Dijkstra.PathExists(_table_to_net(table), int(node)))


def dijkstra_get_path(table: Any, start: int) -> np.ndarray | None:
    """Return the edge indices from a node to the destination (``Dijkstra.GetPath``).

    Parameters
    ----------
    table : array_like
        A ``(nodes, 3)`` routing table.
    start : int
        The start node.

    Returns
    -------
    numpy.ndarray or None
        The ordered int32 edge indices (empty at the destination), or None when unreachable.
    """
    return _path(_opt().Dijkstra.GetPath(_table_to_net(table), int(start)))


def dijkstra_try_get_path(table: Any, start: int) -> tuple[bool, np.ndarray, float]:
    """Walk a routing table from a node, returning the path and its cost (``Dijkstra.TryGetPath``).

    Parameters
    ----------
    table : array_like
        A ``(nodes, 3)`` routing table.
    start : int
        The start node.

    Returns
    -------
    tuple
        ``(found, edge_indices, total_cost)``; an unreachable node gives ``(False, [], inf)``.
    """
    found, path, cost = _opt().Dijkstra.TryGetPath(_table_to_net(table), int(start), None, 0.0)
    return bool(found), to_numpy_ints(path.ToArray()), float(cost)


def shortest_distances(graph: Any, sources: Any, limit: float = np.inf, *, return_predecessors: bool = False,
                       return_nearest_source: bool = False) -> Any:
    """Return the double-precision distance of every node from its nearest source (``Dijkstra.ShortestDistances``).

    The graph is directed, in compressed sparse row form: row ``n`` holds the edges leaving node ``n``. The
    result matches ``scipy.sparse.csgraph.dijkstra(graph, indices=sources, min_only=True, limit=limit)``;
    where SciPy reports -9999 for a missing predecessor or source, this function reports -1.

    Parameters
    ----------
    graph : sparse matrix or tuple
        A square SciPy sparse matrix or array (converted with ``tocsr()``; any other object with ``shape``,
        ``indptr``, ``indices`` and ``data`` is read as CSR), or a ``(row_start, to_node, weight)`` tuple of
        the CSR arrays. Weights must be non-negative; infinity makes
        an edge impassable.
    sources : int or array_like of int
        The source nodes, each at distance zero. Integral floats are accepted; other floats raise ValueError.
    limit : float, optional
        The largest distance computed; farther nodes are reported unreachable. Default infinity.
    return_predecessors : bool, optional
        Also return the predecessor of each node on its shortest path (-1 for a source or an unreachable node).
    return_nearest_source : bool, optional
        Also return the source each shortest path starts from (-1 for an unreachable node).

    Returns
    -------
    numpy.ndarray or tuple
        The float64 distances, infinity where unreachable; with either flag, a tuple ``(distances, ...)``
        holding the requested int32 arrays in the order predecessors, nearest source.

    Raises
    ------
    ValueError
        If the matrix is not square, a source is not integral, or the CSR arrays are inconsistent (raised by
        Numerics).
    TypeError
        If the sources are not numbers.
    """
    dijkstra = _opt().Dijkstra
    import System  # noqa: PLC0415 - the namespace import above starts the runtime

    if isinstance(graph, tuple):
        row_start, to_node, weight = graph
        node_count = len(row_start) - 1
    else:
        rows, cols = graph.shape
        if rows != cols:
            raise ValueError(f"The graph matrix must be square; got shape {graph.shape}.")
        if hasattr(graph, "tocsr"):
            # A CSC (or other) sparse matrix also has indptr/indices/data, which would read as the transpose.
            graph = graph.tocsr()
        row_start, to_node, weight, node_count = graph.indptr, graph.indices, graph.data, rows
    source_array = np.atleast_1d(np.asarray(sources))
    if source_array.dtype.kind == "f":
        if not np.all(np.isfinite(source_array) & (source_array == np.floor(source_array))):
            raise ValueError("The sources must be integral node indices.")
    elif source_array.dtype.kind not in "iu":
        raise TypeError(f"The sources must be integer node indices; got dtype {source_array.dtype}.")
    # The typed array constructor binds once; Array.CreateInstance resolves its overloads on every call.
    ints = System.Array[System.Int32]
    predecessors = ints(node_count) if return_predecessors else None
    nearest = ints(node_count) if return_nearest_source else None
    distances = to_numpy_doubles(dijkstra.ShortestDistances(
        int(node_count), to_net_ints(row_start), to_net_ints(to_node), to_net_doubles(weight),
        to_net_ints(source_array), float(limit), predecessors, nearest))
    if not (return_predecessors or return_nearest_source):
        return distances
    result = [distances]
    if return_predecessors:
        result.append(to_numpy_ints(predecessors))
    if return_nearest_source:
        result.append(to_numpy_ints(nearest))
    return tuple(result)


class Network:
    """A compiled directed network with fixed destinations (``Numerics.Mathematics.Optimization.Network``).

    Parameters
    ----------
    edges : sequence
        :class:`Edge` values or ``(from, to, weight[, index])`` tuples.
    destination_indices : array_like of int
        The destination nodes.

    Attributes
    ----------
    net : Numerics.Mathematics.Optimization.Network
        The wrapped .NET network.
    """

    def __init__(self, edges: Any, destination_indices: Any) -> None:
        self.net = _opt().Network(_edges(edges), to_net_ints(np.atleast_1d(destination_indices)))
        self._table: Any = None

    @property
    def destination_indices(self) -> np.ndarray:
        """The destination nodes (``DestinationIndices``)."""
        return to_numpy_ints(self.net.DestinationIndices)

    @property
    def incoming_edges(self) -> list[list[Edge]]:
        """The incoming edges of each node, in edge order (``IncomingEdges``)."""
        return [[Edge.from_net(e) for e in node] for node in self.net.IncomingEdges]

    @property
    def outgoing_edges(self) -> list[list[Edge]]:
        """The outgoing edges of each node, in edge order (``OutgoingEdges``)."""
        return [[Edge.from_net(e) for e in node] for node in self.net.OutgoingEdges]

    @property
    def node_count(self) -> int:
        """The number of nodes (``NodeCount``)."""
        return int(self.net.NodeCount)

    def _reused_table(self) -> Any:
        """Return the .NET table the reuse overloads fill, allocated once per network."""
        if self._table is None:
            import System  # noqa: PLC0415 - the runtime is running

            self._table = System.Array.CreateInstance(System.Single, self.node_count, 3)
        return self._table

    def solve(self, destinations: Any = None, *, edge_weights: Any = None, reuse: bool = False) -> np.ndarray:
        """Return the routing table toward the destinations (``Network.Solve``).

        Parameters
        ----------
        destinations : int or array_like of int, optional
            Destinations to solve for; None uses the network's. Not allowed with ``edge_weights``.
        edge_weights : array_like, optional
            Custom weights, one per edge in constructor order, toward the network's destinations.
        reuse : bool, optional
            With ``edge_weights``, fill one .NET table kept by the network (the allocation-free
            ``Solve(float[], float[,])`` overload) instead of a new one; default False.

        Returns
        -------
        numpy.ndarray
            The ``(nodes, 3)`` float32 routing table.

        Raises
        ------
        ValueError
            If both ``destinations`` and ``edge_weights`` are given.
        """
        if edge_weights is not None:
            if destinations is not None:
                raise ValueError("Custom edge weights solve toward the network's destinations; omit destinations.")
            if reuse:
                table = self._reused_table()
                self.net.Solve(_singles(edge_weights), table)
                return _table_to_numpy(table)
            return _table_to_numpy(self.net.Solve(_singles(edge_weights)))
        if destinations is None:
            destinations = self.net.DestinationIndices
            return _table_to_numpy(self.net.Solve(destinations))
        return _table_to_numpy(self.net.Solve(_destinations(destinations)))

    def solve_nearest(self, *, edge_weights: Any = None, reuse: bool = False) -> np.ndarray:
        """Return the routing table to the nearest destination in one pass (``Network.SolveNearest``).

        Parameters
        ----------
        edge_weights : array_like, optional
            Custom weights, one per edge; None uses the constructor's.
        reuse : bool, optional
            With ``edge_weights``, fill one .NET table kept by the network; default False.

        Returns
        -------
        numpy.ndarray
            The ``(nodes, 3)`` float32 routing table.
        """
        if edge_weights is None:
            return _table_to_numpy(self.net.SolveNearest())
        table = self._reused_table() if reuse else _table_to_net(np.zeros((self.node_count, 3)))
        self.net.SolveNearest(_singles(edge_weights), table)
        return _table_to_numpy(table)

    def get_path(self, edges_to_remove: Any, start: int, table: Any = None,
                 edge_weights: Any = None) -> np.ndarray | None:
        """Return the shortest path from a node avoiding some edges (``Network.GetPath``).

        Parameters
        ----------
        edges_to_remove : array_like of int
            Edge indices to exclude.
        start : int
            The start node.
        table : array_like, optional
            A routing table solved on this network without exclusions; its route is reused when it
            avoids every excluded edge.
        edge_weights : array_like, optional
            The custom weights ``table`` was solved with.

        Returns
        -------
        numpy.ndarray or None
            The ordered int32 edge indices. Without ``table`` an unreachable start gives None; with a
            table it gives an empty array, as Numerics does.
        """
        removed = to_net_ints(np.asarray(edges_to_remove, dtype=np.int64).ravel())
        if table is None:
            return _path(self.net.GetPath(removed, int(start)))
        if edge_weights is None:
            return _path(self.net.GetPath(removed, int(start), _table_to_net(table)))
        return _path(self.net.GetPath(removed, int(start), _table_to_net(table), _singles(edge_weights)))


@dataclass(frozen=True, slots=True)
class HeapNode:
    """A node of a :class:`BinaryHeap` (``BinaryHeap<T>.Node``).

    Parameters
    ----------
    weight : float
        The priority (stored as float32 by Numerics).
    index : int
        The key the keyed operations address, in [0, heap size).
    value : object, optional
        Any value carried with the node.
    """

    weight: float
    index: int
    value: Any = None


class BinaryHeap:
    """A fixed-capacity keyed binary min-heap (``Numerics.Mathematics.Optimization.BinaryHeap<object>``).

    Parameters
    ----------
    heap_size : int
        The capacity; node indices lie in [0, heap_size).

    Attributes
    ----------
    net : Numerics.Mathematics.Optimization.BinaryHeap
        The wrapped .NET heap of ``object`` values.
    """

    def __init__(self, heap_size: int) -> None:
        import clr  # noqa: PLC0415 - the namespace import starts the runtime
        import System  # noqa: PLC0415

        heap_type = _opt().BinaryHeap[System.Object]
        self.net = heap_type(int(heap_size))
        node = clr.GetClrType(heap_type).GetNestedType("Node")
        self._node_type = node.MakeGenericType(System.Array[System.Type]([clr.GetClrType(System.Object)]))

    def _node(self, node: HeapNode) -> Any:
        """Build the .NET node."""
        import System  # noqa: PLC0415 - the runtime is running

        arguments = System.Array[System.Object]([System.Single(float(node.weight)), System.Int32(int(node.index)),
                                                 node.value])
        return System.Activator.CreateInstance(self._node_type, arguments)

    @property
    def count(self) -> int:
        """The number of nodes in the heap (``Count``)."""
        return int(self.net.Count)

    def __len__(self) -> int:
        """Return the number of nodes in the heap.

        Returns
        -------
        int
            The count.
        """
        return self.count

    def add(self, node: HeapNode) -> None:
        """Add a node (``Add``).

        Parameters
        ----------
        node : HeapNode
            The node.
        """
        self.net.Add(self._node(node))

    def remove_min(self) -> HeapNode:
        """Remove and return the lightest node (``RemoveMin``).

        Returns
        -------
        HeapNode
            The node with the smallest weight.
        """
        node = self.net.RemoveMin()
        return HeapNode(float(node.Weight), int(node.Index), node.Value)

    def decrease_key(self, node: HeapNode) -> None:
        """Lower the weight of the node with the same index, or add it when absent (``DecreaseKey``).

        Parameters
        ----------
        node : HeapNode
            The index to address and the candidate weight; a larger weight is ignored.
        """
        self.net.DecreaseKey(self._node(node))

    def replace(self, node: HeapNode) -> None:
        """Replace the weight and value of the node with the same index (``Replace``).

        Parameters
        ----------
        node : HeapNode
            The index to address and the replacement weight and value; an absent index does nothing.
        """
        self.net.Replace(self._node(node))


__all__ = [
    "BinaryHeap", "Edge", "HeapNode", "Network", "dijkstra_get_path", "dijkstra_path_exists", "dijkstra_solve",
    "dijkstra_solve_nearest", "dijkstra_try_get_path", "shortest_distances",
]

translate_public(globals())
