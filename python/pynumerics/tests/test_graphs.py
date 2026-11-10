"""Shortest paths and the binary heap against scipy and textbook values.

Dijkstra routing tables are checked against ``scipy.sparse.csgraph.shortest_path`` and the heap
against ``heapq``.
"""

from __future__ import annotations

import heapq

import numpy as np
import pytest

from pynumerics import callbacks, graphs
from pynumerics.graphs import BinaryHeap, Edge, HeapNode, Network

csgraph = pytest.importorskip("scipy.sparse.csgraph")
sparse = pytest.importorskip("scipy.sparse")
needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")

# A textbook network (CLRS figure 24.6 with the source as node 0), costs to the destination 4.
TEXTBOOK = [(0, 1, 10.0), (0, 3, 5.0), (1, 2, 1.0), (1, 3, 2.0), (2, 4, 4.0), (3, 1, 3.0), (3, 2, 9.0),
            (3, 4, 2.0), (4, 0, 7.0), (4, 2, 6.0)]


def _random_edges(nodes: int, count: int, seed: int) -> list[tuple[int, int, float]]:
    rng = np.random.default_rng(seed)
    pairs = {(int(a), int(b)) for a, b in rng.integers(0, nodes, size=(count, 2)) if a != b}
    return [(a, b, float(rng.integers(1, 20))) for a, b in sorted(pairs)]


def _scipy_costs(edges: list[tuple[int, int, float]], nodes: int, removed: tuple[int, ...] = ()) -> np.ndarray:
    keep = [e for i, e in enumerate(edges) if i not in removed]
    best: dict[tuple[int, int], float] = {}
    for a, b, w in keep:
        best[(a, b)] = min(w, best.get((a, b), np.inf))
    rows, cols = zip(*best, strict=True) if best else ((), ())
    matrix = sparse.csr_matrix((list(best.values()), (rows, cols)), shape=(nodes, nodes))
    return csgraph.shortest_path(matrix, method="D", directed=True)


def test_textbook_network_costs_and_paths() -> None:
    # Act
    table = graphs.dijkstra_solve(TEXTBOOK, 4)
    found, path, cost = graphs.dijkstra_try_get_path(table, 0)
    # Assert
    np.testing.assert_array_equal(table[:, 2], [7.0, 4.0, 4.0, 2.0, 0.0])
    np.testing.assert_array_equal(table[4], [4.0, -1.0, 0.0])
    assert table.dtype == np.float32
    assert found and cost == 7.0 and path.tolist() == [1, 7]
    np.testing.assert_array_equal(graphs.dijkstra_get_path(table, 0), [1, 7])
    assert graphs.dijkstra_path_exists(table, 0)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_dijkstra_costs_match_scipy_shortest_path(seed) -> None:
    # Arrange
    edges = _random_edges(60, 200, seed)
    nodes = 60
    reference = _scipy_costs(edges, nodes)
    destinations = [5, 17, 42]
    # Act
    single = graphs.dijkstra_solve(edges, 5, nodes)
    merged = graphs.dijkstra_solve(edges, destinations, nodes)
    nearest = graphs.dijkstra_solve_nearest(edges, destinations, nodes)
    # Assert
    np.testing.assert_array_equal(single[:, 2], reference[:, 5].astype(np.float32))
    np.testing.assert_array_equal(merged[:, 2], reference[:, destinations].min(axis=1).astype(np.float32))
    np.testing.assert_array_equal(nearest[:, 2], merged[:, 2])
    for node in range(nodes):
        path = graphs.dijkstra_get_path(single, node)
        if np.isinf(reference[node, 5]):
            assert path is None and not graphs.dijkstra_path_exists(single, node)
            assert graphs.dijkstra_try_get_path(single, node)[0] is False
            np.testing.assert_array_equal(single[node], [-1.0, -1.0, np.inf])
        else:
            assert sum(edges[i][2] for i in path) == reference[node, 5]


def test_edges_accept_dataclasses_tuples_with_indices_and_net_arrays() -> None:
    # Arrange
    edges = [Edge(a, b, w, 100 + i) for i, (a, b, w) in enumerate(TEXTBOOK)]
    tuples = [(a, b, w, 100 + i) for i, (a, b, w) in enumerate(TEXTBOOK)]
    # Act
    net_edges = graphs._edges(edges)
    from_dataclasses = graphs.dijkstra_solve(edges, 4)
    from_tuples = graphs.dijkstra_solve(tuples, 4)
    from_net = graphs.dijkstra_solve(net_edges, 4)
    # Assert
    np.testing.assert_array_equal(from_dataclasses, from_tuples)
    np.testing.assert_array_equal(from_dataclasses, from_net)
    assert from_dataclasses[0, 1] == 101.0
    assert Edge.from_net(net_edges[3]) == edges[3]


def test_network_matches_the_static_solver_and_scipy() -> None:
    # Arrange
    edges = _random_edges(40, 150, 7)
    network = Network(edges, [3, 9])
    weights = np.random.default_rng(4).integers(1, 9, size=len(edges)).astype(np.float32)
    reweighted = [(a, b, float(w)) for (a, b, _), w in zip(edges, weights, strict=True)]
    # Act
    own = network.solve()
    single = network.solve(9)
    custom = network.solve(edge_weights=weights)
    reused = network.solve(edge_weights=weights, reuse=True)
    nearest = network.solve_nearest()
    nearest_custom = network.solve_nearest(edge_weights=weights, reuse=True)
    # Assert
    assert network.node_count == 40 and network.destination_indices.tolist() == [3, 9]
    np.testing.assert_array_equal(own, graphs.dijkstra_solve(edges, [3, 9], 40))
    np.testing.assert_array_equal(single, graphs.dijkstra_solve(edges, 9, 40))
    np.testing.assert_array_equal(custom, graphs.dijkstra_solve(reweighted, [3, 9], 40))
    np.testing.assert_array_equal(reused, custom)
    np.testing.assert_array_equal(nearest[:, 2], own[:, 2])
    np.testing.assert_array_equal(nearest_custom[:, 2], custom[:, 2])
    np.testing.assert_array_equal(own[:, 2], _scipy_costs(edges, 40)[:, [3, 9]].min(axis=1).astype(np.float32))
    assert sum(len(n) for n in network.incoming_edges) == len(edges)
    assert all(e.from_index == node for node, out in enumerate(network.outgoing_edges) for e in out)
    with pytest.raises(ValueError, match="omit destinations"):
        network.solve(3, edge_weights=weights)


def test_network_detours_match_scipy_with_the_edge_removed() -> None:
    # Arrange
    network = Network(TEXTBOOK, [4])
    table = network.solve()
    removed = (7,)
    reference = _scipy_costs(TEXTBOOK, 5, removed)[0, 4]
    # Act
    detour = network.get_path(list(removed), 0)
    reused = network.get_path(list(removed), 0, table)
    weighted = network.get_path(list(removed), 0, table, [w for _, _, w in TEXTBOOK])
    # Assert
    assert sum(TEXTBOOK[i][2] for i in detour) == reference
    np.testing.assert_array_equal(reused, detour)
    np.testing.assert_array_equal(weighted, detour)
    assert network.get_path([], 0, table).tolist() == [1, 7]
    assert network.get_path([7, 6, 4, 0], 0) is None


def test_binary_heap_orders_like_heapq_and_supports_keyed_updates() -> None:
    # Arrange
    rng = np.random.default_rng(3)
    weights = rng.integers(0, 1000, size=50).astype(np.float32)
    heap = BinaryHeap(50)
    reference = []
    for i, w in enumerate(weights):
        heap.add(HeapNode(float(w), i, f"item{i}"))
        heapq.heappush(reference, (float(w), i))
    # Act
    heap.decrease_key(HeapNode(-1.0, 7, "lowered"))
    heap.decrease_key(HeapNode(5000.0, 8, "ignored"))
    heap.replace(HeapNode(9999.0, 9, "raised"))
    reference = [(-1.0 if i == 7 else 9999.0 if i == 9 else w, i) for w, i in reference]
    heapq.heapify(reference)
    first = heap.remove_min()
    # Assert
    assert len(heap) == heap.count == 49
    assert (first.weight, first.index, first.value) == (-1.0, 7, "lowered")
    heapq.heappop(reference)
    popped = [heap.remove_min() for _ in range(49)]
    assert [p.weight for p in popped] == [w for w, _ in sorted(reference)]
    assert popped[-1].value == "raised"
    with pytest.raises(Exception, match="empty"):
        heap.remove_min()


def test_shortest_distances_matches_scipy_dijkstra_min_only() -> None:
    # Arrange: a random sparse graph; SciPy 1.18.1 is the oracle.
    rng = np.random.default_rng(20261108)
    n, m = 2000, 10000
    graph = sparse.csr_matrix((rng.uniform(0.0, 1.0, m) + 1e-3, (rng.integers(0, n, m), rng.integers(0, n, m))),
                              shape=(n, n))
    graph.sum_duplicates()
    sources = rng.choice(n, 5, replace=False)
    expected, expected_pred, expected_src = csgraph.dijkstra(graph, indices=sources, min_only=True,
                                                             return_predecessors=True)
    # Act
    distance, pred, src = graphs.shortest_distances(graph, sources, return_predecessors=True,
                                                    return_nearest_source=True)
    plain = graphs.shortest_distances(graph, sources)
    # Assert: distances agree to 1e-14 relative (bitwise here); ties may pick another predecessor.
    finite = np.isfinite(expected)
    assert np.array_equal(np.isinf(distance), np.isinf(expected))
    np.testing.assert_allclose(distance[finite], expected[finite], rtol=1e-14, atol=0.0)
    assert np.array_equal(plain, distance)
    reached = pred >= 0
    np.testing.assert_array_equal(distance[reached], distance[pred[reached]]
                                  + np.asarray(graph[pred[reached], np.nonzero(reached)[0]]).ravel())
    assert np.array_equal(np.where(expected_src == -9999, -1, expected_src), src)
    assert np.array_equal(np.where(expected_pred == -9999, -1, expected_pred)[~reached], pred[~reached])


def test_shortest_distances_limit_tuple_input_and_errors() -> None:
    # Arrange
    row_start, to_node, weight = [0, 2, 4, 5, 5], [1, 2, 2, 3, 3], [0.1, 0.3, 0.2, 1.2, 0.7]
    matrix = sparse.csr_matrix((weight, to_node, row_start), shape=(4, 4))
    # Act
    limited = graphs.shortest_distances((row_start, to_node, weight), 0, 0.35)
    # Assert
    np.testing.assert_array_equal(limited, csgraph.dijkstra(matrix, indices=[0], min_only=True, limit=0.35))
    assert "shortest_distances" in graphs.__all__
    with pytest.raises(ValueError, match="square"):
        graphs.shortest_distances(sparse.csr_matrix((2, 3)), [0])
    with pytest.raises(ValueError):
        graphs.shortest_distances((row_start, to_node, [0.1, 0.3, -0.2, 1.2, 0.7]), [0])
    with pytest.raises(ValueError):
        graphs.shortest_distances((row_start, to_node, weight), [4])


def test_routing_table_costs_are_float32_of_scipy_double_distances() -> None:
    # Arrange: decimal weights whose float32 sums round differently from their float64 sums.
    rng = np.random.default_rng(7)
    n, m = 300, 1500
    frm, to = rng.integers(0, n, m), rng.integers(0, n, m)
    w = np.round(rng.uniform(0.01, 10.0, m), 2).astype(np.float32)
    edges = [(int(a), int(b), float(c)) for a, b, c in zip(frm, to, w, strict=True)]
    best: dict[tuple[int, int], float] = {}
    for a, b, c in edges:
        best[(b, a)] = min(float(c), best.get((b, a), np.inf))
    rows, cols = zip(*best, strict=True)
    reverse = sparse.csr_matrix((list(best.values()), (rows, cols)), shape=(n, n))
    destinations = [3, 17]
    # Act
    table = graphs.dijkstra_solve_nearest(edges, destinations, n)
    # Assert
    expected = csgraph.dijkstra(reverse, indices=destinations, min_only=True).astype(np.float32)
    assert np.array_equal(table[:, 2].view(np.uint32), expected.view(np.uint32))


def test_shortest_distances_reads_csc_and_other_formats_as_their_matrix() -> None:
    # Arrange: an asymmetric graph, so reading CSC arrays as CSR would give the transpose.
    csr = sparse.csr_matrix(([0.1, 0.3, 0.2, 1.2, 0.7], [1, 2, 2, 3, 3], [0, 2, 4, 5, 5]), shape=(4, 4))
    expected = graphs.shortest_distances(csr, [0])
    # Act
    found = [graphs.shortest_distances(csr.tocsc(), [0]), graphs.shortest_distances(csr.tocoo(), [0]),
             graphs.shortest_distances(sparse.csc_array(csr), [0])]
    # Assert
    np.testing.assert_array_equal(expected, [0.0, 0.1, 0.3, 1.0])
    for distances in found:
        np.testing.assert_array_equal(distances, expected)


def test_shortest_distances_rejects_non_integral_sources() -> None:
    # Arrange
    graph = ([0, 2, 4, 5, 5], [1, 2, 2, 3, 3], [0.1, 0.3, 0.2, 1.2, 0.7])
    # Act and assert: integral floats are node indices; anything else is refused.
    np.testing.assert_array_equal(graphs.shortest_distances(graph, [0.0]), graphs.shortest_distances(graph, [0]))
    with pytest.raises(ValueError, match="integral"):
        graphs.shortest_distances(graph, [1.7])
    with pytest.raises(ValueError, match="integral"):
        graphs.shortest_distances(graph, [np.nan])
    with pytest.raises(TypeError):
        graphs.shortest_distances(graph, ["a"])


def test_negative_weights_raise_value_error() -> None:
    cycle = [(0, 1, 1.0), (1, 0, -2.0), (1, 2, 1.0)]
    with pytest.raises(ValueError, match="must be non-negative"):
        graphs.dijkstra_solve(cycle, 2)
    with pytest.raises(ValueError, match="must be non-negative"):
        graphs.dijkstra_solve_nearest(cycle, [2])
    network = Network(cycle, [2])
    with pytest.raises(ValueError, match="must be non-negative"):
        network.solve_nearest()
    with pytest.raises(ValueError, match="must be non-negative"):
        network.solve(edge_weights=[1.0, -1.0, 1.0])
    np.testing.assert_array_equal(network.solve_nearest(edge_weights=[1.0, 2.0, 1.0])[:, 2], [2.0, 1.0, 0.0])
    with pytest.raises(ValueError, match="must be non-negative"):
        graphs.shortest_distances((np.array([0, 1, 1]), np.array([1]), np.array([-1.0])), [0])
