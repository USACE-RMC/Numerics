# Optimization, shortest paths and dynamic programming

`pynumerics.optimization` wraps the local, global, constrained and multi-objective optimizers of
`Numerics.Mathematics.Optimization`, the exact 0/1 knapsack and discrete dynamic programming.
`pynumerics.graphs` wraps the shortest-path solvers (Dijkstra and `Network`) and the keyed binary heap.

Every optimizer runs the Numerics search in .NET. Python only supplies the objective.

## The optimizer interface

Every optimizer class takes the objective first, then its start point and bounds. Each Numerics
setting is a snake-case property, and you can also pass it to the constructor as a keyword.
`minimize()` and `maximize()` return an `OptimizationResult`:

| Field | Meaning |
| --- | --- |
| `parameters` | best parameters (`BestParameterSet.Values`) |
| `value` | objective at `parameters`, with the sign restored after `maximize()` |
| `iterations`, `function_evaluations` | counts used by the run |
| `status` | an `OptimizationStatus` (`SUCCESS`, `MAXIMUM_ITERATIONS_REACHED`, ...); `success` is a shortcut |
| `hessian` | Hessian at the optimum when `compute_hessian` is on (the default), else None |
| `trace_parameters`, `trace_values` | the best-so-far point and value after every evaluation (copied on first use) |

Settings you leave unset keep the Numerics defaults.

```python
import numpy as np
from pynumerics.optimization import BFGS, NelderMead, Powell, OptimizationStatus

def rosenbrock(x):
    return 100.0 * (x[1] - x[0] ** 2) ** 2 + (1.0 - x[0]) ** 2

lower, upper = [-2.048, -2.048], [2.048, 2.048]
for method in (NelderMead, Powell, BFGS):
    result = method(rosenbrock, [0.0, 0.0], lower, upper).minimize()
    assert result.status is OptimizationStatus.SUCCESS
    print(method.__name__, np.round(result.parameters, 6), result.function_evaluations)

nm = NelderMead(rosenbrock, [0.0, 0.0], lower, upper, max_iterations=2000, simplex_scale=0.1)
result = nm.minimize()
print(result.trace_values[:5], result.trace_parameters.shape, result.hessian)
```

### Local methods

| Class | Method | Constructor |
| --- | --- | --- |
| `NelderMead` | downhill simplex | `(objective, initial_values, lower_bounds, upper_bounds)` |
| `Powell` | conjugate directions | same |
| `BFGS`, `LBFGSB` | quasi-Newton | same, plus `gradient=` |
| `GradientDescent`, `ADAM` | first order | same, plus `alpha=` and `gradient=` |
| `GoldenSection`, `BrentSearch` | one-dimensional | `(objective, lower_bound, upper_bound)` |

```python
import numpy as np
from pynumerics.optimization import LBFGSB, BrentSearch, GoldenSection

def rosenbrock(x):
    return 100.0 * (x[1] - x[0] ** 2) ** 2 + (1.0 - x[0]) ** 2

def gradient(x):
    return np.array([-400.0 * x[0] * (x[1] - x[0] ** 2) - 2.0 * (1.0 - x[0]), 200.0 * (x[1] - x[0] ** 2)])

result = LBFGSB(rosenbrock, [-1.2, 1.0], [-5.0, -5.0], [5.0, 5.0], gradient=gradient, memory_size=5).minimize()
print(result.parameters)

print(GoldenSection(lambda x: (x - 2.0) ** 2 + 1.0, 0.0, 5.0).minimize().parameters)
search = BrentSearch(lambda x: (x - 10.0) ** 2, 0.0, 1.0)
print(search.bracket())          # moves the interval until it brackets the minimum
print(search.minimize().value)
```

### Global methods

`DifferentialEvolution`, `ParticleSwarm`, `ShuffledComplexEvolution` and `SimulatedAnnealing` take
`(objective, lower_bounds, upper_bounds)`. `CMAES` and `DynamicallyDimensionedSearch` also accept an
optional `initial_values`. `MultiStart` and `MLSL` take a start point and a local `method`
(`LocalMethod.BFGS` by default). The searches are seeded (`prng_seed`), so repeated runs give the
same result.

```python
import numpy as np
from pynumerics.optimization import CMAES, MLSL, DifferentialEvolution, LocalMethod, ParticleSwarm

def ackley(x):
    return (-20.0 * np.exp(-0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2)))
            - np.exp(0.5 * (np.cos(2 * np.pi * x[0]) + np.cos(2 * np.pi * x[1]))) + np.e + 20.0)

lower, upper = [-5.0, -5.0], [5.0, 5.0]
print(DifferentialEvolution(ackley, lower, upper, population_size=40).minimize().parameters)
print(ParticleSwarm(ackley, lower, upper, population_size=50).minimize().parameters)
print(CMAES(ackley, lower, upper, initial_values=[1.0, 1.0]).minimize().parameters)

mlsl = MLSL(ackley, [3.0, 3.0], lower, upper, method=LocalMethod.BFGS, sample_size=50)
print(mlsl.minimize().parameters, len(mlsl.local_minimums))
```

## Callbacks, C# objectives and parallel runs

An objective, gradient or constraint can be any of these:

- a Python callable of a float64 array;
- a string of C# source, compiled once with Roslyn (`pynumerics.callbacks.compile_csharp`);
- a .NET delegate.

A Python objective costs about 2.5 microseconds per evaluation to cross from .NET into Python. A
compiled objective runs entirely in .NET.

`DifferentialEvolution` and `ParticleSwarm` can evaluate each generation on several threads with
`use_parallel=True`. This setting stays off for a Python objective, because Python callbacks
serialize on the GIL. Parallel evaluation runs the deferred variant of the search, so its result
differs from the serial run, but it does not depend on the thread count.

```python
from pynumerics import callbacks
from pynumerics.optimization import DifferentialEvolution

if callbacks.roslyn_available():
    source = "x => (x[0] - 1) * (x[0] - 1) + 10 * (x[1] + 2) * (x[1] + 2)"
    de = DifferentialEvolution(source, [-5.0, -5.0], [5.0, 5.0], use_parallel=True)
    print(de.use_parallel, de.minimize().parameters)

python_de = DifferentialEvolution(lambda x: float(x[0] ** 2 + x[1] ** 2), [-5.0, -5.0], [5.0, 5.0],
                                  use_parallel=True)
print(python_de.use_parallel)   # False: a Python objective runs serially
```

## Progress and cancellation

```python
import pynumerics
from pynumerics.optimization import DifferentialEvolution, OptimizationStatus

pynumerics.initialize()   # .NET namespaces such as System import once the runtime is running
import System

source = System.Threading.CancellationTokenSource()
reports = []

def on_progress(report):
    reports.append((report.function_evaluations, report.best_fitness))
    if len(reports) == 2:
        source.Cancel()

de = DifferentialEvolution(lambda x: float(x[0] ** 2 + x[1] ** 2), [-5.0, -5.0], [5.0, 5.0], progress_interval=100)
de.progress = on_progress
de.cancellation_token = source.Token
result = de.minimize()
print(reports, result.status is OptimizationStatus.CANCELLED)
```

## Constrained optimization

`AugmentedLagrange` wraps an inner optimizer and a list of `Constraint` objects. Each constraint has
the form `g(x) relation value`, where the relation is a `ConstraintType`.

```python
from pynumerics.optimization import BFGS, AugmentedLagrange, Constraint, ConstraintType

def objective(x):
    return (x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2

disk = Constraint(lambda x: x[0] ** 2 + x[1] ** 2, 2, 2.0, ConstraintType.LESSER_THAN_OR_EQUAL_TO)
inner = BFGS(objective, [0.5, 0.5], [-5.0, -5.0], [5.0, 5.0])
solver = AugmentedLagrange(objective, inner, [disk])
result = solver.minimize()
print(result.parameters, result.value, solver.mu)
```

## Multi-objective utilities

`ParetoFront` minimizes every objective and takes one row of objective values per point.

```python
import numpy as np
from pynumerics.optimization import ParetoFront

points = np.array([[1.0, 4.0], [2.0, 2.0], [4.0, 1.0], [3.0, 3.0], [5.0, 5.0]])
fronts = ParetoFront.non_dominated_sort(points)       # [array([0, 1, 2]), array([3]), array([4])]
print(fronts, ParetoFront.non_dominated_indices(points))
print(ParetoFront.crowding_distance(points[fronts[0]]))
print(ParetoFront.hypervolume(points, [6.0, 6.0]))    # 20.0
reference = np.array([[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]])
print(ParetoFront.generational_distance(points, reference),
      ParetoFront.inverted_generational_distance(points, reference))
```

## Shortest paths

Edges are `Edge(from_index, to_index, weight, index)` values or `(from, to, weight[, index])`
tuples. Solvers return a float32 routing table with one row per node:

- column 0 is the next node;
- column 1 is the edge to take;
- column 2 is the cost to the destination.

```python
from pynumerics import graphs

edges = [(0, 1, 10.0), (0, 3, 5.0), (1, 2, 1.0), (1, 3, 2.0), (2, 4, 4.0),
         (3, 1, 3.0), (3, 2, 9.0), (3, 4, 2.0), (4, 0, 7.0), (4, 2, 6.0)]
table = graphs.dijkstra_solve(edges, 4)
print(table[:, 2])                                # [7. 4. 4. 2. 0.]
print(graphs.dijkstra_get_path(table, 0))         # edge indices [1 7]
print(graphs.dijkstra_try_get_path(table, 0))

network = graphs.Network(edges, [4])
print(network.solve_nearest()[:, 2])
print(network.solve(edge_weights=[1.0] * len(edges), reuse=True)[:, 2])
print(network.get_path([7], 0))                   # detour that avoids edge 7
```

The table costs accumulate in float64 and are rounded to float32 once, so they equal
`scipy.sparse.csgraph.dijkstra` on the reversed graph, rounded to float32.

`graphs.shortest_distances` returns the float64 distance of every node from the nearest of several
sources on a CSR graph, as `scipy.sparse.csgraph.dijkstra(graph, indices=sources, min_only=True)`
does. It takes a SciPy CSR matrix or a `(row_start, to_node, weight)` tuple, an optional `limit`,
and returns predecessors and nearest sources on request (-1 where SciPy reports -9999). With ten sources on a
random graph of 100,000 nodes and 500,000 edges it takes about 18 ms against 25 ms for SciPy
1.18.1, and it stays ahead at 10,000 and 1,000,000 nodes.

```python
import numpy as np
from pynumerics import graphs

row_start = np.array([0, 2, 4, 5, 5])                  # edges 0->1, 0->2, 1->2, 1->3, 2->3
to_node = np.array([1, 2, 2, 3, 3])
weight = np.array([0.1, 0.3, 0.2, 1.2, 0.7])
distance, predecessor = graphs.shortest_distances((row_start, to_node, weight), [0], return_predecessors=True)
print(distance)                                   # [0.  0.1 0.3 1. ]
print(predecessor)                                # [-1  0  0  2]
```

## Knapsack

```python
from pynumerics import optimization

total, selected = optimization.knapsack([60.0, 100.0, 120.0], [10, 20, 30], 50)
print(total, selected)                            # 220.0 [1 2]
total, selected = optimization.knapsack([60.0, 100.0, 120.0], [[10, 1], [20, 1], [30, 1]], [50, 1])
print(total, selected)                            # one budget per constraint
```

## Dynamic programming

`DynamicProgramming(stages, states, decisions, stage_value, transition)` solves a deterministic
problem. Pass `max_outcomes` and a transition that fills the outcome arrays for a stochastic one.

```python
from pynumerics.optimization import DynamicProgramming

def cost(t, s, a):                                # order a units, hold s units
    return (2.0 * a + 1.0 if a > 0 else 0.0) + 0.5 * s

def step(t, s, a):                                # demand of one unit per stage; -1 is infeasible
    nxt = s + a - 1
    return nxt if 0 <= nxt < 4 else -1

dp = DynamicProgramming(4, 4, 3, cost, step)
dp.solve()
print(dp.value[0], dp.policy)
print(dp.simulate(1))

def transition(t, s, a, next_states, probabilities):
    next_states[0], probabilities[0] = s, 0.5
    next_states[1], probabilities[1] = min(s + a, 3), 0.5
    return 2

mdp = DynamicProgramming(10, 4, 2, lambda t, s, a: float(s - a), transition, max_outcomes=2)
mdp.maximize = True
mdp.discount_factor = 0.9
print(mdp.solve_steady_state(), mdp.steady_state_value, mdp.steady_state_policy)
mdp.solve()
print(mdp.simulate(0, seed=1).states)
```

## Binary heap

```python
from pynumerics.graphs import BinaryHeap, HeapNode

heap = BinaryHeap(10)
heap.add(HeapNode(3.0, 0, "c"))
heap.add(HeapNode(1.0, 1, "a"))
heap.decrease_key(HeapNode(0.5, 0, "c, lowered"))
print(len(heap), heap.remove_min(), heap.remove_min())
```
