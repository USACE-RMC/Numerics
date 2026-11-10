# Data, interpolation, time series and functions

This page covers three modules:

- `pynumerics.data`: paired data, interpolation, searching, line simplification, probability
  ordinates and linear regression (`Numerics.Data`).
- `pynumerics.timeseries`: dated series with arithmetic, smoothing, block series, resampling and
  summary statistics (`Numerics.Data.TimeSeries`).
- `pynumerics.functions`: link functions and univariate transfer functions such as rating curves
  (`Numerics.Functions`).

Every result comes from Numerics through pythonnet. Arrays cross into .NET in one bulk copy
through the array overloads, which give the scalar members' bits: interpolation and extrapolation
(`Bilinear` point pairs included), link and transfer functions, and the dates and values of a time
series. A wrapper's `.net` attribute holds the
wrapped .NET object.

## Interpolation

`Linear`, `CubicSpline` and `Polynomial` share the `Interpolater` members: `interpolate`
(also callable as the object itself), `search`, `x_values`, `y_values` and the search cursor.
Build an interpolator once and reuse it.

```python
import numpy as np
from pynumerics import data

x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
y = np.array([0.0, 0.8, 0.9, 0.1, -0.7])

line = data.Linear(x, y)
print(line([0.5, 2.5]))               # linear, end values held outside the knots
print(line.extrapolate(5.0))          # extends the last segment

spline = data.CubicSpline(x, y)       # natural cubic spline
print(spline(np.linspace(0.0, 4.0, 5)))

poly = data.Polynomial(2, x, y)       # Neville on 3 neighboring knots
print(poly(1.5), poly.error)
```

`Linear` interpolates in transformed spaces. A frequency curve is linear in the Normal quantile
of the probability and the log10 of the value:

```python
from pynumerics import data

p = [0.01, 0.1, 0.5, 0.9, 0.99]
q = [1.0, 10.0, 50.0, 200.0, 1000.0]
curve = data.Linear(p, q, x_transform="normal_z", y_transform="logarithmic")
print(curve([0.05, 0.95]))
```

`Bilinear` interpolates on a rectangular grid, with `y[i, j]` the value at `(x1[i], x2[j])`:

```python
import numpy as np
from pynumerics import data

grid = data.Bilinear([0.0, 1.0, 2.0], [10.0, 20.0], [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
print(grid(0.5, 15.0), grid(np.array([0.5, 1.5]), 20.0))  # point pairs broadcast, then cross once
```

`search` returns the lower bracketing index of a value in ordered data:

```python
from pynumerics import data

print(data.search(2.5, [0.0, 1.0, 2.0, 3.0], method="hunt"))   # 2
```

## Paired data

`OrderedPairedData` holds monotone (x, y) pairs and looks values up in either direction:

```python
from pynumerics import data

rating = data.OrderedPairedData([1.0, 2.0, 4.0, 8.0], [10.0, 40.0, 160.0, 640.0])
print(rating.get_y_from_x(3.0, x_transform="logarithmic", y_transform="logarithmic"))
print(rating.get_x_from_y([20.0, 300.0]))
print(rating.is_valid, rating.trapezoidal_area_under_y())

rating.append(16.0, 2560.0)
print(len(rating), rating[4])
simplified = rating.douglas_peucker_simplify(1.0)
restored = data.OrderedPairedData.from_xml(rating.to_xml())
```

`UncertainOrderedPairedData` holds a distribution of y at each x. `curve_sample` returns the curve
at one quantile of every distribution:

```python
from pynumerics import data, distributions

table = data.UncertainOrderedPairedData(
    [1.0, 2.0, 4.0],
    [distributions.Normal(10.0, 1.0), distributions.Normal(20.0, 2.0), distributions.Normal(40.0, 3.0)])
print(table.curve_sample(0.9).y)      # the 90 percent quantile at each x
print(table.curve_sample().y)         # the means
```

The Ramer-Douglas-Peucker algorithm simplifies a polyline:

```python
from pynumerics import data

kx, ky = data.ramer_douglas_peucker([0, 1, 2, 3, 4], [0.0, 0.1, -0.1, 2.0, 2.1], 0.5)
print(kx, ky)
```

`ProbabilityOrdinates` is the list of exceedance probabilities that frequency curves report:

```python
from pynumerics import data

ordinates = data.ProbabilityOrdinates()          # 25 standard ordinates, 1e-6 to 0.99
print(len(ordinates), ordinates.validate())
custom = data.ProbabilityOrdinates.parse("0.5|0.1|0.01")
print(custom.values)
```

## Time series

A `TimeSeries` has NumPy `datetime64[us]` dates and float values. NaN marks a missing value. A
regular series takes a start date and an interval. A series with explicit dates takes `dates`;
its dates and values cross in one call, and every value keeps its bits.

A regular series steps from the previous date, so month, quarter and year steps from a day after
the 28th settle on the shortest day met: a monthly series from 31 January 2020 runs 31 January,
29 February, 29 March. `shift_all_dates`, `fill_missing_dates` and the resampling methods step the
same way, and a series from a start to an end date stops at the last date on or before the end.
`TimeSeries.anchored_dates` counts each step from the start instead, as pandas
`DateOffset(months=k)` applied to the start; pass its dates as `dates` to build a series on them.

```python
import numpy as np
from pynumerics.timeseries import TimeSeries

rng = np.random.default_rng(1)
daily = TimeSeries(rng.gamma(2.0, 3.0, 3 * 365), interval="one_day", start="2019-01-01")
print(daily.start_date, daily.end_date, len(daily))
print(daily.dates[:3], daily.values[:3])

peaks = TimeSeries([1200.0, 1350.0, 1280.0], interval="irregular",
                   dates=["2015-04-02", "2016-05-19", "2017-03-30"])
print(peaks.values)
```

The arithmetic and missing-data members change the series in place:

```python
import numpy as np
from pynumerics.timeseries import TimeSeries

series = TimeSeries([10.0, np.nan, 30.0, 40.0], start="2024-01-01")
series.interpolate_missing_data(2)
series.multiply(2.0)
series.add(1.0, [0, 3])
print(series.values)
print(series.index_of_date("2024-01-03"), series.index_of_date("2025-01-01"))  # 2 -1
```

`index_of_date` finds the first ordinate on a date (or -1) from one bulk copy of the dates, so it
works on unsorted and irregular series; `index_of` matches the date and the value together.

Block series reduce a series to one value per year, water year, season, month or quarter. Each
maximum keeps the date it occurred:

```python
import numpy as np
from pynumerics.timeseries import TimeSeries

rng = np.random.default_rng(1)
daily = TimeSeries(rng.gamma(2.0, 3.0, 3 * 365), start="2019-01-01")
annual = daily.calendar_year_series("maximum")
water_year = daily.custom_year_series(10, block_function="maximum")
monthly_mean = daily.monthly_series("average")
smoothed = daily.calendar_year_series("maximum", "moving_average", 3)   # 3-day mean maxima
pot = daily.peaks_over_threshold(20.0, min_steps_between_events=2)
print(annual.dates, annual.values)
```

Smoothing, interval changes and resampling return new series:

```python
import numpy as np
from pynumerics.timeseries import TimeSeries

rng = np.random.default_rng(1)
daily = TimeSeries(rng.gamma(2.0, 3.0, 120), start="2019-01-01")
week = daily.moving_average(7)
total = daily.moving_sum(7, min_valid_count=5)
monthly = daily.aggregate_by_interval("one_month", "sum")
hourly = daily.convert_time_interval("one_hour")
synthetic = daily.resample_with_block_bootstrap(365, 30, seed=42)
trend, seasonal, residual = daily.seasonal_decompose(7)
print(daily.summary_statistics()["Mean"], daily.percentiles([0.1, 0.5, 0.9]))
```

`difference` dates each difference at its earlier value, using the series' own dates (gaps and
irregular series included); `label_at_later_date=True` dates it at its later value, as pandas
`Series.diff` does. `convert_time_interval` to a coarser step ends with the last block that holds
data (72 hourly values give 3 daily values, as pandas `resample("D")` does):

```python
import numpy as np
from pynumerics.timeseries import TimeSeries

hourly = TimeSeries(np.arange(72.0), interval="one_hour", start="2020-01-01")
print(hourly.difference().dates[0])                              # 2020-01-01T00:00
print(hourly.difference(label_at_later_date=True).dates[0])      # 2020-01-01T01:00
print(len(hourly.convert_time_interval("one_day")))               # 3
monthly = TimeSeries([1.0, 2.0, 3.0], interval="one_month", start="2020-01-31")
print(monthly.dates)                                              # 31 Jan, 29 Feb, 29 Mar
anchored = TimeSeries.anchored_dates("2020-01-31", 3, "one_month")
print(TimeSeries([1.0, 2.0, 3.0], interval="one_month", dates=anchored).dates)  # 31 Jan, 29 Feb, 31 Mar
```

`to_pandas`, `from_pandas`, `to_xml` and `from_xml` convert a series. `TimeSeries.from_net` wraps
a .NET series without copying it, for example one returned by another Numerics member.

```python
import pandas as pd
from pynumerics.timeseries import TimeSeries

frame = pd.Series([1.0, 2.0, 3.0], index=pd.date_range("2020-01-01", periods=3, freq="D"))
series = TimeSeries.from_pandas(frame, interval="one_day")
print(series.to_pandas())
```

## Link functions

A link maps a parameter to an unbounded space and back. The links are `IdentityLink`, `LogLink`,
`LogitLink`, `ProbitLink`, `ComplementaryLogLogLink`, `FisherZLink` and `YeoJohnsonLink`.
`LinkController` holds one link per parameter.

```python
from pynumerics import functions

logit = functions.LogitLink()
print(logit.link([0.1, 0.5]), logit.inverse_link(0.0), logit.d_link(0.5))

controller = functions.LinkController(["identity", functions.LogLink(), None])
eta = controller.link([1.0, 2.0, 3.0])
print(eta, controller.inverse_link(eta), controller.log_det_jacobian(eta))
```

## Univariate functions

`LinearFunction`, `PowerFunction`, `SegmentedPowerFunction`, `TabularFunction` and
`CompositeFunction` evaluate `y = f(x)` and its inverse. A stochastic function evaluates its
mean curve until `confidence_level` is set in [0, 1].

```python
from pynumerics import data, distributions, functions

rating = functions.PowerFunction(alpha=20.0, beta=1.6, xi=0.5, sigma=0.1)
print(rating([1.0, 2.0, 5.0]))
rating.confidence_level = 0.95
print(rating(2.0), rating.inverse_function(rating(2.0)))

controls = functions.SegmentedPowerFunction(parameters=[1.0, 1.2, 1.6, 3.0, 0.4, 2.1, 0.05])
print(controls(4.0), controls.get_breakpoint(2))

table = functions.TabularFunction(data.UncertainOrderedPairedData(
    [1.0, 2.0, 4.0], [distributions.Normal(10.0, 1.0), distributions.Normal(20.0, 2.0),
                      distributions.Normal(40.0, 3.0)]))
print(table(3.0))

blend = functions.CompositeFunction([functions.LinearFunction(0.0, 1.0), functions.LinearFunction(1.0, 2.0)],
                                    [0.25, 0.75])
print(blend(2.0))

same = functions.function_from_xml(rating.to_xml())
```

`EnsembleFunction` holds parameter sets for one template, such as a posterior sample of a
rating curve:

```python
import numpy as np
from pynumerics import functions

sets = np.array([[20.0, 1.5, 0.5, 0.1], [22.0, 1.6, 0.4, 0.1], [21.0, 1.55, 0.45, 0.1]])
ensemble = functions.EnsembleFunction(functions.PowerFunction(20.0, 1.5, 0.5, 0.1), sets)
print(len(ensemble), ensemble.sample(0.5)(3.0), ensemble.sample_at(0).alpha)
```
