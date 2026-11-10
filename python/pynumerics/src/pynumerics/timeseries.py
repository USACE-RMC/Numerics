"""Time series from ``Numerics.Data``.

:class:`TimeSeries` wraps ``Numerics.Data.TimeSeries``: arithmetic and transforms in place,
missing-data handling, smoothing, block series (annual, water-year, monthly and quarterly maxima
and other block functions), peaks over threshold, interval conversion and aggregation,
bootstrap and k-NN resampling, seasonal decomposition and summary statistics.

Dates are NumPy ``datetime64[us]`` values on the Python side and ``System.DateTime`` (kind
unspecified) on the .NET side; :attr:`TimeSeries.ticks` gives the exact .NET ticks (100 ns).
The dates and values of a series cross in two bulk copies (``IndexTicksToArray`` and
``ValuesToArray``), and a regular series is built from one array (the ``(TimeInterval, DateTime,
IList<double>)`` constructor). A series with explicit dates is built from its ticks and values in
one call (the ``(TimeInterval, IList<long>, IList<double>)`` constructor), so every value keeps its
bits, NaN payloads included.

Regular series step from the previous date (``AddTimeInterval``), so month, quarter and year steps
from a day after the 28th settle on the shortest day met: a monthly series from 31 January 2020
runs 31 January, 29 February, 29 March, as pandas ``DateOffset(months=1)`` applied to the previous
date. :meth:`TimeSeries.anchored_dates` counts each step from the start instead (31 January,
29 February, 31 March, as pandas ``DateOffset(months=k)`` applied to the start).
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np

from pynumerics._util import _data, net_enum
from pynumerics.data import parse_xml
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    to_net_doubles,
    to_net_ints,
    to_net_longs,
    to_numpy_doubles,
    to_numpy_longs,
    to_numpy_matrix,
)

_EPOCH_TICKS = 621355968000000000
"""``DateTime(1970, 1, 1).Ticks``: the .NET tick count of the NumPy epoch."""


class TimeInterval(NetEnum):
    """Time steps of a series (``Numerics.Data.TimeInterval``)."""

    ONE_MINUTE = "OneMinute"
    FIVE_MINUTE = "FiveMinute"
    FIFTEEN_MINUTE = "FifteenMinute"
    THIRTY_MINUTE = "ThirtyMinute"
    ONE_HOUR = "OneHour"
    SIX_HOUR = "SixHour"
    TWELVE_HOUR = "TwelveHour"
    ONE_DAY = "OneDay"
    SEVEN_DAY = "SevenDay"
    ONE_MONTH = "OneMonth"
    ONE_QUARTER = "OneQuarter"
    ONE_YEAR = "OneYear"
    IRREGULAR = "Irregular"


class BlockFunction(NetEnum):
    """Reductions of a block of values (``Numerics.Data.BlockFunctionType``)."""

    MINIMUM = "Minimum"
    MAXIMUM = "Maximum"
    AVERAGE = "Average"
    SUM = "Sum"


class SmoothingFunction(NetEnum):
    """Smoothing applied before blocking (``Numerics.Data.SmoothingFunctionType``)."""

    DIFFERENCE = "Difference"
    MOVING_AVERAGE = "MovingAverage"
    MOVING_SUM = "MovingSum"
    NONE = "None"


class MathFunction(NetEnum):
    """Element-wise operations on a series (``Numerics.Data.MathFunctionType``), for user interfaces."""

    ADD = "Add"
    SUBTRACT = "Subtract"
    MULTIPLY = "Multiply"
    DIVIDE = "Divide"
    LOGARITHM = "Logarithm"
    EXPONENTIATE = "Exponentiate"
    INVERSE = "Inverse"
    REPLACE = "Replace"
    INTERPOLATE = "Interpolate"


class TimeBlockWindow(NetEnum):
    """Block windows of a series (``Numerics.Data.TimeBlockWindow``), for user interfaces."""

    CALENDAR_YEAR = "CalendarYear"
    WATER_YEAR = "WaterYear"
    CUSTOM_YEAR = "CustomYear"
    QUARTER = "Quarter"
    MONTH = "Month"


def to_ticks(dates: Any) -> np.ndarray:
    """Return dates as .NET ``DateTime.Ticks`` (100 ns since 0001-01-01).

    Parameters
    ----------
    dates : array_like
        ``datetime64`` values, ``datetime.datetime`` objects or ISO 8601 strings (a pandas
        ``DatetimeIndex`` included, without a time zone).

    Returns
    -------
    numpy.ndarray
        The int64 ticks; nanosecond input keeps its 100 ns ticks, coarser input is exact.
    """
    array = np.asarray(dates)
    if array.dtype.kind != "M":
        array = np.asarray(dates, dtype="datetime64[us]")
    if np.datetime_data(array.dtype)[0] in ("ns", "ps", "fs", "as"):
        return array.astype("datetime64[ns]").astype(np.int64) // 100 + _EPOCH_TICKS
    return array.astype("datetime64[us]").astype(np.int64) * 10 + _EPOCH_TICKS


def from_ticks(ticks: Any) -> np.ndarray:
    """Return .NET ``DateTime.Ticks`` as ``datetime64[us]`` dates.

    Parameters
    ----------
    ticks : array_like
        The int64 ticks (100 ns since 0001-01-01).

    Returns
    -------
    numpy.ndarray
        The dates; sub-microsecond ticks are floored to the microsecond.
    """
    return ((np.asarray(ticks, dtype=np.int64) - _EPOCH_TICKS) // 10).astype("datetime64[us]")


def _net_date(value: Any) -> Any:
    """Return one date as a ``System.DateTime`` of unspecified kind."""
    import System  # noqa: PLC0415 - the caller started the runtime

    return System.DateTime(int(to_ticks([value])[0]))


def _py_date(net: Any) -> np.datetime64:
    """Return a ``System.DateTime`` as a ``datetime64[us]`` scalar."""
    return from_ticks([int(net.Ticks)])[0]


class TimeSeries:
    """A series of dated values (``Numerics.Data.TimeSeries``).

    Parameters
    ----------
    values : array_like, optional
        The values; NaN marks a missing value.
    interval : TimeInterval or str, optional
        The time step; default one day.
    start : datetime-like, optional
        The first date of a regular series; the later dates step by ``interval``.
    end : datetime-like, optional
        With ``start`` and no ``values``: the last date of a regular series filled with ``fill``.
    dates : array_like, optional
        Explicit dates, one per value, instead of ``start``.
    fill : float, optional
        The value of a series built from ``start`` and ``end``; default NaN.

    Attributes
    ----------
    net : Numerics.Data.TimeSeries
        The wrapped .NET object.

    Notes
    -----
    With neither ``start`` nor ``dates`` the series is empty and :meth:`append` adds ordinates.
    """

    def __init__(self, values: Any = None, interval: TimeInterval | str = TimeInterval.ONE_DAY, start: Any = None,
                 end: Any = None, dates: Any = None, fill: float = float("nan")) -> None:
        data = _data()
        kind = net_enum(data.TimeInterval, interval)
        if dates is not None:
            self.net = self._explicit(kind, dates, values)
        elif start is not None and values is not None:
            self.net = data.TimeSeries(kind, _net_date(start), to_net_doubles(values))
        elif start is not None and end is not None:
            self.net = data.TimeSeries(kind, _net_date(start), _net_date(end), float(fill))
        elif start is None and values is None:
            self.net = data.TimeSeries(kind)
        else:
            raise ValueError("Give values with start or dates, start with end, or nothing for an empty series.")

    @staticmethod
    def _explicit(kind: Any, dates: Any, values: Any) -> Any:
        """Build a .NET series from explicit dates and values: two arrays in one constructor call."""
        ticks = np.ravel(to_ticks(dates))
        numbers = np.asarray(values, dtype=np.float64).ravel()
        if ticks.shape != numbers.shape:
            raise ValueError(f"Expected one value per date; got {numbers.size} values for {ticks.size} dates.")
        return _data().TimeSeries(kind, to_net_longs(ticks), to_net_doubles(numbers))

    @classmethod
    def from_net(cls, net: Any) -> TimeSeries:
        """Wrap an existing .NET time series (for example one from ``TimeSeriesDownload``).

        Parameters
        ----------
        net : Numerics.Data.TimeSeries
            The .NET series; it is shared, not copied.

        Returns
        -------
        TimeSeries
            The wrapper.
        """
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def from_xml(cls, xml: str) -> TimeSeries:
        """Restore a series from its XML (the ``XElement`` constructor).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        TimeSeries
            The series.
        """
        return cls.from_net(_data().TimeSeries(parse_xml(xml)))

    @classmethod
    def from_pandas(cls, series: Any, interval: TimeInterval | str = TimeInterval.IRREGULAR) -> TimeSeries:
        """Build a series from a pandas ``Series`` with a ``DatetimeIndex``.

        Parameters
        ----------
        series : pandas.Series
            The values, indexed by time zone naive dates.
        interval : TimeInterval or str, optional
            The time step to record; default irregular.

        Returns
        -------
        TimeSeries
            The series, with the same dates and values.
        """
        return cls(series.to_numpy(dtype=np.float64), interval=interval, dates=series.index.to_numpy())

    def to_pandas(self) -> Any:
        """Return the series as a pandas ``Series`` indexed by its dates.

        Returns
        -------
        pandas.Series
            The values, with a ``DatetimeIndex``.
        """
        import pandas as pd  # noqa: PLC0415 - pandas is optional

        return pd.Series(self.values, index=pd.DatetimeIndex(self.dates))

    def to_xml(self) -> str:
        """Return the series' ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())

    def copy(self) -> TimeSeries:
        """Return an independent copy (``Clone``).

        Returns
        -------
        TimeSeries
            The copy.
        """
        return TimeSeries.from_net(self.net.Clone())

    def __copy__(self) -> TimeSeries:
        return self.copy()

    def __deepcopy__(self, memo: dict) -> TimeSeries:
        return self.copy()

    # Collection --------------------------------------------------------------------------------

    def __len__(self) -> int:
        return int(self.net.Count)

    def __iter__(self) -> Iterator[tuple[np.datetime64, float]]:
        return iter(zip(list(self.dates), self.values.tolist(), strict=True))

    def __getitem__(self, index: int) -> tuple[np.datetime64, float]:
        ordinate = self.net[int(index)]
        return _py_date(ordinate.Index), float(ordinate.Value)

    def __setitem__(self, index: int, item: tuple[Any, float]) -> None:
        self.net[int(index)] = self._ordinate(*item)

    def _ordinate(self, date: Any, value: float) -> Any:
        """Return a ``SeriesOrdinate<DateTime, double>``."""
        import System  # noqa: PLC0415 - the runtime is running

        return _data().SeriesOrdinate[System.DateTime, System.Double](_net_date(date), float(value))

    @property
    def count(self) -> int:
        """The number of ordinates."""
        return int(self.net.Count)

    @property
    def values(self) -> np.ndarray:
        """The values in collection order, copied in one transfer (``ValuesToArray``)."""
        return to_numpy_doubles(self.net.ValuesToArray())

    @property
    def ticks(self) -> np.ndarray:
        """The dates as int64 .NET ticks, copied in one transfer (``IndexTicksToArray``)."""
        return to_numpy_longs(self.net.IndexTicksToArray())

    @property
    def dates(self) -> np.ndarray:
        """The dates in collection order as ``datetime64[us]``."""
        return from_ticks(self.ticks)

    def append(self, date: Any, value: float) -> None:
        """Append an ordinate (``Add(SeriesOrdinate)``).

        Parameters
        ----------
        date : datetime-like
            The date.
        value : float
            The value (NaN for missing).
        """
        ordinate = self._ordinate(date, value)
        # Bind Add(SeriesOrdinate) explicitly: the group also has Add(double) and IList.Add(object).
        self.net.Add.Overloads[type(ordinate)](ordinate)

    def insert(self, index: int, date: Any, value: float) -> None:
        """Insert an ordinate before a position (``Insert``).

        Parameters
        ----------
        index : int
            The position.
        date : datetime-like
            The date.
        value : float
            The value.
        """
        self.net.Insert(int(index), self._ordinate(date, value))

    def contains(self, date: Any, value: float) -> bool:
        """Return whether an ordinate with this date and value is in the series (``Contains``).

        Parameters
        ----------
        date : datetime-like
            The date.
        value : float
            The value.

        Returns
        -------
        bool
            Whether it is present.
        """
        return bool(self.net.Contains(self._ordinate(date, value)))

    def index_of(self, date: Any, value: float) -> int:
        """Return the position of an ordinate, or -1 (``IndexOf``).

        Parameters
        ----------
        date : datetime-like
            The date.
        value : float
            The value.

        Returns
        -------
        int
            The first position, or -1.
        """
        return int(self.net.IndexOf(self._ordinate(date, value)))

    def index_of_date(self, date: Any) -> int:
        """Return the position of the first ordinate on a date, or -1 when no ordinate has it.

        Unlike :meth:`index_of`, which matches the date and the value (``IndexOf``), this matches
        the date alone, to the .NET tick (100 ns). The dates cross in one transfer
        (``IndexTicksToArray``) and are compared in NumPy, so the series may be unsorted or
        irregular.

        Parameters
        ----------
        date : datetime-like
            The date, as for :func:`to_ticks`.

        Returns
        -------
        int
            The first position whose date equals ``date``, or -1.

        Examples
        --------
        >>> from pynumerics.timeseries import TimeSeries
        >>> series = TimeSeries([1.0, 2.0, 3.0], "one_day", "2020-01-01")  # doctest: +SKIP
        >>> series.index_of_date("2020-01-02"), series.index_of_date("2021-01-01")  # doctest: +SKIP
        (1, -1)
        """
        target = to_ticks([date])[0]
        matches = np.flatnonzero(self.ticks == target)
        return int(matches[0]) if matches.size else -1

    def remove(self, date: Any, value: float) -> bool:
        """Remove the first ordinate with this date and value (``Remove``).

        Parameters
        ----------
        date : datetime-like
            The date.
        value : float
            The value.

        Returns
        -------
        bool
            Whether one was removed.
        """
        return bool(self.net.Remove(self._ordinate(date, value)))

    def remove_at(self, index: int) -> None:
        """Remove the ordinate at a position (``RemoveAt``).

        Parameters
        ----------
        index : int
            The position.
        """
        self.net.RemoveAt(int(index))

    def clear(self) -> None:
        """Remove every ordinate (``Clear``)."""
        self.net.Clear()

    # Properties --------------------------------------------------------------------------------

    @property
    def interval(self) -> TimeInterval:
        """The time step."""
        return TimeInterval(str(self.net.TimeInterval))

    @property
    def start_date(self) -> np.datetime64:
        """The earliest date."""
        return _py_date(self.net.StartDate)

    @property
    def end_date(self) -> np.datetime64:
        """The latest date."""
        return _py_date(self.net.EndDate)

    @property
    def has_missing_values(self) -> bool:
        """Whether any value is NaN."""
        return bool(self.net.HasMissingValues)

    def number_of_missing_values(self) -> int:
        """Return the number of NaN values (``NumberOfMissingValues``).

        Returns
        -------
        int
            The count.
        """
        return int(self.net.NumberOfMissingValues())

    # In-place operations -----------------------------------------------------------------------

    def _apply(self, name: str, indexes: Any, *arguments: Any) -> None:
        """Call an in-place member on every value, or on the values at ``indexes``."""
        member = getattr(self.net, name)
        if indexes is None:
            member(*arguments)
        else:
            member(*arguments, to_net_ints(indexes))

    def add(self, constant: float, indexes: Any = None) -> None:
        """Add a constant to the values in place (``Add(double)``).

        Parameters
        ----------
        constant : float
            The addend.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Add", indexes, float(constant))

    def subtract(self, constant: float, indexes: Any = None) -> None:
        """Subtract a constant from the values in place (``Subtract``).

        Parameters
        ----------
        constant : float
            The subtrahend.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Subtract", indexes, float(constant))

    def multiply(self, constant: float, indexes: Any = None) -> None:
        """Multiply the values by a constant in place (``Multiply``).

        Parameters
        ----------
        constant : float
            The factor.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Multiply", indexes, float(constant))

    def divide(self, constant: float, indexes: Any = None) -> None:
        """Divide the values by a constant in place (``Divide``).

        Parameters
        ----------
        constant : float
            The divisor (not zero).
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Divide", indexes, float(constant))

    def absolute_value(self, indexes: Any = None) -> None:
        """Replace the values by their absolute values in place (``AbsoluteValue``).

        Parameters
        ----------
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("AbsoluteValue", indexes)

    def exponentiate(self, power: float, indexes: Any = None) -> None:
        """Raise the values to a power in place (``Exponentiate``).

        Parameters
        ----------
        power : float
            The exponent.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Exponentiate", indexes, float(power))

    def log_transform(self, base: float = 10.0, indexes: Any = None) -> None:
        """Replace the values by their logarithms in place (``LogTransform``).

        Parameters
        ----------
        base : float, optional
            The logarithm base; default 10.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        if indexes is None:
            self.net.LogTransform(float(base))
        else:
            self.net.LogTransform(to_net_ints(indexes), float(base))

    def inverse(self, indexes: Any = None) -> None:
        """Replace the values by their reciprocals in place (``Inverse``).

        Parameters
        ----------
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        self._apply("Inverse", indexes)

    def standardize(self) -> None:
        """Replace the values by their z-scores ``(x - mean) / sd`` in place (``Standardize``)."""
        self.net.Standardize()

    def replace_missing_data(self, value: float, indexes: Any = None) -> None:
        """Replace NaN values by a constant in place (``ReplaceMissingData``).

        Parameters
        ----------
        value : float
            The replacement.
        indexes : array_like of int, optional
            The positions to change; default all.
        """
        if indexes is None:
            self.net.ReplaceMissingData(float(value))
        else:
            self.net.ReplaceMissingData(to_net_ints(indexes), float(value))

    def interpolate_missing_data(self, max_number_of_missing: int, indexes: Any = None) -> None:
        """Fill gaps of NaN values by linear interpolation in place (``InterpolateMissingData``).

        Parameters
        ----------
        max_number_of_missing : int
            The longest gap filled; longer gaps stay missing.
        indexes : array_like of int, optional
            The positions to consider; default all.
        """
        if indexes is None:
            self.net.InterpolateMissingData(int(max_number_of_missing))
        else:
            self.net.InterpolateMissingData(int(max_number_of_missing), to_net_ints(indexes))

    def sort_by_time(self, descending: bool = False) -> None:
        """Sort the ordinates by date in place (``SortByTime``).

        Parameters
        ----------
        descending : bool, optional
            Whether the latest date comes first; default False.
        """
        self._sort("SortByTime", descending)

    def sort_by_value(self, descending: bool = False) -> None:
        """Sort the ordinates by value in place (``SortByValue``).

        Parameters
        ----------
        descending : bool, optional
            Whether the largest value comes first; default False.
        """
        self._sort("SortByValue", descending)

    def _sort(self, name: str, descending: bool) -> None:
        """Call a sort member with a ``ListSortDirection``."""
        import clr  # noqa: PLC0415 - the runtime is running
        import System  # noqa: PLC0415

        direction = clr.GetClrType(_data().TimeSeries).GetMethod(name).GetParameters()[0].ParameterType
        getattr(self.net, name)(System.Enum.ToObject(direction, 1 if descending else 0))

    def shift_all_dates(self, new_start_date: Any) -> None:
        """Move the first date and restep the others by the interval in place (``ShiftAllDates``).

        Parameters
        ----------
        new_start_date : datetime-like
            The new first date; an irregular series moves only its first date.
        """
        self.net.ShiftAllDates(_net_date(new_start_date))

    # Derived series ----------------------------------------------------------------------------

    def cumulative_sum(self) -> TimeSeries:
        """Return the running sum (``CumulativeSum``).

        Returns
        -------
        TimeSeries
            The cumulative series.
        """
        return TimeSeries.from_net(self.net.CumulativeSum())

    def difference(self, lag: int = 1, differences: int = 1, label_at_later_date: bool = False) -> TimeSeries:
        """Return the lagged differences ``x[t + lag] - x[t]``, applied ``differences`` times (``Difference``).

        Parameters
        ----------
        lag : int, optional
            The lag; default 1.
        differences : int, optional
            The order of differencing; default 1.
        label_at_later_date : bool, optional
            If True, date each difference at its later value, as pandas ``Series.diff`` does
            (``Difference(int, int, bool)``); default False, the earlier value.

        Returns
        -------
        TimeSeries
            The differenced series, shorter by ``lag * differences``. By default each difference
            carries the date of its earlier value, so the result starts at the start date; with
            ``label_at_later_date`` it starts ``lag * differences`` ordinates later, as
            ``s.diff(lag).dropna()``. The dates are the series' own, so a series with gaps or an
            irregular interval keeps its real dates.
        """
        if not label_at_later_date:
            return TimeSeries.from_net(self.net.Difference(int(lag), int(differences)))
        import System  # noqa: PLC0415 - the runtime is started

        member = self.net.Difference.Overloads[System.Int32, System.Int32, System.Boolean]
        return TimeSeries.from_net(member(int(lag), int(differences), True))

    def moving_average(self, period: int, min_valid_count: int | None = None) -> TimeSeries:
        """Return the trailing moving average over ``period`` steps (``MovingAverage``).

        Parameters
        ----------
        period : int
            The window length, less than the series length.
        min_valid_count : int, optional
            The fewest non-missing values a window needs, else NaN; default ``period``.

        Returns
        -------
        TimeSeries
            One value per full window, dated at the window's last date.
        """
        return TimeSeries.from_net(self.net.MovingAverage(int(period), self._optional(min_valid_count)))

    def moving_sum(self, period: int, min_valid_count: int | None = None) -> TimeSeries:
        """Return the trailing moving sum over ``period`` steps (``MovingSum``).

        Parameters
        ----------
        period : int
            The window length, less than the series length.
        min_valid_count : int, optional
            The fewest non-missing values a window needs, else NaN; default ``period``.

        Returns
        -------
        TimeSeries
            One value per full window, dated at the window's last date.
        """
        return TimeSeries.from_net(self.net.MovingSum(int(period), self._optional(min_valid_count)))

    @staticmethod
    def _optional(value: int | None) -> Any:
        """Return a ``Nullable<int>`` argument."""
        import System  # noqa: PLC0415 - the runtime is running

        return None if value is None else System.Nullable[System.Int32](int(value))

    def smoothed_series(self, smoothing_function: SmoothingFunction | str, period: int = 1) -> TimeSeries:
        """Return the series smoothed by a moving average, moving sum or difference (``SmoothedSeries``).

        Parameters
        ----------
        smoothing_function : SmoothingFunction or str
            The smoothing.
        period : int, optional
            The window length or lag; default 1 (a copy for the moving functions).

        Returns
        -------
        TimeSeries
            The smoothed series.
        """
        return TimeSeries.from_net(self.net.SmoothedSeries(net_enum(_data().SmoothingFunctionType,
                                                                    smoothing_function), int(period)))

    def shift_dates_by_day(self, number_of_days: int) -> TimeSeries:
        """Return the series with every date moved by whole days (``ShiftDatesByDay``).

        Parameters
        ----------
        number_of_days : int
            The shift; negative moves earlier.

        Returns
        -------
        TimeSeries
            The shifted series.
        """
        return TimeSeries.from_net(self.net.ShiftDatesByDay(int(number_of_days)))

    def shift_dates_by_month(self, number_of_months: int) -> TimeSeries:
        """Return the series with every date moved by whole months (``ShiftDatesByMonth``).

        Parameters
        ----------
        number_of_months : int
            The shift; negative moves earlier.

        Returns
        -------
        TimeSeries
            The shifted series.
        """
        return TimeSeries.from_net(self.net.ShiftDatesByMonth(int(number_of_months)))

    def shift_dates_by_year(self, number_of_years: int) -> TimeSeries:
        """Return the series with every date moved by whole years (``ShiftDatesByYear``).

        Parameters
        ----------
        number_of_years : int
            The shift; negative moves earlier.

        Returns
        -------
        TimeSeries
            The shifted series.
        """
        return TimeSeries.from_net(self.net.ShiftDatesByYear(int(number_of_years)))

    def clip(self, start_date: Any, end_date: Any) -> TimeSeries:
        """Return the ordinates dated within ``[start_date, end_date]`` (``ClipTimeSeries``).

        Parameters
        ----------
        start_date, end_date : datetime-like
            The window, inside the series' own dates.

        Returns
        -------
        TimeSeries
            The clipped series.
        """
        return TimeSeries.from_net(self.net.ClipTimeSeries(_net_date(start_date), _net_date(end_date)))

    def fill_missing_dates(self, start_date: Any, end_date: Any, value: float = 0.0) -> TimeSeries:
        """Return the series on every step from start to end, missing steps given a value (``FillMissingDates``).

        Parameters
        ----------
        start_date, end_date : datetime-like
            The range.
        value : float, optional
            The value of the added steps; default 0.

        Returns
        -------
        TimeSeries
            The completed series.
        """
        return TimeSeries.from_net(_data().TimeSeries.FillMissingDates(self.net, _net_date(start_date),
                                                                       _net_date(end_date), float(value)))

    def convert_time_interval(self, interval: TimeInterval | str, average: bool = True) -> TimeSeries | None:
        """Return the series at another time step (``ConvertTimeInterval``).

        A finer step interpolates linearly (``average``) or splits each value evenly; a coarser
        step averages (``average``) or sums consecutive blocks, ending with the last block that
        holds data, so 72 hourly values give 3 daily values. Both steps must have a fixed length (see
        :meth:`time_interval_in_hours`).

        Parameters
        ----------
        interval : TimeInterval or str
            The new step.
        average : bool, optional
            Whether values are rates (averaged) rather than volumes (summed); default True.

        Returns
        -------
        TimeSeries or None
            The converted series (a copy at the same step), or None when Numerics cannot convert.
        """
        net = self.net.ConvertTimeInterval(net_enum(_data().TimeInterval, interval), bool(average))
        return None if net is None else TimeSeries.from_net(net)

    def aggregate_by_interval(self, interval: TimeInterval | str,
                              block_function: BlockFunction | str = BlockFunction.AVERAGE,
                              min_valid_count: int = 1) -> TimeSeries:
        """Return one value per calendar-aligned block of a coarser step (``AggregateByInterval``).

        Blocks start at midnight, the first of a month, quarter or year (seven-day blocks at the
        first date); missing values are skipped and empty blocks are NaN.

        Parameters
        ----------
        interval : TimeInterval or str
            The block length (not irregular).
        block_function : BlockFunction or str, optional
            The reduction; default average.
        min_valid_count : int, optional
            The fewest non-missing values a block needs, else NaN; default 1.

        Returns
        -------
        TimeSeries
            One ordinate per block from the first to the last, dated at the block start.
        """
        data = _data()
        return TimeSeries.from_net(self.net.AggregateByInterval(
            net_enum(data.TimeInterval, interval), net_enum(data.BlockFunctionType, block_function),
            int(min_valid_count)))

    # Block series ------------------------------------------------------------------------------

    def _block(self, name: str, block_function: Any, smoothing_function: Any, period: int, *leading: int) -> Any:
        """Call a block-series member and wrap the result."""
        data = _data()
        return TimeSeries.from_net(getattr(self.net, name)(
            *leading, net_enum(data.BlockFunctionType, block_function),
            net_enum(data.SmoothingFunctionType, smoothing_function), int(period)))

    def calendar_year_series(self, block_function: BlockFunction | str = BlockFunction.MAXIMUM,
                             smoothing_function: SmoothingFunction | str = SmoothingFunction.NONE,
                             period: int = 1) -> TimeSeries:
        """Return one value per calendar year (``CalendarYearSeries``), for example the annual maxima.

        Parameters
        ----------
        block_function : BlockFunction or str, optional
            The reduction; default maximum (dated at the maximum).
        smoothing_function : SmoothingFunction or str, optional
            Smoothing applied first; default none.
        period : int, optional
            The smoothing window; default 1.

        Returns
        -------
        TimeSeries
            An irregular series, one ordinate per year with data.
        """
        return self._block("CalendarYearSeries", block_function, smoothing_function, period)

    def custom_year_series(self, start_month: int = 10, end_month: int | None = None,
                           block_function: BlockFunction | str = BlockFunction.MAXIMUM,
                           smoothing_function: SmoothingFunction | str = SmoothingFunction.NONE,
                           period: int = 1) -> TimeSeries:
        """Return one value per year starting in a given month, or per seasonal window (``CustomYearSeries``).

        With only ``start_month`` the year runs from that month (10, the default, is the US water
        year, labelled by the year it ends in); with ``end_month`` too, only the months from
        ``start_month`` to ``end_month`` (wrapping over the new year) are blocked.

        Parameters
        ----------
        start_month : int, optional
            The first month, 1 to 12; default 10.
        end_month : int, optional
            The last month of a seasonal window, 1 to 12.
        block_function : BlockFunction or str, optional
            The reduction; default maximum.
        smoothing_function : SmoothingFunction or str, optional
            Smoothing applied first; default none.
        period : int, optional
            The smoothing window; default 1.

        Returns
        -------
        TimeSeries
            An irregular series, one ordinate per year with data.
        """
        months = (int(start_month),) if end_month is None else (int(start_month), int(end_month))
        return self._block("CustomYearSeries", block_function, smoothing_function, period, *months)

    def monthly_series(self, block_function: BlockFunction | str = BlockFunction.MAXIMUM,
                       smoothing_function: SmoothingFunction | str = SmoothingFunction.NONE,
                       period: int = 1) -> TimeSeries:
        """Return one value per calendar month (``MonthlySeries``).

        Parameters
        ----------
        block_function : BlockFunction or str, optional
            The reduction; default maximum.
        smoothing_function : SmoothingFunction or str, optional
            Smoothing applied first; default none.
        period : int, optional
            The smoothing window; default 1.

        Returns
        -------
        TimeSeries
            An irregular series, one ordinate per month with data.
        """
        return self._block("MonthlySeries", block_function, smoothing_function, period)

    def quarterly_series(self, block_function: BlockFunction | str = BlockFunction.MAXIMUM,
                         smoothing_function: SmoothingFunction | str = SmoothingFunction.NONE,
                         period: int = 1) -> TimeSeries:
        """Return one value per calendar quarter (``QuarterlySeries``).

        Parameters
        ----------
        block_function : BlockFunction or str, optional
            The reduction; default maximum.
        smoothing_function : SmoothingFunction or str, optional
            Smoothing applied first; default none.
        period : int, optional
            The smoothing window; default 1.

        Returns
        -------
        TimeSeries
            An irregular series, one ordinate per quarter with data.
        """
        return self._block("QuarterlySeries", block_function, smoothing_function, period)

    def peaks_over_threshold(self, threshold: float, min_steps_between_events: int = 1,
                             smoothing_function: SmoothingFunction | str = SmoothingFunction.NONE,
                             period: int = 1) -> TimeSeries:
        """Return the peak of each cluster of values above a threshold (``PeaksOverThresholdSeries``).

        A cluster ends once the values have stayed at or below the threshold for more than
        ``min_steps_between_events`` steps after its peak.

        Parameters
        ----------
        threshold : float
            The threshold.
        min_steps_between_events : int, optional
            The independence separation in steps; default 1.
        smoothing_function : SmoothingFunction or str, optional
            Smoothing applied first; default none.
        period : int, optional
            The smoothing window; default 1.

        Returns
        -------
        TimeSeries
            An irregular series of the cluster peaks.
        """
        return TimeSeries.from_net(self.net.PeaksOverThresholdSeries(
            float(threshold), int(min_steps_between_events),
            net_enum(_data().SmoothingFunctionType, smoothing_function), int(period)))

    # Resampling and decomposition --------------------------------------------------------------

    def resample_with_knn(self, time_steps: int, k: int, seed: int = 12345) -> TimeSeries:
        """Return a synthetic series by k-nearest-neighbor conditional bootstrap (``ResampleWithKNN``).

        Parameters
        ----------
        time_steps : int
            The length of the synthetic series.
        k : int
            The number of neighbors.
        seed : int, optional
            The Mersenne Twister seed; default 12345.

        Returns
        -------
        TimeSeries
            The synthetic series, from the start date at this interval.
        """
        return TimeSeries.from_net(self.net.ResampleWithKNN(int(time_steps), int(k), int(seed)))

    def resample_with_block_bootstrap(self, time_steps: int, block_size: int, seed: int = 12345) -> TimeSeries:
        """Return a synthetic series by moving-block bootstrap (``ResampleWithBlockBootstrap``).

        Parameters
        ----------
        time_steps : int
            The length of the synthetic series.
        block_size : int
            The block length.
        seed : int, optional
            The Mersenne Twister seed; default 12345.

        Returns
        -------
        TimeSeries
            The synthetic series, from the start date at this interval.
        """
        return TimeSeries.from_net(self.net.ResampleWithBlockBootstrap(int(time_steps), int(block_size), int(seed)))

    def seasonal_decompose(self, period: int) -> tuple[TimeSeries, np.ndarray, TimeSeries]:
        """Split the series into trend, seasonal and residual parts (``SeasonalDecompose``).

        The trend is the trailing moving average over ``period``; the seasonal part keeps the
        harmonics of ``period`` in the FFT of the detrended series.

        Parameters
        ----------
        period : int
            The season length in steps (at least 2, at most half the series).

        Returns
        -------
        tuple of (TimeSeries, numpy.ndarray, TimeSeries)
            The trend, the seasonal component (one value per ordinate) and the residual (where
            the trend is defined).
        """
        result = self.net.SeasonalDecompose(int(period))
        return (TimeSeries.from_net(result.Item1), to_numpy_doubles(result.Item2),
                TimeSeries.from_net(result.Item3))

    # Statistics --------------------------------------------------------------------------------

    def min_value(self) -> float:
        """Return the smallest non-missing value (``MinValue``).

        Returns
        -------
        float
            The minimum.
        """
        return float(self.net.MinValue())

    def max_value(self) -> float:
        """Return the largest non-missing value (``MaxValue``).

        Returns
        -------
        float
            The maximum.
        """
        return float(self.net.MaxValue())

    def mean_value(self) -> float:
        """Return the mean of the non-missing values (``MeanValue``).

        Returns
        -------
        float
            The mean.
        """
        return float(self.net.MeanValue())

    def standard_deviation(self) -> float:
        """Return the sample standard deviation of the non-missing values (``StandardDeviation``).

        Returns
        -------
        float
            The standard deviation (n - 1 denominator).
        """
        return float(self.net.StandardDeviation())

    def percentiles(self, k: Any) -> np.ndarray:
        """Return percentiles of the non-missing values (``Percentiles``).

        Parameters
        ----------
        k : array_like
            Probabilities in [0, 1].

        Returns
        -------
        numpy.ndarray
            One percentile per probability (``Statistics.Percentile`` on the sorted values).
        """
        return to_numpy_doubles(self.net.Percentiles(to_net_doubles(k)))

    def summary_percentiles(self) -> np.ndarray:
        """Return the 5, 25, 50, 75 and 95 percent percentiles (``SummaryPercentiles``).

        Returns
        -------
        numpy.ndarray
            The five percentiles.
        """
        return to_numpy_doubles(self.net.SummaryPercentiles())

    def duration(self) -> np.ndarray:
        """Return the flow-duration curve (``Duration``).

        Returns
        -------
        numpy.ndarray
            Shape ``(n, 2)``: the Weibull exceedance percent and the value, largest value first.
        """
        return to_numpy_matrix(self.net.Duration())

    def monthly_percentiles(self, k: Any) -> np.ndarray:
        """Return percentiles of the values of each calendar month (``MonthlyPercentiles``).

        Parameters
        ----------
        k : array_like
            Probabilities in [0, 1].

        Returns
        -------
        numpy.ndarray
            Shape ``(12, len(k))``, January first; a month without data is zeros.
        """
        return to_numpy_matrix(self.net.MonthlyPercentiles(to_net_doubles(k)))

    def monthly_summary_statistics(self) -> np.ndarray:
        """Return summary statistics of each calendar month (``MonthlySummaryStatistics``).

        Returns
        -------
        numpy.ndarray
            Shape ``(12, 8)``, January first: minimum, 5, 25, 50, 75 and 95 percent percentiles,
            maximum and mean of the non-missing values.
        """
        return to_numpy_matrix(self.net.MonthlySummaryStatistics())

    def monthly_frequency(self) -> np.ndarray:
        """Return the number of ordinates in each calendar month (``MonthlyFrequency``).

        Returns
        -------
        numpy.ndarray
            Twelve counts, January first.
        """
        return to_numpy_doubles(self.net.MonthlyFrequency())

    def summary_statistics(self) -> dict[str, float]:
        """Return the record length, missing count, moments and percentiles (``SummaryStatistics``).

        Returns
        -------
        dict of str to float
            Keys ``"Record Length"``, ``"Missing Values"``, ``"Minimum"``, ``"Maximum"``,
            ``"Mean"``, ``"Std Dev"``, ``"Skewness"``, ``"Kurtosis"`` and ``"1%"`` to ``"99%"``.
        """
        return {str(pair.Key): float(pair.Value) for pair in self.net.SummaryStatistics()}

    def summary_hypothesis_test(self, split_location: int = -1) -> dict[str, float]:
        """Return the p-values of normality, independence, trend and homogeneity tests (``SummaryHypothesisTest``).

        Parameters
        ----------
        split_location : int, optional
            The index splitting the two samples of the two-sample tests; default the middle.

        Returns
        -------
        dict of str to float
            Test name to p-value (Jarque-Bera, Ljung-Box, Wald-Wolfowitz, Mann-Whitney,
            Mann-Kendall, Welch t and F tests).
        """
        return {str(pair.Key): float(pair.Value) for pair in self.net.SummaryHypothesisTest(int(split_location))}

    # Interval arithmetic -----------------------------------------------------------------------

    @staticmethod
    def add_time_interval(date: Any, interval: TimeInterval | str) -> np.datetime64:
        """Return a date one step later (``TimeSeries.AddTimeInterval``).

        Parameters
        ----------
        date : datetime-like
            The date.
        interval : TimeInterval or str
            The step (months, quarters and years are calendar steps).

        Returns
        -------
        numpy.datetime64
            The later date.
        """
        data = _data()
        return _py_date(data.TimeSeries.AddTimeInterval(_net_date(date), net_enum(data.TimeInterval, interval)))

    @staticmethod
    def anchored_dates(start: Any, count: int, interval: TimeInterval | str) -> np.ndarray:
        """Return ``count`` dates that count each step from the start (``TimeSeries.AnchoredDates``).

        Parameters
        ----------
        start : datetime-like
            The first date.
        count : int
            The number of dates, the start included.
        interval : TimeInterval or str
            The step (not irregular).

        Returns
        -------
        numpy.ndarray
            The ``datetime64[us]`` dates. Month, quarter and year dates keep the start's day of
            month where the month has it, as pandas ``Timestamp(start) + DateOffset(months=k)``:
            from 31 January 2020 the monthly dates are 31 January, 29 February, 31 March, 30 April.
            Fixed steps, and calendar steps from day 28 or earlier, give the dates of a regular
            series. Pass the dates as ``dates`` to :class:`TimeSeries` to build a series on them.

        Raises
        ------
        ValueError
            If ``interval`` is irregular or ``count`` is negative (.NET ``ArgumentException``).
        """
        data = _data()
        dates = data.TimeSeries.AnchoredDates(_net_date(start), int(count), net_enum(data.TimeInterval, interval))
        return from_ticks(np.fromiter((d.Ticks for d in dates), dtype=np.int64, count=len(dates)))

    @staticmethod
    def subtract_time_interval(date: Any, interval: TimeInterval | str) -> np.datetime64:
        """Return a date one step earlier (``TimeSeries.SubtractTimeInterval``).

        Parameters
        ----------
        date : datetime-like
            The date.
        interval : TimeInterval or str
            The step.

        Returns
        -------
        numpy.datetime64
            The earlier date.
        """
        data = _data()
        return _py_date(data.TimeSeries.SubtractTimeInterval(_net_date(date), net_enum(data.TimeInterval, interval)))

    @staticmethod
    def time_interval_in_hours(interval: TimeInterval | str) -> float:
        """Return the nominal length of a step in hours (``TimeSeries.TimeIntervalInHours``).

        Parameters
        ----------
        interval : TimeInterval or str
            The step.

        Returns
        -------
        float
            The hours of a fixed step up to seven days; NaN for months, quarters, years and
            irregular steps, which have no fixed length.
        """
        data = _data()
        return float(data.TimeSeries.TimeIntervalInHours(net_enum(data.TimeInterval, interval)))


__all__ = [
    "BlockFunction", "MathFunction", "SmoothingFunction", "TimeBlockWindow", "TimeInterval", "TimeSeries",
    "from_ticks", "to_ticks",
]

translate_public(globals())
