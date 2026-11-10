"""Paired data, interpolation, line simplification and regression from ``Numerics.Data``.

:class:`Ordinate`, :class:`OrderedPairedData`, :class:`UncertainOrdinate` and
:class:`UncertainOrderedPairedData` wrap the paired-data types; :class:`Linear`,
:class:`CubicSpline` and :class:`Polynomial` (all :class:`Interpolater`) and :class:`Bilinear`
wrap the interpolators; :func:`search` wraps ``Search``; :func:`ramer_douglas_peucker` and
:func:`perpendicular_distance` wrap ``LineSimplification``; :class:`ProbabilityOrdinates` and
:class:`LinearRegression` wrap their namesakes; :func:`time_series_index_ticks` and
:func:`time_series_values` read a ``Numerics.Data.TimeSeries`` in bulk (the full wrapper is
:class:`pynumerics.timeseries.TimeSeries`). Lookups on arrays cross once into the
``IList<double>`` overloads (``GetYFromX(IList, Transform, Transform, ExtrapolationSides)``,
``GetXFromY(...)``, ``Interpolater.Interpolate(IList)``), so the interpolation arithmetic is
Numerics' own.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Any

import numpy as np

from pynumerics import distributions
from pynumerics._util import (
    _data,
    bind,
    bind_list,
    evaluate,
    namespace,
    net_enum,
)
from pynumerics.enums import Extrapolation, NetEnum, SortOrder, Transform
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _numerics_matrix,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_longs,
    to_numpy_matrix,
)

_XELEMENT: list[Any] = []


def parse_xml(xml: str) -> Any:
    """Parse XML text into a ``System.Xml.Linq.XElement`` (``XElement.Parse``).

    The element is what the ``XElement`` constructors and ``FromXElement`` members of Numerics
    take; the text is what the wrappers' ``to_xml`` methods return.

    Parameters
    ----------
    xml : str
        One XML element.

    Returns
    -------
    System.Xml.Linq.XElement
        The parsed element.
    """
    if not _XELEMENT:
        _data()
        import clr  # noqa: PLC0415 - the runtime is running
        import System  # noqa: PLC0415

        # On .NET Core and later, System.Xml.Linq only forwards XElement to System.Private.Xml.Linq;
        # the .NET Framework defines it in System.Xml.Linq itself.
        name = "System.Private.Xml.Linq" if System.Type.GetType("System.Xml.Linq.XElement, System.Private.Xml.Linq") \
            else "System.Xml.Linq"
        clr.AddReference(name)
        from System.Xml.Linq import XElement  # noqa: PLC0415

        _XELEMENT.append(XElement.Parse)
    return _XELEMENT[0](str(xml))


def _sort_order(net: Any) -> SortOrder:
    """Return the Python member of a .NET ``SortOrder``."""
    return SortOrder(str(net))


class Ordinate:
    """One (x, y) point (``Numerics.Data.Ordinate``).

    Parameters
    ----------
    x, y : float
        The coordinates.

    Attributes
    ----------
    net : Numerics.Data.Ordinate
        The wrapped .NET value.
    """

    def __init__(self, x: float, y: float) -> None:
        self.net = _data().Ordinate(float(x), float(y))

    @classmethod
    def _from_net(cls, net: Any) -> Ordinate:
        """Wrap an existing .NET ordinate."""
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def from_xml(cls, xml: str) -> Ordinate:
        """Restore an ordinate from its XML (the ``XElement`` constructor).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        Ordinate
            The ordinate.
        """
        return cls._from_net(_data().Ordinate(parse_xml(xml)))

    def __repr__(self) -> str:
        return f"Ordinate({self.x!r}, {self.y!r})"

    @property
    def x(self) -> float:
        """The x value."""
        return float(self.net.X)

    @property
    def y(self) -> float:
        """The y value."""
        return float(self.net.Y)

    @property
    def is_valid(self) -> bool:
        """Whether both values are finite numbers."""
        return bool(self.net.IsValid)

    def transform(self, x_transform: Transform | str = Transform.NONE,
                  y_transform: Transform | str = Transform.NONE) -> Ordinate:
        """Return the ordinate in transformed spaces (``Transform``).

        Parameters
        ----------
        x_transform, y_transform : Transform or str, optional
            ``"logarithmic"`` is log10, ``"normal_z"`` the standard Normal quantile; default none.

        Returns
        -------
        Ordinate
            The transformed ordinate.
        """
        data = _data()
        return Ordinate._from_net(self.net.Transform(net_enum(data.Transform, x_transform),
                                                     net_enum(data.Transform, y_transform)))

    def ordinate_valid(self, other: Ordinate, strict_x: bool = True, strict_y: bool = True,
                       x_order: SortOrder | str = SortOrder.ASCENDING, y_order: SortOrder | str = SortOrder.ASCENDING,
                       other_is_next: bool = True) -> bool:
        """Return whether this ordinate and a neighbor satisfy an order (``OrdinateValid``).

        Parameters
        ----------
        other : Ordinate
            The neighbor.
        strict_x, strict_y : bool, optional
            Whether equal values break the order; default True.
        x_order, y_order : SortOrder or str, optional
            The required orders; default ascending.
        other_is_next : bool, optional
            Whether ``other`` follows this ordinate; default True.

        Returns
        -------
        bool
            Whether the pair is valid.
        """
        data = _data()
        return bool(self.net.OrdinateValid(other.net, bool(strict_x), bool(strict_y),
                                           net_enum(data.SortOrder, x_order), net_enum(data.SortOrder, y_order),
                                           bool(other_is_next)))

    def ordinate_errors(self, other: Ordinate | None = None, strict_x: bool = True, strict_y: bool = True,
                        x_order: SortOrder | str = SortOrder.ASCENDING,
                        y_order: SortOrder | str = SortOrder.ASCENDING, other_is_next: bool = True) -> list[str]:
        """Return the validation messages of this ordinate, or of it and a neighbor (``OrdinateErrors``).

        Parameters
        ----------
        other : Ordinate, optional
            The neighbor; without it only this ordinate's own values are checked.
        strict_x, strict_y : bool, optional
            Whether equal values break the order; default True.
        x_order, y_order : SortOrder or str, optional
            The required orders; default ascending.
        other_is_next : bool, optional
            Whether ``other`` follows this ordinate; default True.

        Returns
        -------
        list of str
            The messages (empty when valid).
        """
        if other is None:
            return [str(m) for m in self.net.OrdinateErrors()]
        data = _data()
        return [str(m) for m in self.net.OrdinateErrors(other.net, bool(strict_x), bool(strict_y),
                                                        net_enum(data.SortOrder, x_order),
                                                        net_enum(data.SortOrder, y_order), bool(other_is_next))]

    def to_xml(self) -> str:
        """Return the ordinate's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())


class OrderedPairedData:
    """Monotone (x, y) pairs with transformed linear lookups in both directions.

    Parameters
    ----------
    x, y : array_like
        Paired values of equal length.
    strict_x : bool, optional
        Whether x must be strictly monotone; default True.
    x_order : SortOrder or str, optional
        Order of x; default ascending.
    strict_y : bool, optional
        Whether y must be strictly monotone; default True.
    y_order : SortOrder or str, optional
        Order of y; default ascending.

    Attributes
    ----------
    net : Numerics.Data.OrderedPairedData
        The wrapped .NET object.
    """

    def __init__(self, x: Any, y: Any, strict_x: bool = True, x_order: SortOrder | str = SortOrder.ASCENDING,
                 strict_y: bool = True, y_order: SortOrder | str = SortOrder.ASCENDING) -> None:
        data = _data()
        self.net = data.OrderedPairedData(to_net_doubles(x), to_net_doubles(y), bool(strict_x),
                                          net_enum(data.SortOrder, x_order), bool(strict_y),
                                          net_enum(data.SortOrder, y_order))

    def _lookup(self, name: str) -> tuple[Any, Any]:
        """Return the scalar and list overloads of ``GetYFromX`` or ``GetXFromY`` with extrapolation sides.

        Both are bound once per wrapped object: pythonnet resolves the four-overload group on every
        call otherwise, which costs more than a short lookup. The bound overloads are the members the
        group resolves to for these arguments.
        """
        cache = self.__dict__.setdefault("_bound", {})
        entry = cache.get(name)
        if entry is None or entry[0] is not self.net:
            data = _data()
            import System  # noqa: PLC0415 - the runtime is running (self.net exists)
            from System.Collections.Generic import IList  # noqa: PLC0415

            group = getattr(self.net, name)
            tail = (data.Transform, data.Transform, data.ExtrapolationSides)
            entry = cache[name] = (self.net, group.Overloads[(System.Double, *tail)],
                                   group.Overloads[(IList[System.Double], *tail)])
        return entry[1], entry[2]

    def __len__(self) -> int:
        return int(self.net.Count)

    @property
    def x(self) -> np.ndarray:
        """The x values."""
        return np.array([float(o.X) for o in self.net], dtype=np.float64)

    @property
    def y(self) -> np.ndarray:
        """The y values."""
        return np.array([float(o.Y) for o in self.net], dtype=np.float64)

    @property
    def is_valid(self) -> bool:
        """Whether the pairs satisfy the configured order and strictness."""
        return bool(self.net.IsValid)

    @property
    def errors(self) -> list[str]:
        """The validation errors (``GetErrors``)."""
        return [str(message) for message in self.net.GetErrors()]

    def get_y_from_x(self, x: Any, x_transform: Transform | str = Transform.NONE,
                     y_transform: Transform | str = Transform.NONE,
                     extrapolation: Extrapolation | str = Extrapolation.NONE) -> Any:
        """Interpolate y at x (``GetYFromX``).

        Parameters
        ----------
        x : float or array_like
            Lookup values.
        x_transform, y_transform : Transform or str, optional
            Interpolation spaces; default none.
        extrapolation : Extrapolation or str, optional
            Sides on which to extend linearly; default none (endpoint hold).

        Returns
        -------
        float or numpy.ndarray
            Interpolated y, shaped like ``x``.
        """
        data = _data()
        xt, yt = net_enum(data.Transform, x_transform), net_enum(data.Transform, y_transform)
        ext = net_enum(data.ExtrapolationSides, extrapolation)
        scalar, vector = self._lookup("GetYFromX")
        return evaluate(x, lambda v: scalar(v, xt, yt, ext), lambda arr: vector(arr, xt, yt, ext))

    def get_x_from_y(self, y: Any, x_transform: Transform | str = Transform.NONE,
                     y_transform: Transform | str = Transform.NONE,
                     extrapolation: Extrapolation | str = Extrapolation.NONE) -> Any:
        """Interpolate x at y (``GetXFromY``).

        Parameters
        ----------
        y : float or array_like
            Lookup values.
        x_transform, y_transform : Transform or str, optional
            Interpolation spaces; default none.
        extrapolation : Extrapolation or str, optional
            Sides on which to extend linearly; default none (endpoint hold).

        Returns
        -------
        float or numpy.ndarray
            Interpolated x, shaped like ``y``.
        """
        data = _data()
        xt, yt = net_enum(data.Transform, x_transform), net_enum(data.Transform, y_transform)
        ext = net_enum(data.ExtrapolationSides, extrapolation)
        scalar, vector = self._lookup("GetXFromY")
        return evaluate(y, lambda v: scalar(v, xt, yt, ext), lambda arr: vector(arr, xt, yt, ext))

    @classmethod
    def _from_net(cls, net: Any) -> OrderedPairedData:
        """Wrap an existing .NET paired data object."""
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def from_ordinates(cls, ordinates: Sequence[Ordinate], strict_x: bool = True,
                       x_order: SortOrder | str = SortOrder.ASCENDING, strict_y: bool = True,
                       y_order: SortOrder | str = SortOrder.ASCENDING) -> OrderedPairedData:
        """Build paired data from ordinates (the ``IList<Ordinate>`` constructor).

        Parameters
        ----------
        ordinates : sequence of Ordinate
            The points, in order.
        strict_x : bool, optional
            Whether x must be strictly monotone; default True.
        x_order : SortOrder or str, optional
            Order of x; default ascending.
        strict_y : bool, optional
            Whether y must be strictly monotone; default True.
        y_order : SortOrder or str, optional
            Order of y; default ascending.

        Returns
        -------
        OrderedPairedData
            The paired data.
        """
        data = _data()
        import System  # noqa: PLC0415 - the runtime is running

        array = System.Array[data.Ordinate]([o.net for o in ordinates])
        return cls._from_net(data.OrderedPairedData(array, bool(strict_x), net_enum(data.SortOrder, x_order),
                                                    bool(strict_y), net_enum(data.SortOrder, y_order)))

    @classmethod
    def from_xml(cls, xml: str) -> OrderedPairedData:
        """Restore paired data from its XML (the ``XElement`` constructor).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        OrderedPairedData
            The paired data.
        """
        return cls._from_net(_data().OrderedPairedData(parse_xml(xml)))

    def to_xml(self) -> str:
        """Return the paired data's ``SaveToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.SaveToXElement().ToString())

    def __iter__(self) -> Iterator[Ordinate]:
        for net in self.net:
            yield Ordinate._from_net(net)

    def __getitem__(self, index: int) -> Ordinate:
        return Ordinate._from_net(self.net[int(index)])

    def __setitem__(self, index: int, ordinate: Ordinate) -> None:
        self.net[int(index)] = ordinate.net

    def __contains__(self, ordinate: Ordinate) -> bool:
        return bool(self.net.Contains(ordinate.net))

    @property
    def count(self) -> int:
        """The number of ordinates."""
        return int(self.net.Count)

    @property
    def strict_x(self) -> bool:
        """Whether x must be strictly monotone."""
        return bool(self.net.StrictX)

    @strict_x.setter
    def strict_x(self, value: bool) -> None:
        self.net.StrictX = bool(value)

    @property
    def strict_y(self) -> bool:
        """Whether y must be strictly monotone."""
        return bool(self.net.StrictY)

    @strict_y.setter
    def strict_y(self, value: bool) -> None:
        self.net.StrictY = bool(value)

    @property
    def order_x(self) -> SortOrder:
        """The required order of x."""
        return _sort_order(self.net.OrderX)

    @order_x.setter
    def order_x(self, value: SortOrder | str) -> None:
        self.net.OrderX = net_enum(_data().SortOrder, value)

    @property
    def order_y(self) -> SortOrder:
        """The required order of y."""
        return _sort_order(self.net.OrderY)

    @order_y.setter
    def order_y(self, value: SortOrder | str) -> None:
        self.net.OrderY = net_enum(_data().SortOrder, value)

    @property
    def use_smart_search(self) -> bool:
        """Whether lookups use bisection and hunting (True) or a sequential scan."""
        return bool(self.net.UseSmartSearch)

    @use_smart_search.setter
    def use_smart_search(self, value: bool) -> None:
        self.net.UseSmartSearch = bool(value)

    @property
    def x_search_start(self) -> int:
        """The cursor the next x search starts from."""
        return int(self.net.XSearchStart)

    @x_search_start.setter
    def x_search_start(self, value: int) -> None:
        self.net.XSearchStart = int(value)

    @property
    def y_search_start(self) -> int:
        """The cursor the next y search starts from."""
        return int(self.net.YSearchStart)

    @y_search_start.setter
    def y_search_start(self, value: int) -> None:
        self.net.YSearchStart = int(value)

    def append(self, x: float, y: float) -> None:
        """Append an ordinate (``Add``); the validity is updated.

        Parameters
        ----------
        x, y : float
            The point.
        """
        self.net.Add(_data().Ordinate(float(x), float(y)))

    def insert(self, index: int, x: float, y: float) -> None:
        """Insert an ordinate before a position (``Insert``).

        Parameters
        ----------
        index : int
            The position.
        x, y : float
            The point.
        """
        self.net.Insert(int(index), _data().Ordinate(float(x), float(y)))

    def index_of(self, x: float, y: float) -> int:
        """Return the position of an ordinate, or -1 (``IndexOf``).

        Parameters
        ----------
        x, y : float
            The point.

        Returns
        -------
        int
            The first position holding the point, or -1.
        """
        return int(self.net.IndexOf(float(x), float(y)))

    def remove(self, x: float, y: float) -> bool:
        """Remove the first ordinate equal to a point (``Remove``).

        Parameters
        ----------
        x, y : float
            The point.

        Returns
        -------
        bool
            Whether an ordinate was removed.
        """
        return bool(self.net.Remove(_data().Ordinate(float(x), float(y))))

    def remove_at(self, index: int) -> None:
        """Remove the ordinate at a position (``RemoveAt``).

        Parameters
        ----------
        index : int
            The position.
        """
        self.net.RemoveAt(int(index))

    def remove_range(self, index: int, count: int) -> None:
        """Remove a run of ordinates (``RemoveRange``).

        Parameters
        ----------
        index : int
            The first position.
        count : int
            How many to remove.
        """
        self.net.RemoveRange(int(index), int(count))

    def clear(self) -> None:
        """Remove every ordinate (``Clear``)."""
        self.net.Clear()

    def copy(self) -> OrderedPairedData:
        """Return an independent copy (``Clone``).

        Returns
        -------
        OrderedPairedData
            The copy, validated.
        """
        return OrderedPairedData._from_net(self.net.Clone())

    def copy_without_validation(self) -> OrderedPairedData:
        """Return an independent copy that skips revalidation (``CloneWithoutValidation``).

        Returns
        -------
        OrderedPairedData
            The copy.
        """
        return OrderedPairedData._from_net(self.net.CloneWithoutValidation())

    def invert(self) -> OrderedPairedData:
        """Return the paired data with x and y swapped (``Invert``).

        Returns
        -------
        OrderedPairedData
            The inverted paired data.
        """
        return OrderedPairedData._from_net(self.net.Invert())

    def trapezoidal_area_under_y(self) -> float:
        """Return the trapezoidal area between the curve and the x axis (``TrapezoidalAreaUnderY``).

        Returns
        -------
        float
            ``sum (x[i+1] - x[i]) (y[i] + y[i+1]) / 2``.
        """
        return float(self.net.TrapezoidalAreaUnderY())

    def trapezoidal_area_under_x(self) -> float:
        """Return the trapezoidal area between the curve and the y axis (``TrapezoidalAreaUnderX``).

        Returns
        -------
        float
            ``sum (y[i+1] - y[i]) (x[i] + x[i+1]) / 2``.
        """
        return float(self.net.TrapezoidalAreaUnderX())

    def search_x(self, x: float) -> int:
        """Return the lower bracketing index of x with the configured search (``SearchX``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The index ``i`` with x between ``x[i]`` and ``x[i + 1]``.
        """
        return int(self.net.SearchX(float(x)))

    def search_y(self, y: float) -> int:
        """Return the lower bracketing index of y with the configured search (``SearchY``).

        Parameters
        ----------
        y : float
            The value.

        Returns
        -------
        int
            The index ``i`` with y between ``y[i]`` and ``y[i + 1]``.
        """
        return int(self.net.SearchY(float(y)))

    def binary_search_x(self, x: float) -> int:
        """Return the index of x by binary search (``BinarySearchX``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The index of an equal x, or the bitwise complement of the insertion point (as
            ``List.BinarySearch``).
        """
        return int(self.net.BinarySearchX(float(x)))

    def binary_search_y(self, y: float) -> int:
        """Return the index of y by binary search (``BinarySearchY``).

        Parameters
        ----------
        y : float
            The value.

        Returns
        -------
        int
            The index of an equal y, or the bitwise complement of the insertion point.
        """
        return int(self.net.BinarySearchY(float(y)))

    def sequential_search_x(self, x: float) -> int:
        """Return the lower bracketing index of x by a sequential scan (``SequentialSearchX``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.SequentialSearchX(float(x)))

    def sequential_search_y(self, y: float) -> int:
        """Return the lower bracketing index of y by a sequential scan (``SequentialSearchY``).

        Parameters
        ----------
        y : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.SequentialSearchY(float(y)))

    def bisection_search_x(self, x: float) -> int:
        """Return the lower bracketing index of x by bisection (``BisectionSearchX``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.BisectionSearchX(float(x)))

    def bisection_search_y(self, y: float) -> int:
        """Return the lower bracketing index of y by bisection (``BisectionSearchY``).

        Parameters
        ----------
        y : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.BisectionSearchY(float(y)))

    def hunt_search_x(self, x: float) -> int:
        """Return the lower bracketing index of x by hunting from the cursor (``HuntSearchX``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.HuntSearchX(float(x)))

    def hunt_search_y(self, y: float) -> int:
        """Return the lower bracketing index of y by hunting from the cursor (``HuntSearchY``).

        Parameters
        ----------
        y : float
            The value.

        Returns
        -------
        int
            The lower bracketing index.
        """
        return int(self.net.HuntSearchY(float(y)))

    def douglas_peucker_simplify(self, tolerance: float) -> OrderedPairedData:
        """Return the curve simplified by the Douglas-Peucker algorithm (``DouglasPeuckerSimplify``).

        Parameters
        ----------
        tolerance : float
            The largest perpendicular distance a removed point may have from the simplified curve.

        Returns
        -------
        OrderedPairedData
            The kept ordinates, end points included.
        """
        return OrderedPairedData._from_net(self.net.DouglasPeuckerSimplify(float(tolerance)))

    def visvalingam_whyatt_simplify(self, number_to_keep: int) -> OrderedPairedData:
        """Return the curve simplified by the Visvalingam-Whyatt algorithm (``VisvaligamWhyattSimplify``).

        Points with the smallest effective triangle area are removed one at a time.

        Parameters
        ----------
        number_to_keep : int
            How many ordinates remain.

        Returns
        -------
        OrderedPairedData
            The kept ordinates.
        """
        return OrderedPairedData._from_net(self.net.VisvaligamWhyattSimplify(int(number_to_keep)))

    def lang_simplify(self, tolerance: float, look_ahead: int) -> OrderedPairedData:
        """Return the curve simplified by the Lang algorithm (``LangSimplify``).

        Parameters
        ----------
        tolerance : float
            The largest perpendicular distance a removed point may have from its chord.
        look_ahead : int
            The number of points each search region spans.

        Returns
        -------
        OrderedPairedData
            The kept ordinates (a copy when ``look_ahead <= 1`` or ``tolerance <= 0``).
        """
        return OrderedPairedData._from_net(self.net.LangSimplify(float(tolerance), int(look_ahead)))


class Interpolater:
    """Shared members of the one-dimensional interpolators (``Numerics.Data.Interpolater``).

    Construct :class:`Linear`, :class:`CubicSpline` or :class:`Polynomial`. Build once and
    reuse: the interpolator keeps a search cursor between calls.

    Attributes
    ----------
    net : Numerics.Data.Interpolater
        The wrapped .NET object.
    """

    net: Any

    def _bind(self) -> None:
        """Bind the overloads of ``Interpolate`` once per wrapped object."""
        self._interpolate = bind(self.net.Interpolate)
        self._interpolate_list = bind_list(self.net.Interpolate)

    def __len__(self) -> int:
        return int(self.net.Count)

    @property
    def count(self) -> int:
        """The number of knots."""
        return int(self.net.Count)

    @property
    def x_values(self) -> np.ndarray:
        """The knot abscissas."""
        return to_numpy_doubles(self.net.XValues)

    @property
    def y_values(self) -> np.ndarray:
        """The knot ordinates."""
        return to_numpy_doubles(self.net.YValues)

    @property
    def sort_order(self) -> SortOrder:
        """The order of the abscissas."""
        return _sort_order(self.net.SortOrder)

    @property
    def use_smart_search(self) -> bool:
        """Whether :meth:`search` hunts or bisects (True, default) or scans sequentially."""
        return bool(self.net.UseSmartSearch)

    @use_smart_search.setter
    def use_smart_search(self, value: bool) -> None:
        self.net.UseSmartSearch = bool(value)

    @property
    def search_start(self) -> int:
        """The cursor the next search starts from."""
        return int(self.net.SearchStart)

    @search_start.setter
    def search_start(self, value: int) -> None:
        self.net.SearchStart = int(value)

    def interpolate(self, x: Any) -> Any:
        """Interpolate at x (``Interpolate``).

        Parameters
        ----------
        x : float or array_like
            Lookup values.

        Returns
        -------
        float or numpy.ndarray
            Interpolated values, shaped like ``x``.
        """
        return evaluate(x, self._interpolate, self._interpolate_list)

    def __call__(self, x: Any) -> Any:
        return self.interpolate(x)

    def base_interpolate(self, x: float, index: int) -> float:
        """Evaluate the interpolant of the knot interval starting at an index (``BaseInterpolate``).

        Parameters
        ----------
        x : float
            The lookup value.
        index : int
            The lower bracketing knot, as :meth:`search` returns it.

        Returns
        -------
        float
            The interpolated value.
        """
        return float(self.net.BaseInterpolate(float(x), int(index)))

    def search(self, x: float) -> int:
        """Return the lower bracketing knot of x and move the cursor (``Search``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The index ``i`` with x between ``x_values[i]`` and ``x_values[i + 1]``.
        """
        return int(self.net.Search(float(x)))

    def sequential_search(self, x: float) -> int:
        """Return the lower bracketing knot of x by a sequential scan from the cursor (``SequentialSearch``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing knot.
        """
        return int(self.net.SequentialSearch(float(x)))

    def bisection_search(self, x: float) -> int:
        """Return the lower bracketing knot of x by bisection (``BisectionSearch``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing knot.
        """
        return int(self.net.BisectionSearch(float(x)))

    def hunt_search(self, x: float) -> int:
        """Return the lower bracketing knot of x by hunting from the cursor (``HuntSearch``).

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        int
            The lower bracketing knot.
        """
        return int(self.net.HuntSearch(float(x)))


class Linear(Interpolater):
    """Linear interpolation in optionally transformed spaces (``Numerics.Data.Linear``).

    :meth:`interpolate` holds the end values outside the knots; :meth:`extrapolate` extends the
    end segments.

    Parameters
    ----------
    x, y : array_like
        Knots; x sorted in ``sort_order``.
    sort_order : SortOrder or str, optional
        Order of x; default ascending.
    x_transform, y_transform : Transform or str, optional
        Interpolation spaces; default none.

    Attributes
    ----------
    net : Numerics.Data.Linear
        The wrapped .NET object.
    """

    def __init__(self, x: Any, y: Any, sort_order: SortOrder | str = SortOrder.ASCENDING,
                 x_transform: Transform | str = Transform.NONE, y_transform: Transform | str = Transform.NONE) -> None:
        data = _data()
        self.net = data.Linear(to_net_doubles(x), to_net_doubles(y), net_enum(data.SortOrder, sort_order))
        self.net.XTransform = net_enum(data.Transform, x_transform)
        self.net.YTransform = net_enum(data.Transform, y_transform)
        self._bind()

    def _bind(self) -> None:
        """Bind the overloads of ``Interpolate`` and ``Extrapolate`` once per wrapped object."""
        super()._bind()
        self._extrapolate = bind(self.net.Extrapolate)
        self._extrapolate_list = bind_list(self.net.Extrapolate)

    @property
    def x_transform(self) -> Transform:
        """The interpolation space of x."""
        return Transform(str(self.net.XTransform))

    @x_transform.setter
    def x_transform(self, value: Transform | str) -> None:
        self.net.XTransform = net_enum(_data().Transform, value)

    @property
    def y_transform(self) -> Transform:
        """The interpolation space of y."""
        return Transform(str(self.net.YTransform))

    @y_transform.setter
    def y_transform(self, value: Transform | str) -> None:
        self.net.YTransform = net_enum(_data().Transform, value)

    def extrapolate(self, x: Any) -> Any:
        """Evaluate the end segments' straight lines beyond the knots (``Extrapolate``).

        Inside the knots this equals :meth:`interpolate`; outside it extends the first or
        last segment in the transformed spaces instead of holding the end values. An array
        crosses once through ``Extrapolate(IList<double>)``, with the scalar call's bits.

        Parameters
        ----------
        x : float or array_like
            Lookup values.

        Returns
        -------
        float or numpy.ndarray
            Values, shaped like ``x``.
        """
        return evaluate(x, self._extrapolate, self._extrapolate_list)

    def shallow_clone(self) -> Linear:
        """Return an interpolator over the same knot lists with its own search cursor (``ShallowClone``).

        Returns
        -------
        Linear
            The clone; it shares the knots, so it is cheap for one interpolator per thread.
        """
        clone = Linear.__new__(Linear)
        clone.net = self.net.ShallowClone()
        clone._bind()
        return clone

    def __deepcopy__(self, memo: dict) -> Linear:
        """Return an independent interpolator over the same knots (its own .NET object and search cursor)."""
        return Linear(self.x_values, self.y_values, self.sort_order, self.x_transform, self.y_transform)


class CubicSpline(Interpolater):
    """Natural cubic spline interpolation (``Numerics.Data.CubicSpline``).

    The second derivatives vanish at both end knots (Press et al., Numerical Recipes, 3rd ed.,
    section 3.3); outside the knots the end cubic pieces are extended.

    Parameters
    ----------
    x, y : array_like
        Knots; x sorted in ``sort_order``.
    sort_order : SortOrder or str, optional
        Order of x; default ascending.

    Attributes
    ----------
    net : Numerics.Data.CubicSpline
        The wrapped .NET object.
    """

    def __init__(self, x: Any, y: Any, sort_order: SortOrder | str = SortOrder.ASCENDING) -> None:
        data = _data()
        self.net = data.CubicSpline(to_net_doubles(x), to_net_doubles(y), net_enum(data.SortOrder, sort_order))
        self._bind()

    def set_second_derivatives(self) -> None:
        """Recompute the knot second derivatives of the natural spline (``SetSecondDerivatives``)."""
        self.net.SetSecondDerivatives()


class Polynomial(Interpolater):
    """Polynomial interpolation through ``order + 1`` neighboring knots (``Numerics.Data.Polynomial``).

    Neville's algorithm on the knots starting at the bracketing index (Press et al., Numerical
    Recipes, 3rd ed., section 3.2); :attr:`error` is the estimate of the last evaluation.

    Parameters
    ----------
    order : int
        The polynomial degree; less than the number of knots.
    x, y : array_like
        Knots; x sorted in ``sort_order``.
    sort_order : SortOrder or str, optional
        Order of x; default ascending.

    Attributes
    ----------
    net : Numerics.Data.Polynomial
        The wrapped .NET object.
    """

    def __init__(self, order: int, x: Any, y: Any, sort_order: SortOrder | str = SortOrder.ASCENDING) -> None:
        data = _data()
        self.net = data.Polynomial(int(order), to_net_doubles(x), to_net_doubles(y),
                                   net_enum(data.SortOrder, sort_order))
        self._bind()

    @property
    def order(self) -> int:
        """The polynomial degree."""
        return int(self.net.Order)

    @order.setter
    def order(self, value: int) -> None:
        self.net.Order = int(value)

    @property
    def error(self) -> float:
        """The error estimate of the last interpolation: its last Neville correction."""
        return float(self.net.Error)


def _pair_overload(member: Any) -> Any:
    """Return the ``(IList<double>, IList<double>)`` overload of a .NET method group."""
    import System  # noqa: PLC0415 - the runtime is running
    from System.Collections.Generic import IList  # noqa: PLC0415

    return member.Overloads[IList[System.Double], IList[System.Double]]


class Bilinear:
    """Bilinear interpolation on a rectangular grid (``Numerics.Data.Bilinear``).

    Outside the grid the nearest edge or corner value is held. Arrays of points cross once
    through ``Interpolate(IList<double>, IList<double>)``, with the scalar call's bits.

    Parameters
    ----------
    x1, x2 : array_like
        Grid coordinates, sorted in ``sort_order``.
    y : array_like
        Shape ``(len(x1), len(x2))``: ``y[i, j]`` is the value at ``(x1[i], x2[j])``.
    sort_order : SortOrder or str, optional
        Order of both coordinates; default ascending.
    x1_transform, x2_transform, y_transform : Transform or str, optional
        Interpolation spaces; default none.

    Attributes
    ----------
    net : Numerics.Data.Bilinear
        The wrapped .NET object.
    """

    def __init__(self, x1: Any, x2: Any, y: Any, sort_order: SortOrder | str = SortOrder.ASCENDING,
                 x1_transform: Transform | str = Transform.NONE, x2_transform: Transform | str = Transform.NONE,
                 y_transform: Transform | str = Transform.NONE) -> None:
        data = _data()
        self.net = data.Bilinear(to_net_doubles(x1), to_net_doubles(x2), to_net_matrix(y),
                                 net_enum(data.SortOrder, sort_order))
        self.net.X1Transform = net_enum(data.Transform, x1_transform)
        self.net.X2Transform = net_enum(data.Transform, x2_transform)
        self.net.YTransform = net_enum(data.Transform, y_transform)
        self._interpolate = bind(self.net.Interpolate, 2)
        self._interpolate_list = _pair_overload(self.net.Interpolate)

    @property
    def x1_values(self) -> np.ndarray:
        """The first grid coordinates."""
        return to_numpy_doubles(self.net.X1Values)

    @property
    def x2_values(self) -> np.ndarray:
        """The second grid coordinates."""
        return to_numpy_doubles(self.net.X2Values)

    @property
    def y_values(self) -> np.ndarray:
        """The grid values, shape ``(len(x1), len(x2))``."""
        return to_numpy_matrix(self.net.YValues)

    @property
    def sort_order(self) -> SortOrder:
        """The order of the coordinates."""
        return _sort_order(self.net.SortOrder)

    @property
    def use_smart_search(self) -> bool:
        """Whether the searches hunt or bisect (True, default) or scan sequentially."""
        return bool(self.net.UseSmartSearch)

    @use_smart_search.setter
    def use_smart_search(self, value: bool) -> None:
        self.net.UseSmartSearch = bool(value)

    @property
    def x1_transform(self) -> Transform:
        """The interpolation space of x1."""
        return Transform(str(self.net.X1Transform))

    @x1_transform.setter
    def x1_transform(self, value: Transform | str) -> None:
        self.net.X1Transform = net_enum(_data().Transform, value)

    @property
    def x2_transform(self) -> Transform:
        """The interpolation space of x2."""
        return Transform(str(self.net.X2Transform))

    @x2_transform.setter
    def x2_transform(self, value: Transform | str) -> None:
        self.net.X2Transform = net_enum(_data().Transform, value)

    @property
    def y_transform(self) -> Transform:
        """The interpolation space of y."""
        return Transform(str(self.net.YTransform))

    @y_transform.setter
    def y_transform(self, value: Transform | str) -> None:
        self.net.YTransform = net_enum(_data().Transform, value)

    def interpolate(self, x1: Any, x2: Any) -> Any:
        """Interpolate at points (``Interpolate``).

        Parameters
        ----------
        x1, x2 : float or array_like
            Coordinates, broadcast against each other.

        Returns
        -------
        float or numpy.ndarray
            Interpolated values, shaped like the broadcast coordinates.
        """
        a, b = np.broadcast_arrays(np.asarray(x1, dtype=np.float64), np.asarray(x2, dtype=np.float64))
        member = self._interpolate
        if a.ndim == 0:
            return float(member(float(a), float(b)))
        values = self._interpolate_list(to_net_doubles(a.ravel()), to_net_doubles(b.ravel()))
        return to_numpy_doubles(values).reshape(a.shape)

    def __call__(self, x1: Any, x2: Any) -> Any:
        return self.interpolate(x1, x2)


class SearchMethod(NetEnum):
    """Bracketing search methods (``Numerics.Data.Search.Method``)."""

    SEQUENTIAL = "Sequential"
    BISECTION = "Bisection"
    HUNT = "Hunt"


def search(x: float, values: Any, method: SearchMethod | str = SearchMethod.BISECTION, start: int = 0,
           order: SortOrder | str = SortOrder.ASCENDING) -> int:
    """Return the lower bracketing index of x in ordered values (``Search.Sequential``, ``Bisection`` or ``Hunt``).

    Parameters
    ----------
    x : float
        The value.
    values : array_like or OrderedPairedData
        Ordered values, or paired data searched on x (in its own order).
    method : SearchMethod or str, optional
        The search; default bisection.
    start : int, optional
        The index the sequential and hunt searches start from; default 0.
    order : SortOrder or str, optional
        Order of ``values``; default ascending (ignored for paired data).

    Returns
    -------
    int
        The index ``i`` with x between ``values[i]`` and ``values[i + 1]``: -1 below the first
        value and ``len(values)`` above the last.
    """
    data = _data()
    member = getattr(data.Search, SearchMethod(method).value)
    if isinstance(values, OrderedPairedData):
        return int(member(float(x), values.net, int(start)))
    return int(member(float(x), to_net_doubles(values), int(start), net_enum(data.SortOrder, order)))


def _predictors(x: Any) -> Any:
    """Return a predictor ``Matrix``: one row per observation, a 1-D input as a single column."""
    array = np.asarray(x, dtype=np.float64)
    return _numerics_matrix(array.reshape(-1, 1) if array.ndim == 1 else array)


class LinearRegression:
    """Ordinary least-squares regression ``y = X b + e`` (``Numerics.Data.LinearRegression``).

    Numerics fits the model by singular value decomposition at construction; the estimates,
    standard errors, t statistics and goodness of fit match R's ``lm``.

    Parameters
    ----------
    x : array_like
        Predictors: one row per observation and one column per predictor; a 1-D array is a single
        predictor.
    y : array_like
        Responses, one per row of ``x`` (at least three).
    has_intercept : bool, optional
        Whether to estimate an intercept, the first parameter; default True.

    Attributes
    ----------
    net : Numerics.Data.LinearRegression
        The wrapped .NET object.
    """

    def __init__(self, x: Any, y: Any, has_intercept: bool = True) -> None:
        linear_algebra = namespace("Numerics.Mathematics.LinearAlgebra")
        self.net = _data().LinearRegression(_predictors(x), linear_algebra.Vector(to_net_doubles(y)),
                                            bool(has_intercept))

    @property
    def has_intercept(self) -> bool:
        """Whether the model has an intercept."""
        return bool(self.net.HasIntercept)

    @property
    def parameters(self) -> np.ndarray:
        """The estimates, the intercept first when the model has one."""
        return to_numpy_doubles(self.net.Parameters)

    @property
    def parameter_names(self) -> list[str]:
        """The parameter names (``"Intercept"``, then the predictor names)."""
        return [str(name) for name in self.net.ParameterNames]

    @property
    def parameter_standard_errors(self) -> np.ndarray:
        """The standard errors of the estimates."""
        return to_numpy_doubles(self.net.ParameterStandardErrors)

    @property
    def parameter_t_stats(self) -> np.ndarray:
        """The t statistics of the estimates."""
        return to_numpy_doubles(self.net.ParameterTStats)

    @property
    def covariance(self) -> np.ndarray:
        """The unscaled covariance ``(X^T X)^-1``; times ``standard_error**2`` it is the estimates' covariance."""
        return to_numpy_matrix(self.net.Covariance.ToArray())

    @property
    def residuals(self) -> np.ndarray:
        """The residuals ``y - X b``."""
        return to_numpy_doubles(self.net.Residuals)

    @property
    def standard_error(self) -> float:
        """The residual standard error ``sqrt(SSE / (n - p))``."""
        return float(self.net.StandardError)

    @property
    def sample_size(self) -> int:
        """The number of observations."""
        return int(self.net.SampleSize)

    @property
    def degrees_of_freedom(self) -> int:
        """The residual degrees of freedom ``n - p``."""
        return int(self.net.DegreesOfFreedom)

    @property
    def r_squared(self) -> float:
        """The coefficient of determination."""
        return float(self.net.RSquared)

    @property
    def adjusted_r_squared(self) -> float:
        """The adjusted coefficient of determination."""
        return float(self.net.AdjRSquared)

    @property
    def x(self) -> np.ndarray:
        """The design matrix, with the leading column of ones when the model has an intercept."""
        return to_numpy_matrix(self.net.X.ToArray())

    @property
    def y(self) -> np.ndarray:
        """The responses."""
        return to_numpy_doubles(self.net.Y.ToArray())

    def predict(self, x: Any) -> np.ndarray:
        """Return the fitted mean at new predictor rows (``Predict``).

        Parameters
        ----------
        x : array_like
            Predictor rows, without the intercept column (or with it, as a leading column of ones).

        Returns
        -------
        numpy.ndarray
            One prediction per row.
        """
        return to_numpy_doubles(self.net.Predict(_predictors(x)))

    def prediction_intervals(self, x: Any, alpha: float = 0.1) -> np.ndarray:
        """Return two-sided prediction intervals at new predictor rows (``PredictionIntervals``).

        The Student-t interval of R's ``predict(..., interval = "prediction")``, with level
        ``1 - alpha``.

        Parameters
        ----------
        x : array_like
            Predictor rows, as for :meth:`predict`.
        alpha : float, optional
            One minus the interval level; default 0.1 (90 percent intervals).

        Returns
        -------
        numpy.ndarray
            Shape ``(rows, 3)``: the lower bound, the upper bound and the prediction.
        """
        return to_numpy_matrix(self.net.PredictionIntervals(_predictors(x), float(alpha)))

    def summary(self) -> list[str]:
        """Return the R-style text summary (``Summary``).

        Returns
        -------
        list of str
            The summary lines.
        """
        return [str(line) for line in self.net.Summary()]


def _net_distribution(distribution: Any) -> Any:
    """Return the .NET distribution behind a pynumerics wrapper, or a .NET distribution as it is."""
    return getattr(distribution, "net", distribution)


def _wrap_distribution(net: Any) -> Any:
    """Wrap a .NET distribution, keeping None as None."""
    return None if net is None else distributions.wrap(net)


class UncertainOrdinate:
    """An x value with a distribution of y (``Numerics.Data.UncertainOrdinate``).

    Parameters
    ----------
    x : float
        The x value.
    distribution : UnivariateDistribution
        The distribution of y at x.

    Attributes
    ----------
    net : Numerics.Data.UncertainOrdinate
        The wrapped .NET value.
    """

    def __init__(self, x: float, distribution: Any) -> None:
        self.net = _data().UncertainOrdinate(float(x), _net_distribution(distribution))

    @classmethod
    def _from_net(cls, net: Any) -> UncertainOrdinate:
        """Wrap an existing .NET uncertain ordinate."""
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def from_xml(cls, xml: str, distribution_type: str | None = None) -> UncertainOrdinate:
        """Restore an uncertain ordinate from its XML (the ``XElement`` constructors).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.
        distribution_type : str, optional
            A ``UnivariateDistributionType`` name: read the legacy layout, whose distribution
            parameters are attributes of the element (``X="1" Mu="10" Sigma="2"``); default the
            current layout, with a child ``Distribution`` element.

        Returns
        -------
        UncertainOrdinate
            The uncertain ordinate.
        """
        data = _data()
        element = parse_xml(xml)
        if distribution_type is None:
            return cls._from_net(data.UncertainOrdinate(element))
        kind = net_enum(namespace("Numerics.Distributions").UnivariateDistributionType, distribution_type)
        return cls._from_net(data.UncertainOrdinate(element, kind))

    @property
    def x(self) -> float:
        """The x value."""
        return float(self.net.X)

    @property
    def y(self) -> Any:
        """The distribution of y (None when unset)."""
        return _wrap_distribution(self.net.Y)

    @property
    def is_valid(self) -> bool:
        """Whether x is finite and the distribution's parameters are valid."""
        return bool(self.net.IsValid)

    def get_ordinate(self, probability: float | None = None) -> Ordinate:
        """Return the point at a quantile of y, or at its mean (``GetOrdinate``).

        Parameters
        ----------
        probability : float, optional
            The non-exceedance probability of y; default the mean.

        Returns
        -------
        Ordinate
            ``(x, y.inverse_cdf(probability))`` or ``(x, y.mean)``.
        """
        net = self.net.GetOrdinate() if probability is None else self.net.GetOrdinate(float(probability))
        return Ordinate._from_net(net)

    def _comparison(self, other: UncertainOrdinate, strict_x: bool, strict_y: bool, x_order: Any, y_order: Any,
                    other_is_next: bool, allow_different_types: bool, tail_probability: float | None) -> tuple:
        """Return the .NET arguments of ``OrdinateValid`` and ``OrdinateErrors``."""
        data = _data()
        arguments = (other.net, bool(strict_x), bool(strict_y), net_enum(data.SortOrder, x_order),
                     net_enum(data.SortOrder, y_order), bool(other_is_next), bool(allow_different_types))
        return arguments if tail_probability is None else (*arguments, float(tail_probability))

    def ordinate_valid(self, other: UncertainOrdinate, strict_x: bool = True, strict_y: bool = True,
                       x_order: SortOrder | str = SortOrder.ASCENDING, y_order: SortOrder | str = SortOrder.ASCENDING,
                       other_is_next: bool = True, allow_different_types: bool = False,
                       tail_probability: float | None = None) -> bool:
        """Return whether this ordinate and a neighbor satisfy an order (``OrdinateValid``).

        The y order is checked on the distributions' quantiles (at the tail probability when given).

        Parameters
        ----------
        other : UncertainOrdinate
            The neighbor.
        strict_x, strict_y : bool, optional
            Whether equal values break the order; default True.
        x_order, y_order : SortOrder or str, optional
            The required orders; default ascending.
        other_is_next : bool, optional
            Whether ``other`` follows this ordinate; default True.
        allow_different_types : bool, optional
            Whether the two distributions may be of different families; default False.
        tail_probability : float, optional
            The tail probability the y order is checked at; default the Numerics default.

        Returns
        -------
        bool
            Whether the pair is valid.
        """
        return bool(self.net.OrdinateValid(*self._comparison(other, strict_x, strict_y, x_order, y_order,
                                                             other_is_next, allow_different_types, tail_probability)))

    def ordinate_errors(self, other: UncertainOrdinate | None = None, strict_x: bool = True, strict_y: bool = True,
                        x_order: SortOrder | str = SortOrder.ASCENDING,
                        y_order: SortOrder | str = SortOrder.ASCENDING, other_is_next: bool = True,
                        allow_different_types: bool = False, tail_probability: float | None = None) -> list[str]:
        """Return the validation messages of this ordinate, or of it and a neighbor (``OrdinateErrors``).

        Parameters
        ----------
        other : UncertainOrdinate, optional
            The neighbor; without it only this ordinate is checked.
        strict_x, strict_y : bool, optional
            Whether equal values break the order; default True.
        x_order, y_order : SortOrder or str, optional
            The required orders; default ascending.
        other_is_next : bool, optional
            Whether ``other`` follows this ordinate; default True.
        allow_different_types : bool, optional
            Whether the two distributions may be of different families; default False.
        tail_probability : float, optional
            The tail probability the y order is checked at; default the Numerics default.

        Returns
        -------
        list of str
            The messages (empty when valid).
        """
        if other is None:
            return [str(m) for m in self.net.OrdinateErrors()]
        return [str(m) for m in self.net.OrdinateErrors(*self._comparison(
            other, strict_x, strict_y, x_order, y_order, other_is_next, allow_different_types, tail_probability))]

    def to_xml(self) -> str:
        """Return the uncertain ordinate's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())


class UncertainOrderedPairedData:
    """Ordered x values with a distribution of y at each (``Numerics.Data.UncertainOrderedPairedData``).

    Parameters
    ----------
    x : array_like
        The x values.
    distributions : sequence of UnivariateDistribution
        One distribution of y per x value.
    strict_x : bool, optional
        Whether x must be strictly monotone; default True.
    x_order : SortOrder or str, optional
        Order of x; default ascending.
    strict_y : bool, optional
        Whether the y quantiles must be strictly monotone; default True.
    y_order : SortOrder or str, optional
        Order of y; default ascending.
    distribution_type : str, optional
        The ``UnivariateDistributionType`` of the y distributions; default the type of the first.

    Attributes
    ----------
    net : Numerics.Data.UncertainOrderedPairedData
        The wrapped .NET object.
    """

    def __init__(self, x: Any, distributions: Sequence[Any], strict_x: bool = True,
                 x_order: SortOrder | str = SortOrder.ASCENDING, strict_y: bool = True,
                 y_order: SortOrder | str = SortOrder.ASCENDING, distribution_type: str | None = None) -> None:
        data = _data()
        nd = namespace("Numerics.Distributions")
        import System  # noqa: PLC0415 - the runtime is running

        nets = [_net_distribution(d) for d in distributions]
        kind = nets[0].Type if distribution_type is None else net_enum(nd.UnivariateDistributionType,
                                                                       distribution_type)
        self.net = data.UncertainOrderedPairedData(
            to_net_doubles(x), System.Array[nd.UnivariateDistributionBase](nets), bool(strict_x),
            net_enum(data.SortOrder, x_order), bool(strict_y), net_enum(data.SortOrder, y_order), kind)

    @classmethod
    def _from_net(cls, net: Any) -> UncertainOrderedPairedData:
        """Wrap an existing .NET object."""
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def from_xml(cls, xml: str) -> UncertainOrderedPairedData:
        """Restore the paired data from its XML (the ``XElement`` constructor).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        UncertainOrderedPairedData
            The paired data.
        """
        return cls._from_net(_data().UncertainOrderedPairedData(parse_xml(xml)))

    def to_xml(self) -> str:
        """Return the paired data's ``SaveToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.SaveToXElement().ToString())

    def __len__(self) -> int:
        return int(self.net.Count)

    def __iter__(self) -> Iterator[UncertainOrdinate]:
        for net in self.net:
            yield UncertainOrdinate._from_net(net)

    def __getitem__(self, index: int) -> UncertainOrdinate:
        return UncertainOrdinate._from_net(self.net[int(index)])

    def __setitem__(self, index: int, ordinate: UncertainOrdinate) -> None:
        self.net[int(index)] = ordinate.net

    def __contains__(self, ordinate: UncertainOrdinate) -> bool:
        return bool(self.net.Contains(ordinate.net))

    @property
    def count(self) -> int:
        """The number of ordinates."""
        return int(self.net.Count)

    @property
    def x(self) -> np.ndarray:
        """The x values."""
        return np.array([float(o.X) for o in self.net], dtype=np.float64)

    @property
    def distributions(self) -> list[Any]:
        """The distribution of y at each x value."""
        return [_wrap_distribution(o.Y) for o in self.net]

    @property
    def distribution(self) -> str:
        """The ``UnivariateDistributionType`` name of the y distributions."""
        return str(self.net.Distribution)

    @property
    def is_valid(self) -> bool:
        """Whether the ordinates satisfy the configured orders."""
        return bool(self.net.IsValid)

    @property
    def errors(self) -> list[str]:
        """The validation errors (``GetErrors``)."""
        return [str(message) for message in self.net.GetErrors()]

    @property
    def strict_x(self) -> bool:
        """Whether x must be strictly monotone."""
        return bool(self.net.StrictX)

    @strict_x.setter
    def strict_x(self, value: bool) -> None:
        self.net.StrictX = bool(value)

    @property
    def strict_y(self) -> bool:
        """Whether the y quantiles must be strictly monotone."""
        return bool(self.net.StrictY)

    @strict_y.setter
    def strict_y(self, value: bool) -> None:
        self.net.StrictY = bool(value)

    @property
    def order_x(self) -> SortOrder:
        """The required order of x."""
        return _sort_order(self.net.OrderX)

    @order_x.setter
    def order_x(self, value: SortOrder | str) -> None:
        self.net.OrderX = net_enum(_data().SortOrder, value)

    @property
    def order_y(self) -> SortOrder:
        """The required order of y."""
        return _sort_order(self.net.OrderY)

    @order_y.setter
    def order_y(self, value: SortOrder | str) -> None:
        self.net.OrderY = net_enum(_data().SortOrder, value)

    @property
    def allow_different_distribution_types(self) -> bool:
        """Whether the y distributions may be of different families."""
        return bool(self.net.AllowDifferentDistributionTypes)

    @allow_different_distribution_types.setter
    def allow_different_distribution_types(self, value: bool) -> None:
        self.net.AllowDifferentDistributionTypes = bool(value)

    @property
    def validation_tail_probability(self) -> float:
        """The tail probability the y order is validated at (NaN: the Numerics default)."""
        return float(self.net.ValidationTailProbability)

    @validation_tail_probability.setter
    def validation_tail_probability(self, value: float) -> None:
        self.net.ValidationTailProbability = float(value)

    def validate(self) -> None:
        """Revalidate every ordinate (``Validate``)."""
        self.net.Validate()

    def curve_sample(self, probability: float | None = None) -> OrderedPairedData:
        """Return the curve at a quantile of every distribution, or at their means (``CurveSample``).

        Parameters
        ----------
        probability : float, optional
            The non-exceedance probability, clipped to [0, 1]; default the means.

        Returns
        -------
        OrderedPairedData
            The sampled curve with this object's orders.
        """
        net = self.net.CurveSample() if probability is None else self.net.CurveSample(float(probability))
        return OrderedPairedData._from_net(net)

    def copy(self) -> UncertainOrderedPairedData:
        """Return an independent copy (``Clone``).

        Returns
        -------
        UncertainOrderedPairedData
            The copy.
        """
        return UncertainOrderedPairedData._from_net(self.net.Clone())

    def append(self, x: float, distribution: Any) -> None:
        """Append an uncertain ordinate (``Add``).

        Parameters
        ----------
        x : float
            The x value.
        distribution : UnivariateDistribution
            The distribution of y.
        """
        self.net.Add(_data().UncertainOrdinate(float(x), _net_distribution(distribution)))

    def _ordinates(self, x: Any, distributions: Sequence[Any]) -> Any:
        """Return a ``UncertainOrdinate[]`` of paired values and distributions."""
        data = _data()
        import System  # noqa: PLC0415 - the runtime is running

        values = np.asarray(x, dtype=np.float64).ravel().tolist()
        return System.Array[data.UncertainOrdinate]([
            data.UncertainOrdinate(v, _net_distribution(d)) for v, d in zip(values, distributions, strict=True)])

    def extend(self, x: Any, distributions: Sequence[Any]) -> None:
        """Append several uncertain ordinates (``AddRange``).

        Parameters
        ----------
        x : array_like
            The x values.
        distributions : sequence of UnivariateDistribution
            One distribution of y per x value.
        """
        self.net.AddRange(self._ordinates(x, distributions))

    def insert(self, index: int, x: float, distribution: Any) -> None:
        """Insert an uncertain ordinate before a position (``Insert``).

        Parameters
        ----------
        index : int
            The position.
        x : float
            The x value.
        distribution : UnivariateDistribution
            The distribution of y.
        """
        self.net.Insert(int(index), _data().UncertainOrdinate(float(x), _net_distribution(distribution)))

    def insert_range(self, index: int, x: Any, distributions: Sequence[Any]) -> None:
        """Insert several uncertain ordinates before a position (``InsertRange``).

        Parameters
        ----------
        index : int
            The position.
        x : array_like
            The x values.
        distributions : sequence of UnivariateDistribution
            One distribution of y per x value.
        """
        self.net.InsertRange(int(index), self._ordinates(x, distributions))

    def index_of(self, ordinate: UncertainOrdinate) -> int:
        """Return the position of an uncertain ordinate, or -1 (``IndexOf``).

        Parameters
        ----------
        ordinate : UncertainOrdinate
            The ordinate (equal x and an equal distribution).

        Returns
        -------
        int
            The first position holding it, or -1.
        """
        return int(self.net.IndexOf(ordinate.net))

    def remove(self, ordinate: UncertainOrdinate) -> bool:
        """Remove the first equal uncertain ordinate (``Remove``).

        Parameters
        ----------
        ordinate : UncertainOrdinate
            The ordinate.

        Returns
        -------
        bool
            Whether one was removed.
        """
        return bool(self.net.Remove(ordinate.net))

    def remove_at(self, index: int) -> None:
        """Remove the uncertain ordinate at a position (``RemoveAt``).

        Parameters
        ----------
        index : int
            The position.
        """
        self.net.RemoveAt(int(index))

    def remove_range(self, index: Any, count: int | None = None) -> None:
        """Remove a run of uncertain ordinates, or the ordinates at several rows (``RemoveRange``).

        Parameters
        ----------
        index : int or sequence of int
            The first position of the run, or the rows to remove (any order).
        count : int, optional
            How many to remove from ``index``; omit when ``index`` lists rows.
        """
        if count is None:
            self.net.RemoveRange(to_net_ints(index))
        else:
            self.net.RemoveRange(int(index), int(count))

    def clear(self) -> None:
        """Remove every uncertain ordinate (``Clear``)."""
        self.net.Clear()


class ProbabilityOrdinates:
    """A list of exceedance probabilities for frequency curves (``Numerics.Data.ProbabilityOrdinates``).

    Parameters
    ----------
    probabilities : array_like, optional
        The probabilities; default the 25 standard ordinates from 1e-6 to 0.99.

    Attributes
    ----------
    net : Numerics.Data.ProbabilityOrdinates
        The wrapped .NET object.
    """

    DEFAULT_DELIMITER = "|"
    """The delimiter of :meth:`to_delimited_string` and :meth:`parse` (``DefaultDelimiter``)."""

    def __init__(self, probabilities: Any = None) -> None:
        data = _data()
        self.net = data.ProbabilityOrdinates() if probabilities is None else \
            data.ProbabilityOrdinates(to_net_doubles(probabilities))

    @classmethod
    def parse(cls, value: str, delimiter: str = DEFAULT_DELIMITER) -> ProbabilityOrdinates:
        """Read probabilities from delimited text (``Parse``); unreadable tokens are skipped.

        Parameters
        ----------
        value : str
            The text, for example ``"0.01|0.1|0.5"``.
        delimiter : str, optional
            The separator; default ``"|"``.

        Returns
        -------
        ProbabilityOrdinates
            The probabilities.
        """
        instance = cls.__new__(cls)
        instance.net = _data().ProbabilityOrdinates.Parse(str(value), str(delimiter))
        return instance

    def __len__(self) -> int:
        return int(self.net.Count)

    def __getitem__(self, index: int) -> float:
        return float(self.net[int(index)])

    def __setitem__(self, index: int, value: float) -> None:
        self.net[int(index)] = float(value)

    @property
    def values(self) -> np.ndarray:
        """The probabilities, copied in one transfer."""
        return to_numpy_doubles(self.net)

    def to_delimited_string(self, delimiter: str = DEFAULT_DELIMITER) -> str:
        """Return the probabilities as delimited text in round-trip precision (``ToDelimitedString``).

        Parameters
        ----------
        delimiter : str, optional
            The separator; default ``"|"``.

        Returns
        -------
        str
            The text.
        """
        return str(self.net.ToDelimitedString(str(delimiter)))

    def from_delimited_string(self, value: str, delimiter: str = DEFAULT_DELIMITER) -> None:
        """Replace the probabilities with those read from delimited text (``FromDelimitedString``).

        Parameters
        ----------
        value : str
            The text.
        delimiter : str, optional
            The separator; default ``"|"``.
        """
        self.net.FromDelimitedString(str(value), str(delimiter))

    def reset_to_defaults(self) -> None:
        """Replace the probabilities with the 25 standard ordinates (``ResetToDefaults``)."""
        self.net.ResetToDefaults()

    def validate(self) -> tuple[bool, list[str]]:
        """Check that the list is non-empty, in (0, 1) and strictly increasing (``Validate``).

        Returns
        -------
        tuple of (bool, list of str)
            Whether the list is valid, and the messages.
        """
        result = self.net.Validate()
        return bool(result.Item1), [str(m) for m in result.Item2]

    def append(self, value: float) -> None:
        """Append a probability (``Add``).

        Parameters
        ----------
        value : float
            The probability.
        """
        self.net.Add(float(value))

    def extend(self, values: Any) -> None:
        """Append several probabilities (``AddRange``).

        Parameters
        ----------
        values : array_like
            The probabilities.
        """
        self.net.AddRange(to_net_doubles(values))

    def insert(self, index: int, value: float) -> None:
        """Insert a probability before a position (``Insert``).

        Parameters
        ----------
        index : int
            The position.
        value : float
            The probability.
        """
        self.net.Insert(int(index), float(value))

    def insert_range(self, index: int, values: Any) -> None:
        """Insert several probabilities before a position (``InsertRange``).

        Parameters
        ----------
        index : int
            The position.
        values : array_like
            The probabilities.
        """
        self.net.InsertRange(int(index), to_net_doubles(values))

    def remove(self, value: float) -> bool:
        """Remove the first equal probability (``Remove``).

        Parameters
        ----------
        value : float
            The probability.

        Returns
        -------
        bool
            Whether one was removed.
        """
        return bool(self.net.Remove(float(value)))

    def remove_at(self, index: int) -> None:
        """Remove the probability at a position (``RemoveAt``).

        Parameters
        ----------
        index : int
            The position.
        """
        self.net.RemoveAt(int(index))

    def remove_range(self, index: int, count: int) -> None:
        """Remove a run of probabilities (``RemoveRange``).

        Parameters
        ----------
        index : int
            The first position.
        count : int
            How many to remove.
        """
        self.net.RemoveRange(int(index), int(count))

    def clear(self) -> None:
        """Remove every probability (``Clear``)."""
        self.net.Clear()


def _ordinate_list(x: Any, y: Any) -> Any:
    """Return a ``List<Ordinate>`` of paired arrays, built in .NET from one unvalidated paired data object."""
    data = _data()
    from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

    unordered = getattr(data.SortOrder, "None")
    pairs = data.OrderedPairedData(to_net_doubles(x), to_net_doubles(y), False, unordered, False, unordered)
    return List[data.Ordinate](pairs)


def ramer_douglas_peucker(x: Any, y: Any, epsilon: float) -> tuple[np.ndarray, np.ndarray]:
    """Simplify a polyline by the Ramer-Douglas-Peucker algorithm (``LineSimplification.RamerDouglasPeucker``).

    The points cross once, as two arrays; the kept points come back one ordinate at a time.

    Parameters
    ----------
    x, y : array_like
        The polyline vertices, in order (at least two).
    epsilon : float
        The largest perpendicular distance a removed vertex may have from the simplified line.

    Returns
    -------
    tuple of numpy.ndarray
        The x and y of the kept vertices, end points included.
    """
    data = _data()
    from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

    output = List[data.Ordinate]()
    result = data.LineSimplification.RamerDouglasPeucker(_ordinate_list(x, y), float(epsilon), output)
    kept = result[1] if isinstance(result, tuple) else output
    return (np.array([float(o.X) for o in kept], dtype=np.float64),
            np.array([float(o.Y) for o in kept], dtype=np.float64))


def perpendicular_distance(point: Any, line_start: Any, line_end: Any) -> float:
    """Return the distance of a point from the line through two points (``LineSimplification.PerpendicularDistance``).

    Parameters
    ----------
    point, line_start, line_end : Ordinate or tuple of float
        The point and two points on the line, as ordinates or ``(x, y)`` pairs.

    Returns
    -------
    float
        The perpendicular distance (the distance to ``line_start`` when the two line points coincide).
    """
    data = _data()

    def net(p: Any) -> Any:
        return p.net if isinstance(p, Ordinate) else data.Ordinate(float(p[0]), float(p[1]))

    return float(data.LineSimplification.PerpendicularDistance(net(point), net(line_start), net(line_end)))


def time_series_index_ticks(series: Any) -> np.ndarray:
    """Return the dates of a time series as ``DateTime.Ticks`` (``TimeSeries.IndexTicksToArray``).

    One ``long[]`` crosses instead of one ``DateTime`` per ordinate. A tick is 100 nanoseconds
    since 0001-01-01 00:00; the ``DateTime.Kind`` is not included. For example,
    ``(ticks - 621355968000000000) // 10`` are microseconds since 1970-01-01, which
    ``.astype("datetime64[us]")`` turns into NumPy dates.

    Parameters
    ----------
    series : Numerics.Data.TimeSeries
        The .NET time series (or an object holding it as ``.net``).

    Returns
    -------
    numpy.ndarray
        The int64 ticks of each ordinate date, in collection order.
    """
    net = getattr(series, "net", series)
    return to_numpy_longs(net.IndexTicksToArray())


def time_series_values(series: Any) -> np.ndarray:
    """Return the values of a time series (``Series.ValuesToArray``).

    Parameters
    ----------
    series : Numerics.Data.TimeSeries
        The .NET time series (or an object holding it as ``.net``).

    Returns
    -------
    numpy.ndarray
        The float64 value of each ordinate, in collection order, aligned with
        :func:`time_series_index_ticks`; missing values are NaN.
    """
    return to_numpy_doubles(getattr(series, "net", series).ValuesToArray())


__all__ = [
    "Bilinear", "CubicSpline", "Interpolater", "Linear", "LinearRegression", "Ordinate", "OrderedPairedData",
    "Polynomial", "ProbabilityOrdinates", "SearchMethod", "UncertainOrderedPairedData", "UncertainOrdinate",
    "parse_xml", "perpendicular_distance", "ramer_douglas_peucker", "search", "time_series_index_ticks",
    "time_series_values",
]

translate_public(globals())
