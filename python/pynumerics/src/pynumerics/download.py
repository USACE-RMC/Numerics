"""Downloads of public hydrologic and climate records (``Numerics.Data.TimeSeriesDownload``).

The four providers are the U.S. Geological Survey Water Data API (:func:`from_usgs`), NOAA NCEI
GHCN-Daily station files (:func:`from_ghcn`), the Water Survey of Canada (:func:`from_chmn`) and the
Australian Bureau of Meteorology Water Data Online service (:func:`from_abom`). Each downloader
fetches the provider's response and parses it with the matching parser, which is also public
(:func:`parse_usgs`, :func:`parse_ghcn`, :func:`parse_chmn`, :func:`parse_abom`), so a saved response
gives the same series offline.

Every function returns a :class:`DownloadedSeries`: the :class:`~pynumerics.timeseries.TimeSeries`,
the raw text (USGS downloads only), the optional qualifiers, one per value, and whether the dates
are UTC. Requests, retries, units, gap filling and parsing all run in Numerics; nothing is parsed in
Python. The dates and values cross in two bulk copies, and the qualifiers in one.

``timeout`` cancels a download after that many seconds (``CancellationTokenSource.CancelAfter``) and
raises :class:`TimeoutError`. A provider timeout also raises :class:`TimeoutError`, a network failure
:class:`ConnectionError`, an invalid argument :class:`ValueError`, and a provider error or "No data
found" :class:`RuntimeError`, each with the Numerics message.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from pynumerics._util import namespace, net_enum
from pynumerics.enums import NetEnum
from pynumerics.errors import net_errors, translate_public
from pynumerics.timeseries import TimeSeries, _net_date

_SEPARATOR = "\x1f"
"""Joins the qualifiers into one string for the transfer; no provider qualifier contains it."""


class TimeSeriesType(NetEnum):
    """The kinds of series a provider offers (``TimeSeriesDownload.TimeSeriesType``)."""

    DAILY_DISCHARGE = "DailyDischarge"
    DAILY_STAGE = "DailyStage"
    DAILY_PRECIPITATION = "DailyPrecipitation"
    DAILY_SNOW = "DailySnow"
    INSTANTANEOUS_DISCHARGE = "InstantaneousDischarge"
    INSTANTANEOUS_STAGE = "InstantaneousStage"
    PEAK_DISCHARGE = "PeakDischarge"
    PEAK_STAGE = "PeakStage"
    MEASURED_DISCHARGE = "MeasuredDischarge"
    MEASURED_STAGE = "MeasuredStage"


class DepthUnit(NetEnum):
    """Depth units for precipitation and snowfall (``TimeSeriesDownload.DepthUnit``)."""

    MILLIMETERS = "Millimeters"
    CENTIMETERS = "Centimeters"
    INCHES = "Inches"


class DischargeUnit(NetEnum):
    """Discharge units (``TimeSeriesDownload.DischargeUnit``)."""

    CUBIC_FEET_PER_SECOND = "CubicFeetPerSecond"
    CUBIC_METERS_PER_SECOND = "CubicMetersPerSecond"


class HeightUnit(NetEnum):
    """Stage (gage height) units (``TimeSeriesDownload.HeightUnit``)."""

    FEET = "Feet"
    METERS = "Meters"


class UsgsStatistic(NetEnum):
    """The daily statistic of a USGS daily series (``TimeSeriesDownload.UsgsStatistic``)."""

    MEAN = "Mean"
    MAXIMUM = "Maximum"
    MINIMUM = "Minimum"


@dataclass(frozen=True)
class DownloadedSeries:
    """A downloaded or parsed series with its raw text and qualifiers.

    Parameters
    ----------
    series : TimeSeries
        The series, wrapping the ``Numerics.Data.TimeSeries`` Numerics returned.
    raw_text : str or None
        The text the service returned, after ``#`` comment lines naming the source, each request and the
        retrieval time (USGS downloads only); None otherwise.
    qualifiers : tuple of str or None
        With ``qualifiers=True``: one entry per value, in series order (USGS qualifier field, GHCN day
        flags, Canadian data symbol, BOM quality code; empty for a day filled as missing). None otherwise.
    utc : bool
        True when the dates are UTC (USGS instantaneous values and field measurements, Canadian
        real-time values); otherwise the dates are calendar dates or the provider's local time.
    """

    series: TimeSeries
    raw_text: str | None
    qualifiers: tuple[str, ...] | None
    utc: bool

    @property
    def dates(self) -> np.ndarray:
        """The dates as ``datetime64[us]`` (see :attr:`utc` for their time zone)."""
        return self.series.dates

    @property
    def values(self) -> np.ndarray:
        """The values; NaN marks a missing value."""
        return self.series.values


def _download() -> Any:
    """Return ``Numerics.Data.TimeSeriesDownload`` after starting the runtime."""
    return namespace("Numerics.Data").TimeSeriesDownload


def _enum(name: str, value: Any) -> Any:
    """Return a member of a ``TimeSeriesDownload`` nested enum."""
    return net_enum(getattr(_download(), name), value)


def _date(value: Any) -> Any:
    """Return a date as a ``System.DateTime``, or None."""
    return None if value is None else _net_date(value)


def _qualifier_list(qualifiers: bool) -> Any:
    """Return a new .NET ``List<string>`` for the qualifier overloads, or None."""
    if not qualifiers:
        return None
    from System import String  # noqa: PLC0415 - the runtime is running
    from System.Collections.Generic import List  # noqa: PLC0415

    return List[String]()


def _result(net: Any, raw_text: Any = None, qualifiers: Any = None) -> DownloadedSeries:
    """Wrap a Numerics series, copying the qualifiers in one transfer."""
    import System  # noqa: PLC0415 - the runtime is running

    codes = None
    if qualifiers is not None:
        joined = str(System.String.Join(_SEPARATOR, qualifiers))
        codes = tuple(joined.split(_SEPARATOR)) if int(qualifiers.Count) > 0 else ()
    utc = int(net.Count) > 0 and net[0].Index.Kind == System.DateTimeKind.Utc
    return DownloadedSeries(TimeSeries.from_net(net), None if raw_text is None else str(raw_text), codes, bool(utc))


def _run(start: Any, timeout: float | None) -> Any:
    """Run a download task to completion and return its result.

    ``start`` takes a ``CancellationToken`` and returns the task. With a ``timeout`` the token is
    cancelled after that many seconds. The task is awaited with ``GetAwaiter().GetResult()``, so a
    failure surfaces as the .NET exception itself, not an ``AggregateException``; the wait polls, so
    Ctrl+C cancels the download.
    """
    import System  # noqa: PLC0415 - the runtime is running
    from System.Net.Http import HttpRequestException  # noqa: PLC0415
    from System.Threading import CancellationTokenSource  # noqa: PLC0415

    if timeout is not None and not timeout > 0:
        raise ValueError("timeout must be positive.")
    source = CancellationTokenSource()
    try:
        if timeout is not None:
            source.CancelAfter(System.TimeSpan.FromSeconds(float(timeout)))
        task = start(source.Token)
        try:
            while not task.IsCompleted:
                time.sleep(0.01)
        except KeyboardInterrupt:
            source.Cancel()
            raise
        try:
            return task.GetAwaiter().GetResult()
        except System.OperationCanceledException as error:
            if source.IsCancellationRequested:
                raise TimeoutError(f"The download did not complete within {timeout} seconds.") from error
            raise
        except System.TimeoutException as error:
            raise TimeoutError(str(error.Message)) from error
        except (HttpRequestException, System.IO.IOException) as error:
            raise ConnectionError(str(error.Message)) from error
    finally:
        source.Dispose()


def from_usgs(site: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE, *, start: Any = None,
              end: Any = None, statistic: UsgsStatistic | str = UsgsStatistic.MEAN, api_key: str | None = None,
              qualifiers: bool = False, timeout: float | None = None) -> DownloadedSeries:
    """Download a series from the USGS Water Data API (``TimeSeriesDownload.FromUSGS``).

    Parameters
    ----------
    site : str
        The USGS site number, 8 to 15 digits (for example ``"01646500"``).
    series_type : TimeSeriesType or str, optional
        Daily, instantaneous, peak or measured discharge or stage; default daily discharge.
    start, end : datetime-like, optional
        The inclusive date range; default the period of record. Calendar dates for daily values and
        peaks, UTC for instantaneous values and field measurements.
    statistic : UsgsStatistic or str, optional
        The daily statistic for daily types; default the mean.
    api_key : str, optional
        A USGS Water Data API key for this download; default :func:`get_usgs_api_key`.
    qualifiers : bool, optional
        True to return the qualifier field of each value as well.
    timeout : float, optional
        Seconds before the download is cancelled; default no limit beyond the per-request limits.

    Returns
    -------
    DownloadedSeries
        The series and the raw text. Daily values form a regular one-day series with NaN for missing
        days; instantaneous values and field measurements are irregular and UTC; peaks are irregular, and
        a peak with a time of day is dated in the site's local standard time.

    Raises
    ------
    ValueError
        For a site number that is not 8 to 15 digits, an unsupported type or ``start`` after ``end``.
    TimeoutError
        When ``timeout`` or a request's time limit is exceeded.
    ConnectionError
        When the service cannot be reached.
    RuntimeError
        When the service returns an error or no data (``"No data found ..."``).
    """
    net = _download()
    with net_errors():
        options = net.UsgsDownloadOptions()
        options.StartDate = _date(start)
        options.EndDate = _date(end)
        options.Statistic = _enum("UsgsStatistic", statistic)
        options.ApiKey = api_key
        kind = _enum("TimeSeriesType", series_type)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            result = _run(lambda token: net.FromUSGS(site, kind, options, token), timeout)
        else:
            result = _run(lambda token: net.FromUSGS(site, kind, options, codes, token), timeout)
        return _result(result.Item1, result.Item2, codes)


def from_ghcn(site: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_PRECIPITATION,
              unit: DepthUnit | str = DepthUnit.INCHES, *, qualifiers: bool = False,
              timeout: float | None = None) -> DownloadedSeries:
    """Download daily precipitation or snowfall from a GHCN-Daily station file (``TimeSeriesDownload.FromGHCN``).

    Parameters
    ----------
    site : str
        The 11-character GHCN station code (for example ``"USW00094728"``).
    series_type : TimeSeriesType or str, optional
        Daily precipitation (default) or daily snow.
    unit : DepthUnit or str, optional
        The depth unit; default inches.
    qualifiers : bool, optional
        True to return the three flag characters of each day as well.
    timeout : float, optional
        Seconds before the download is cancelled.

    Returns
    -------
    DownloadedSeries
        A regular one-day series: missing days are NaN and days after today are left out.

    Raises
    ------
    ValueError
        For an invalid station code or type.
    TimeoutError
        When ``timeout`` or the request's time limit is exceeded.
    ConnectionError
        When the service cannot be reached or the station file does not exist.
    RuntimeError
        When the service returns an empty file, an HTML page or a rate-limit response.
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        depth = _enum("DepthUnit", unit)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = _run(lambda token: net.FromGHCN(site, kind, depth, token), timeout)
        else:
            series = _run(lambda token: net.FromGHCN(site, kind, depth, codes, token), timeout)
        return _result(series, None, codes)


def from_chmn(station: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE,
              discharge_unit: DischargeUnit | str = DischargeUnit.CUBIC_METERS_PER_SECOND,
              height_unit: HeightUnit | str = HeightUnit.METERS, *, start: Any = None, end: Any = None,
              qualifiers: bool = False, timeout: float | None = None) -> DownloadedSeries:
    """Download a series from the Water Survey of Canada (``TimeSeriesDownload.FromCHMN``).

    Parameters
    ----------
    station : str
        The 7-character station id (for example ``"08MF005"``).
    series_type : TimeSeriesType or str, optional
        Daily, instantaneous (real-time) or peak discharge or stage; default daily discharge.
    discharge_unit : DischargeUnit or str, optional
        The discharge unit; default cubic meters per second.
    height_unit : HeightUnit or str, optional
        The stage unit; default meters.
    start, end : datetime-like, optional
        The inclusive date range; default the period of record (the last 18 months for real-time data).
    qualifiers : bool, optional
        True to return the data symbol of each value as well (for example ``"B"`` for ice conditions).
    timeout : float, optional
        Seconds before the download is cancelled.

    Returns
    -------
    DownloadedSeries
        The series; real-time dates are UTC.

    Raises
    ------
    ValueError
        For a station id that is not 7 characters or an unsupported type.
    TimeoutError
        When ``timeout`` or the request's time limit is exceeded.
    ConnectionError
        When the service cannot be reached.
    RuntimeError
        When the service returns an error, or an empty, header-only or HTML response.
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        flow = _enum("DischargeUnit", discharge_unit)
        height = _enum("HeightUnit", height_unit)
        first, last = _date(start), _date(end)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = _run(lambda token: net.FromCHMN(station, kind, flow, height, first, last, token), timeout)
        else:
            series = _run(lambda token: net.FromCHMN(station, kind, flow, height, first, last, codes, token), timeout)
        return _result(series, None, codes)


def from_abom(station: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE,
              discharge_unit: DischargeUnit | str = DischargeUnit.CUBIC_METERS_PER_SECOND,
              height_unit: HeightUnit | str = HeightUnit.METERS, depth_unit: DepthUnit | str = DepthUnit.MILLIMETERS, *,
              start: Any = None, end: Any = None, qualifiers: bool = False,
              timeout: float | None = None) -> DownloadedSeries:
    """Download a series from the BOM Water Data Online service (``TimeSeriesDownload.FromABOM``).

    Parameters
    ----------
    station : str
        The BOM station number (for example ``"410730"``).
    series_type : TimeSeriesType or str, optional
        Daily or instantaneous discharge or stage, or daily precipitation; default daily discharge.
    discharge_unit : DischargeUnit or str, optional
        The discharge unit; default cubic meters per second.
    height_unit : HeightUnit or str, optional
        The stage unit; default meters.
    depth_unit : DepthUnit or str, optional
        The rainfall unit; default millimeters.
    start, end : datetime-like, optional
        The date range; default the period of record to today.
    qualifiers : bool, optional
        True to return the KiWIS quality code of each value as well.
    timeout : float, optional
        Seconds before the download is cancelled.

    Returns
    -------
    DownloadedSeries
        The series; daily types use the midnight-to-midnight (``.24HR``) series.

    Raises
    ------
    ValueError
        For an invalid station number or an unsupported type.
    TimeoutError
        When ``timeout`` is exceeded.
    RuntimeError
        When the service returns an error or has no matching series.
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        flow = _enum("DischargeUnit", discharge_unit)
        height = _enum("HeightUnit", height_unit)
        depth = _enum("DepthUnit", depth_unit)
        first, last = _date(start), _date(end)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = _run(lambda token: net.FromABOM(station, kind, flow, height, depth, first, last, token), timeout)
        else:
            series = _run(lambda token: net.FromABOM(station, kind, flow, height, depth, first, last, codes, token),
                          timeout)
        return _result(series, None, codes)


def parse_usgs(text: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE,
               peak_standard_offset_hours: float | None = None, *, qualifiers: bool = False) -> DownloadedSeries:
    """Parse USGS Water Data API CSV text offline (``TimeSeriesDownload.ParseUSGS``).

    Parameters
    ----------
    text : str
        The raw text :func:`from_usgs` returns, or a Water Data API CSV response; ``#`` lines are skipped.
    series_type : TimeSeriesType or str, optional
        The type the text holds; default daily discharge.
    peak_standard_offset_hours : float, optional
        For peaks with a time of day: the standard UTC offset of the site's time zone in hours (for
        example -5 for EST).
    qualifiers : bool, optional
        True to return the qualifier field of each value as well.

    Returns
    -------
    DownloadedSeries
        The series :func:`from_usgs` would return for the same text.

    Raises
    ------
    ValueError
        For precipitation or snow types, or a timed peak without ``peak_standard_offset_hours``.
    """
    net = _download()
    import System  # noqa: PLC0415 - the runtime is running

    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        hours = peak_standard_offset_hours
        offset = None if hours is None else System.TimeSpan.FromHours(float(hours))
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = net.ParseUSGS(text, kind, offset)
        else:
            series = net.ParseUSGS(text, kind, offset, codes)
        return _result(series, None, codes)


def parse_ghcn(text: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_PRECIPITATION,
               unit: DepthUnit | str = DepthUnit.INCHES, *, qualifiers: bool = False) -> DownloadedSeries:
    """Parse a GHCN-Daily station file offline (``TimeSeriesDownload.ParseGHCN``).

    Parameters
    ----------
    text : str
        The ``.dly`` station file text.
    series_type : TimeSeriesType or str, optional
        Daily precipitation (default) or daily snow.
    unit : DepthUnit or str, optional
        The depth unit; default inches. Precipitation is stored in tenths of millimeters and snowfall in
        millimeters.
    qualifiers : bool, optional
        True to return the three flag characters of each day as well.

    Returns
    -------
    DownloadedSeries
        A regular one-day series with NaN for missing days; days after today are left out.

    Raises
    ------
    ValueError
        For a type other than daily precipitation or snow.
    RuntimeError
        For a record too short to parse (a file cut off mid-line).
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        depth = _enum("DepthUnit", unit)
        codes = _qualifier_list(qualifiers)
        series = net.ParseGHCN(text, kind, depth) if codes is None else net.ParseGHCN(text, kind, depth, codes)
        return _result(series, None, codes)


def parse_chmn(text: str, station: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE,
               discharge_unit: DischargeUnit | str = DischargeUnit.CUBIC_METERS_PER_SECOND,
               height_unit: HeightUnit | str = HeightUnit.METERS, *, qualifiers: bool = False) -> DownloadedSeries:
    """Parse a Water Survey of Canada CSV response offline (``TimeSeriesDownload.ParseCHMN``).

    Parameters
    ----------
    text : str
        The CSV text of a daily_data, real_time_data or peak_data response.
    station : str
        The 7-character station id; rows of other stations are skipped.
    series_type : TimeSeriesType or str, optional
        Daily, instantaneous (real-time) or peak discharge or stage; default daily discharge.
    discharge_unit : DischargeUnit or str, optional
        The discharge unit; default cubic meters per second.
    height_unit : HeightUnit or str, optional
        The stage unit; default meters.
    qualifiers : bool, optional
        True to return the data symbol of each value as well.

    Returns
    -------
    DownloadedSeries
        The series :func:`from_chmn` would return for the same text.

    Raises
    ------
    ValueError
        For a station id that is not 7 characters or an unsupported type.
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        flow = _enum("DischargeUnit", discharge_unit)
        height = _enum("HeightUnit", height_unit)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = net.ParseCHMN(text, station, kind, flow, height)
        else:
            series = net.ParseCHMN(text, station, kind, flow, height, codes)
        return _result(series, None, codes)


def parse_abom(text: str, series_type: TimeSeriesType | str = TimeSeriesType.DAILY_DISCHARGE,
               discharge_unit: DischargeUnit | str = DischargeUnit.CUBIC_METERS_PER_SECOND,
               height_unit: HeightUnit | str = HeightUnit.METERS, depth_unit: DepthUnit | str = DepthUnit.MILLIMETERS,
               *, qualifiers: bool = False) -> DownloadedSeries:
    """Parse a BOM KiWIS getTimeseriesValues JSON response offline (``TimeSeriesDownload.ParseABOM``).

    Parameters
    ----------
    text : str
        The JSON response.
    series_type : TimeSeriesType or str, optional
        Daily or instantaneous discharge or stage, or daily precipitation; default daily discharge.
    discharge_unit : DischargeUnit or str, optional
        The discharge unit; default cubic meters per second.
    height_unit : HeightUnit or str, optional
        The stage unit; default meters.
    depth_unit : DepthUnit or str, optional
        The rainfall unit; default millimeters.
    qualifiers : bool, optional
        True to return the quality code of each value as well (empty when the response has no Quality
        Code column).

    Returns
    -------
    DownloadedSeries
        The series :func:`from_abom` would return for the same response.

    Raises
    ------
    ValueError
        For an unsupported type.
    RuntimeError
        For a KiWIS error response or malformed JSON.
    """
    net = _download()
    with net_errors():
        kind = _enum("TimeSeriesType", series_type)
        flow = _enum("DischargeUnit", discharge_unit)
        height = _enum("HeightUnit", height_unit)
        depth = _enum("DepthUnit", depth_unit)
        codes = _qualifier_list(qualifiers)
        if codes is None:
            series = net.ParseABOM(text, kind, flow, height, depth)
        else:
            series = net.ParseABOM(text, kind, flow, height, depth, codes)
        return _result(series, None, codes)


def is_connected_to_internet(timeout: float | None = None) -> bool:
    """Return whether at least one provider endpoint answers (``TimeSeriesDownload.IsConnectedToInternet``).

    Parameters
    ----------
    timeout : float, optional
        Seconds to wait; the probe has its own two-second limit per endpoint.

    Returns
    -------
    bool
        True when any provider endpoint completes an HTTP request. Downloads do not need this check.
    """
    net = _download()
    with net_errors():
        return bool(_run(lambda _token: net.IsConnectedToInternet(), timeout))


def get_usgs_api_key() -> str | None:
    """Return the default USGS Water Data API key (``TimeSeriesDownload.UsgsApiKey``).

    Returns
    -------
    str or None
        The key sent in the ``X-Api-Key`` header when a download gives none, or None.
    """
    key = _download().UsgsApiKey
    return None if key is None else str(key)


def set_usgs_api_key(key: str | None) -> None:
    """Set the default USGS Water Data API key (``TimeSeriesDownload.UsgsApiKey``).

    Parameters
    ----------
    key : str or None
        A free key from https://api.waterdata.usgs.gov/signup/, or None for anonymous use. It raises the
        hourly request limit and is never written to the raw text.
    """
    _download().UsgsApiKey = key


__all__ = [
    "DepthUnit", "DischargeUnit", "DownloadedSeries", "HeightUnit", "TimeSeriesType", "UsgsStatistic",
    "from_abom", "from_chmn", "from_ghcn", "from_usgs", "get_usgs_api_key", "is_connected_to_internet",
    "parse_abom", "parse_chmn", "parse_ghcn", "parse_usgs", "set_usgs_api_key",
]

translate_public(globals())
