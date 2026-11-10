"""TimeSeriesDownload wrapper: offline parsing of recorded provider responses, with the maintained clients as oracles.

The fixtures in ``tests/fixtures/download`` are small cuts of real responses (USGS Water Data API, NOAA NCEI
GHCN-Daily, Water Survey of Canada, BOM Water Data Online). The oracles parse the same records with the maintained
Python clients: dataretrieval (USGS, served the recorded GeoJSON from a localhost server), wetterdienst (its GHCN
daily unit factors applied to the NCEI by-station CSV of the same days) and kiwis-pie (BOM, served the recorded JSON).
The USGS peak oracle is the official NWIS peak-flow file. Canada has no maintained Python client for these services,
so the csv module is its oracle. Every value comes from Numerics.

Tests that reach the real services run only with ``PYNUMERICS_NETWORK_TESTS=1`` and make one small request each.
"""

from __future__ import annotations

import csv
import http.server
import io
import json
import os
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from pynumerics import download
from pynumerics.download import DownloadedSeries, TimeSeriesType

FIXTURES = Path(__file__).parent / "fixtures" / "download"
NETWORK = pytest.mark.skipif(os.environ.get("PYNUMERICS_NETWORK_TESTS") != "1",
                             reason="live provider tests run only with PYNUMERICS_NETWORK_TESTS=1")


def _text(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _same(actual: np.ndarray, expected: np.ndarray) -> None:
    """Assert equal values bit for bit, NaN where the oracle is missing."""
    actual, expected = np.asarray(actual, dtype=np.float64), np.asarray(expected, dtype=np.float64)
    assert actual.shape == expected.shape
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
    np.testing.assert_array_equal(actual[~np.isnan(actual)], expected[~np.isnan(expected)])


@pytest.fixture(scope="module")
def usgs_server():
    """A localhost stand-in for the Water Data API that serves the recorded GeoJSON to dataretrieval."""
    routes = {"/schema": _text("dataretrieval-daily-schema.json"),
              "/daily/items": _text("dataretrieval-daily-05082500-ice.json"),
              "/continuous/items": _text("dataretrieval-continuous-01646500.json")}

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *args):  # noqa: D401 - silence the request log
            pass

        def do_GET(self):  # noqa: N802 - http.server API
            path = self.path.split("?")[0]
            body = next((text for key, text in routes.items() if path.endswith(key)), None)
            data = (body or '{"type":"FeatureCollection","features":[],"numberReturned":0,"links":[]}').encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json" if path.endswith("/schema") else "application/geo+json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def _dataretrieval(base_url: str, getter: str, **query):
    dataretrieval = pytest.importorskip("dataretrieval")
    from dataretrieval import waterdata  # noqa: PLC0415
    from dataretrieval.waterdata.configuration import WaterdataConfiguration  # noqa: PLC0415

    with dataretrieval.configure(WaterdataConfiguration(base_url=base_url)):
        frame, _ = getattr(waterdata, getter)(**query)
    return frame


# USGS -----------------------------------------------------------------------------------------


def test_parse_usgs_daily_matches_dataretrieval(usgs_server) -> None:
    # Arrange
    frame = _dataretrieval(usgs_server, "get_daily", monitoring_location_id="USGS-05082500", parameter_code="00060",
                           statistic_id="00003", skip_geometry=True)
    frame = frame.sort_values("time")
    # Act
    result = download.parse_usgs(_text("usgs-daily-05082500-ice.csv"), "daily_discharge", qualifiers=True)
    # Assert
    days = np.asarray(frame["time"], dtype="datetime64[D]").astype("datetime64[us]")
    np.testing.assert_array_equal(result.dates, days)
    _same(result.values, frame["value"].astype(float).to_numpy())
    assert result.series.interval.value == "OneDay" and not result.utc and result.raw_text is None
    pairs = zip(frame["time"], frame["qualifier"], strict=True)
    ice = {str(t)[:10] for t, q in pairs if q is not None and "ICE" in str(q)}
    ours = {str(d)[:10] for d, q in zip(result.dates, result.qualifiers, strict=True) if "ICE" in q}
    assert ours == ice and len(ice) == 50


def test_parse_usgs_instantaneous_matches_dataretrieval_in_utc(usgs_server) -> None:
    # Arrange
    frame = _dataretrieval(usgs_server, "get_continuous", monitoring_location_id="USGS-01646500",
                           parameter_code="00060",
                           time="2024-03-10T06:00:00Z/2024-03-10T12:00:00Z")
    frame = frame.sort_values("time")
    expected = np.asarray(frame["time"].dt.tz_convert("UTC").dt.tz_localize(None), dtype="datetime64[us]")
    # Act
    result = download.parse_usgs(_text("usgs-instantaneous-01646500.csv"), TimeSeriesType.INSTANTANEOUS_DISCHARGE)
    # Assert
    np.testing.assert_array_equal(result.dates, expected)
    _same(result.values, frame["value"].astype(float).to_numpy())
    assert result.utc and result.qualifiers is None


def test_parse_usgs_peaks_match_nwis_peak_file() -> None:
    # Arrange
    rows = list(csv.DictReader(io.StringIO(_text("usgs-peaks-01646500.expected.csv"))))
    dates = np.asarray([r["index"] for r in rows], dtype="datetime64[us]")
    values = np.asarray([np.frombuffer(bytes.fromhex(r["bits"]), dtype=">f8")[0] for r in rows])
    # Act
    result = download.parse_usgs(_text("usgs-peaks-01646500.csv"), "peak_discharge", peak_standard_offset_hours=-5)
    # Assert
    np.testing.assert_array_equal(result.dates, dates)
    _same(result.values, values)
    assert result.series.interval.value == "Irregular"


def test_parse_usgs_invalid_inputs_raise_value_error() -> None:
    # Act / Assert
    with pytest.raises(ValueError, match="daily precipitation or daily snow"):
        download.parse_usgs("x,y,time,value\n", "daily_snow")
    with pytest.raises(ValueError, match="standard UTC offset"):
        download.parse_usgs(_text("usgs-peaks-01646500.csv"), "peak_stage")
    with pytest.raises(ValueError):
        download.parse_usgs("", "not_a_type")


# GHCN -----------------------------------------------------------------------------------------


def _wetterdienst_daily(element: str) -> tuple[np.ndarray, np.ndarray]:
    """Apply wetterdienst's GHCN daily unit factors to the NCEI by-station CSV of the fixture's days."""
    pl = pytest.importorskip("polars")
    pytest.importorskip("wetterdienst")
    from wetterdienst.provider.noaa.ghcn.api import NoaaGhcnValues  # noqa: PLC0415

    rows = list(csv.DictReader(io.StringIO(_text("ncei-access-USC00042319.csv"))))
    frame = pl.DataFrame({"parameter": [element.lower()] * len(rows),
                          "value": [r[element].strip() or None for r in rows]})
    values = NoaaGhcnValues._apply_daily_factors(frame)["value"].to_numpy().astype(np.float64)  # noqa: SLF001
    return np.asarray([r["DATE"] for r in rows], dtype="datetime64[us]"), values


@pytest.mark.parametrize(("element", "series_type", "max_ulp"),
                         [("PRCP", "daily_precipitation", 1), ("SNOW", "daily_snow", 0)])
def test_parse_ghcn_matches_wetterdienst(element: str, series_type: str, max_ulp: int) -> None:
    # Arrange
    dates, expected = _wetterdienst_daily(element)
    # Act
    result = download.parse_ghcn(_text("ghcn-USC00042319.dly"), series_type, "millimeters")
    # Assert: every day the CSV has
    index = {d: i for i, d in enumerate(result.dates)}
    ours = result.values[[index[d] for d in dates]]
    np.testing.assert_array_equal(np.isnan(ours), np.isnan(expected))
    if max_ulp:  # wetterdienst multiplies by 1/10; Numerics divides by 10, the correctly rounded tenth
        np.testing.assert_array_max_ulp(ours[~np.isnan(ours)], expected[~np.isnan(expected)], maxulp=max_ulp)
    else:
        _same(ours, expected)
    if element == "PRCP":  # December 1984 has no PRCP record: 31 filled days
        december = (result.dates >= np.datetime64("1984-12-01")) & (result.dates <= np.datetime64("1984-12-31"))
        assert december.sum() == 31 and np.isnan(result.values[december]).all()


def test_parse_ghcn_snow_is_millimeters_and_flags_are_returned() -> None:
    # Act
    snow = download.parse_ghcn(_text("ghcn-USW00094728.dly"), "daily_snow", "inches", qualifiers=True)
    # Assert: the 30 December 2000 storm at Central Park is 305 mm, 12 inches
    day = int(np.flatnonzero(snow.dates == np.datetime64("2000-12-30"))[0])
    assert snow.values[day] == 305 / 25.4
    assert snow.qualifiers[int(np.flatnonzero(snow.dates == np.datetime64("2001-01-23"))[0])] == " I0"
    assert len(snow.qualifiers) == len(snow.values)


def test_parse_ghcn_cut_record_raises() -> None:
    # Arrange
    lines = _text("ghcn-USW00094728.dly").split("\n")
    # Act / Assert
    with pytest.raises(RuntimeError, match="line 2"):
        download.parse_ghcn(lines[0] + "\n" + lines[1][:150])


# Canada and Australia ---------------------------------------------------------------------------


def test_parse_chmn_daily_matches_csv_module() -> None:
    # Arrange
    rows = [r for r in csv.reader(io.StringIO(_text("chmn-05BH004-daily.csv").lstrip("﻿")))][1:]
    observed = {np.datetime64(r[1], "us"): (float(r[3]) if r[3] else np.nan, r[4]) for r in rows if r}
    # Act
    result = download.parse_chmn(_text("chmn-05BH004-daily.csv"), "05BH004", qualifiers=True)
    # Assert
    for date, value, symbol in zip(result.dates, result.values, result.qualifiers, strict=True):
        if date in observed:
            assert value == observed[date][0] and symbol == observed[date][1]
        else:
            assert np.isnan(value) and symbol == ""
    assert len(result.values) == 55 and not result.utc


def test_parse_chmn_realtime_is_utc() -> None:
    # Act
    result = download.parse_chmn(_text("chmn-08MF005-realtime.csv"), "08MF005", "instantaneous_discharge")
    # Assert
    assert result.utc and result.dates[0] == np.datetime64("2026-11-01T07:00:00")


@pytest.mark.filterwarnings("ignore:The inplace keyword:Warning")  # inside kiwis-pie
def test_parse_abom_matches_kiwis_pie(monkeypatch) -> None:
    # Arrange
    pytest.importorskip("kiwis_pie")
    import requests  # noqa: PLC0415
    from kiwis_pie import KIWIS  # noqa: PLC0415

    text = _text("bom-410730-daily-quality.json")

    class Response:
        status_code, content, url = 200, text.encode(), "http://bom.invalid/"

        def raise_for_status(self) -> None:
            pass

        def json(self):
            return json.loads(text)

    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: Response())
    frame = KIWIS("http://bom.invalid/").get_timeseries_values(
        ts_id="1573010", to="1990-07-06", **{"from": "1990-05-28"},
        return_fields=["Timestamp", "Value", "Quality Code"], keep_tz=False)
    # Act
    result = download.parse_abom(text, "daily_discharge", qualifiers=True)
    # Assert
    _same(result.values, frame["Value"].astype(float).to_numpy())
    assert list(result.qualifiers) == [str(int(q)) for q in frame["Quality Code"]]
    days = np.asarray(frame.index.normalize().tz_localize(None), dtype="datetime64[us]")
    np.testing.assert_array_equal(result.dates, days)


def test_parse_abom_instantaneous_keeps_wall_clock() -> None:
    # Arrange
    data = json.loads(_text("bom-410730-instantaneous.json"))[0]["data"]
    # Act
    result = download.parse_abom(_text("bom-410730-instantaneous.json"), "instantaneous_discharge")
    # Assert
    np.testing.assert_array_equal(result.dates, np.asarray([d[0][:19] for d in data], dtype="datetime64[us]"))
    _same(result.values, np.asarray([d[1] for d in data], dtype=np.float64))
    assert not result.utc


# Errors, timeout and speed ---------------------------------------------------------------------


def test_timeout_cancels_and_raises_timeout_error() -> None:
    # Arrange
    from pynumerics._util import namespace  # noqa: PLC0415

    namespace("Numerics")
    from System import TimeSpan  # noqa: PLC0415
    from System.Threading.Tasks import Task  # noqa: PLC0415

    started = time.perf_counter()
    # Act / Assert
    with pytest.raises(TimeoutError, match="0.2 seconds"):
        download._run(lambda token: Task.Delay(TimeSpan.FromSeconds(30), token), 0.2)  # noqa: SLF001
    assert time.perf_counter() - started < 5
    with pytest.raises(ValueError, match="positive"):
        download._run(lambda token: Task.Delay(1, token), 0)  # noqa: SLF001


def test_downloaded_series_fields() -> None:
    # Act
    result = download.parse_usgs(_text("usgs-daily-05082500-ice.csv"))
    # Assert
    assert isinstance(result, DownloadedSeries)
    np.testing.assert_array_equal(result.values, result.series.values)
    np.testing.assert_array_equal(result.dates, result.series.dates)
    with pytest.raises(AttributeError):
        result.utc = True  # frozen


def test_parse_usgs_bulk_speed() -> None:
    # Arrange: 40,000 daily records with qualifiers, about a century of record
    days = np.arange(np.datetime64("1900-01-01"), np.datetime64("1900-01-01") + np.timedelta64(40000, "D"))
    text = "x,y,time,value,approval_status,qualifier\n" + "".join(
        f",,{d},{i % 997 + 0.5},Approved,{'ICE' if i % 50 == 0 else ''}\n" for i, d in enumerate(days))
    download.parse_usgs(text[:2000])  # warm up
    # Act
    started = time.perf_counter()
    result = download.parse_usgs(text, qualifiers=True)
    elapsed = time.perf_counter() - started
    # Assert: the parse runs in Numerics and the dates, values and qualifiers cross in bulk
    assert len(result.values) == 40000 and result.qualifiers.count("ICE") == 800
    assert elapsed < 1.0, f"{elapsed:.3f} s"


# Live (opt-in) ------------------------------------------------------------------------------------


@NETWORK
def test_live_usgs_peaks() -> None:
    result = download.from_usgs("01646500", "peak_discharge", timeout=120)
    assert len(result.values) > 90 and result.raw_text.startswith("# U.S. Geological Survey")


@NETWORK
def test_live_chmn_window() -> None:
    result = download.from_chmn("08MF005", start="2021-11-01", end="2021-11-30", qualifiers=True, timeout=60)
    assert len(result.values) == 30 and len(result.qualifiers) == 30


@NETWORK
def test_live_abom_window() -> None:
    result = download.from_abom("410730", start="2021-11-01", end="2021-11-30", qualifiers=True, timeout=60)
    assert len(result.values) >= 29 and all(q != "" for q in result.qualifiers)


@NETWORK
def test_live_connectivity() -> None:
    assert download.is_connected_to_internet(timeout=30)
