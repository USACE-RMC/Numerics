# Downloading time series

`pynumerics.download` wraps `Numerics.Data.TimeSeriesDownload`, which reads public records from
four providers:

| Function | Provider | Series |
|----------|----------|--------|
| `from_usgs` | U.S. Geological Survey Water Data API | daily, instantaneous, peak and measured discharge or stage |
| `from_ghcn` | NOAA NCEI GHCN-Daily station files | daily precipitation and snowfall |
| `from_chmn` | Water Survey of Canada | daily, real-time and peak discharge or stage |
| `from_abom` | Bureau of Meteorology Water Data Online | daily and instantaneous discharge or stage, daily rainfall |

Requests, retries, units, gap filling and parsing all run in Numerics; the wrapper only passes the
arguments and copies the result back in bulk. Each downloader fetches the provider's response and
parses it with a public parser (`parse_usgs`, `parse_ghcn`, `parse_chmn`, `parse_abom`), so a saved
response gives the same series offline. Every function returns a `DownloadedSeries`:

- `series`: the `pynumerics.timeseries.TimeSeries` (block series, smoothing and the other
  time-series members apply directly);
- `dates` and `values`: NumPy arrays (`datetime64[us]` and float64, NaN for a missing value);
- `raw_text`: the text the service sent (USGS downloads only);
- `qualifiers`: with `qualifiers=True`, one string per value;
- `utc`: True when the dates are UTC rather than calendar dates or provider local time.

## Downloading

The examples on this page run against the live services only when the environment variable
`PYNUMERICS_NETWORK_TESTS` is `1`; otherwise they print nothing.

```python
import os
from pynumerics import download

if os.environ.get("PYNUMERICS_NETWORK_TESTS") == "1":
    peaks = download.from_usgs("01646500", "peak_discharge", timeout=120)
    print(len(peaks.values), peaks.values.max())

    daily = download.from_usgs("01646500", start="2020-10-01", end="2021-09-30", statistic="maximum")
    annual = daily.series.custom_year_series(start_month=10)    # water-year maxima, in Numerics
    print(annual.values)

    rain = download.from_ghcn("USW00094728", "daily_precipitation", "millimeters")
    flow = download.from_chmn("08MF005", start="2021-11-01", end="2021-11-30", qualifiers=True)
    print(rain.dates[-1], flow.qualifiers[:5])
```

- USGS site numbers have 8 to 15 digits. Daily values form a regular one-day series with NaN for
  missing days; instantaneous values and field measurements are irregular and UTC; peaks are an
  irregular series, and a peak with a time of day is dated in the site's local standard time, as in
  the official peak-flow file.
- GHCN precipitation is stored in tenths of millimeters and snowfall in millimeters; both are
  returned in the requested unit. Days after today are left out.
- Canadian real-time values are UTC. BOM daily series are the midnight-to-midnight (`.24HR`) series.
- The USGS service limits anonymous use; `set_usgs_api_key` sets a free key from
  https://api.waterdata.usgs.gov/signup/ for every later download.

## Parsing saved responses

The parsers take the provider's text. A USGS response is CSV; the raw text `from_usgs` returns can
be saved and parsed again:

```python
from pynumerics import download

text = (
    "x,y,time,value,approval_status,qualifier\n"
    ",,2026-02-12,1250,Provisional,\n"
    ",,2026-02-13,1190,Provisional,ESTIMATED\n"
    ",,2026-02-15,,Provisional,ICE\n"
)
result = download.parse_usgs(text, "daily_discharge", qualifiers=True)
print(result.dates)        # 2026-02-14 is filled as a missing day
print(result.values)       # [1250. 1190.   nan   nan]
print(result.qualifiers)   # ('', 'ESTIMATED', '', 'ICE')
```

A GHCN-Daily station file has one fixed-width record per station, month and element: the station,
year, month and element, then 31 days of a five-character value and three flag characters:

```python
from pynumerics import download

days = [25, 0, -9999] + [0] * 28
record = "USW00094728202401PRCP" + "".join(f"{v:5d}   " for v in days)
rain = download.parse_ghcn(record, "daily_precipitation", "millimeters")
print(rain.values[:3])     # [2.5 0.  nan]: tenths of a millimeter, -9999 is missing
```

Water Survey of Canada CSV and BOM KiWIS JSON parse the same way:

```python
from pynumerics import download

csv_text = (
    "ID,Date,Parameter/Paramètre,Value/Valeur,Symbol/Symbole\n"
    "05BH004,1911-01-01,discharge/débit,19.5,B\n"
    "05BH004,1911-01-02,discharge/débit,18.1,B\n"
)
canada = download.parse_chmn(csv_text, "05BH004", qualifiers=True)
print(canada.values, canada.qualifiers)

json_text = ('[{"ts_id":"1573010","columns":"Timestamp,Value,Quality Code","data":'
             '[["1990-05-30T00:00:00.000+10:00",1.813,10],["1990-05-31T00:00:00.000+10:00",5.906,90]]}]')
bom = download.parse_abom(json_text, "daily_discharge", "cubic_feet_per_second", qualifiers=True)
print(bom.values, bom.qualifiers)
```

## Qualifiers

| Provider | Qualifier |
|----------|-----------|
| USGS | the `qualifier` field, for example `ICE` or `ESTIMATED,ICE` |
| GHCN | the day's measurement, quality and source flags, for example `T 1` for a trace |
| Canada | the data symbol, for example `B` (ice conditions) or `E` (estimated) |
| BOM | the KiWIS quality code, for example `10`, `90` or `-1` |

A day filled as missing has an empty qualifier. The series is the same with or without them.

## Errors and time limits

`timeout` cancels a download after that many seconds and raises `TimeoutError`. Each request also
has its own limit inside Numerics, and transient failures (server errors, rate limits with a short
`Retry-After`, connections that drop mid-body, USGS pages cut off mid-record) are retried there.

| Failure | Python exception |
|---------|------------------|
| invalid site number, type or unit | `ValueError` |
| `timeout` or a request time limit exceeded | `TimeoutError` |
| the service cannot be reached, or a station file does not exist | `ConnectionError` |
| a service error, a spent rate limit, or "No data found" | `RuntimeError` |

```python
from pynumerics import download

try:
    download.parse_usgs("x,y,time,value\n", "daily_snow")
except ValueError as error:
    print(error)
```
