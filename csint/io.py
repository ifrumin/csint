"""Readers: CSV, HDF5, ABF (pyabf), ANDI netCDF (.cdf), Agilent .D/.ch (rainbow).

Every reader returns a list[Trace] (a file may hold several channels).
Time axes are converted to MINUTES.
"""
from __future__ import annotations

import os

import numpy as np

from .trace import Trace


def load(path: str, **kw) -> list[Trace]:
    """Sniff by extension and dispatch."""
    p = path.lower()
    if os.path.isdir(path) and p.endswith(".d"):
        return read_agilent_d(path)
    if p.endswith(".ch"):
        return read_agilent_ch(path)
    if p.endswith(".abf"):
        return read_abf(path, **kw)
    if p.endswith(".cdf") or p.endswith(".nc"):
        return read_andi_cdf(path)
    if p.endswith((".h5", ".hdf5")):
        return read_hdf5(path, **kw)
    if p.endswith((".csv", ".txt", ".tsv")):
        return read_csv(path, **kw)
    raise ValueError(f"unrecognized chromatogram file: {path}")


# ---------------------------------------------------------------- CSV / text
def read_csv(path: str, time_unit: str = "auto", **kw) -> list[Trace]:
    """Two+ columns: time, signal[, signal2 ...]. Header optional.

    time_unit: 'min', 's', or 'auto' (heuristic: if median dt > 0.05 the axis
    is assumed to be seconds only when total span > 500 — otherwise ambiguous
    and 'auto' raises; pass explicitly for odd data).
    """
    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None:
        df = pd.read_csv(path, sep=None, engine="python", comment="#", **kw)
        df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").dropna()
        if df.shape[1] < 2:
            raise ValueError("CSV needs at least time + one signal column")
        t = df.iloc[:, 0].to_numpy(float)
        cols = [str(c) for c in df.columns]
        data = [df.iloc[:, c].to_numpy(float) for c in range(1, df.shape[1])]
        name0 = cols[0]
    else:
        # numpy fallback (no pandas): comma/whitespace-delimited numeric CSV
        with open(path) as fh:
            first = fh.readline()
        delim = "," if "," in first else None
        has_header = any(ch.isalpha() for ch in first)
        cols = [c.strip() for c in (first.split(delim) if delim else first.split())] \
            if has_header else []
        arr = np.genfromtxt(path, delimiter=delim,
                            skip_header=1 if has_header else 0, comments="#")
        arr = arr[~np.isnan(arr).any(axis=1)]
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError("CSV needs at least time + one signal column")
        t = arr[:, 0]
        data = [arr[:, c] for c in range(1, arr.shape[1])]
        if len(cols) != arr.shape[1]:
            cols = [f"col{i}" for i in range(arr.shape[1])]
        name0 = cols[0]

    unit = _resolve_time_unit(t, time_unit, name0)
    if unit == "s":
        t = t / 60.0
    return [
        Trace(t, y, {"source": path, "channel": cols[i + 1], "format": "csv"})
        for i, y in enumerate(data)
    ]


def _resolve_time_unit(t: np.ndarray, requested: str, colname) -> str:
    if requested in ("min", "s"):
        return requested
    name = str(colname).lower()
    if "min" in name:
        return "min"
    if "sec" in name or name.strip() in ("s", "t/s", "time (s)", "time_s"):
        return "s"
    span = t[-1] - t[0]
    if span > 500:  # >500 "minutes" is an implausible GC run; assume seconds
        return "s"
    if span < 500 and np.median(np.diff(t)) < 0.05:
        return "min"
    raise ValueError(
        "ambiguous time unit — pass time_unit='min' or 's' explicitly"
    )


# ------------------------------------------------------------------- HDF5
def read_hdf5(path: str, time_key: str | None = None,
              signal_keys: list[str] | None = None) -> list[Trace]:
    """Generic HDF5: find a 1-D time dataset + same-length 1-D signals.

    If keys are not given, search for datasets named like time/t/minutes and
    pair every other same-length 1-D float dataset with it.
    """
    import h5py

    with h5py.File(path, "r") as f:
        ds = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1 and obj.size > 8 \
                    and np.issubdtype(obj.dtype, np.number):
                ds[name] = obj[()]

        f.visititems(visit)

    if not ds:
        raise ValueError("no 1-D numeric datasets found in HDF5 file")

    if time_key is None:
        cands = [k for k in ds if k.split("/")[-1].lower()
                 in ("t", "time", "time_min", "minutes", "time_s", "t_s")]
        if len(cands) != 1:
            raise ValueError(
                f"cannot infer time dataset (candidates: {cands}); pass time_key"
            )
        time_key = cands[0]
    t = np.asarray(ds[time_key], float)
    unit = "s" if time_key.lower().endswith(("_s", "time_s")) or (t[-1] - t[0]) > 500 else "min"
    if unit == "s":
        t = t / 60.0

    if signal_keys is None:
        signal_keys = [k for k in ds if k != time_key and len(ds[k]) == len(t)]
        if not signal_keys:
            raise ValueError("no signal datasets matching time length")
    return [
        Trace(t, np.asarray(ds[k], float),
              {"source": path, "channel": k, "format": "hdf5"})
        for k in signal_keys
    ]


# -------------------------------------------------------------------- ABF
def read_abf(path: str, channels: list[int] | None = None) -> list[Trace]:
    """Axon Binary Format via pyabf. Sweeps are concatenated if multiple."""
    import pyabf

    abf = pyabf.ABF(path)
    chans = channels if channels is not None else list(abf.channelList)
    traces = []
    for ch in chans:
        ys, ts = [], []
        for sw in abf.sweepList:
            abf.setSweep(sw, channel=ch, absoluteTime=True)
            ts.append(abf.sweepX.copy())
            ys.append(abf.sweepY.copy())
        t = np.concatenate(ts) / 60.0  # pyabf gives seconds
        y = np.concatenate(ys).astype(float)
        traces.append(
            Trace(t, y, {
                "source": path, "channel": f"ch{ch}:{abf.adcNames[ch]}",
                "units": abf.adcUnits[ch], "rate_hz": abf.dataRate,
                "format": "abf",
            })
        )
    return traces


# ------------------------------------------------------- ANDI netCDF (.cdf)
def read_andi_cdf(path: str) -> list[Trace]:
    """ANDI/AIA chromatography netCDF (ASTM E1947).

    Uses ordinate_values + either actual_run_time_length/actual_delay_time or
    actual_sampling_interval to rebuild the time axis (seconds → minutes).
    """
    import netCDF4

    nc = netCDF4.Dataset(path, "r")
    try:
        if "ordinate_values" not in nc.variables:
            raise ValueError("not an ANDI chromatography file (no ordinate_values)")
        y = np.asarray(nc.variables["ordinate_values"][:], float)
        n = len(y)

        delay = _cdf_scalar(nc, "actual_delay_time", 0.0)
        interval = _cdf_scalar(nc, "actual_sampling_interval", None)
        if interval is None:
            runlen = _cdf_scalar(nc, "actual_run_time_length", None)
            if runlen is None:
                raise ValueError("cannot reconstruct time axis (no interval/run length)")
            interval = (runlen - delay) / max(n - 1, 1)
        t = (delay + interval * np.arange(n)) / 60.0  # ANDI times are seconds

        detector = getattr(nc, "detector_name", "")
        units = getattr(nc, "detector_unit", "")
        sample = getattr(nc, "sample_name", "")
        return [Trace(t, y, {
            "source": path, "channel": str(detector) or "andi",
            "units": str(units), "sample": str(sample), "format": "andi-cdf",
        })]
    finally:
        nc.close()


def _cdf_scalar(nc, name, default):
    if name in nc.variables:
        v = np.asarray(nc.variables[name][:]).ravel()
        if v.size:
            return float(v[0])
    return default


# --------------------------------------------------------- Agilent .D / .ch
def read_agilent_d(dirpath: str) -> list[Trace]:
    """All signal channels of an Agilent .D directory via rainbow."""
    import rainbow as rb

    datadir = rb.read(dirpath)
    traces = []
    for f in datadir.datafiles:
        traces.extend(_rainbow_to_traces(f, dirpath))
    if not traces:
        raise ValueError(f"no readable signals in {dirpath}")
    return traces


def read_agilent_ch(path: str) -> list[Trace]:
    import rainbow.agilent.chemstation as cs

    f = cs.parse_ch(path)
    if f is None:
        raise ValueError(f"rainbow could not parse {path}")
    return _rainbow_to_traces(f, path)


def _rainbow_to_traces(f, source: str) -> list[Trace]:
    t = np.asarray(f.xlabels, float)  # rainbow: minutes
    data = np.asarray(f.data, float)
    if data.ndim == 1:
        data = data[:, None]
    ylabels = list(getattr(f, "ylabels", range(data.shape[1])))
    out = []
    for i in range(data.shape[1]):
        out.append(Trace(t, data[:, i], {
            "source": source, "file": getattr(f, "name", ""),
            "channel": str(ylabels[i]) if i < len(ylabels) else str(i),
            "units": getattr(f, "metadata", {}).get("unit", ""),
            "format": "agilent",
        }))
    return out
