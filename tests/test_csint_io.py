"""Round-trip tests: write real files in each format, read back, compare.

.ch/.D cannot be fabricated faithfully — excluded here; tested when a real
.D folder is provided.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from csint import load                      # noqa: E402
from csint.synth import standard_suite      # noqa: E402

TOL = 1e-6


def _reference():
    c = standard_suite()[0]
    return c.t, c.y


def test_csv(tmp):
    t, y = _reference()
    p = os.path.join(tmp, "a.csv")
    np.savetxt(p, np.column_stack([t, y]), delimiter=",",
               header="time_min,signal", comments="")
    tr = load(p)[0]
    assert np.allclose(tr.t, t, atol=1e-9) and np.allclose(tr.y, y, atol=1e-4)
    return "csv ok"


def test_hdf5(tmp):
    import h5py
    t, y = _reference()
    p = os.path.join(tmp, "a.h5")
    with h5py.File(p, "w") as f:
        f["time_min"] = t
        f["fid"] = y
    tr = load(p)[0]
    assert np.allclose(tr.t, t) and np.allclose(tr.y, y)
    return "hdf5 ok"


def test_andi_cdf(tmp):
    import netCDF4
    t, y = _reference()
    p = os.path.join(tmp, "a.cdf")
    nc = netCDF4.Dataset(p, "w", format="NETCDF3_CLASSIC")
    nc.createDimension("point_number", len(y))
    v = nc.createVariable("ordinate_values", "f4", ("point_number",))
    v[:] = y
    nc.createDimension("scalar", 1)
    d = nc.createVariable("actual_delay_time", "f4", ("scalar",)); d[:] = t[0] * 60
    i = nc.createVariable("actual_sampling_interval", "f4", ("scalar",))
    i[:] = float(np.median(np.diff(t))) * 60
    nc.detector_name = "FID"
    nc.close()
    tr = load(p)[0]
    assert np.allclose(tr.t, t, atol=1e-4), np.max(np.abs(tr.t - t))
    assert np.allclose(tr.y, y, atol=1e-2)  # f4 storage
    return "andi-cdf ok"


def test_abf(tmp):
    from pyabf.abfWriter import writeABF1
    t, y = _reference()
    rate = 1.0 / (float(np.median(np.diff(t))) * 60.0)
    p = os.path.join(tmp, "a.abf")
    writeABF1(np.asarray(y, np.float32)[np.newaxis, :], p, rate)
    tr = load(p)[0]
    assert len(tr.y) == len(y)
    assert np.allclose(tr.y, y, atol=0.05 * max(1.0, np.abs(y).max()))
    dt_match = abs(tr.dt - float(np.median(np.diff(t)))) < 1e-6
    assert dt_match, (tr.dt, float(np.median(np.diff(t))))
    return "abf ok"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        for fn in (test_csv, test_hdf5, test_andi_cdf, test_abf):
            try:
                print(fn(tmp))
            except Exception as e:  # noqa: BLE001
                print(f"{fn.__name__} FAILED: {type(e).__name__}: {e}")
                raise SystemExit(1)
    print("ALL IO ROUNDTRIPS PASSED")
