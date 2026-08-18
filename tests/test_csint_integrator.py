"""Regression tests for the integrator, report, and CLI (no reference data
needed — synthetic truth only). Run: python3 tests/test_integrator.py"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from csint.trace import Trace                        # noqa: E402
from csint.synth import standard_suite               # noqa: E402
from csint.integrator import Params, integrate       # noqa: E402
from csint.report import peak_table, table_rows      # noqa: E402


def case(name):
    return [c for c in standard_suite() if c.name == name][0]


def test_clean_areas():
    c = case("clean_separated")
    res = integrate(Trace(c.t, c.y, {}), Params(slope_sensitivity=5, peak_width=0.05))
    assert len(res.peaks) == 4, len(res.peaks)
    for pk, tp in zip(res.peaks, c.peaks):
        err = abs(pk.area_min - tp.area) / tp.area
        assert err < 0.001, (pk.rt, err)
        assert pk.code == "BB"
    return "clean: 4 peaks, all <0.1%"


def test_area_sum():
    c = case("noisy_offset")
    res = integrate(Trace(c.t, c.y, {}),
                    Params(slope_sensitivity=400, peak_width=0.05, height_reject=6,
                           timed_events=[(1.5, "AreaSum", True), (4.5, "AreaSum", False)]))
    plus = [pk for pk in res.peaks if pk.code.endswith("+")]
    assert len(plus) == 1
    assert abs(plus[0].area_min - 60.0) < 3.0, plus[0].area_min
    return "area-sum: one '+' peak ~60"


def test_negative_dip():
    c = case("negative_dip")
    res = integrate(Trace(c.t, c.y, {}),
                    Params(slope_sensitivity=300, peak_width=0.05, height_reject=4,
                           negative_peaks=True))
    codes = [pk.code for pk in res.peaks]
    assert sum(1 for x in codes if "N" in x) == 1, codes
    assert len(res.peaks) == 3, codes
    return "negative dip: 3 peaks incl. one N"


def test_skims():
    c = case("rider_on_tail")
    res = integrate(Trace(c.t, c.y, {}),
                    Params(slope_sensitivity=300, peak_width=0.07, height_reject=4,
                           tail_skim_ratio=3))
    skimmed = [pk for pk in res.peaks if pk.code.endswith("T")]
    parent = max(res.peaks, key=lambda pk: pk.height)
    assert len(skimmed) >= 1
    assert abs(parent.area_min - 500) / 500 < 0.02, parent.area_min
    return "skims: parent within 2%, child T-coded"


def test_integration_off_window():
    c = case("clean_separated")
    res = integrate(Trace(c.t, c.y, {}),
                    Params(slope_sensitivity=5, peak_width=0.05,
                           timed_events=[(3.0, "Integration", False),
                                         (5.0, "Integration", True)]))
    assert len(res.peaks) == 3, len(res.peaks)   # peak at 4.0 suppressed
    return "integration-off: window peak suppressed"


def test_report_and_cli():
    c = case("clean_separated")
    res = integrate(Trace(c.t, c.y, {}), Params(slope_sensitivity=5, peak_width=0.05))
    txt = peak_table(res)
    assert "Totals" in txt and txt.count("\n") >= 6
    rows = table_rows(res)
    assert len(rows) == 4 and abs(sum(r["area_percent"] for r in rows) - 100) < 1e-6
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "a.csv")
        np.savetxt(p, np.column_stack([c.t, c.y]), delimiter=",",
                   header="time_min,signal", comments="")
        from csint.__main__ import main
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            main([p, "--ss", "5", "--pw", "0.05",
                  "--json", os.path.join(tmp, "o.json")])
        out = buf.getvalue()
        assert "Totals" in out
        assert os.path.exists(os.path.join(tmp, "o.json"))
    return "report+CLI: table renders, JSON written"


if __name__ == "__main__":
    fails = 0
    for fn in (test_clean_areas, test_area_sum, test_negative_dip,
               test_skims, test_integration_off_window, test_report_and_cli):
        try:
            print("PASS:", fn())
        except AssertionError as e:
            print(f"FAIL: {fn.__name__}: {e}")
            fails += 1
    raise SystemExit(fails)
