"""Tier-1 validation: integrate the synthetic suite, compare to exact truth.

Honest reporting: every truth peak is MATCHED (with area/RT errors), MISSED,
or a detected peak is EXTRA. No averaging away of failures.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from csint.trace import Trace                    # noqa: E402
from csint.synth import standard_suite           # noqa: E402
from csint.integrator import Params, integrate   # noqa: E402

RT_TOL = 0.06  # min, matching tolerance truth<->detected

# Per-case parameters (explicit, like a method would be).
# slope_sensitivity in units/min. Noise sd=2 at 20 Hz gives derivative noise
# of roughly sd*sqrt(2)/dt ~ 100 units/min on filter1; filters 2/3 average
# it down. Values chosen per-case like an operator would; autointegrate later.
CASE_PARAMS = {
    "clean_separated":    Params(slope_sensitivity=5,   peak_width=0.05),
    "noisy_offset":       Params(slope_sensitivity=400, peak_width=0.05, height_reject=6),
    "drift":              Params(slope_sensitivity=300, peak_width=0.05, height_reject=4),
    "fused_pair":         Params(slope_sensitivity=300, peak_width=0.07, height_reject=4),
    "rider_on_tail":      Params(slope_sensitivity=300, peak_width=0.07, height_reject=4,
                                 tail_skim_ratio=3),
    "small_and_shoulder": Params(slope_sensitivity=60,  peak_width=0.05, height_reject=1),
    "negative_dip":       Params(slope_sensitivity=300, peak_width=0.05, height_reject=4,
                                 negative_peaks=True),
}


def run(make_plots=True):
    cases = standard_suite()
    all_rows = []
    for c in cases:
        pr = CASE_PARAMS[c.name]
        res = integrate(Trace(c.t, c.y, {"name": c.name}), pr)
        det = list(res.peaks)
        used = set()
        rows = []
        for tp in c.peaks:
            best, bestd = None, RT_TOL
            for j, dp in enumerate(det):
                if j in used:
                    continue
                d = abs(dp.rt - tp.apex_time)
                if d < bestd:
                    best, bestd = j, d
            if best is None:
                rows.append((c.name, f"truth@{tp.apex_time:.3f}", "MISS",
                             np.nan, np.nan, tp.area, ""))
            else:
                used.add(best)
                dp = det[best]
                # negative truth: detected magnitude + 'N' code required
                if tp.area < 0:
                    stat = "match" if "N" in dp.code else "SIGN?"
                    aerr = 100.0 * (dp.area_min - abs(tp.area)) / abs(tp.area)
                else:
                    stat = "match"
                    aerr = 100.0 * (dp.area_min - tp.area) / tp.area
                rows.append((c.name, f"truth@{tp.apex_time:.3f}", stat,
                             aerr, (dp.rt - tp.apex_time) * 60.0, tp.area, dp.code))
        for j, dp in enumerate(det):
            if j not in used:
                rows.append((c.name, f"det@{dp.rt:.3f}", "EXTRA",
                             np.nan, np.nan, dp.area_min, dp.code))
        all_rows.extend(rows)

        if make_plots:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(11, 3.2))
            ax.plot(c.t, c.y, lw=0.6, color="#1f77b4")
            for dp in det:
                ax.plot([dp.start_t, dp.end_t], [dp.baseline_y0, dp.baseline_y1],
                        color="#d62728", lw=1.0)
                ax.annotate(f"{dp.code}\n{dp.area_min:.1f}", (dp.rt, 0),
                            xytext=(dp.rt, ax.get_ylim()[1] * 0.02),
                            fontsize=6, ha="center", color="#d62728")
                ax.axvline(dp.start_t, color="#999", lw=0.3)
                ax.axvline(dp.end_t, color="#999", lw=0.3)
            for tp in c.peaks:
                ax.axvline(tp.apex_time, color="green", lw=0.4, alpha=0.4)
            ax.set_title(f"{c.name} — red=detected baseline/bounds, green=truth apex", fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(os.path.dirname(__file__), f"tier1_{c.name}.png"), dpi=110)
            plt.close(fig)

    # report
    print(f"{'case':20s} {'peak':14s} {'stat':6s} {'area err %':>10s} {'RT err s':>9s} {'truth A':>8s} code")
    misses = extras = 0
    aerrs = []
    for r in all_rows:
        name, pk, stat, aerr, rterr, ta, code = r
        if stat == "MISS":
            misses += 1
        if stat == "EXTRA":
            extras += 1
        if stat == "match":
            aerrs.append(abs(aerr))
        print(f"{name:20s} {pk:14s} {stat:6s} "
              f"{aerr:10.2f} {rterr:9.2f} {ta:8.2f} {code}")
    print(f"\nmatched {len(aerrs)}, missed {misses}, extra {extras}")
    if aerrs:
        print(f"|area err| median {np.median(aerrs):.2f}%  worst {np.max(aerrs):.2f}%")
    return misses, extras, aerrs


if __name__ == "__main__":
    run()
