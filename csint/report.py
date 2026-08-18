"""ChemStation-style peak table + overlay report."""
from __future__ import annotations

import numpy as np

from .integrator import IntegrationResult


def peak_table(res: IntegrationResult) -> str:
    """Text table in ChemStation report style. Area in units*s."""
    lines = [
        "Peak RetTime Type  Width      Area      Height     Area  ",
        "  #  [min]         [min]   [unit*s]     [unit]       %   ",
        "---- ------- ---- ------- ---------- ---------- ---------",
    ]
    for i, (p, ap) in enumerate(zip(res.peaks, res.area_percent()), 1):
        lines.append(
            f"{i:4d} {p.rt:7.3f} {p.code:<4s} {p.width_ah:7.4f} "
            f"{p.area:10.4g} {p.height:10.4g} {ap:9.5f}"
        )
    tot = sum(p.area for p in res.peaks)
    toth = sum(p.height for p in res.peaks)
    lines.append("")
    lines.append(f"Totals :                 {tot:10.4g} {toth:10.4g}")
    return "\n".join(lines)


def table_rows(res: IntegrationResult) -> list[dict]:
    return [
        dict(peak=i, rt_min=round(p.rt, 4), type=p.code,
             width_min=round(p.width_ah, 4), area=round(p.area, 3),
             height=round(p.height, 3), area_percent=round(ap, 4),
             start_min=round(p.start_t, 4), end_min=round(p.end_t, 4))
        for i, (p, ap) in enumerate(zip(res.peaks, res.area_percent()), 1)
    ]


def overlay(res: IntegrationResult, path: str, title: str = "") -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tr = res.trace
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.plot(tr.t, tr.y, lw=0.7, color="#1f77b4", zorder=2)
    ymax = float(np.max(tr.y))
    for i, p in enumerate(res.peaks, 1):
        ax.plot([p.start_t, p.end_t], [p.baseline_y0, p.baseline_y1],
                color="#d62728", lw=1.2, zorder=3)
        for x, yv in ((p.start_t, p.baseline_y0), (p.end_t, p.baseline_y1)):
            ax.plot([x, x], [yv, min(yv + 0.03 * ymax, ymax)], color="#888",
                    lw=0.5, zorder=1)
        ax.annotate(f"{i}\n{p.rt:.3f}\n{p.code}",
                    (p.rt, p.height + max(p.baseline_y0, p.baseline_y1)),
                    textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=7, color="#d62728", zorder=4)
    ax.set_xlabel("min")
    ax.set_ylabel(res.trace.meta.get("units", "signal"))
    ax.set_title(title or f"{tr.meta.get('source','trace')} — csint overlay",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
