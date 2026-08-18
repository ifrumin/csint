"""Synthetic chromatograms with exactly known ground truth.

Peaks are exponentially modified Gaussians (EMG) — the standard model for
tailing chromatographic peaks. The analytic area of an EMG with amplitude
parameter A (area) is exactly A, so truth areas are exact by construction,
not fitted.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.special import erfcx


@dataclass
class TruthPeak:
    rt_min: float          # location parameter mu of the EMG (min)
    area: float            # exact analytic area (signal-units * min)
    sigma_min: float       # Gaussian sigma (min)
    tau_min: float         # exponential tail time constant (min); 0 => Gaussian
    apex_time: float = 0.0     # numerically located apex time (filled in)
    apex_height: float = 0.0   # numerically located apex height (filled in)


@dataclass
class SyntheticCase:
    name: str
    t: np.ndarray
    y: np.ndarray
    peaks: list[TruthPeak]
    baseline_desc: str
    noise_sd: float
    meta: dict = field(default_factory=dict)

    def truth_dict(self) -> dict:
        return {
            "name": self.name,
            "baseline": self.baseline_desc,
            "noise_sd": self.noise_sd,
            "peaks": [asdict(p) for p in self.peaks],
            "meta": self.meta,
        }


def emg(t: np.ndarray, area: float, mu: float, sigma: float, tau: float) -> np.ndarray:
    """EMG with exact analytic area `area`. tau=0 degenerates to a Gaussian.

    Numerically stable form using erfcx (scaled complementary error function).
    """
    if tau <= 1e-12:
        return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((t - mu) ** 2) / (2 * sigma ** 2))
    # Two-branch numerically stable EMG.
    # v >= 0 (front/near apex): h = A/(2tau) * erfcx(v) * exp(-(t-mu)^2/(2 sigma^2))
    # v < 0  (far tail): h = A/(2tau) * exp(u) * erfc(v), u = sigma^2/(2 tau^2) - (t-mu)/tau
    from scipy.special import erfc

    v = (sigma / tau - (t - mu) / sigma) / np.sqrt(2)
    out = np.empty_like(np.asarray(t, dtype=float))
    m = v >= 0
    out[m] = erfcx(v[m]) * np.exp(-((t[m] - mu) ** 2) / (2 * sigma ** 2))
    u = (sigma ** 2) / (2 * tau ** 2) - (t[~m] - mu) / tau
    out[~m] = np.exp(u) * erfc(v[~m])
    out *= area / (2 * tau)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("EMG produced non-finite values")
    return out


def _fill_apex(t, peaks_y, peaks):
    for y, p in zip(peaks_y, peaks):
        i = int(np.argmax(np.abs(y)))   # supports negative peaks
        p.apex_time = float(t[i])
        p.apex_height = float(y[i])


def build_case(
    name: str,
    peaks: list[TruthPeak],
    run_min: float = 10.0,
    rate_hz: float = 20.0,
    noise_sd: float = 0.0,
    baseline: str = "flat0",
    seed: int = 0,
) -> SyntheticCase:
    """baseline: 'flat0' | 'offset:<c>' | 'drift:<a>' (a units/min) |
    'exp:<h>,<k>' (h*exp(-k*t)) | combinations joined by '+'.
    """
    n = int(run_min * 60 * rate_hz) + 1
    t = np.linspace(0.0, run_min, n)
    comps = []
    for part in baseline.split("+"):
        part = part.strip()
        if part == "flat0" or not part:
            comps.append(np.zeros_like(t))
        elif part.startswith("offset:"):
            comps.append(np.full_like(t, float(part.split(":")[1])))
        elif part.startswith("drift:"):
            comps.append(float(part.split(":")[1]) * t)
        elif part.startswith("exp:"):
            h, k = (float(x) for x in part.split(":")[1].split(","))
            comps.append(h * np.exp(-k * t))
        else:
            raise ValueError(f"unknown baseline component {part!r}")
    base = np.sum(comps, axis=0)

    peaks_y = [emg(t, p.area, p.rt_min, p.sigma_min, p.tau_min) for p in peaks]
    _fill_apex(t, peaks_y, peaks)
    y = base + (np.sum(peaks_y, axis=0) if peaks_y else 0.0)

    rng = np.random.default_rng(seed)
    if noise_sd > 0:
        y = y + rng.normal(0.0, noise_sd, size=n)

    return SyntheticCase(name, t, y, peaks, baseline, noise_sd,
                         {"rate_hz": rate_hz, "run_min": run_min, "seed": seed})


def standard_suite() -> list[SyntheticCase]:
    """Validation tier-1 suite. Areas in unit·min, heights implied."""
    cases = []

    # 1. Clean well-separated peaks, flat baseline — sanity.
    cases.append(build_case(
        "clean_separated",
        [TruthPeak(2.0, 10.0, 0.020, 0.010),
         TruthPeak(4.0, 50.0, 0.030, 0.015),
         TruthPeak(6.5, 5.0, 0.035, 0.020),
         TruthPeak(8.5, 100.0, 0.040, 0.030)],
        noise_sd=0.0))

    # 2. Same with realistic noise + offset.
    cases.append(build_case(
        "noisy_offset",
        [TruthPeak(2.0, 10.0, 0.020, 0.010),
         TruthPeak(4.0, 50.0, 0.030, 0.015),
         TruthPeak(6.5, 5.0, 0.035, 0.020),
         TruthPeak(8.5, 100.0, 0.040, 0.030)],
        noise_sd=2.0, baseline="offset:50", seed=1))

    # 3. Linear drift — tests baseline tracking.
    cases.append(build_case(
        "drift",
        [TruthPeak(2.5, 20.0, 0.025, 0.012),
         TruthPeak(5.0, 20.0, 0.030, 0.015),
         TruthPeak(7.5, 20.0, 0.035, 0.020)],
        noise_sd=1.0, baseline="offset:20+drift:15", seed=2))

    # 4. Fused pair (valley, no baseline between) — drop/valley logic.
    cases.append(build_case(
        "fused_pair",
        [TruthPeak(5.00, 60.0, 0.040, 0.020),
         TruthPeak(5.18, 30.0, 0.040, 0.020)],
        noise_sd=1.0, seed=3))

    # 5. Rider on a tailing parent — skim logic.
    cases.append(build_case(
        "rider_on_tail",
        [TruthPeak(3.0, 500.0, 0.050, 0.250),   # heavy-tailed parent
         TruthPeak(3.8, 8.0, 0.025, 0.012),     # rider on the tail
         TruthPeak(4.6, 8.0, 0.025, 0.012)],    # second rider
        noise_sd=1.0, seed=4))

    # 6. Small peaks near rejection limits + one shoulder.
    cases.append(build_case(
        "small_and_shoulder",
        [TruthPeak(2.0, 1.0, 0.020, 0.010),
         TruthPeak(5.00, 40.0, 0.035, 0.015),
         TruthPeak(5.07, 12.0, 0.035, 0.015),   # shoulder (no valley)
         TruthPeak(8.0, 0.5, 0.025, 0.010)],
        noise_sd=0.5, seed=5))

    # 7. Negative dip between positive peaks (septum purge / pressure dip).
    cases.append(build_case(
        "negative_dip",
        [TruthPeak(2.5, 30.0, 0.030, 0.015),
         TruthPeak(5.0, -10.0, 0.030, 0.015),   # negative peak
         TruthPeak(7.5, 30.0, 0.030, 0.015)],
        noise_sd=1.0, baseline="offset:50", seed=6))

    return cases


def save_suite(cases: list[SyntheticCase], outdir: str) -> None:
    import os

    os.makedirs(outdir, exist_ok=True)
    for c in cases:
        np.savetxt(
            os.path.join(outdir, f"{c.name}.csv"),
            np.column_stack([c.t, c.y]),
            delimiter=",", header="time_min,signal", comments="",
        )
        with open(os.path.join(outdir, f"{c.name}.truth.json"), "w") as fh:
            json.dump(c.truth_dict(), fh, indent=1)
