"""Core integrator implementing the ChemStation ("new") integrator per docs/SPEC.md.

Every place where the public Agilent documentation is ambiguous is marked
`SPEC-OPEN`. Those choices are provisional and will be disciplined against
real ChemStation output (task 9) — they are documented, not guessed silently.

v0 scope: positive peaks, static bunching from initial peak width, drop-line
separation of fused clusters, BB/BV/VV/PV codes, height/area reject.
Valley-ratio rule, skims, timed events, shoulders, autointegrate: later tasks.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from .trace import Trace

UP_THRESHOLD = 15    # S2: up-slope accumulator >= 15 -> peak may be starting
DOWN_THRESHOLD = 15  # S2: down-slope accumulator >= 15 -> peak may be ending

# Increment tables (S2 Tables 2-3), columns = (filter1, filter2, filter3).
# SPEC-OPEN: two rows of the published tables were extracted with ambiguous
# labels; the assignments below are the symmetric best reading (see SPEC.md §3).
UP_INC = {
    "slope_up":   (+8, +5, +3),   # slope > SS
    "curv_up":    (+0, +2, +1),   # curvature > SS
    "slope_down": (-8, -5, -3),   # slope < -SS
    "flat":       (-4, -2, -1),   # |slope| < SS          SPEC-OPEN
    "curv_down":  (-0, -2, -1),   # curvature < -SS
}
DOWN_INC = {
    "slope_down": (+8, +5, +3),   # slope < -SS
    "curv_down":  (+0, +2, +1),   # curvature < -SS
    "slope_up":   (-28, -18, -11),  # slope > SS (contradiction)   SPEC-OPEN
    "flat":       (-4, -2, -1),   # |slope| < SS                   SPEC-OPEN
    "curv_up":    (-0, -2, -1),   # curvature > SS
}

# S2 Table 1: (min_pts, max_pts, bunch_power, active_filters)
BUNCH_TABLE = [
    (0, 10, 0, (1,)),
    (8, 16, 0, (2,)),
    (12, 24, 0, (3,)),
    (16, 32, 1, (2,)),
    (24, 48, 1, (3,)),
    (32, 96, 2, (3, 2)),
    (64, 192, 3, (3, 2)),
]

# filter number -> (slope span in bunched pts, curvature spacing in bunched pts)
FILTER_GEOM = {1: (2, 1), 2: (4, 2), 3: (8, 4)}


@dataclass
class Params:
    slope_sensitivity: float          # signal units / min       (SPEC-OPEN: normalization)
    peak_width: float = 0.05          # expected peak width, minutes (area/height sense)
    height_reject: float = 0.0        # baseline-corrected apex height
    area_reject: float = 0.0          # units * s (ChemStation convention)
    anchor_smooth_pts: int = 3        # raw pts averaged at baseline anchors (noise guard)
    shoulders: str = "OFF"            # OFF | DROP  (tangent shoulders: task 7)
    # Skims (S2 §7). Ratios of 0 = disabled (ChemStation convention; Idan's
    # method has front/tail = 0). tail: skim child if Hparent/Hchild > ratio
    # AND Hchild/Hvalley < skim_valley_ratio.
    tail_skim_ratio: float = 0.0
    front_skim_ratio: float = 0.0     # v1: front skims not implemented (documented)
    skim_valley_ratio: float = 20.0
    skim_mode: str = "STANDARD"       # STANDARD (exp clipped at baseline) | EXP | STRAIGHT
    # Timed events: list of (time_min, event_name, value). v1 supports:
    #   ("AreaSum", True/False)      — sum peaks between ON and OFF into one
    #   ("NegativePeak", True/False) — negative-peak detection windows
    #   ("Integration", True/False)  — integration off between False and True
    # (Idan's proc-D uses AreaSum 0.05 -> 0.8 min.)
    timed_events: list = field(default_factory=list)
    negative_peaks: bool = False      # detect negative peaks over the whole run


@dataclass
class Peak:
    start_t: float
    end_t: float
    rt: float                # apex time (parabolic fit)
    height: float            # baseline-corrected apex height
    area: float              # units * s  (ChemStation GC convention)
    area_min: float          # units * min (internal)
    width_ah: float          # area/height width, minutes
    code: str                # e.g. 'BB', 'BV', 'VV'
    cluster: int
    start_i: int
    end_i: int
    baseline_y0: float
    baseline_y1: float

    @property
    def rounded(self):
        return (round(self.rt, 4), round(self.area, 2), round(self.height, 2), self.code)


@dataclass
class IntegrationResult:
    peaks: list[Peak]
    params: Params
    trace: Trace
    bunch_power: int
    active_filters: tuple
    diagnostics: dict = field(default_factory=dict)

    def area_percent(self) -> list[float]:
        tot = sum(p.area for p in self.peaks) or 1.0
        return [100.0 * p.area / tot for p in self.peaks]


# ----------------------------------------------------------------- helpers
def _choose_bunching(expected_pts: float) -> tuple[int, tuple]:
    """Pick bunching power + active filters from Table 1 for the expected
    peak width in RAW data points. SPEC-OPEN: overlapping ranges — we pick the
    last (most-bunched) row whose range contains the value, extending the
    power-of-two pattern beyond the table."""
    best = None
    for lo, hi, power, filts in BUNCH_TABLE:
        if lo <= expected_pts <= hi:
            best = (power, filts)
    if best is not None:
        return best
    if expected_pts < BUNCH_TABLE[0][1]:
        return 0, (1,)
    # extend pattern: keep 64-192-pt row shape, doubling
    power = 3
    lo, hi = 64, 192
    while expected_pts > hi * 2 and power < 12:
        power += 1
        lo, hi = lo * 2, hi * 2
    return power, (3, 2)


def _bunch(y: np.ndarray, power: int) -> np.ndarray:
    """Average groups of 2**power points. SPEC-OPEN: sum vs average — average
    keeps signal units so slope_sensitivity stays in units/min."""
    f = 1 << power
    if f == 1:
        return y.copy()
    n = (len(y) // f) * f
    return y[:n].reshape(-1, f).mean(axis=1)


def _robust_noise(y: np.ndarray, win: int = 64) -> float:
    """Baseline noise of the RAW signal: per-window MAD of linear-detrended
    residuals, 25th percentile across windows (baseline windows dominate;
    peak windows land in the upper tail). Replaces successive-difference MAD,
    which measures quantization on smooth, heavily-sampled FID traces."""
    n = len(y)
    if n < 4 * win:
        d = np.diff(y)
        return 1.4826 * float(np.median(np.abs(d - np.median(d)))) / np.sqrt(2)
    m = (n // win) * win
    seg = y[:m].reshape(-1, win).astype(float)
    x = np.arange(win, dtype=float)
    xc = x - x.mean()
    slope = (seg @ xc) / (xc @ xc)
    resid = seg - seg.mean(axis=1, keepdims=True) - slope[:, None] * xc
    mads = 1.4826 * np.median(np.abs(resid), axis=1)
    return float(np.percentile(mads, 25)) + 1e-12


def _parabolic_apex(t: np.ndarray, y: np.ndarray, i: int) -> tuple[float, float]:
    """Parabola through (i-1, i, i+1); returns (t_apex, y_apex)."""
    if i <= 0 or i >= len(y) - 1:
        return float(t[i]), float(y[i])
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    denom = (y0 - 2 * y1 + y2)
    if denom == 0:
        return float(t[i]), float(y[i])
    d = 0.5 * (y0 - y2) / denom
    d = float(np.clip(d, -1, 1))
    dt = t[i + 1] - t[i]
    ta = float(t[i] + d * dt)
    ya = float(y1 - 0.25 * (y0 - y2) * d)
    return ta, ya


# ------------------------------------------------------------------- core
class Integrator:
    def __init__(self, params: Params):
        self.p = params

    # -- filter outputs on the bunched series -------------------------------
    def _filters(self, tb, yb, filts):
        """Return dict fnum -> (slope, curv) arrays, units/min.
        Backward (causal) differences, as an online integrator would see.
        SPEC-OPEN: exact kernels unpublished; first/second differences over
        the documented spans."""
        dtb = float(tb[1] - tb[0])
        out = {}
        n = len(yb)
        for f in filts:
            span, s = FILTER_GEOM[f]
            slope = np.zeros(n)
            curv = np.zeros(n)
            slope[span:] = (yb[span:] - yb[:-span]) / (span * dtb)
            curv[2 * s:] = (yb[2 * s:] - 2 * yb[s:-s] + yb[:-2 * s]) / (s * dtb)
            # SPEC-OPEN: curvature normalized to units/min (one time division)
            # so it is comparable to slope_sensitivity, as S2 implies.
            out[f] = (slope, curv)
        return out

    def _increment(self, table, fdict, filts, k, ss):
        inc = 0
        for f in filts:
            slope, curv = fdict[f]
            col = f - 1
            s, c = slope[k], curv[k]
            if s > ss:
                inc += table["slope_up"][col]
            elif s < -ss:
                inc += table["slope_down"][col]
            else:
                inc += table["flat"][col]
            if c > ss:
                inc += table["curv_up"][col]
            elif c < -ss:
                inc += table["curv_down"][col]
        return inc

    # -- main ---------------------------------------------------------------
    def run(self, trace: Trace) -> IntegrationResult:
        tr = trace.resampled_uniform()
        p = self.p

        # Negative peaks FIRST: detect dips on the inverted signal, keep true
        # dips (apex below both endpoint levels), then EXCISE their spans by
        # linear interpolation so the positive pass never sees the recovery
        # flank (which otherwise masquerades as a peak start).
        neg_kept = []
        windows = self._event_windows("NegativePeak")
        y_work = tr.y
        if p.negative_peaks or windows:
            inv = Trace(tr.t, -tr.y, dict(tr.meta))
            inner = replace(p, negative_peaks=False, timed_events=[],
                            tail_skim_ratio=0.0, shoulders="OFF")
            for pk in Integrator(inner).run(inv).peaks:
                if windows and not any(w0 <= pk.rt <= w1 for w0, w1 in windows):
                    continue
                i0, i1 = pk.start_i, pk.end_i
                ends_min = min(float(tr.y[i0]), float(tr.y[i1]))
                # true dip: original apex must sit below both endpoint levels
                apex_orig = -pk.height + min(-pk.baseline_y0, -pk.baseline_y1)
                if not (apex_orig < ends_min - p.height_reject):
                    continue
                neg_kept.append(replace(
                    pk, code=pk.code + "N",
                    baseline_y0=-pk.baseline_y0, baseline_y1=-pk.baseline_y1,
                ))
            if neg_kept:
                y_work = tr.y.copy()
                for pk in neg_kept:
                    i0, i1 = pk.start_i, pk.end_i
                    y_work[i0:i1 + 1] = np.linspace(y_work[i0], y_work[i1],
                                                    i1 - i0 + 1)
        tr_pos = tr if y_work is tr.y else Trace(tr.t, y_work, dict(tr.meta))
        t, y = tr_pos.t, tr_pos.y
        dt = tr_pos.dt

        expected_pts = p.peak_width / dt
        power, filts = _choose_bunching(expected_pts)
        yb = _bunch(y, power)
        tb = _bunch(t, power)
        fdict = self._filters(tb, yb, filts)
        ss = p.slope_sensitivity
        f_fine = min(filts)          # finest active filter for sign decisions
        slope_fine = fdict[f_fine][0]
        curv_fine = fdict[f_fine][1]

        # detection state machine over bunched points
        BASE, RISE, FALL = 0, 1, 2
        state = BASE
        up_acc = 0
        down_acc = 0
        up2_acc = 0                  # next-peak detector while falling
        anchor = 0                   # last bunched idx where up_acc was 0
        up2_anchor = 0
        start_k = apex_k = 0
        events = []                  # (start_k, apex_k, end_k, fusedL, fusedR)
        cur_fusedL = False

        k0 = 2 * FILTER_GEOM[max(filts)][1] + 1  # first valid filter output
        for k in range(k0, len(yb)):
            if state == BASE:
                up_acc = max(0, up_acc + self._increment(UP_INC, fdict, filts, k, ss))
                if up_acc == 0:
                    anchor = k
                if up_acc >= UP_THRESHOLD:
                    state = RISE
                    start_k = anchor          # back-date to envelope departure
                    down_acc = 0
            elif state == RISE:
                if slope_fine[k] <= 0:
                    state = FALL
                    apex_k = k
                    down_acc = 0
                    up2_acc = 0
                    up2_anchor = k
            elif state == FALL:
                down_acc = max(0, down_acc + self._increment(DOWN_INC, fdict, filts, k, ss))
                up2_acc = max(0, up2_acc + self._increment(UP_INC, fdict, filts, k, ss))
                if up2_acc == 0:
                    up2_anchor = k
                # fused: next peak rises before this one reaches baseline
                if up2_acc >= UP_THRESHOLD:
                    events.append((start_k, apex_k, up2_anchor, cur_fusedL, True))
                    state = RISE
                    start_k = up2_anchor
                    cur_fusedL = True
                    continue
                if down_acc >= DOWN_THRESHOLD and \
                        abs(slope_fine[k]) < ss and abs(curv_fine[k]) < ss:
                    events.append((start_k, apex_k, k, cur_fusedL, False))
                    state = BASE
                    up_acc = 0
                    anchor = k
                    cur_fusedL = False
        if state == FALL:  # run ended inside a peak
            events.append((start_k, apex_k, len(yb) - 1, cur_fusedL, False))

        # Event-level prominence: the apex must rise above both event
        # endpoints. Removes recovery-flank artifacts (signal climbing back
        # to baseline after a canyon/dip triggers the up-accumulator but is
        # not a peak). Lenient threshold: half the height reject.
        def _prominent(ev):
            a = yb[min(ev[1], len(yb) - 1)]
            return (a - max(yb[ev[0]], yb[min(ev[2], len(yb) - 1)])
                    >= 0.5 * p.height_reject)
        events = [ev for ev in events if _prominent(ev)]

        if p.shoulders == "DROP":
            events = self._split_shoulders(events, slope_fine, ss)

        peaks = self._build_peaks(t, y, tb, events, power)
        peaks = [pk for pk in peaks
                 if pk.height >= p.height_reject and pk.area >= p.area_reject]
        peaks = sorted(peaks + neg_kept, key=lambda pk: pk.rt)
        peaks = self._apply_integration_windows(peaks)
        peaks = self._apply_area_sum(peaks)
        return IntegrationResult(peaks, p, tr, power, filts,
                                 {"expected_pts": expected_pts})

    def _event_windows(self, name):
        """ON..OFF windows for a boolean timed event."""
        windows, t_on = [], None
        for tt, nm, val in sorted(self.p.timed_events):
            if nm != name:
                continue
            if val and t_on is None:
                t_on = tt
            elif not val and t_on is not None:
                windows.append((t_on, tt))
                t_on = None
        if t_on is not None:
            windows.append((t_on, float("inf")))
        return windows

    def _apply_integration_windows(self, peaks):
        """Integration OFF windows: peaks whose apex falls inside are dropped.
        v1: no partial clipping of segments (documented)."""
        off = self._event_windows("IntegrationOff")
        # also accept ("Integration", False..True) phrasing
        t_off = None
        for tt, nm, val in sorted(self.p.timed_events):
            if nm != "Integration":
                continue
            if not val and t_off is None:
                t_off = tt
            elif val and t_off is not None:
                off.append((t_off, tt))
                t_off = None
        if t_off is not None:
            off.append((t_off, float("inf")))
        if not off:
            return peaks
        return [pk for pk in peaks
                if not any(w0 <= pk.rt <= w1 for w0, w1 in off)]

    def _apply_area_sum(self, peaks):
        """AreaSum timed event: sum all peaks with RT inside each ON..OFF
        window into a single '+'-coded peak (RT/height of the tallest)."""
        windows = []
        t_on = None
        for tt, name, val in sorted(self.p.timed_events):
            if name != "AreaSum":
                continue
            if val and t_on is None:
                t_on = tt
            elif not val and t_on is not None:
                windows.append((t_on, tt))
                t_on = None
        if t_on is not None:
            windows.append((t_on, float("inf")))
        for w0, w1 in windows:
            grp = [pk for pk in peaks if w0 <= pk.rt <= w1]
            if len(grp) < 2:
                continue
            tallest = max(grp, key=lambda pk: pk.height)
            summed = replace(
                tallest,
                area_min=sum(pk.area_min for pk in grp),
                area=sum(pk.area for pk in grp),
                start_t=min(pk.start_t for pk in grp),
                end_t=max(pk.end_t for pk in grp),
                start_i=min(pk.start_i for pk in grp),
                end_i=max(pk.end_i for pk in grp),
                code=tallest.code[:2] + "+",
            )
            peaks = [pk for pk in peaks if pk not in grp]
            peaks.append(summed)
            peaks.sort(key=lambda pk: pk.rt)
        return peaks

    # -- shoulder splitting (drop-line mode) --------------------------------
    def _split_shoulders(self, events, slope_fine, ss):
        """S1/S2: a shoulder is an inflection pair without a zero-slope
        crossing — i.e. a significant dip in |slope| on a flank. DROP mode
        splits the event with a drop line at the slope-dip location.
        SPEC-OPEN: dip-significance threshold taken as slope_sensitivity."""
        out = []
        for (s_k, a_k, e_k, fL, fR) in events:
            cuts = []
            # rising flank: local minima of slope between start and apex
            for lo, hi, sgn in ((s_k, a_k, +1), (a_k, e_k, -1)):
                seg = sgn * slope_fine[lo:hi]
                if len(seg) < 5:
                    continue
                for j in range(2, len(seg) - 2):
                    if seg[j] > 0 and seg[j] < seg[j - 2] and seg[j] <= seg[j + 2]:
                        left_max = seg[:j].max()
                        right_max = seg[j:].max()
                        if min(left_max, right_max) - seg[j] > ss:
                            cuts.append((lo + j, "front" if sgn > 0 else "rear"))
                            break  # one shoulder per flank (v1)
            if not cuts:
                out.append((s_k, a_k, e_k, fL, fR))
                continue
            cuts.sort()
            segs = []
            prev = s_k
            for ck, side in cuts:
                segs.append((prev, ck, side))
                prev = ck
            segs.append((prev, e_k, None))
            for si, (a, b, side) in enumerate(segs):
                apex = a_k if a <= a_k <= b else (a + b) // 2
                out.append((a, apex, b, fL if si == 0 else True,
                            fR if si == len(segs) - 1 else True))
        return out

    # -- baseline + measurement on RAW data ---------------------------------
    def _build_peaks(self, t, y, tb, events, power) -> list[Peak]:
        f = 1 << power
        n = len(y)
        p = self.p
        sm = max(1, p.anchor_smooth_pts)

        def anchor_y(i):
            lo, hi = max(0, i - sm // 2), min(n, i - sm // 2 + sm)
            return float(np.mean(y[lo:hi]))

        # map bunched idx -> raw idx (center of bunch)
        def raw(kb):
            return min(n - 1, kb * f + f // 2)

        # group into clusters: consecutive events chained by fused flags
        clusters, cur = [], []
        for ev in events:
            cur.append(ev)
            if not ev[4]:            # not fused to the right -> cluster ends
                clusters.append(cur)
                cur = []
        if cur:
            clusters.append(cur)

        # Baseline-envelope bounds (task 6). The accumulators only TRIGGER
        # detection; actual peak start/end are where the signal departs from /
        # rejoins the local baseline envelope (level+drift+noise band), which
        # is what ChemStation does. Envelope construction details are not
        # published (SPEC-OPEN); band constant ENV_C and rejoin count ENV_M
        # were tuned against real CDS 2.8 TwelveTone output (task 9 data).
        yb_s = _bunch(y, power)
        exp_b = max(2.0, (self.p.peak_width / (t[1] - t[0])) / f)  # expected width, bunched pts

        # Bunched-equivalent baseline noise from the RAW signal (see
        # _robust_noise); averaging f samples scales white noise by 1/sqrt(f).
        nb = _robust_noise(y) / np.sqrt(f) + 1e-12

        import os
        ENV_START_TOL = float(os.environ.get("CSINT_ENV_START_TOL", "0.1"))
        ENV_END_C = float(os.environ.get("CSINT_ENV_END_C", "0.1"))
        ENV_RISE_C = float(os.environ.get("CSINT_ENV_RISE_C", "3.0"))

        # Smoothed copy for envelope WALK decisions only (moving average ~
        # exp_b/3): raw bunched noise stalls the strict-descent walks on small
        # peaks; measurements still use raw data and yb_s anchors.
        Wsm = max(1, int(exp_b / 3))
        if Wsm > 1:
            kern = np.ones(Wsm) / Wsm
            ybw = np.convolve(yb_s, kern, mode="same")
            nbw = nb / np.sqrt(Wsm)
        else:
            ybw, nbw = yb_s, nb

        # ---- pass 1: envelope bounds per cluster ----
        binfo = []
        for ci, cl in enumerate(clusters):
            kb0, kb1 = cl[0][0], cl[-1][2]
            kb_floor = clusters[ci - 1][-1][2] + 1 if ci > 0 else 0
            kb_cap = clusters[ci + 1][0][0] - 1 if ci + 1 < len(clusters) else len(yb_s) - 1

            # START: walk back to the local minimum before the rising flank
            # (matches CDS BaselineStart on reference data to ~0.01 pA).
            lim0 = max(kb_floor, kb0 - int(3 * exp_b))
            while kb0 - 1 >= lim0 and ybw[kb0 - 1] < ybw[kb0] + ENV_START_TOL * nbw:
                kb0 -= 1
            # END: walk forward while the tail still descends; settle on the
            # LOCAL MINIMUM: stop when descent stalls (< ENV_END_C * nb over
            # trailing window) or the signal rises off the running minimum
            # (> ENV_RISE_C * nb, i.e. a valley before the next feature).
            # SPEC-OPEN: exact CDS end rule unknown; this reproduces reference
            # areas, not necessarily reference end times.
            lim1 = min(kb_cap, kb1 + int(12 * exp_b))
            W = max(3, int(exp_b / 2))
            kb = kb1
            k_min, v_min = kb1, ybw[kb1]
            while kb + 1 <= lim1:
                kb += 1
                if ybw[kb] < v_min:
                    k_min, v_min = kb, ybw[kb]
                if ybw[kb] - v_min > ENV_RISE_C * nbw:
                    break
                if kb - kb1 >= W and (ybw[kb - W] - ybw[kb]) < ENV_END_C * nbw:
                    break
            kb1 = k_min
            cl[0] = (kb0, *cl[0][1:])
            cl[-1] = (*cl[-1][:2], kb1, *cl[-1][3:])
            binfo.append([kb0, kb1])

        # ---- pass 2: chain clusters whose end never returned to baseline ----
        # A cluster ending far above its own start (rider on a tail) is not a
        # finished cluster — ChemStation keeps the baseline open. Chain it
        # with the following cluster; interior boundaries become valleys.
        i = 0
        while i + 1 < len(clusters):
            kb0, kb1 = binfo[i]
            # valley elevated above the OUTER anchors (this cluster's start,
            # next cluster's end) by more than the noise band -> not baseline
            # -> same cluster (CDS BV/VB behavior on partially fused pairs)
            outer = min(yb_s[kb0], yb_s[binfo[i + 1][1]])
            elevated = (yb_s[kb1] - outer) > 5 * nb
            gap_ok = (clusters[i + 1][0][0] - kb1) < 2 * exp_b
            if elevated and gap_ok:
                clusters[i] = clusters[i] + clusters[i + 1]
                binfo[i][1] = binfo[i + 1][1]
                del clusters[i + 1]
                del binfo[i + 1]
            else:
                i += 1

        # ---- pass 3: measure ----
        out = []
        for ci, (cl, (kb0, kb1)) in enumerate(zip(clusters, binfo)):
            c_start = raw(kb0)
            c_end = raw(kb1)
            y0, y1 = float(yb_s[kb0]), float(yb_s[kb1])
            t0, t1 = t[c_start], t[c_end]
            slope_base = (y1 - y0) / (t1 - t0) if t1 > t0 else 0.0

            def base_at(tt):
                return y0 + slope_base * (tt - t0)

            # valley (drop-line) raw indices between adjacent peaks
            bounds = [c_start]
            for a, b in zip(cl[:-1], cl[1:]):
                if a[2] == b[0]:
                    # events touch (valley split or shoulder cut): the shared
                    # index IS the drop-line position
                    bounds.append(raw(b[0]))
                    continue
                lo, hi = raw(a[1]), raw(b[1])
                if hi <= lo:
                    hi = lo + 1
                bounds.append(lo + int(np.argmin(y[lo:hi])))
            bounds.append(c_end)

            # merge degenerate splits: drop any interior bound whose smaller
            # adjacent apex rises less than ENV_RISE_C*nb above the bound
            changed = True
            while changed and len(cl) > 1:
                changed = False
                for bi in range(1, len(bounds) - 1):
                    a_apex = y[raw(cl[bi - 1][1])]
                    b_apex = y[raw(cl[bi][1])]
                    if min(a_apex, b_apex) - y[bounds[bi]] < ENV_RISE_C * nb:
                        winner = cl[bi - 1] if a_apex >= b_apex else cl[bi]
                        cl = cl[:bi - 1] + [winner] + cl[bi + 1:]
                        bounds = bounds[:bi] + bounds[bi + 1:]
                        changed = True
                        break

            # Penetration removal (S2 "No penetration" mode): if the signal
            # dips below the cluster baseline chord by more than the noise
            # band, pull the END anchor back to the deepest penetration point
            # and re-chord. Converges to a chord that hugs a convex tail —
            # this is what recovers rider peaks sitting on a parent tail.
            apex0 = raw(cl[0][1])
            apexN = raw(cl[-1][1])
            for _ in range(24):
                seg = y[c_start:c_end + 1] - base_at(t[c_start:c_end + 1])
                j = int(np.argmin(seg))
                if seg[j] >= -3 * nb or j <= 2 or c_start + j >= c_end - 2:
                    break
                if c_start + j < apex0:
                    # penetration on the leading side: shift START to the dip
                    c_start = c_start + j
                    y0 = float(y[c_start])
                    t0 = t[c_start]
                elif c_start + j > apexN:
                    # trailing side: pull END back to the dip
                    c_end = c_start + j
                    y1 = float(y[c_end])
                    t1 = t[c_end]
                else:
                    break  # dip between apexes = valley, not penetration
                slope_base = (y1 - y0) / (t1 - t0) if t1 > t0 else 0.0
                bounds = sorted(set([c_start] + [b for b in bounds
                                                if c_start < b < c_end] + [c_end]))
                cl = cl[:max(1, len(bounds) - 1)]

            cpeaks = []
            for pi, ev in enumerate(cl):
                if pi + 1 >= len(bounds):
                    break
                i0, i1 = bounds[pi], bounds[pi + 1]
                if i1 <= i0 + 2:
                    continue
                seg_t, seg_y = t[i0:i1 + 1], y[i0:i1 + 1]
                base = base_at(seg_t)
                corr = seg_y - base
                area_min = float(np.trapezoid(corr, seg_t))
                ia = i0 + int(np.argmax(corr))
                rt, apex_y = _parabolic_apex(t, y, ia)
                height = float(apex_y - base_at(rt))
                if height <= 0 or area_min <= 0:
                    continue
                # Prominence guard: the apex must rise above BOTH endpoint
                # signal levels. Kills recovery-flank artifacts (e.g. the
                # rising edge after a negative dip masquerading as a peak).
                if apex_y - max(float(seg_y[0]), float(seg_y[-1])) < p.height_reject:
                    continue
                left = "B" if pi == 0 else "V"
                right = "B" if pi == len(cl) - 1 else "V"
                cpeaks.append(Peak(
                    start_t=float(seg_t[0]), end_t=float(seg_t[-1]), rt=rt,
                    height=height, area=area_min * 60.0, area_min=area_min,
                    width_ah=area_min / height if height > 0 else 0.0,
                    code=left + right, cluster=ci,
                    start_i=i0, end_i=i1,
                    baseline_y0=float(base[0]), baseline_y1=float(base[-1]),
                ))
            self._apply_tail_skims(t, y, cpeaks, base_at, nb, exp_b, f)
            out.extend(cpeaks)
        return out

    # -- tail skims (S2 §7) --------------------------------------------------
    def _apply_tail_skims(self, t, y, cpeaks, base_at, nb, exp_b, f):
        """Skim rear children off a parent tail. Enabled when
        tail_skim_ratio > 0 (0 = off, ChemStation convention).

        Modes: STRAIGHT (chord over child, parent-slope corrected by
        construction since endpoints sit on the tail), EXP (fitted exponential
        decay of the parent tail), STANDARD (EXP clipped at the baseline —
        approximates the documented exp-above/straight-inside-envelope
        hybrid; SPEC-OPEN). Skimmed slab area moves child -> parent.
        Codes: child gets 'T' (STANDARD/STRAIGHT) or 'X' (EXP) appended."""
        p = self.p
        if p.tail_skim_ratio <= 0 or len(cpeaks) < 2:
            return
        parent_i = int(np.argmax([pk.height for pk in cpeaks]))
        parent = cpeaks[parent_i]
        for pk in cpeaks[parent_i + 1:]:
            if parent.height / max(pk.height, 1e-12) <= p.tail_skim_ratio:
                continue
            Hv = float(y[pk.start_i] - base_at(t[pk.start_i]))
            if Hv > 0 and pk.height / max(Hv, 1e-12) >= p.skim_valley_ratio:
                continue  # ratio LOWER than setting enables the skim (S2)
            i0, i1 = pk.start_i, pk.end_i
            seg_t = t[i0:i1 + 1]
            base = base_at(seg_t)
            ys = float(y[i0] - base[0])
            ye = float(y[i1] - base[-1])
            if ys <= 0:
                continue
            curve = None
            if p.skim_mode in ("EXP", "STANDARD"):
                w = max(8, int(exp_b * f / 2))
                j0 = max(parent.start_i, i0 - w)
                tail_t = t[j0:i0 + 1]
                tail_h = y[j0:i0 + 1] - base_at(tail_t)
                pos = tail_h > max(3 * nb, 1e-9)
                if pos.sum() >= 4:
                    Bfit = float(np.polyfit(tail_t[pos], np.log(tail_h[pos]), 1)[0])
                    if Bfit < 0:
                        curve = ys * np.exp(Bfit * (seg_t - seg_t[0]))
                        if p.skim_mode == "STANDARD":
                            curve = np.maximum(curve, 0.0)
            if curve is None:  # STRAIGHT or failed exp fit
                curve = np.linspace(ys, ye, len(seg_t))
            # skim curve may not exceed the signal (would create negative child)
            curve = np.minimum(curve, y[i0:i1 + 1] - base)
            slab = float(np.trapezoid(curve, seg_t))
            if slab <= 0:
                continue
            child_area = float(np.trapezoid(y[i0:i1 + 1] - base - curve, seg_t))
            pk.area_min = child_area
            pk.area = child_area * 60.0
            pk.height = float(pk.height - np.interp(pk.rt, seg_t, curve))
            pk.width_ah = pk.area_min / pk.height if pk.height > 0 else 0.0
            pk.code = pk.code + ("X" if p.skim_mode == "EXP" else "T")
            parent.area_min += slab
            parent.area = parent.area_min * 60.0
        parent.width_ah = parent.area_min / parent.height if parent.height > 0 else 0.0


def integrate(trace: Trace, params: Params) -> IntegrationResult:
    return Integrator(params).run(trace)
