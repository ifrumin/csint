"""Unified chromatogram trace model."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Trace:
    """A single-channel chromatogram.

    t : time axis in MINUTES (ChemStation convention), strictly increasing.
    y : signal in native units (counts, pA, mAU ...), same length as t.
    meta : free-form provenance (source file, detector, sampling rate ...).
    """

    t: np.ndarray
    y: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.t = np.asarray(self.t, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        if self.t.ndim != 1 or self.y.ndim != 1 or len(self.t) != len(self.y):
            raise ValueError("t and y must be 1-D arrays of equal length")
        if len(self.t) < 8:
            raise ValueError("trace too short")
        dt = np.diff(self.t)
        if np.any(dt <= 0):
            raise ValueError("time axis must be strictly increasing")

    @property
    def dt(self) -> float:
        """Median sampling interval (minutes)."""
        return float(np.median(np.diff(self.t)))

    @property
    def rate_hz(self) -> float:
        return 1.0 / (self.dt * 60.0)

    def is_uniform(self, rtol: float = 1e-3) -> bool:
        d = np.diff(self.t)
        return bool(np.max(np.abs(d - self.dt)) <= rtol * self.dt)

    def resampled_uniform(self) -> "Trace":
        """Return a uniformly sampled copy (linear interp) if needed."""
        if self.is_uniform():
            return self
        n = len(self.t)
        tu = np.linspace(self.t[0], self.t[-1], n)
        yu = np.interp(tu, self.t, self.y)
        meta = dict(self.meta, resampled=True)
        return Trace(tu, yu, meta)

    def crop(self, t0: float | None = None, t1: float | None = None) -> "Trace":
        m = np.ones(len(self.t), dtype=bool)
        if t0 is not None:
            m &= self.t >= t0
        if t1 is not None:
            m &= self.t <= t1
        return Trace(self.t[m], self.y[m], dict(self.meta))
