"""csint — ChemStation-compatible chromatographic peak integrator (WIP).

Status: readers + synthetic-truth suite implemented; core integrator under
construction against docs/SPEC.md. Nothing here is validated against real
ChemStation output yet.
"""
from .trace import Trace
from .io import load

__all__ = ["Trace", "load"]
__version__ = "0.0.1"
