# csint — ChemStation-compatible peak integrator (library in this repo)

`csint/` is a validated reimplementation of the OpenLab CDS 2.8 TwelveTone
(ChemStation) integrator. USE THIS for any peak integration in this repo —
do not write ad-hoc integration code and do not extend peak_analyzer.py with
new integration logic.

## Validation status (2026-08-10, vs real CDS 2.8 output from GC-HSS1)

- Multi-peak (HITEMP MIX A, 6 inj): 18/18 peaks, |ΔArea| ≤ 0.06%,
  peak codes 18/18 (BB/BV/VB), ΔRT ≤ 0.11 s, blank 0/0.
- Single-peak (pentadione, 6 inj): 5/6 ≤ 0.01%; one −2.3% (CDS late-end
  rule, open item in docs/CSINT_SPEC.md).
- Synthetic exact-truth suite: 22/23; regression tests 6/6
  (`python tests/test_csint_integrator.py`).

Known limits: small peaks near noise ±15–25%; shoulders/front-skims/
autointegrate not implemented (our methods don't use them).

## Use

```python
from csint import load
from csint.integrator import Params, integrate
from csint.report import peak_table, table_rows, overlay

tr = load(path)[0]        # .csv/.h5/.abf/ANDI .cdf/Agilent .ch/.D
res = integrate(tr, Params(
    slope_sensitivity=1000,   # = proc-D.pmx SlopeSensitivity (pA/min)
    peak_width=0.1,           # min
    height_reject=1.7, area_reject=1.0,
    timed_events=[(0.05, "AreaSum", True), (0.8, "AreaSum", False)],
))
print(peak_table(res))    # ChemStation-style table; res.peaks for objects
```

CLI: `python -m csint FILE --ss 1000 --pw 0.1 --hr 1.7 --ar 1 --png o.png`

OpenLab .dx files are zips — read the .CH inside. Parameters map 1:1 to the
CDS processing method (proc-D.pmx values above are the lab default).
Skims: `tail_skim_ratio` (0=off). Negative peaks: `negative_peaks=True`.

## Rules for future sessions

- Algorithm changes require the gates: `python tests/test_csint_integrator.py`
  and `python csint_validation/run_tier1.py` must pass, and any change to
  detection/baseline behavior must be re-diffed against the stored CDS
  reference sets (C:\CDSProjects\Test\Results — pentadione + HITEMP MIX A,
  processed with proc-D) before claiming correctness.
- Ambiguities vs Agilent's spec are marked SPEC-OPEN in code and
  docs/CSINT_SPEC.md — resolve them with reference data, never by guessing.
