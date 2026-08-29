# csint — ChemStation-compatible GC peak integrator

Reimplementation of the Agilent ChemStation / OpenLab CDS "TwelveTone"
integrator, built from Agilent's published algorithm documentation and
multiple file formats as input, validated against real OpenLab CDS 2.8 output from an Agilent 8900 HSS-GC-FID
 (Headapace delivered via valve&loop sampler - Agilent 7697A)).

## Status (2026-08-10)

Validated:
- Tier-2 (vs real CDS 2.8 TwelveTone, pentadione set, identical raw data and
  parameters): 5/6 injections |ΔArea| ≤ 0.01%, ΔRT ≤ 0.08 s, ΔHeight ≤ 0.5%,
  peak codes match. 1/6 at −2.3% (CDS chose a later peak end; open item).
- Tier-1 (synthetic exact-truth suite, 7 cases): 22/23 peaks, 0 false
  positives; includes skimmed riders (parent −0.3%), fused pairs (±2%),
  drift, negative dip (BBN, +3.6%).

Known limits (documented, not hidden):
- Small peaks near the noise floor: area ±15–25%.
- Shoulder detection: not functional for tailing-peak shoulders (parked;
  Idan's method runs shoulders OFF).
- Front skims, autointegrate, baseline-hold/now/next-valley timed events:
  not yet implemented.
- inj02-style late-end rule and several SPEC-OPEN constants: see docs/SPEC.md.

## Usage

```
python -m csint FILE --ss 1000 --pw 0.1 --hr 1.7 --ar 1 \
       --png overlay.png --json peaks.json
```

FILE FORMATS: .csv, .h5, .abf, ANDI .cdf, Agilent .ch, or a .D directory.
Parameters mirror ChemStation: --ss slope sensitivity (units/min),
--pw initial peak width (min), --hr height reject, --ar area reject (units·s).

Python API:

```python
from csint import load
from csint.integrator import Params, integrate
from csint.report import peak_table

tr = load("run.dx-extracted.CH")[0]      # or .csv/.abf/.cdf/...
res = integrate(tr, Params(slope_sensitivity=1000, peak_width=0.1,
                           height_reject=1.7, area_reject=1.0))
print(peak_table(res))
```

Skims: `tail_skim_ratio=3` (0=off, ChemStation convention), modes
STANDARD/EXP/STRAIGHT. Timed events: `timed_events=[(0.05,"AreaSum",True),
(0.8,"AreaSum",False)]`.

Note: OpenLab .dx files are zip containers — the .CH inside is what csint
reads (`unzip run.dx`), or pass the .D directory for classic ChemStation data.

## Layout

- `csint/` — package (io, integrator, synth, report, CLI)
- `docs/SPEC.md` — algorithm spec digest + all SPEC-OPEN items
- `validation/` — tier-1 synthetic suite + tier-2 comparison artifacts
- `tests/` — IO round-trip tests

## Validation protocol

Tier-1: `python validation/run_tier1.py` — synthetic chromatograms with
analytically exact areas. Tier-2: comparison against CDS-produced peak tables
(.rx InjectionACAML) on identical raw data. No accuracy claim is made from
tier-1 alone; "matches ChemStation" claims come only from tier-2.
