# ChemStation integrator — implementation spec digest

Compiled 2026-08-10 from Agilent public documentation. Target: reproduce the
OpenLab ChemStation ("new") integrator behavior. Items marked **[OPEN]** are
uncertain and must be resolved by testing against real ChemStation output or a
better source — they are NOT to be silently guessed.

## Sources

- S1: Agilent "Understanding Your ChemStation" G2070-91126
  https://www.agilent.com/cs/library/usermanuals/public/G2070-91126_Understanding.pdf
- S2: Agilent OpenLab ChemStation Data Analysis Reference CS-LTS 01.11
  https://www.agilent.com/cs/library/usermanuals/public/CS-LTS_01.11_Reference_en.pdf
- S3: Agilent e-seminar "Integrating My Results in ChemStation" (2008)
  https://www.agilent.com/cs/library/eseminars/Public/Integrating%20My%20Results%20in%20ChemStation_121008.pdf

## 1. Sampling / bunching (S2, Table 1)

Data points are bunched in powers of two so the expected peak width stays in
the working range of the recognition filters:

| Expected peak width (data pts) | Filter(s) used | Bunching |
|---|---|---|
| 0–10   | First         | none |
| 8–16   | Second        | none |
| 12–24  | Third         | none |
| 16–32  | Second        | once (2x) |
| 24–48  | Third         | once (2x) |
| 32–96  | Third, second | twice (4x) |
| 64–192 | Third, second | three times (8x) |

Ranges overlap intentionally (hysteresis). Bunched value = sum/average of the
bunch **[OPEN: sum vs average — affects slope-sensitivity scaling; resolve
empirically]**. Pattern extends by powers of 2 for wider peaks.

## 2. Peak recognition filters (S2)

First and second derivative over bunched points:
- Filter 1: slope from 2 contiguous points; curvature from 3 contiguous points.
- Filter 2: slope from 4 contiguous points; curvature from 3 non-contiguous points.
- Filter 3: slope from 8 contiguous points; curvature from 3 non-contiguous points.

**[OPEN]** Exact convolution coefficients are not published. Implement as
first/second central differences over the stated spans; verify scaling against
real ChemStation results.

Filter outputs are compared against **slope sensitivity** (user parameter).
Units: signal units per unit time **[OPEN: exact normalization; S3 gives only
practical values ~50 → 20 for typical FID work]**.

## 3. Up/down-slope accumulators (S2, Tables 2–3)

A peak START candidate is declared when the up-slope accumulator reaches ≥ 15.
A peak END candidate is declared when the down-slope accumulator reaches ≥ 15.

Up-slope accumulator increments:

| Condition | F1 | F2 | F3 |
|---|---|---|---|
| slope > SS            | +8 | +5 | +3 |
| curvature > SS        | +0 | +2 | +1 |
| slope < −SS           | −8 | −5 | −3 |
| |slope| < SS (flat)   | −4 | −2 | −1 |
| curvature < −SS       | −0 | −2 | −1 |

Down-slope accumulator increments:

| Condition | F1 | F2 | F3 |
|---|---|---|---|
| slope < −SS           | +8 | +5 | +3 |
| curvature < −SS       | +0 | +2 | +1 |
| slope > SS            | −11 | −7 | −4 |
| |slope| < SS (flat)   | −28 | −18 | −11 |
| curvature > SS        | −0 | −2 | −1 |

**[OPEN]** The rows rendered here as "|slope| < SS" were extracted as
"Slope < |Slope Sensitivity|" / "Slope > |Slope Sensitivity|" — label/sign
ambiguity from PDF extraction. The −28/−18/−11 row in particular needs
verification (a strong reset when signal flattens is plausible, but which
condition triggers it is uncertain). Verify against real data before trusting
merged-peak behavior.

## 4. Cardinal points (S1, S2)

- Start: up-slope accumulator hits threshold; start point is set at/back-dated
  to where the envelope was left **[OPEN: exact back-dating rule]**.
- Apex: parabolic fit through the highest data points.
- End: down-slope accumulator hits threshold and slope/curvature return within
  limits.
- Shoulders (when enabled): from curvature (2nd derivative); an inflection
  without a zero-slope crossing → shoulder. Modes: drop-line shoulders or
  tangent shoulders (S3).

## 5. Peak width tracking (S2)

- GC: width = area/height.
- LC/CE: width = 0.3×(right inflection − left inflection) + 0.7×(area/height).
- Running update: new = 0.75×existing + 0.25×current-peak width.
- Width can be reset/fixed by timed events (Fixed/Auto Peak Width).

## 6. Baseline allocation (S1, S2)

- Baseline is re-allocated during the run at a frequency set by peak width;
  tracks drift via a baseline envelope.
- Unresolved cluster separation: **peak valley ratio** parameter.
  With peak heights H1, H2 (baseline-corrected) and valley height Hv:
  ratio = H2/Hv if H1 ≥ H2 else H1/Hv. If ratio < user value → drop line;
  else → valley baseline.
- Baseline penetration: "No penetration" mode shifts peak start/end until no
  penetration remains; "Advanced" mode re-baselines the cluster.
- Default baseline mode: baseline-to-baseline across a fused cluster with drop
  lines at valleys, unless valley baselines forced by events.

## 7. Solvent peak + skims (S2, S3)

Solvent peak detected by slope threshold (units mV/s); trailing riders are
auto-skimmed. Skim modes:
- Standard: hybrid — exponential fit while signal well above baseline,
  straight line inside the baseline envelope.
- Exponential (old): skim curve Hb(t) = H0·exp(−B·(t−t0)) + A·t + C
  (H0 height above baseline at skim start t0; B decay; A parent baseline
  slope; C offset). Curve passes under each child; area under curve goes to
  parent.
- New exponential: one exponential approximates the parent edge; children
  after the first are separated by drop lines from the first child's end.
- Straight: line start→end of child, corrected for parent slope.

Skim criteria (peak qualifies as child/rider):
- Tail peak skim height ratio: Hp/Hc > threshold → skim (typ. 3).
- Front peak skim height ratio: Hp/Hc > threshold → skim (typ. 6).
- Skim valley ratio: Hc/Hv < threshold → skim (typ. 20).

## 8. Peak separation codes (S2)

Char 1 = start, char 2 = end: B baseline, P penetration, V valley, H forced
horizontal, F forced point, M manual, U unassigned.
Char 3 flag: A aborted, D distorted, blank normal.
Char 4 type: S solvent, N negative, + area-summed, T tangent (standard),
X old-exponential skim, E new-exponential skim, R recalculated parent,
f/b front/rear shoulder tangent, F/B front/rear shoulder drop-line,
m/n/t/x manual variants, U unassigned.

## 9. Timed events (S2)

Initial events: slope sensitivity, peak width, area reject, height reject,
shoulders on/off.
Timed: Area reject, Area sum on/off, Area sum slice, Auto/Fixed peak width,
Baseline at valleys, Baseline backwards, Baseline hold, Baseline next valley,
Baseline now, Detect shoulders, Height reject, Integration on/off, Max area,
Max height, Negative peak on/off, Set (low) baseline from range, Slope
sensitivity, Solvent peak, Split peak, Tail tangent skim on/off, Tangent skim
mode, Unassigned peaks, Update peak height, Use baseline from range.

## 10. Autointegrate (S2)

1. Examine first and last 1% of the chromatogram; noise = 3×sd of linear
   regression ÷ sqrt(percent points used); also take slope.
2. Temporary height reject & slope sensitivity from noise; temporary peak
   width = 0.5% (LC) or 0.3–0.2% (GC) of run time.
3. Area reject 0; trial integration; adjust and repeat until ≥5 peaks or
   height reject reaches 0; max 10 trials.
4. Peak width refined from detected peaks (biased to early peaks), using only
   peaks with symmetry 0.8–1.3.
5. Refine height reject & slope sensitivity from inter-peak baseline; area
   reject = 90% of the min area of the most symmetric peak.
6. Re-integrate with final values.

## 11. Known gaps / classic-vs-new

- User's ChemStation revision unknown → which integrator generated their
  reports is unresolved. The classic (pre-A.06) integrator differs in skim
  and shoulder handling. **[OPEN — ask Idan; blocked-on-user]**
- Slope sensitivity normalization vs bunching **[OPEN]**.
- Exact filter kernels **[OPEN]**.
- Start-point back-dating rule **[OPEN]**.
All OPEN items are testable against paired .D + report data once provided.

## 12. Reference-tuning findings (2026-08-10, CDS 2.8 TwelveTone, GC-HSS1 FID)

Method parameters recovered from proc-D.pmx (TransformationChains/DefaultChain,
UTF-16 XAML with HTML-escaped inner XML): SS=1000, PW=0.1 min, AreaReject=1,
HeightReject=1.7, ShouldersMode=OFF, AreaPercentReject=0, AreaSum 0.05→0.8 min,
TangentSkimMode=Standard, Front/TailSkimHeightRatio=0 (=off), SkimValleyRatio=20,
PeakToValleyRatio=500, BaselineCorrection=Advanced.

Empirically established against 6 reference injections (pentadione, 100 Hz FID):

- Peak START = local minimum of the (smoothed) bunched trace immediately
  before the rising flank. Matches CDS BeginTime/BaselineStart to ~0.005 min /
  ~0.05 pA. Implemented as strict-descent walk-back (tol 0.1·noise).
- Peak END = local minimum found by walking the descending tail until it
  stalls or rises (valley). Matches CDS EndTime on 5/6 injections to
  ~0.005 min. One injection (inj02): CDS continued 0.42 min further to a
  later, lower minimum — rule for choosing between minima UNRESOLVED
  [OPEN: inj02-late-end].
- Walk decisions must run on a smoothed copy (moving avg ~exp_width/3);
  raw-noise stalls the walks on small peaks. Measurements stay on raw data.
- Noise estimator: per-window (64 pt) MAD of linear-detrended raw residuals,
  25th percentile across windows. Measured on inj01: 0.008 pA — the FID
  baseline genuinely is that quiet; earlier "quantization artifact" suspicion
  was wrong in magnitude. Estimator validated on noisy synthetic (sd=2.0)
  without gate regression. [RESOLVED 2026-08-10]
- Penetration handling: shift START toward dips before the first apex, END
  toward dips after the last apex (both directions needed; end-only collapses
  clusters when the chord sits high).
- CDS results storage: peak tables live in the .rx (zip) → Base/InjectionACAML
  <Peak> elements, val attributes, Area in pA·s; written only after
  Data Analysis processes AND saves. Unprocessed results carry
  TransformationChainState=NoMethodProvided.
- Verified unusable reference sets (FID data byte-level zero): Dec-2025
  testmix series, Jan-02 ETH-TIG series; salicylaldehyde 174902 aborted
  (6 KB .CH, zero delta_time).
