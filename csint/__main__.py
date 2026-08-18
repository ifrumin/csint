"""CLI: python -m csint FILE [options]

Integrates a chromatogram (.csv/.h5/.abf/.cdf/.ch/.D) with ChemStation-style
parameters and prints the peak table; optionally writes an overlay PNG and a
CSV table.
"""
import argparse
import json
import sys

from . import load
from .integrator import Params, integrate
from .report import peak_table, table_rows, overlay


def main(argv=None):
    ap = argparse.ArgumentParser(prog="csint")
    ap.add_argument("file")
    ap.add_argument("--ss", type=float, required=True,
                    help="slope sensitivity (signal units/min)")
    ap.add_argument("--pw", type=float, default=0.05,
                    help="initial peak width, minutes")
    ap.add_argument("--hr", type=float, default=0.0, help="height reject")
    ap.add_argument("--ar", type=float, default=0.0, help="area reject (units*s)")
    ap.add_argument("--shoulders", choices=["OFF", "DROP"], default="OFF")
    ap.add_argument("--channel", type=int, default=0,
                    help="channel index if the file holds several")
    ap.add_argument("--png", help="write overlay image here")
    ap.add_argument("--json", help="write peak table as JSON here")
    a = ap.parse_args(argv)

    traces = load(a.file)
    tr = traces[a.channel]
    res = integrate(tr, Params(slope_sensitivity=a.ss, peak_width=a.pw,
                               height_reject=a.hr, area_reject=a.ar,
                               shoulders=a.shoulders))
    print(f"# {a.file} ch={a.channel} ({tr.meta.get('channel','')}) "
          f"{len(tr.t)} pts {tr.rate_hz:.0f} Hz")
    print(peak_table(res))
    if a.png:
        overlay(res, a.png)
        print(f"overlay -> {a.png}", file=sys.stderr)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(table_rows(res), fh, indent=1)
        print(f"table -> {a.json}", file=sys.stderr)


if __name__ == "__main__":
    main()
