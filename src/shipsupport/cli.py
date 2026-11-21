"""Command line: shipsupport run CASE.yaml [-o OUT] [--methods rigid fem] [--no-plots]"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .analysis import METHODS, run_case
from .loader import load_case
from .model import InputError
from .solvers import ContactError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="shipsupport",
                                     description="Dry-docking block reaction analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="analyse a case file")
    run.add_argument("case", type=Path, help="YAML case file")
    run.add_argument("-o", "--out", type=Path, default=None,
                     help="output directory (default: <case name>_results next to the case)")
    run.add_argument("--methods", nargs="+", choices=METHODS, default=None,
                     help="override the methods in the case file")
    run.add_argument("--no-plots", action="store_true", help="skip PNG plots")
    run.add_argument("--strict", action="store_true",
                     help="exit with code 2 if any check fails")
    args = parser.parse_args(argv)

    from . import report

    # the summary uses characters such as Σ, δ and →; never fail on a narrow console encoding
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")

    try:
        case = load_case(args.case)
        result = run_case(case, args.methods)
    except (InputError, ContactError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    summary = report.format_summary(result)
    print(summary)

    out = args.out or args.case.with_name(f"{args.case.stem}_results")
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.txt").write_text(summary, encoding="utf-8")
    written = [out / "summary.txt", *report.write_csv(result, out)]
    if not args.no_plots:
        try:
            written.append(report.plot_results(result, out / "results.png"))
            written.append(report.plot_block_forces(result, out / "block_forces.png"))
            for extra in (report.plot_welding_data(result, out / "welding_data.png"),
                          report.plot_spacing(result, out / "spacing_design.png")):
                if extra:
                    written.append(extra)
        except ImportError:
            print("note: matplotlib not installed, plots skipped", file=sys.stderr)
    print("\nOutput files:")
    for p in written:
        print(f"  {p}")

    failed = any(s.passed is False for s in result.strength.values()) or \
        (result.spacing is not None and not result.spacing.fits)
    return 2 if args.strict and failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
