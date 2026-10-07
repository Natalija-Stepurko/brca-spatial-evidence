"""`bse run <stage>`: the pipeline, one stage per module. Stages arrive one pull request at a time;
until a stage exists, asking for it says so and exits non-zero."""
import argparse
import importlib
import sys

STAGES = [
    ("data", "download every input with sha256 provenance; copy study 2 candidate tables at a pinned commit"),
    ("compartment", "spot labels, normalisation, tumour log-fold-change and fraction, both nulls; P1-P5"),
    ("tiles", "one tile per spot at 0.5 um/px; frozen encoder features cached"),
    ("histology", "PCA + ridge per gene, every section held out once, nulls and baselines; P6-P8"),
    ("report", "verdicts and figures for the project page"),
]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bse", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="run one stage")
    run.add_argument("stage", choices=[s for s, _ in STAGES])
    sub.add_parser("stages", help="list the stages")
    a = ap.parse_args(argv)
    if a.cmd == "stages":
        for s, d in STAGES:
            print(f"{s:<12} {d}")
        return 0
    try:
        mod = importlib.import_module(f"bse.stages.{a.stage}")
    except ModuleNotFoundError:
        print(f"stage {a.stage!r} is not implemented yet; see the open pull requests", file=sys.stderr)
        return 2
    return mod.run() or 0


if __name__ == "__main__":
    sys.exit(main())
