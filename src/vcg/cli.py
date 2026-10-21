"""Command line entry point: ``python -m vcg.cli <step>``.

Steps
-----
prepare    raw collections -> data/processed/*.parquet
evaluate   processed tables -> results/tables/records.parquet (one row per
           study/sex/endpoint/dose group/configuration/draw)
summarise  records -> results/tables/*.csv, results/summary.json, results/figures/*.png
all        the three steps in order
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd

from . import config as C
from . import evaluation, plots, prepare

log = logging.getLogger("vcg")


def step_prepare(args: argparse.Namespace) -> None:
    prepare.build(C.RAW_DIR, C.PROCESSED_DIR, workers=args.workers)


def step_evaluate(args: argparse.Namespace) -> None:
    data = prepare.load_processed(C.PROCESSED_DIR)
    t0 = time.time()
    records = evaluation.run(data, workers=args.workers, n_draws=args.draws, seed=args.seed, limit=args.limit)
    C.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    records.to_parquet(C.TABLES_DIR / "records.parquet", index=False)
    log.info("%d records written in %.0f s", len(records), time.time() - t0)


def step_summarise(args: argparse.Namespace) -> None:
    records = pd.read_parquet(C.TABLES_DIR / "records.parquet")
    data = prepare.load_processed(C.PROCESSED_DIR)
    tables = evaluation.summarise(records)
    for name, table in tables.items():
        table.to_csv(C.TABLES_DIR / f"{name}.csv", index=False, float_format="%.4f")
    figures = plots.make_all(tables, C.FIGURES_DIR)

    lvl = "L3 +laboratory"
    by_level = tables["agreement_by_level"]
    head = by_level[(by_level["mode"] == "vcg") & (by_level["level"] == lvl)].iloc[0]
    ep = evaluation.endpoint_level(records)
    units = ep[(ep["band"] == "") & (ep["mode"] == "vcg") & (ep["level"] == lvl)]
    summary = {
        "headline_level": lvl,
        "agreement": round(float(head["agreement"]), 4),
        "false_positive": round(float(head["false_positive"]), 4),
        "false_negative": round(float(head["false_negative"]), 4),
        "fp_given_conc_nonsig": round(float(head["fp_given_conc_nonsig"]), 4),
        "fn_given_conc_sig": round(float(head["fn_given_conc_sig"]), 4),
        "sensitivity": round(float(head["sensitivity"]), 4),
        "specificity": round(float(head["specificity"]), 4),
        "kappa": round(float(head["kappa"]), 4),
        "trivial_agreement_never_significant": round(float(head["trivial_agreement"]), 4),
        "concurrent_significant_rate": round(float(head["conc_sig_rate"]), 4),
        "vcg_significant_rate": round(float(head["vcg_sig_rate"]), 4),
        "n_endpoint_units": int(head["n_endpoints"]),
        "n_studies": int(units["study"].nunique()),
        "n_study_sex_units": int(units[["study", "sex"]].drop_duplicates().shape[0]),
        "n_records": int(len(records)),
        "n_draws": int(records["draw"].max() + 1),
        "studies_in_processed_data": int(data["studies"].shape[0]),
        "by_level": {r["level"]: {"agreement": round(float(r["agreement"]), 4), "kappa": round(float(r["kappa"]), 4),
                                  "coverage_units": int(r["n_endpoints"])}
                     for _, r in by_level[by_level["mode"] == "vcg"].iterrows()},
        "partial_by_level": {r["level"]: {"agreement": round(float(r["agreement"]), 4), "kappa": round(float(r["kappa"]), 4)}
                             for _, r in by_level[by_level["mode"] == "partial"].iterrows()},
        "figures": [f.relative_to(C.ROOT).as_posix() for f in figures],
    }
    (C.RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    log.info("summary: %s", json.dumps(summary, indent=2))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="vcg", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "evaluate", "summarise", "all"])
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--draws", type=int, default=C.N_DRAWS, help="virtual control groups per configuration")
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--limit", type=int, default=None, help="only evaluate the first N study/sex units (smoke test)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    steps = {"prepare": step_prepare, "evaluate": step_evaluate, "summarise": step_summarise}
    for name in (["prepare", "evaluate", "summarise"] if args.step == "all" else [args.step]):
        log.info("== %s", name)
        steps[name](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
