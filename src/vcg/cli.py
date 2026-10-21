"""Command line entry point: ``python -m vcg.cli <step>``.

Steps
-----
prepare    raw collections -> data/processed/*.parquet
evaluate   processed tables -> results/tables/records.parquet (one row per
           study/sex/endpoint/dose group/configuration/draw)
summarise  records -> results/tables/*.csv, results/summary.json, results/figures/*.png
methods    processed tables -> results/tables/method_endpoints.parquet: the
           pre-specified method comparison (docs/analysis_plan.md), one row per
           study/sex/endpoint/method/mode with the draws averaged
posthoc    post-hoc control analysis, not in the plan: the kept concurrent half
           alone -> results/tables/method_endpoints_posthoc.parquet
compare    method_endpoints -> results/tables/methods_*.csv, results/methods_summary.json,
           results/figures/methods_comparison.png
all        all six steps in order
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
from . import evaluation, methods, plots, prepare

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


def step_methods(args: argparse.Namespace) -> None:
    data = prepare.load_processed(C.PROCESSED_DIR)
    t0 = time.time()
    df = methods.run(data, workers=args.workers, n_draws=args.draws, seed=args.seed, limit=args.limit)
    C.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    df.to_parquet(C.TABLES_DIR / "method_endpoints.parquet", index=False)
    log.info("%d endpoint rows written in %.0f s", len(df), time.time() - t0)


def step_posthoc(args: argparse.Namespace) -> None:
    """Post-hoc control analysis, not in the plan: the kept concurrent half alone."""
    data = prepare.load_processed(C.PROCESSED_DIR)
    df = methods.run(data, workers=args.workers, n_draws=args.draws, seed=args.seed, limit=args.limit, posthoc=True)
    df.to_parquet(C.TABLES_DIR / "method_endpoints_posthoc.parquet", index=False)
    log.info("%d post-hoc endpoint rows written", len(df))


def step_compare(args: argparse.Namespace) -> None:
    df = pd.read_parquet(C.TABLES_DIR / "method_endpoints.parquet")
    posthoc_path = C.TABLES_DIR / "method_endpoints_posthoc.parquet"
    if posthoc_path.exists():
        cats = ["set", "study", "sex", "species", "design", "domain", "endpoint", "method", "mode"]
        df = pd.concat([df, pd.read_parquet(posthoc_path)], ignore_index=True)
        df[cats] = df[cats].astype(str).astype("category")
    tables = methods.compare(df, seed=args.seed)
    for name, table in tables.items():
        table.to_csv(C.TABLES_DIR / f"{name}.csv", index=False, float_format="%.4f")
    variance = methods.between_study_variance(prepare.load_processed(C.PROCESSED_DIR))
    variance.to_csv(C.TABLES_DIR / "between_study_variance.csv", index=False, float_format="%.5g")
    by_domain = variance.groupby("domain").agg(strata_endpoints=("icc", "size"), median_icc=("icc", "median"),
                                               q25_icc=("icc", lambda x: x.quantile(0.25)),
                                               q75_icc=("icc", lambda x: x.quantile(0.75))).reset_index()
    by_domain["implied_fp_individual"] = methods.implied_false_positive_rate(by_domain["median_icc"])
    by_domain["implied_fp_whole_group"] = methods.implied_false_positive_rate(by_domain["median_icc"],
                                                                              sampling="whole_group")
    by_domain.to_csv(C.TABLES_DIR / "between_study_variance_by_domain.csv", index=False, float_format="%.4f")
    fig = plots.method_comparison(tables["methods_vs_baseline"], C.FIGURES_DIR / "methods_comparison.png")
    vs = tables["methods_vs_baseline"]
    summary = {
        "cutoff_year": C.CUTOFF_YEAR,
        "matching_level": C.METHOD_LEVEL,
        "units_by_set": df.drop_duplicates(["set", "study", "sex"]).groupby("set", observed=True).size().to_dict(),
        "studies_by_set": df.drop_duplicates(["set", "study"]).groupby("set", observed=True).size().to_dict(),
        "verdicts_holdout": {f"{r.method}/{r.mode}/{r.domain}": f"{r.verdict} (kappa {r.kappa_diff:+.3f}, "
                                                                f"95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f})"
                             for r in vs[vs["set"] == "holdout"].itertuples()},
        "figure": fig.relative_to(C.ROOT).as_posix(),
    }
    (C.RESULTS_DIR / "methods_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("method comparison written: %d paired comparisons", len(vs))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="vcg", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "evaluate", "summarise", "methods", "posthoc", "compare", "all"])
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--draws", type=int, default=C.N_DRAWS, help="virtual control groups per configuration")
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--limit", type=int, default=None, help="only evaluate the first N study/sex units (smoke test)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    steps = {"prepare": step_prepare, "evaluate": step_evaluate, "summarise": step_summarise,
             "methods": step_methods, "posthoc": step_posthoc, "compare": step_compare}
    for name in (list(steps) if args.step == "all" else [args.step]):
        log.info("== %s", name)
        steps[name](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
