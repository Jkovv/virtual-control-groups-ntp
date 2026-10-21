"""Figures for the results section.  All functions take the summary tables
produced by :func:`vcg.evaluation.summarise` and write PNG files."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import config as C  # noqa: E402

DOMAIN_LABELS = {
    "terminal_bw": "Terminal body weight",
    "inlife_bw": "In-life body weight",
    "organ_weight": "Organ weights",
    "clin_chem": "Clinical chemistry",
    "hematology": "Hematology",
    "nonneoplastic lesion": "Non-neoplastic lesions",
    "neoplasm": "Neoplasms (poly-3)",
}
DOMAIN_ORDER = list(DOMAIN_LABELS)
LEVEL_ORDER = list(C.MATCH_LEVELS)
LEVEL_SHORT = {"L0 design": "design", "L1 +route": "+route", "L2 +vehicle": "+vehicle",
               "L3 +laboratory": "+lab", "L4 +period (3y)": "+period"}

# neutral palette, readable in print
COLORS = {"agree": "#3b6ea5", "fp": "#d1495b", "fn": "#edae49", "partial": "#8fb8de", "grey": "#9a9a9a"}

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def _pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def headline(tables: dict[str, pd.DataFrame], out: Path, level: str = "L3 +laboratory") -> Path:
    """Agreement with the concurrent-control conclusion by endpoint domain, full
    versus partial VCG, with the FP/FN decomposition of the disagreements."""
    t = tables["agreement_by_domain_level"]
    t = t[t["level"] == level]
    full = t[t["mode"] == "vcg"].set_index("domain").reindex(DOMAIN_ORDER).dropna(subset=["agreement"])
    part = t[t["mode"] == "partial"].set_index("domain").reindex(full.index)
    overall = tables["agreement_by_level"]
    ov_full = overall[(overall["mode"] == "vcg") & (overall["level"] == level)].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.35, 1]})
    y = np.arange(len(full))
    ax = axes[0]
    ax.barh(y + 0.2, full["agreement"], height=0.38, color=COLORS["agree"], label="Virtual control group")
    ax.barh(y - 0.2, part["agreement"], height=0.38, color=COLORS["partial"], label="Partial (half concurrent)")
    ax.scatter(full["trivial_agreement"], y, marker="|", s=160, color="k", zorder=3,
               label='Baseline: always "not significant"')
    for yi, (a, n, kp) in enumerate(zip(full["agreement"], full["n_endpoints"], full["kappa"])):
        ax.text(1.0, yi + 0.2, f"{_pct(a)}  n={int(n)}  kappa={kp:.2f}", va="center", fontsize=7.5)
    ax.set_yticks(y, [DOMAIN_LABELS.get(d, d) for d in full.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.5)
    ax.axvline(ov_full["agreement"], color="k", lw=0.8, ls="--")
    ax.set_xlabel("Share of endpoints with the same significance conclusion (alpha = 0.05)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=3, fontsize=7.5, frameon=False)
    ax.set_ylim(len(full) - 0.3, -0.7)
    ax.text(ov_full["agreement"] + 0.02, -0.55, f"all endpoints {_pct(ov_full['agreement'])}, kappa {ov_full['kappa']:.2f}",
            fontsize=7.5, va="center")
    ax.set_title(f"Agreement with concurrent control, matching {LEVEL_SHORT.get(level, level)}", fontsize=9.5, loc="left")

    ax = axes[1]
    ax.barh(y, full["false_positive"], height=0.6, color=COLORS["fp"], label="VCG significant, concurrent not")
    ax.barh(y, full["false_negative"], height=0.6, left=full["false_positive"], color=COLORS["fn"],
            label="Concurrent significant, VCG not")
    ax.set_yticks(y, [])
    ax.invert_yaxis()
    ax.set_xlabel("Share of endpoints (disagreements)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1, fontsize=7.5, frameon=False)
    ax.set_title("Where the disagreements come from", fontsize=9.5, loc="left")
    fig.suptitle("Virtual control groups built from NTP historical controls versus the real concurrent control",
                 fontsize=11, x=0.01, ha="left", y=1.02)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def strictness(tables: dict[str, pd.DataFrame], out: Path) -> Path:
    """Agreement, FP, FN and coverage as the matching gets stricter."""
    t = tables["agreement_by_domain_level"]
    cov = tables["coverage_by_level"].set_index("level").reindex(LEVEL_ORDER)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    x = np.arange(len(LEVEL_ORDER))
    for dom in DOMAIN_ORDER:
        d = t[(t["mode"] == "vcg") & (t["domain"] == dom)].set_index("level").reindex(LEVEL_ORDER)
        if d["agreement"].notna().sum() == 0:
            continue
        axes[0].plot(x, d["agreement"], marker="o", ms=4, label=DOMAIN_LABELS.get(dom, dom))
        axes[1].plot(x, d["false_positive"], marker="o", ms=4)
        axes[1].plot(x, d["false_negative"], marker="s", ms=3, ls=":", color=axes[1].lines[-1].get_color())
    axes[2].bar(x, cov["coverage"], color=COLORS["grey"])
    for ax, title in zip(axes, ["Agreement", "False positive (o) / false negative (dotted)", "Coverage of endpoints"]):
        ax.set_xticks(x, [LEVEL_SHORT[l] for l in LEVEL_ORDER], rotation=20)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_ylim(0, 1.02)
    axes[0].legend(fontsize=7, frameon=False)
    axes[2].set_title("Coverage: share of study/sex/endpoint units with a pool", fontsize=9, loc="left")
    fig.suptitle("Effect of matching strictness (each level adds one criterion)", fontsize=11, x=0.01, ha="left", y=1.04)
    fig.subplots_adjust(wspace=0.3)
    fig.savefig(out)
    plt.close(fig)
    return out


def drift(tables: dict[str, pd.DataFrame], out: Path) -> Path:
    """Agreement versus the time gap between the historical and the current study."""
    t = tables["agreement_by_drift_domain"]
    t = t[t["mode"] == "vcg"]
    bands = [f"{lo}-{hi}" for lo, hi in C.DRIFT_BANDS]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    x = np.arange(len(bands))
    for dom in DOMAIN_ORDER:
        d = t[t["domain"] == dom].set_index("band").reindex(bands)
        if d["agreement"].notna().sum() < 2:
            continue
        axes[0].plot(x, d["agreement"], marker="o", ms=4, label=DOMAIN_LABELS.get(dom, dom))
        axes[1].plot(x, d["false_positive"], marker="o", ms=4)
    all_ = tables["agreement_by_drift"]
    all_ = all_[all_["mode"] == "vcg"].set_index("band").reindex(bands)
    axes[0].plot(x, all_["agreement"], color="k", lw=2, label="all endpoints")
    axes[1].plot(x, all_["false_positive"], color="k", lw=2)
    for ax, title in zip(axes, ["Agreement", "False positive rate"]):
        ax.set_xticks(x, [b.replace("-", " to ") + " y" for b in bands])
        ax.set_xlabel("Years between historical control study and current study")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_ylim(0, 1.02)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"Drift: pools restricted to a time-gap band (matching {LEVEL_SHORT[C.DRIFT_LEVEL]})",
                 fontsize=11, x=0.01, ha="left")
    fig.savefig(out)
    plt.close(fig)
    return out


def endpoints(tables: dict[str, pd.DataFrame], out: Path, level: str = "L3 +laboratory",
              min_units: int = 30, top: int = 12) -> Path:
    """Most and least suitable endpoints, ranked by agreement."""
    t = tables["agreement_by_endpoint"]
    t = t[(t["level"] == level) & (t["n_endpoints"] >= min_units)].sort_values("agreement")
    if t.empty:
        t = tables["agreement_by_endpoint"]
        t = t[t["level"] == level].sort_values("agreement")
    low, high = t.head(top), t.tail(top)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    fig.subplots_adjust(wspace=1.0)
    for ax, d, title in zip(axes, [high.iloc[::-1], low], ["Highest agreement", "Lowest agreement"]):
        y = np.arange(len(d))
        ax.barh(y, d["agreement"], color=COLORS["agree"])
        ax.barh(y, d["false_positive"], left=d["agreement"], color=COLORS["fp"])
        ax.barh(y, d["false_negative"], left=d["agreement"] + d["false_positive"], color=COLORS["fn"])
        labels = [f"{DOMAIN_LABELS.get(dom, dom).split(' (')[0]}: {e}"[:44] for dom, e in zip(d["domain"], d["endpoint"])]
        ax.set_yticks(y, labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_title(f"{title} (n >= {min_units} study/sex units)", fontsize=10, loc="left")
        for yi, (a, n) in enumerate(zip(d["agreement"], d["n_endpoints"])):
            ax.text(0.01, yi, f"{_pct(a)} n={int(n)}", va="center", fontsize=6.5, color="white")
    fig.suptitle(f"Endpoint suitability for virtual controls (matching {LEVEL_SHORT.get(level, level)}); "
                 "blue agreement, red VCG-only significant, yellow concurrent-only significant",
                 fontsize=9.5, x=0.01, ha="left", y=1.02)
    fig.savefig(out)
    plt.close(fig)
    return out


def make_all(tables: dict[str, pd.DataFrame], fig_dir: Path = C.FIGURES_DIR) -> list[Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    return [headline(tables, fig_dir / "headline_agreement.png"),
            strictness(tables, fig_dir / "agreement_vs_strictness.png"),
            drift(tables, fig_dir / "agreement_vs_time_gap.png"),
            endpoints(tables, fig_dir / "endpoint_suitability.png")]
