"""
analysis.py
Runs Parts 2-4 and writes all outputs to outputs/.
Usage: python analysis.py
"""

import os
import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

DB_PATH = "cell_counts.db"
OUTPUT_DIR = "outputs"
CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

FLAT_QUERY = """
    SELECT
        sa.sample_id                        AS sample,
        sub.subject_id                      AS subject,
        sub.project_id                      AS project,
        sub.condition,
        sub.age,
        sub.sex,
        sub.treatment,
        sa.sample_type,
        sa.time_from_treatment_start,
        COALESCE(sa.response, '')           AS response,
        cc.b_cell,
        cc.cd8_t_cell,
        cc.cd4_t_cell,
        cc.nk_cell,
        cc.monocyte
    FROM cell_counts cc
    JOIN samples  sa  ON sa.sample_id   = cc.sample_id
    JOIN subjects sub ON sub.subject_id = sa.subject_id
"""


def load_data() -> pd.DataFrame:
    if not os.path.exists(DB_PATH):
        raise FileNotFoundError(f"'{DB_PATH}' not found. Run load_data.py first.")
    with sqlite3.connect(DB_PATH) as conn:
        df = pd.read_sql(FLAT_QUERY, conn)
    df["time_from_treatment_start"] = df["time_from_treatment_start"].astype(int)
    return df


def compute_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["total_count"] = df[CELL_TYPES].sum(axis=1)
    freq = df.melt(
        id_vars=[c for c in df.columns if c not in CELL_TYPES],
        value_vars=CELL_TYPES,
        var_name="cell_type",
        value_name="cell_count",
    )
    freq["percentage"] = (freq["cell_count"] / freq["total_count"] * 100).round(4)
    return freq


def run_part2(df: pd.DataFrame) -> pd.DataFrame:
    print("Part 2: computing cell-type frequencies...")
    freq = compute_frequencies(df)
    out = os.path.join(OUTPUT_DIR, "frequency_summary.csv")
    freq.to_csv(out, index=False)
    print(f"  {len(freq)} rows written to {out}")
    return freq


def run_part3(freq: pd.DataFrame) -> pd.DataFrame:
    print("Part 3: statistical analysis (melanoma / miraclib / PBMC)...")
    subset = freq[
        (freq["condition"] == "melanoma")
        & (freq["treatment"] == "miraclib")
        & (freq["sample_type"] == "PBMC")
        & (freq["response"].isin(["yes", "no"]))
    ].copy()

    results = []
    for ct in CELL_TYPES:
        ct_data = subset[subset["cell_type"] == ct]
        yes_vals = ct_data.loc[ct_data["response"] == "yes", "percentage"].values
        no_vals  = ct_data.loc[ct_data["response"] == "no",  "percentage"].values
        stat, p = mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
        results.append({
            "cell_type":              ct,
            "n_responders":           len(yes_vals),
            "n_non_responders":       len(no_vals),
            "median_responders":      round(float(pd.Series(yes_vals).median()), 4),
            "median_non_responders":  round(float(pd.Series(no_vals).median()), 4),
            "mann_whitney_u":         round(stat, 4),
            "p_value":                round(p, 6),
            "significant_p05":        p < 0.05,
        })

    results_df = pd.DataFrame(results)
    out = os.path.join(OUTPUT_DIR, "statistical_results.csv")
    results_df.to_csv(out, index=False)
    print(f"  Results written to {out}")

    sig = results_df.loc[results_df["significant_p05"], "cell_type"].tolist()
    print(f"  Significant populations (p < 0.05): {', '.join(sig) if sig else 'none'}")

    _make_boxplot(subset)
    return results_df


def _make_boxplot(subset: pd.DataFrame) -> None:
    colors = {"yes": "#4878CF", "no": "#D65F5F"}
    fig, axes = plt.subplots(1, len(CELL_TYPES), figsize=(16, 5))

    for ax, ct in zip(axes, CELL_TYPES):
        ct_data = subset[subset["cell_type"] == ct]
        groups = [
            ct_data.loc[ct_data["response"] == "yes", "percentage"].values,
            ct_data.loc[ct_data["response"] == "no",  "percentage"].values,
        ]
        bp = ax.boxplot(
            groups,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=1.5),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker="o", markersize=3, alpha=0.5),
        )
        for patch, key in zip(bp["boxes"], ["yes", "no"]):
            patch.set_facecolor(colors[key])
            patch.set_alpha(0.7)
        ax.set_title(ct.replace("_", " ").title(), fontsize=10)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Resp.", "Non-resp."], fontsize=8)
        ax.set_ylabel("Percentage (%)" if ct == CELL_TYPES[0] else "", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)

    fig.legend(
        handles=[
            mpatches.Patch(facecolor=colors["yes"], alpha=0.7, label="Responder"),
            mpatches.Patch(facecolor=colors["no"],  alpha=0.7, label="Non-responder"),
        ],
        loc="upper right", fontsize=9, framealpha=0.8,
    )
    fig.suptitle(
        "Cell-type percentages: Responders vs Non-responders (Melanoma, miraclib, PBMC)",
        fontsize=11,
    )
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "responder_boxplot.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Boxplot saved to {out}")


def run_part4(df: pd.DataFrame) -> None:
    print("Part 4: baseline subset queries...")
    sub = df[
        (df["condition"] == "melanoma")
        & (df["sample_type"] == "PBMC")
        & (df["time_from_treatment_start"] == 0)
        & (df["treatment"] == "miraclib")
    ].copy()

    # Q1: samples per project
    q1 = sub.groupby("project")["sample"].nunique().reset_index()
    q1.columns = ["project", "n_samples"]
    q1.to_csv(os.path.join(OUTPUT_DIR, "q1_samples_per_project.csv"), index=False)

    # Q2: unique subjects by response
    q2 = (
        sub[sub["response"].isin(["yes", "no"])]
        .groupby("response")["subject"].nunique()
        .reset_index()
    )
    q2.columns = ["response", "n_subjects"]
    q2.to_csv(os.path.join(OUTPUT_DIR, "q2_responder_counts.csv"), index=False)

    # Q3: unique subjects by sex
    q3 = sub.groupby("sex")["subject"].nunique().reset_index()
    q3.columns = ["sex", "n_subjects"]
    q3.to_csv(os.path.join(OUTPUT_DIR, "q3_sex_counts.csv"), index=False)

    # Q4: mean b_cell, melanoma males, time=0, miraclib, responders
    mean_bcell = round(
        df[
            (df["condition"] == "melanoma")
            & (df["sample_type"] == "PBMC")
            & (df["time_from_treatment_start"] == 0)
            & (df["treatment"] == "miraclib")
            & (df["sex"] == "M")
            & (df["response"] == "yes")
        ]["b_cell"].mean(),
        2,
    )

    print(f"  Samples per project:\n{q1.to_string(index=False)}")
    print(f"  Subjects by response:\n{q2.to_string(index=False)}")
    print(f"  Subjects by sex:\n{q3.to_string(index=False)}")
    print(f"  Mean B-cell (melanoma males, time=0, miraclib, response=yes): {mean_bcell}")

    with open(os.path.join(OUTPUT_DIR, "subset_summary.txt"), "w") as f:
        f.write(f"Q1 - samples per project:\n{q1.to_string(index=False)}\n\n")
        f.write(f"Q2 - subjects by response:\n{q2.to_string(index=False)}\n\n")
        f.write(f"Q3 - subjects by sex:\n{q3.to_string(index=False)}\n\n")
        f.write(f"Q4 - mean B-cell (melanoma males, time=0, miraclib, response=yes): {mean_bcell}\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_data()
    freq = run_part2(df)
    run_part3(freq)
    run_part4(df)
    print("Done. All outputs written to outputs/")


if __name__ == "__main__":
    main()
