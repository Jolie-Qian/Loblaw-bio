"""
dashboard.py
Streamlit dashboard for the cell-count clinical trial analysis.
Usage: streamlit run dashboard.py
"""

import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from scipy.stats import mannwhitneyu

DB_PATH = "cell_counts.db"
CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
CELL_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}

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

st.set_page_config(page_title="Cell Count Analysis", layout="wide")


@st.cache_data
def load_raw() -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as conn:
        df = pd.read_sql(FLAT_QUERY, conn)
    df["time_from_treatment_start"] = df["time_from_treatment_start"].astype(int)
    return df


@st.cache_data
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


df = load_raw()
freq = compute_frequencies(df)

st.sidebar.title("Navigation")
section = st.sidebar.radio(
    "Section",
    ["Overview", "Part 2: Cell-Type Frequencies", "Part 3: Statistical Analysis", "Part 4: Subset Query"],
)


if section == "Overview":
    st.title("Cell Count Analysis Dashboard")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Samples",  f"{df['sample'].nunique():,}")
    col2.metric("Total Subjects", f"{df['subject'].nunique():,}")
    col3.metric("Projects",       df["project"].nunique())
    col4.metric("Conditions",     df["condition"].nunique())

    st.subheader("Dataset Breakdown")
    c1, c2 = st.columns(2)
    with c1:
        for label, col in [("Samples by Condition", "condition"), ("Samples by Treatment", "treatment")]:
            st.markdown(f"**{label}**")
            counts = df.groupby(col)["sample"].nunique().reset_index()
            counts.columns = [col.title(), "Samples"]
            st.dataframe(counts, width="stretch", hide_index=True)
    with c2:
        for label, col in [("Samples by Sample Type", "sample_type"), ("Samples by Project", "project")]:
            st.markdown(f"**{label}**")
            counts = df.groupby(col)["sample"].nunique().reset_index()
            counts.columns = [col.replace("_", " ").title(), "Samples"]
            st.dataframe(counts, width="stretch", hide_index=True)


elif section == "Part 2: Cell-Type Frequencies":
    st.title("Part 2: Cell-Type Frequencies")
    st.markdown("Relative frequency of each immune cell population per sample.")

    col1, col2, col3 = st.columns(3)
    sel_cond  = col1.selectbox("Condition",   ["All"] + sorted(freq["condition"].unique().tolist()))
    sel_treat = col2.selectbox("Treatment",   ["All"] + sorted(freq["treatment"].unique().tolist()))
    sel_stype = col3.selectbox("Sample Type", ["All"] + sorted(freq["sample_type"].unique().tolist()))

    filt = freq.copy()
    if sel_cond  != "All": filt = filt[filt["condition"]   == sel_cond]
    if sel_treat != "All": filt = filt[filt["treatment"]   == sel_treat]
    if sel_stype != "All": filt = filt[filt["sample_type"] == sel_stype]

    st.markdown(f"Showing **{len(filt):,}** rows ({filt['sample'].nunique():,} samples)")
    st.dataframe(
        filt[["sample", "subject", "project", "condition", "treatment",
              "response", "sample_type", "time_from_treatment_start",
              "cell_type", "cell_count", "total_count", "percentage"]].head(500),
        width="stretch",
        hide_index=True,
    )
    st.caption("Capped at 500 rows. Full data in outputs/frequency_summary.csv.")

    st.subheader("Median percentage by cell type")
    median_pct = (
        filt.groupby("cell_type")["percentage"]
        .median().reset_index()
        .rename(columns={"cell_type": "Cell Type", "percentage": "Median %"})
    )
    median_pct["Cell Type"] = median_pct["Cell Type"].map(CELL_LABELS)
    st.dataframe(median_pct, width="stretch", hide_index=True)


elif section == "Part 3: Statistical Analysis":
    st.title("Part 3: Statistical Analysis")
    st.markdown("Responders vs non-responders — melanoma patients on **miraclib**, PBMC samples only.")

    subset = freq[
        (freq["condition"] == "melanoma")
        & (freq["treatment"] == "miraclib")
        & (freq["sample_type"] == "PBMC")
        & (freq["response"].isin(["yes", "no"]))
    ].copy()

    st.markdown(
        f"Cohort: **{subset['sample'].nunique():,}** samples, "
        f"**{subset['subject'].nunique():,}** subjects"
    )

    results = []
    for ct in CELL_TYPES:
        ct_data  = subset[subset["cell_type"] == ct]
        yes_vals = ct_data.loc[ct_data["response"] == "yes", "percentage"].values
        no_vals  = ct_data.loc[ct_data["response"] == "no",  "percentage"].values
        stat, p  = mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
        results.append({
            "Cell Type":             CELL_LABELS[ct],
            "N Responders":          len(yes_vals),
            "N Non-responders":      len(no_vals),
            "Median Resp. (%)":      round(float(pd.Series(yes_vals).median()), 4),
            "Median Non-resp. (%)":  round(float(pd.Series(no_vals).median()), 4),
            "Mann-Whitney U":        round(stat, 2),
            "p-value":               round(p, 6),
            "Significant (p<0.05)":  p < 0.05,
        })

    results_df = pd.DataFrame(results)

    def highlight_sig(row):
        color = "background-color: #d4edda" if row["Significant (p<0.05)"] else ""
        return [color] * len(row)

    st.subheader("Mann-Whitney U Test Results")
    st.dataframe(results_df.style.apply(highlight_sig, axis=1), width="stretch", hide_index=True)

    sig = results_df.loc[results_df["Significant (p<0.05)"], "Cell Type"].tolist()
    if sig:
        st.success(f"Significantly different (p < 0.05): {', '.join(sig)}")
    else:
        st.info("No cell types reached significance at p < 0.05.")

    st.subheader("Boxplot: Responders vs Non-responders")
    colors = {"yes": "#4878CF", "no": "#D65F5F"}
    fig, axes = plt.subplots(1, len(CELL_TYPES), figsize=(16, 5))
    for ax, ct in zip(axes, CELL_TYPES):
        ct_data = subset[subset["cell_type"] == ct]
        groups  = [
            ct_data.loc[ct_data["response"] == "yes", "percentage"].values,
            ct_data.loc[ct_data["response"] == "no",  "percentage"].values,
        ]
        bp = ax.boxplot(
            groups, patch_artist=True,
            medianprops=dict(color="black", linewidth=1.5),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker="o", markersize=3, alpha=0.5),
        )
        for patch, key in zip(bp["boxes"], ["yes", "no"]):
            patch.set_facecolor(colors[key])
            patch.set_alpha(0.7)
        ax.set_title(CELL_LABELS[ct], fontsize=10)
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
    st.pyplot(fig)
    plt.close()


elif section == "Part 4: Subset Query":
    st.title("Part 4: Subset Query")
    st.markdown("Baseline cohort: **melanoma / PBMC / time=0 / miraclib**")

    sub = df[
        (df["condition"] == "melanoma")
        & (df["sample_type"] == "PBMC")
        & (df["time_from_treatment_start"] == 0)
        & (df["treatment"] == "miraclib")
    ].copy()

    st.subheader("Q1: Samples per project")
    q1 = sub.groupby("project")["sample"].nunique().reset_index()
    q1.columns = ["Project", "Samples"]
    st.dataframe(q1, width="stretch", hide_index=True)

    st.subheader("Q2: Subjects by response")
    q2 = (
        sub[sub["response"].isin(["yes", "no"])]
        .groupby("response")["subject"].nunique().reset_index()
    )
    q2.columns = ["Response", "Unique Subjects"]
    q2["Response"] = q2["Response"].map({"yes": "Responders", "no": "Non-responders"})
    st.dataframe(q2, width="stretch", hide_index=True)

    st.subheader("Q3: Subjects by sex")
    q3 = sub.groupby("sex")["subject"].nunique().reset_index()
    q3.columns = ["Sex", "Unique Subjects"]
    q3["Sex"] = q3["Sex"].map({"M": "Male", "F": "Female"})
    st.dataframe(q3, width="stretch", hide_index=True)

    st.subheader("Q4: Mean B-cell count — melanoma males, time=0, miraclib, response=yes")
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
    st.metric("Mean B-cell count", f"{mean_bcell:,.2f}")

    with st.expander("Note on quintazide"):
        st.markdown(
            "The exam prompt mentions **quintazide**. This analysis targets **miraclib** "
            "as specified. Quintazide is present in the dataset and can be analysed with "
            "the same pipeline by changing the treatment filter."
        )
