# Loblaw Bio - Cell Count Analysis

## Setup

```bash
pip install -r requirements.txt
```

## Run the pipeline

```bash
python load_data.py
python analysis.py
```

## Launch dashboard

```bash
streamlit run dashboard.py
```

Dashboard: https://loblaw-bio.streamlit.app

## Database schema

Four normalized tables: projects -> subjects -> samples -> cell_counts

- projects: one row per project
- subjects: one row per patient, stores condition, age, sex, treatment
- samples: one row per biological sample, stores sample_type, time_from_treatment_start, response (response is kept at sample level because in a clinical trial it is assessed per visit, not a fixed patient attribute)
- cell_counts: one row per sample, stores the five immune cell counts (kept separate from samples so additional assay types can be added without altering sample metadata)

Indexes on condition/treatment, sample_type/time/response, and subject_id cover the main query paths. This layout scales cleanly to hundreds of projects and millions of samples because subject metadata is stored once and cell_counts can be partitioned by project if needed.

## Code structure

- load_data.py: initializes SQLite schema and loads cell-count.csv
- analysis.py: runs Parts 2-4, writes all outputs to outputs/
- dashboard.py: Streamlit dashboard with four sections matching Parts 1-4
- requirements.txt: Python dependencies
- Makefile: make install / make pipeline / make dashboard
