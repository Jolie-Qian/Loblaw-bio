# Loblaw Bio - Cell Count Analysis

## Overview

This project analyzes immune cell population data from a clinical trial dataset and provides:

* SQLite database creation and data loading
* Relative frequency analysis of immune cell populations
* Statistical comparison of responders versus non-responders
* Clinical subset analysis
* Interactive Streamlit dashboard

The project is designed to be fully reproducible in GitHub Codespaces.

---

## Setup

Install all required dependencies:

```bash
make setup
```

Alternatively:

```bash
pip install -r requirements.txt
```

---

## Reproducing Results

Run the complete analysis pipeline:

```bash
make pipeline
```

This command will:

1. Create and initialize the SQLite database
2. Load all records from `cell-count.csv`
3. Generate the Part 2 relative frequency summary table
4. Perform the Part 3 statistical analysis
5. Generate all required visualizations
6. Produce the Part 4 subset analysis outputs

All generated outputs are written to the `outputs/` directory.

The pipeline can also be executed manually:

```bash
python load_data.py
python analysis.py
```

---

## Dashboard

Launch the interactive dashboard:

```bash
make dashboard
```

or

```bash
streamlit run dashboard.py
```

Public dashboard:

https://loblaw-bio.streamlit.app

---

## Database Schema

The SQLite database follows a normalized relational design consisting of four tables:

```text
projects
    │
    └── subjects
            │
            └── samples
                    │
                    └── cell_counts
```

### projects

Stores project-level information.

| Column     | Description               |
| ---------- | ------------------------- |
| project_id | Unique project identifier |

### subjects

Stores patient-level metadata.

| Column     | Description               |
| ---------- | ------------------------- |
| subject_id | Unique patient identifier |
| project_id | Associated project        |
| condition  | Disease indication        |
| age        | Patient age               |
| sex        | Biological sex            |
| treatment  | Treatment assignment      |

### samples

Stores biological sample metadata.

| Column                    | Description                                |
| ------------------------- | ------------------------------------------ |
| sample_id                 | Unique sample identifier                   |
| subject_id                | Associated patient                         |
| sample_type               | Sample source                              |
| time_from_treatment_start | Timepoint relative to treatment initiation |
| response                  | Clinical response status                   |

### cell_counts

Stores immune cell population measurements.

| Column     | Description       |
| ---------- | ----------------- |
| sample_id  | Associated sample |
| b_cell     | B-cell count      |
| cd8_t_cell | CD8 T-cell count  |
| cd4_t_cell | CD4 T-cell count  |
| nk_cell    | NK-cell count     |
| monocyte   | Monocyte count    |

### Design Rationale

The schema separates project-level, patient-level, sample-level, and measurement-level data into independent tables.

This design was chosen because it:

* Minimizes data duplication
* Preserves data integrity through foreign-key relationships
* Simplifies future maintenance
* Supports efficient analytical queries

Response status is stored at the sample level rather than the subject level because in a clinical trial, response is formally assessed per visit, not a fixed property of the patient. A subject could theoretically have different response evaluations at different timepoints; the schema does not assume otherwise.

The schema scales naturally to hundreds of projects, thousands of patients, millions of samples, and additional assay types. New measurement tables can be added without modifying existing clinical metadata tables.

To support the required analyses, indexes are created on:

* condition and treatment
* sample_type, response, and time_from_treatment_start
* subject foreign-key relationships

These indexes optimize the primary query paths used throughout Parts 2-4.

---

## Code Structure

### load_data.py

Responsible for:

* Initializing the SQLite schema
* Creating tables and indexes
* Loading data from `cell-count.csv`
* Populating the database

### analysis.py

Responsible for:

* Generating the Part 2 frequency summary table
* Performing the Part 3 responder versus non-responder analysis
* Running statistical significance testing
* Creating required plots
* Executing the Part 4 subset analysis
* Writing all outputs to the `outputs/` directory

### dashboard.py

Provides an interactive Streamlit dashboard containing:

* Database overview
* Relative frequency summaries
* Statistical comparison results
* Clinical subset analysis results
* Interactive visualizations

### requirements.txt

Contains all Python dependencies required to run the project.

### Makefile

Provides standardized commands:

```bash
make setup
make pipeline
make dashboard
```

These commands allow the project to run consistently in GitHub Codespaces and local environments.

---

## Output Files

Running the pipeline generates the following output files:

```text
outputs/
├── frequency_summary.csv
├── statistical_results.csv
├── responder_boxplot.png
├── subset_summary.txt
├── q1_samples_per_project.csv
├── q2_responder_counts.csv
└── q3_sex_counts.csv
```

---

## Technologies Used

* Python
* SQLite
* Pandas
* SciPy
* Matplotlib
* Streamlit

---

## Running in GitHub Codespaces

The project is designed to run without modification in GitHub Codespaces.

To reproduce all results:

```bash
make setup
make pipeline
make dashboard
```

No additional configuration is required.
