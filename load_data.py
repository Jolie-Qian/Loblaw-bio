"""
load_data.py
Establishes the 4-table normalized SQLite schema and populates it from cell-count.csv.
Usage: python load_data.py
"""

import csv
import os
import sqlite3

CSV_PATH = "cell-count.csv"
DB_PATH = "cell_counts.db"

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    project_id  TEXT PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL REFERENCES projects(project_id),
    condition   TEXT NOT NULL,
    age         INTEGER,
    sex         TEXT NOT NULL CHECK(sex IN ('M', 'F')),
    treatment   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                 TEXT NOT NULL,
    time_from_treatment_start   INTEGER NOT NULL,
    response                    TEXT CHECK(response IN ('yes','no') OR response IS NULL)
);

CREATE TABLE IF NOT EXISTS cell_counts (
    sample_id   TEXT PRIMARY KEY REFERENCES samples(sample_id),
    b_cell      INTEGER NOT NULL DEFAULT 0,
    cd8_t_cell  INTEGER NOT NULL DEFAULT 0,
    cd4_t_cell  INTEGER NOT NULL DEFAULT 0,
    nk_cell     INTEGER NOT NULL DEFAULT 0,
    monocyte    INTEGER NOT NULL DEFAULT 0
);

-- Optimization indexes for upstream analysis queries
CREATE INDEX IF NOT EXISTS idx_subjects_cohort ON subjects(condition, treatment);
CREATE INDEX IF NOT EXISTS idx_samples_lookup ON samples(sample_type, time_from_treatment_start, response);
CREATE INDEX IF NOT EXISTS idx_samples_subject_fk ON samples(subject_id);
"""


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()
    print("Schema ready.")


def load_csv(conn: sqlite3.Connection, csv_path: str) -> None:
    cursor = conn.cursor()

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cursor.execute(
                "INSERT OR IGNORE INTO projects (project_id) VALUES (?)",
                (row["project"],)
            )

            cursor.execute(
                """INSERT OR IGNORE INTO subjects
                   (subject_id, project_id, condition, age, sex, treatment)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    row["subject"],
                    row["project"],
                    row["condition"],
                    int(row["age"]) if row["age"] else None,
                    row["sex"],
                    row["treatment"]
                )
            )

            raw_resp = row["response"].strip()
            response = raw_resp if raw_resp in ("yes", "no") else None

            cursor.execute(
                """INSERT OR IGNORE INTO samples
                   (sample_id, subject_id, sample_type, time_from_treatment_start, response)
                   VALUES (?, ?, ?, ?, ?)""",
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                    response
                )
            )

            cursor.execute(
                """INSERT OR IGNORE INTO cell_counts
                   (sample_id, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    row["sample"],
                    int(row["b_cell"]),
                    int(row["cd8_t_cell"]),
                    int(row["cd4_t_cell"]),
                    int(row["nk_cell"]),
                    int(row["monocyte"])
                )
            )

    conn.commit()

def main() -> None:
    if not os.path.exists(CSV_PATH):
        raise FileNotFoundError(f"Source data file '{CSV_PATH}' missing.")

    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    print(f"Initializing SQLite database at: {DB_PATH}")
    with sqlite3.connect(DB_PATH) as conn:
        init_db(conn)
        print("Loading clinical trial CSV data...")
        load_csv(conn, CSV_PATH)

        cursor = conn.cursor()
        total_samples = cursor.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        total_subjects = cursor.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
        print(f"Ingestion complete: {total_samples} samples, {total_subjects} subjects loaded successfully.")


if __name__ == "__main__":
    main()
