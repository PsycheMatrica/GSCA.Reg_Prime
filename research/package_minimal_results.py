#!/usr/bin/env python3
"""Package only the completed minimal GSCA experiment (Python standard library).

Requires gh with GH_TOKEN and a checkout containing the numerical parent commit.
No fitting or simulation runs here. Artifact files and seed records are retained.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import zipfile


BRANCH = "research/indicator-auxiliary-minimal-20261106"
BASELINE = "9d586446970a47ebd420516a0451143ba608fff9"
NUMERICAL_SOURCES = (
    "GSCA_Reg.m", "GSCA_Indicator_OOB.m", "ALS_GSCA_Reg.m", "GC_Boot.m",
    "Run_Example_GSCA_Reg.m", "ACSI_Comp_Replica.csv",
    "research/Run_Minimal_Indicator_Test.m",
    "research/run_minimal_indicator_example.m",
    "research/gsca_correction_dgp.m",
    "research/GSCA_Reg_OriginalFixture.m",
    "research/test_gsca_original_flow.m",
    "research/test_gsca_minimal_update.m",
    "research/MINIMAL_INDICATOR_UPDATE.md",
    ".github/workflows/gsca-minimal.yml",
)
PACKAGE_SOURCES = (
    "research/package_minimal_results.py", ".github/workflows/package-minimal.yml",
)
OPTIONAL_REPORT = "research/GSCA_Minimal_Test_Report_KO.md"
SIM_FILES = (
    "minimal_simulation.mat", "summary.csv", "replicates.csv",
    "summary.json", "replicates.json", "config.json",
)
EXAMPLE_FILES = ("minimal_example.mat", "paths.csv", "paths.json", "config.json")
CELLS = [(lam, beta) for lam in (0, 0.1) for beta in (0, 0.1, 0.3)]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def command(*args):
    result = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout


def api(repo, endpoint):
    return json.loads(command("gh", "api", f"repos/{repo}/{endpoint}"))


def read_json(path):
    def invalid(value):
        raise ValueError(f"Nonstandard JSON number {value} in {path}")
    return json.loads(path.read_text(encoding="utf-8-sig"), parse_constant=invalid)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_csv(path, rows):
    require(bool(rows), f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def row_list(value):
    return value if isinstance(value, list) else [value]


def check_csv(path, rows):
    """Check the parallel CSV and JSON records, including failure text and seeds."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        csv_rows = list(reader)
        require(set(reader.fieldnames or []) == set(rows[0]), f"CSV schema mismatch: {path}")
    require(len(csv_rows) == len(rows), f"CSV row count mismatch: {path}")
    for index, (left, right) in enumerate(zip(csv_rows, rows), 1):
        for key, expected in right.items():
            value = left[key]
            context = f"{path.name} row {index} field {key}"
            if expected is None:
                require(value.strip().lower() in ("", "nan", "na", "inf", "+inf", "-inf"), context)
            elif isinstance(expected, bool):
                require(value.lower() in (("1", "true") if expected else ("0", "false")), context)
            elif isinstance(expected, (int, float)):
                # MATLAB writetable may use fewer digits than jsonencode.
                require(math.isclose(float(value), expected, rel_tol=1e-5, abs_tol=1e-10), context)
            else:
                require(value == expected, context)


def check_simulation(folder, n):
    config = read_json(folder / "config.json")
    for key, expected in {"N": n, "Jp": 3, "Lambdas": [0, 0.1], "Betas": [0, 0.1, 0.3],
                          "Replications": 100, "Boot": 199, "Seed": 20261106}.items():
        require(config.get(key) == expected, f"Unexpected N={n} config {key}")
    summary = row_list(read_json(folder / "summary.json"))
    records = row_list(read_json(folder / "replicates.json"))
    require(len(summary) == 6 and len(records) == 600, f"Incomplete N={n} outputs")
    require([(r["Lambda"], r["Beta"]) for r in summary] == CELLS, f"Incorrect N={n} conditions")
    for condition, (lam, beta) in enumerate(CELLS):
        group = records[condition * 100:(condition + 1) * 100]
        s = summary[condition]
        require(s["N"] == n and s["Jp"] == 3, "Summary model mismatch")
        for rep, row in enumerate(group, 1):
            dataset_index = ((6 if n == 200 else 0) + condition) * 10000 + rep
            for key, expected in {"N": n, "Jp": 3, "Lambda": lam, "Beta": beta, "Rep": rep,
                                  "DatasetIndex": dataset_index, "DataSeed": 20261106 + dataset_index,
                                  "BootstrapSeed": 21261106 + dataset_index}.items():
                require(row.get(key) == expected, f"N={n} condition={condition + 1} rep={rep} invalid {key}")
            require(isinstance(row["Success"], bool), "Success must be a logical field")
            if row["Success"]:
                for key in ("PathEstimate", "PathBootstrapSE", "PathCI025", "PathCI975",
                            "IndicatorDeltaMean", "IndicatorLower95", "IndicatorP",
                            "LegacyComponentDeltaMean", "IndicatorEPE"):
                    require(row[key] is not None and math.isfinite(row[key]), f"Nonfinite successful result: {key}")
                require(row["IndicatorReject"] == (row["IndicatorLower95"] > 0),
                        "Rejection differs from original lower-bound rule")
            else:
                require(bool(row["FailureReason"]), "Failed dataset lacks its reason")
        completed = sum(row["Success"] for row in group)
        require(s["Attempted"] == 100 and s["Completed"] == completed
                and s["Failed"] == 100 - completed, "Summary counts disagree with records")
    check_csv(folder / "summary.csv", summary)
    check_csv(folder / "replicates.csv", records)
    return summary, records, config


def check_example(folder):
    config = read_json(folder / "config.json")
    for key, expected in {"Lambda": 0, "Boot": 199, "Seed": 20261106}.items():
        require(config.get(key) == expected, f"Unexpected ACSI config {key}")
    paths = row_list(read_json(folder / "paths.json"))
    expected = {(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6), (5, 6)}
    require(len(paths) == 9 and {(r["Source"], r["Target"]) for r in paths} == expected,
            "Incomplete or incorrect ACSI paths")
    for row in paths:
        require(all(isinstance(v, (int, float)) and math.isfinite(v) for v in row.values()),
                "Nonfinite ACSI result")
    check_csv(folder / "paths.csv", paths)
    return paths, config


def result_folder(artifact, name, files):
    candidates = [p for p in artifact.rglob(name) if p.is_dir()]
    require(len(candidates) == 1, f"Expected one {name} folder in {artifact}")
    folder = candidates[0]
    for filename in files:
        path = folder / filename
        require(path.is_file() and path.stat().st_size > 0, f"Missing or empty result file {path}")
    require({p.name for p in folder.iterdir()} == set(files), f"Unexpected result files in {folder}")
    return folder


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default="PsycheMatrica/GSCA.Reg_Prime")
    parser.add_argument("--run-id", type=int, default=42339690205)
    parser.add_argument("--artifact-root", type=Path, default=Path("minimal_source_artifacts"))
    parser.add_argument("--output-dir", type=Path, default=Path("minimal_package"))
    args = parser.parse_args()
    require(bool(os.environ.get("GH_TOKEN")), "GH_TOKEN is required for artifact download")
    root = Path(__file__).resolve().parents[1]
    os.chdir(root)
    for folder in (args.artifact_root, args.output_dir):
        require(not folder.exists(), f"Output already exists; use a fresh destination: {folder}")
    run = api(args.repo, f"actions/runs/{args.run_id}")
    require(run["status"] == "completed" and run["conclusion"] == "success",
            "The numerical parent run must be completed with conclusion=success")
    require(run["head_branch"] == BRANCH, "Unexpected numerical parent branch")
    require(run["path"] == ".github/workflows/gsca-minimal.yml", "Unexpected parent workflow")
    parent_sha = run["head_sha"]
    for relative in NUMERICAL_SOURCES:
        path = root / relative
        require(path.is_file(), f"Missing required source {relative}")
        require(path.read_bytes() == command("git", "show", f"{parent_sha}:{relative}"),
                f"Source differs from the successful numerical run: {relative}")
    for relative in ("ALS_GSCA_Reg.m", "GC_Boot.m"):
        require((root / relative).read_bytes() == command("git", "show", f"{BASELINE}:{relative}"),
                f"Original solver/bootstrap helper changed: {relative}")
    artifacts = api(args.repo, f"actions/runs/{args.run_id}/artifacts?per_page=100")
    require(artifacts["total_count"] <= 100, "Unexpected artifact pagination")
    selected = {}
    for name in ("gsca-minimal-validation", "gsca-minimal-n80", "gsca-minimal-n200"):
        matches = [a for a in artifacts["artifacts"] if a["name"] == name and not a["expired"]]
        require(len(matches) == 1, f"Missing or ambiguous artifact {name}")
        selected[name] = matches[0]
        destination = args.artifact_root / name
        destination.mkdir(parents=True)
        command("gh", "run", "download", str(args.run_id), "--repo", args.repo,
                "--name", name, "--dir", str(destination))
    folders = {
        "N80": result_folder(args.artifact_root / "gsca-minimal-n80", "N80", SIM_FILES),
        "N200": result_folder(args.artifact_root / "gsca-minimal-n200", "N200", SIM_FILES),
        "ACSI_example": result_folder(args.artifact_root / "gsca-minimal-validation", "ACSI_example", EXAMPLE_FILES),
    }
    summary80, records80, config80 = check_simulation(folders["N80"], 80)
    summary200, records200, config200 = check_simulation(folders["N200"], 200)
    paths, example_config = check_example(folders["ACSI_example"])
    summary = summary80 + summary200
    records = records80 + records200
    bundle = args.output_dir / "GSCA_Minimal_Indicator_Test"
    bundle.mkdir(parents=True)
    for relative in NUMERICAL_SOURCES + PACKAGE_SOURCES:
        source = root / relative
        require(source.is_file(), f"Missing packaging source {relative}")
        target = bundle / "code" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    if (root / OPTIONAL_REPORT).is_file():
        shutil.copyfile(root / OPTIONAL_REPORT, bundle / "GSCA_Minimal_Test_Report_KO.md")
    for name, folder in folders.items():
        shutil.copytree(folder, bundle / "results" / "raw" / name)
    combined = bundle / "results" / "combined"
    combined.mkdir(parents=True)
    write_json(combined / "summary.json", summary)
    write_csv(combined / "summary.csv", summary)
    write_json(combined / "replicates.json", records)
    write_csv(combined / "replicates.csv", records)
    write_json(combined / "ACSI_paths.json", paths)
    write_csv(combined / "ACSI_paths.csv", paths)
    write_json(combined / "configs.json", {"N80": config80, "N200": config200, "ACSI_example": example_config})
    provenance = {
        "Repository": args.repo, "NumericalRunID": args.run_id, "NumericalRunURL": run["html_url"],
        "NumericalCommit": parent_sha, "OriginalBaselineCommit": BASELINE,
        "PackagingCommit": command("git", "rev-parse", "HEAD").decode().strip(),
        "NumericalRunStatus": run["status"], "NumericalRunConclusion": run["conclusion"],
        "NumericalRunCreatedAt": run["created_at"], "NumericalRunUpdatedAt": run["updated_at"],
        "MATLABRelease": "R2026a", "RequiredToolbox": "Statistics and Machine Learning Toolbox",
        "Attempted": len(records), "Completed": sum(r["Success"] for r in records),
        "Failed": sum(not r["Success"] for r in records), "SimulationCells": len(summary),
        "ACSIPaths": len(paths), "BootstrapDrawsPerDataset": 199,
        "Artifacts": [{k: a.get(k) for k in ("id", "name", "size_in_bytes", "digest", "created_at")}
                      for a in selected.values()],
        "KoreanReportIncluded": (bundle / "GSCA_Minimal_Test_Report_KO.md").is_file(),
    }
    write_json(bundle / "provenance.json", provenance)
    (bundle / "README.txt").write_text(
        "Minimal GSCA indicator and auxiliary-regression evaluation update\n\n"
        "The successful numerical parent run is identified in provenance.json.\n"
        "Original ALS/CV/fitting, bootstrap draws and parameter summaries are retained.\n"
        "Only in-bag auxiliary replacement and fixed raw target-indicator evaluation change.\n"
        "No variance correction, nested bootstrap or predictive-null calibration is included.\n\n"
        "results/raw preserves every uploaded MAT, CSV and JSON file, including seeds,\n"
        "failure reasons, original TABLE.B matrices and the complete ACSI GSCA output.\n"
        "results/combined contains 12 simulation summaries, 1200 attempted-dataset records\n"
        "and the 9-path ACSI table. Summaries condition on successful datasets and report\n"
        "all failures. A completed experiment may include failed datasets; none are redrawn.\n"
        "Rejection is the original TABLE.B column 5 lower bound > 0. Rates are not labeled\n"
        "Type I error or power: beta=0 is not a calibrated predictive-null boundary.\n"
        "LegacyComponentDeltaMean has a different estimand and units and is descriptive.\n\n"
        "To rerun, open MATLAB in code/, addpath('research'), then run:\n"
        "test_gsca_minimal_update; test_gsca_original_flow;\n"
        "Run_Minimal_Indicator_Test(80); Run_Minimal_Indicator_Test(200);\n"
        "run_minimal_indicator_example;\n"
        "Requires MATLAB and Statistics and Machine Learning Toolbox.\n"
        "Research gsca_correction_dgp.m is only the shared controlled data generator.\n"
        "No previous correction/nested experiment implementation or results are packaged.\n",
        encoding="utf-8")
    manifest = [{"path": p.relative_to(bundle).as_posix(), "bytes": p.stat().st_size,
                 "sha256": sha256(p)} for p in sorted(bundle.rglob("*")) if p.is_file()]
    write_json(bundle / "file_manifest.json", manifest)
    zip_path = args.output_dir / "GSCA_Minimal_Indicator_Test_Package.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(bundle.rglob("*")):
            if path.is_file():
                info = zipfile.ZipInfo(path.relative_to(args.output_dir).as_posix(), (2026, 11, 6, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, path.read_bytes(), compresslevel=9)
    data = zip_path.read_bytes()
    metadata = {
        "file": zip_path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
        "git_blob_sha": hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest(),
        "run_id": args.run_id, "numerical_commit": parent_sha,
        "attempted": provenance["Attempted"], "completed": provenance["Completed"], "failed": provenance["Failed"],
    }
    write_json(args.output_dir / "package_metadata.json", metadata)
    (args.output_dir / (zip_path.name + ".sha256")).write_text(
        f"{metadata['sha256']}  {zip_path.name}\n", encoding="ascii")
    (args.output_dir / (zip_path.name + ".gitblobsha")).write_text(
        f"{metadata['git_blob_sha']}  {zip_path.name}\n", encoding="ascii")
    print("MINIMAL_SUMMARY_BEGIN", flush=True)
    print(json.dumps(summary, separators=(",", ":"), allow_nan=False), flush=True)
    print("MINIMAL_SUMMARY_END", flush=True)
    print("MINIMAL_ACSI_BEGIN", flush=True)
    print(json.dumps(paths, separators=(",", ":"), allow_nan=False), flush=True)
    print("MINIMAL_ACSI_END", flush=True)
    print("MINIMAL_PACKAGE_META " + json.dumps(metadata, separators=(",", ":")), flush=True)
    encoded = base64.b64encode(data).decode("ascii")
    print("MINIMAL_PACKAGE_BEGIN", flush=True)
    for offset in range(0, len(encoded), 8192):
        print(encoded[offset:offset + 8192], flush=True)
    print("MINIMAL_PACKAGE_END", flush=True)


if __name__ == "__main__":
    main()
