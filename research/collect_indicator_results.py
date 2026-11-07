#!/usr/bin/env python3
"""Collect the nine fixed-indicator GSCA studies and draw descriptive figures.

Usage: python research/collect_indicator_results.py --artifacts DIR --output DIR
       python research/collect_indicator_results.py --check-shard-fixture DIR
Requires Python 3.9+ and matplotlib; all data handling uses the standard library.
Missing experiments, conflicting duplicate artifacts, or broken pairing fail
explicitly. Smoke experiments are ignored. No statistical theorem is inferred
from observed simulation improvements or near-boundary rejection rates.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path


PRIMARY = (
    "primary_n80_jp1", "primary_n80_jp3",
    "primary_n200_jp1", "primary_n200_jp3",
)
NESTED = (
    "nested_n80_sample_boundary",
    "nested_n80_sample_positive",
    "nested_n80_oracle_boundary",
    "nested_n200_hetero_sample_boundary",
    "nested_n80_sample_boundary_inner199",
)
EXPECTED = PRIMARY + NESTED
METHODS = (
    ("Raw", "Raw percentile", "#7a7a7a"),
    ("Shifted", "Shifted percentile", "#e69f00"),
    ("OuterRaw", "Nested raw", "#0072b2"),
    ("OuterCorrected", "Nested corrected", "#009e73"),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def load_json(path):
    with path.open(encoding="utf-8-sig") as stream:
        return json.load(stream)


def table_rows(value, description):
    rows = value if isinstance(value, list) else [value]
    require(rows and all(isinstance(row, dict) for row in rows),
            f"{description}: expected a MATLAB struct or nonempty struct array")
    return rows


def scalar_list(value):
    return value if isinstance(value, list) else [value]


def numeric_cell(value):
    if value is None or value.strip().lower() in ("", "nan", "null"):
        return None
    lower = value.strip().lower()
    if lower in ("true", "false"):
        return lower == "true"
    try:
        number = float(value)
    except ValueError:
        return value
    if not math.isfinite(number):
        return None
    return int(number) if number.is_integer() else number


def load_replicates(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames is not None, f"Missing CSV header: {path}")
        return [{key: numeric_cell(value) for key, value in row.items()} for row in reader]


def experiment_name(directory, config):
    # Config.OutputDir may contain an absolute Linux or Windows runner path.
    output_parts = str(config.get("OutputDir", "")).replace("\\", "/").split("/")
    candidates = {part for part in directory.parts if part in EXPECTED}
    candidates.update(part for part in output_parts if part in EXPECTED)
    require(len(candidates) <= 1, f"Conflicting experiment identity in {directory}")
    return next(iter(candidates)) if candidates else None


def fingerprint(study):
    config = {key: value for key, value in study["config"].items() if key != "OutputDir"}
    payload = {"config": config, "summary": study["summary"],
               "truth": study["truth"], "replicates": study["replicates"]}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def group_key(row):
    return row.get("Condition"), row.get("Reference")


def same_value(left, right):
    # MATLAB's CSV and JSON writers can round the same float differently.
    if finite(left) and finite(right):
        return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)
    return left == right


def validate_study(name, study):
    cfg, summaries, truths, reps = (study[key] for key in ("config", "summary", "truth", "replicates"))
    needed = ("Ns", "Jps", "Lambdas", "DGPs", "Effects", "ReferenceModes", "Replications",
              "Boot", "Outer", "InnerBoot", "Alpha", "Seed")
    require(all(key in cfg for key in needed), f"{name}: missing required configuration fields")
    groups = math.prod(len(scalar_list(cfg[key])) for key in needed[:6])
    require(len(summaries) == groups, f"{name}: expected {groups} summary rows, found {len(summaries)}")
    require(len(reps) == groups * cfg["Replications"], f"{name}: incomplete replicates.csv")
    require(len(truths) * len(scalar_list(cfg["ReferenceModes"])) == groups,
            f"{name}: truth/summary counts disagree")
    require(len(summaries) == (24 if name in PRIMARY else 1),
            f"{name}: does not match the prespecified main experiment layout")
    require(len({group_key(row) for row in summaries}) == len(summaries), f"{name}: duplicate summary keys")
    require(len({row.get("Condition") for row in truths}) == len(truths), f"{name}: duplicate truth keys")
    require(len({(group_key(row), row.get("Rep")) for row in reps}) == len(reps),
            f"{name}: duplicate replicate keys")
    for row in summaries:
        selected = [rep for rep in reps if group_key(rep) == group_key(row)]
        require(len(selected) == cfg["Replications"] == row.get("Attempted"),
                f"{name}: replicate and summary counts disagree for {group_key(row)}")
        require(sum(bool(rep.get("Valid")) for rep in selected) == row.get("Completed"),
                f"{name}: valid-replicate count disagrees with summary")
        for key in ("N", "Jp", "Lambda", "DGP", "Effect", "Beta", "Reference"):
            require(all(same_value(rep.get(key), row.get(key)) for rep in selected),
                    f"{name}: inconsistent replicate field {key}")


def discover(artifacts):
    candidates, ignored, duplicates = {}, [], []
    for summary_path in sorted(artifacts.rglob("summary.json")):
        directory = summary_path.parent
        config_path = directory / "config.json"
        if not config_path.exists():
            if any(part in EXPECTED for part in directory.parts):
                raise ValueError(f"Main experiment is missing config.json: {directory}")
            ignored.append(str(directory))
            continue
        config = load_json(config_path)
        require(isinstance(config, dict), f"Configuration is not an object: {config_path}")
        if "Replications" in config:
            config.setdefault("ReplicationStart", 1)
            if config.get("TotalReplications") in (None, []):
                config["TotalReplications"] = config["Replications"]
        name = experiment_name(directory, config)
        if name is None:
            ignored.append(str(directory))
            continue
        for filename in ("truth.json", "replicates.csv"):
            require((directory / filename).is_file(), f"{name}: missing {filename} in {directory}")
        study = {
            "config": config,
            "summary": table_rows(load_json(summary_path), str(summary_path)),
            "truth": table_rows(load_json(directory / "truth.json"), f"{name}/truth.json"),
            "replicates": load_replicates(directory / "replicates.csv"),
            "source": str(directory),
        }
        validate_study(name, study)
        start = config.get("ReplicationStart", 1)
        identity = (start, config["Replications"])
        copies = candidates.setdefault(name, {})
        if identity in copies:
            require(fingerprint(study) == fingerprint(copies[identity]),
                    f"Conflicting duplicate results for {name}: {copies[identity]['source']} and {directory}")
            duplicates.append({"Experiment": name, "IgnoredIdenticalCopy": str(directory)})
        else:
            copies[identity] = study
    missing = [name for name in EXPECTED if name not in candidates]
    require(not missing, "Missing required main experiments: " + ", ".join(missing))
    studies = {}
    for name in EXPECTED:
        parts = list(candidates[name].values())
        if name in PRIMARY:
            require(len(parts) == 1, f"{name}: primary experiment must be unsharded")
            study = parts[0]
            require(study["config"].get("ReplicationStart", 1) == 1,
                    f"{name}: primary experiment starts after replicate 1")
            require(study["config"]["TotalReplications"] == study["config"]["Replications"],
                    f"{name}: primary experiment is incomplete")
            studies[name] = study
        else:
            total = 200 if name in NESTED[:2] else 100
            if len(parts) == 1 and parts[0]["config"]["Replications"] == total:
                study = parts[0]
                require(study["config"].get("ReplicationStart", 1) == 1,
                        f"{name}: full experiment starts after replicate 1")
                studies[name] = study
            else:
                studies[name] = merge_nested_parts(name, parts, total)
    return studies, ignored, duplicates


def recompute_summary(template, truth, cfg, reps):
    """Mirror v4 summarize() using every saved row, including failures."""
    result = {key: None for key in template}
    result.update({key: value for key, value in truth.items() if key in result})
    result["Reference"] = reps[0]["Reference"]
    good = [row for row in reps if row.get("Valid")]
    result.update({
        "Attempted": len(reps), "Completed": len(good), "DatasetFailures": len(reps) - len(good),
        "BootstrapFitFailures": sum(row["FitFailures"] for row in reps if finite(row.get("FitFailures"))),
        "OuterFailedFits": sum(row["OuterFailures"] for row in reps if finite(row.get("OuterFailures"))),
        "OuterFailedDatasets": sum(not row.get("OuterValid") for row in reps) if cfg["Outer"] > 0 else 0,
        "OuterRawCompleted": 0, "OuterCorrectedCompleted": 0,
    })
    if not good:
        return result
    for field in ("Raw", "Corrected", "Correction"):
        require(all(finite(row.get(field)) for row in good), f"Valid replicate has missing {field}")
    raw = [row["Raw"] for row in good]
    corrected = [row["Corrected"] for row in good]
    correction = [row["Correction"] for row in good]
    result.update({"MeanRaw": statistics.fmean(raw), "MeanCorrected": statistics.fmean(corrected),
                   "MeanCorrection": statistics.fmean(correction)})
    if len(good) > 1:
        result["EmpiricalSDRaw"] = statistics.stdev(raw)
        result["EmpiricalSDCorrected"] = statistics.stdev(corrected)
        result["MeanRawMCSE"] = result["EmpiricalSDRaw"] / math.sqrt(len(good))
        result["MeanCorrectedMCSE"] = result["EmpiricalSDCorrected"] / math.sqrt(len(good))
        result["CorrectionMCSE"] = statistics.stdev(correction) / math.sqrt(len(good))
    if truth["TruthValid"]:
        target, target_mcse = truth["TruthMean"], truth["TargetMCSE"]
        require(finite(target) and finite(target_mcse), "Valid truth has missing mean/MCSE")
        for suffix, values in (("Raw", raw), ("Corrected", corrected)):
            result["Bias" + suffix] = statistics.fmean(value - target for value in values)
            result["RMSE" + suffix] = math.sqrt(statistics.fmean((value - target) ** 2 for value in values))
            if len(good) > 1:
                result["Bias" + suffix + "MCSE"] = math.sqrt(result["Mean" + suffix + "MCSE"] ** 2 + target_mcse ** 2)
        for suffix, _, _ in METHODS:
            lower = [row["Lower" + suffix] for row in good if finite(row.get("Lower" + suffix))]
            result["LowerCoverage" + suffix] = (sum(value <= target for value in lower) / len(lower)) if lower else None
    for suffix in ("Legacy", "Raw", "Shifted", "OuterRaw", "OuterCorrected"):
        values = [row["p" + suffix] for row in good if finite(row.get("p" + suffix))]
        stats = wilson(sum(value < cfg["Alpha"] for value in values), len(values))
        result["Reject" + suffix] = stats["Rate"]
        result["MCSE" + suffix] = stats["MCSE"]
        result["Reject" + suffix + "Lo"], result["Reject" + suffix + "Hi"] = stats["Wilson95"]
        if suffix.startswith("Outer"):
            result[suffix + "Completed"] = len(values)
    return result


def merge_nested_parts(name, parts, total):
    require(len(parts) == 4, f"{name}: require all four shards; found {len(parts)}")
    parts.sort(key=lambda part: part["config"].get("ReplicationStart", 1))
    first = parts[0]
    shard_size = total // 4
    shared_exclusions = {"Replications", "ReplicationStart", "OutputDir"}
    shared = lambda cfg: {key: value for key, value in cfg.items() if key not in shared_exclusions}
    for index, part in enumerate(parts):
        cfg = part["config"]
        expected_start = 1 + index * shard_size
        require(cfg.get("ReplicationStart") == expected_start and cfg["Replications"] == shard_size,
                f"{name}: missing/overlapping shard at expected start {expected_start}")
        require(cfg.get("TotalReplications") == total, f"{name}: shard TotalReplications differs from {total}")
        require(shared(cfg) == shared(first["config"]), f"{name}: shard configurations differ beyond start/count/output")
        require(part["truth"] == first["truth"], f"{name}: shard calibration or independent truth differs")
        require(sorted(row["Rep"] for row in part["replicates"]) == list(range(expected_start, expected_start + shard_size)),
                f"{name}: shard has incomplete or overlapping global replicate indices")
        template = part["summary"][0]
        for key in ("Condition", "N", "Jp", "Lambda", "DGP", "Effect", "Beta", "Reference", "CalibrationMethod"):
            require(template.get(key) == first["summary"][0].get(key), f"{name}: shard metadata differs in {key}")
        for row in part["replicates"]:
            if finite(row.get("Beta")):
                require(row["DataSeed"] == cfg["Seed"] + 20000000 + row["Condition"] * 20000 + row["Rep"],
                        f"{name}: data seed does not use the global replicate index")
                require(row["BootstrapSeed"] == cfg["Seed"] + 30000000 + row["Condition"] * 20000 + row["Rep"],
                        f"{name}: bootstrap seed does not use the global replicate index")
        recomputed = recompute_summary(template, part["truth"][0], cfg, part["replicates"])
        differences = {key: [recomputed.get(key), value] for key, value in template.items()
                       if not same_value(recomputed.get(key), value)}
        require(not differences, f"{name}: Python recomputation differs from a saved shard summary: {differences}")
    reps = sorted([row for part in parts for row in part["replicates"]], key=lambda row: row["Rep"])
    require([row["Rep"] for row in reps] == list(range(1, total + 1)), f"{name}: merged global Rep indices are incomplete")
    config = dict(first["config"])
    config.update({"Replications": total, "ReplicationStart": 1,
                   "OutputDir": "results_indicator/" + name})
    merged = {"config": config, "truth": first["truth"], "replicates": reps,
              "source": [part["source"] for part in parts],
              "shard_configs": [part["config"] for part in parts],
              "summary": [recompute_summary(first["summary"][0], first["truth"][0], config, reps)]}
    validate_study(name, merged)
    return merged


def check_shard_fixture(directory):
    """Validate the production merge against a MATLAB unsplit R=4 fixture."""
    def read_part(label):
        path = directory / label
        for filename in ("summary.json", "truth.json", "config.json", "replicates.csv"):
            require((path / filename).is_file(), f"Shard fixture missing {path / filename}")
        study = {
            "config": load_json(path / "config.json"),
            "summary": table_rows(load_json(path / "summary.json"), str(path / "summary.json")),
            "truth": table_rows(load_json(path / "truth.json"), str(path / "truth.json")),
            "replicates": load_replicates(path / "replicates.csv"),
            "source": str(path),
        }
        validate_study("shard_fixture", study)
        return study
    full = read_part("full")
    require(full["config"]["Replications"] == full["config"].get("TotalReplications") == 4 and
            full["config"].get("ReplicationStart") == 1,
            "Unsplit shard fixture must contain global replications 1 through 4")
    parts = [read_part("part" + str(index)) for index in range(1, 5)]
    merged = merge_nested_parts("shard_fixture", parts, 4)
    require(full["truth"] == merged["truth"], "Fixture merged and unsplit independent truth differ")
    require(set(full["summary"][0]) == set(merged["summary"][0]), "Fixture summary schemas differ")
    differences = {key: [merged["summary"][0].get(key), value] for key, value in full["summary"][0].items()
                   if not same_value(merged["summary"][0].get(key), value)}
    require(not differences, f"Fixture merged summary differs from MATLAB unsplit summary: {differences}")
    full_reps = sorted(full["replicates"], key=lambda row: row["Rep"])
    require(len(full_reps) == len(merged["replicates"]) == 4, "Fixture replicate counts differ")
    for expected, actual in zip(full_reps, merged["replicates"]):
        require(set(expected) == set(actual), "Fixture replicate schemas differ")
        differences = {key: [actual.get(key), value] for key, value in expected.items()
                       if not same_value(actual.get(key), value)}
        require(not differences, f"Fixture global replicate {expected['Rep']} differs: {differences}")
    print(json.dumps({"ShardFixture": str(directory), "Status": "passed",
                      "MatchedGlobalReplicates": 4, "SummaryFieldsCompared": len(full["summary"][0]),
                      "TruthIdentical": True, "UsesProductionMerge": True,
                      "RelativeTolerance": 1e-12, "AbsoluteTolerance": 1e-14}))


def improvement_counts(rows):
    eligible = [row for row in rows if row.get("TruthValid") and
                all(finite(row.get(key)) for key in ("BiasRaw", "BiasCorrected", "RMSERaw", "RMSECorrected"))]
    bias = lambda row: abs(row["BiasCorrected"]) < abs(row["BiasRaw"])
    rmse = lambda row: row["RMSECorrected"] < row["RMSERaw"]
    keys = ("Experiment", "N", "Jp", "Lambda", "DGP", "Effect", "Reference", "TruthMean", "TargetMCSE",
            "BiasRaw", "BiasCorrected", "RMSERaw", "RMSECorrected")
    return {
        "Rows": len(rows), "EligibleFiniteTruthRows": len(eligible),
        "ExcludedRows": len(rows) - len(eligible),
        "AbsoluteBiasImproved": sum(bias(row) for row in eligible),
        "RMSEImproved": sum(rmse(row) for row in eligible),
        "BothImproved": sum(bias(row) and rmse(row) for row in eligible),
        "AbsoluteBiasTies": sum(abs(row["BiasCorrected"]) == abs(row["BiasRaw"]) for row in eligible),
        "RMSETies": sum(row["RMSECorrected"] == row["RMSERaw"] for row in eligible),
        "RowsWithoutBothImprovements": [{key: row.get(key) for key in keys}
                                      for row in eligible if not (bias(row) and rmse(row))],
    }


def wilson(successes, total):
    if not total:
        return {"Rejections": 0, "Completed": 0, "Rate": None, "MCSE": None, "Wilson95": [None, None]}
    p = successes / total
    z = 1.95996398454005
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return {"Rejections": successes, "Completed": total, "Rate": p,
            "MCSE": math.sqrt(p * (1 - p) / total),
            "Wilson95": [max(0.0, center - half), min(1.0, center + half)]}


def coverage_sensitivity(lower_values, truth_mean, target_mcse):
    """Vary the target used to score saved one-sided lower bounds.

    The resulting range is sensitivity to an approximate Monte Carlo target
    interval, not a confidence interval for the unknown coverage probability.
    """
    values = [value for value in lower_values if finite(value)]
    valid_truth = finite(truth_mean)
    valid_mcse = finite(target_mcse) and target_mcse >= 0
    targets = [truth_mean, None, None]
    if valid_truth and valid_mcse:
        targets[1:] = [truth_mean - 1.96 * target_mcse,
                       truth_mean + 1.96 * target_mcse]
    coverages = [sum(value <= target for value in values) / len(values)
                 if values and finite(target) else None for target in targets]
    return {
        "Completed": len(values),
        "TruthMean": truth_mean if valid_truth else None,
        "TargetMCSE": target_mcse if valid_mcse else None,
        "TargetMinus1p96MCSE": targets[1],
        "TargetPlus1p96MCSE": targets[2],
        "CoverageAtTruthMean": coverages[0],
        "CoverageAtTargetMinus1p96MCSE": coverages[1],
        "CoverageAtTargetPlus1p96MCSE": coverages[2],
        "CoverageSensitivityRange": [coverages[1], coverages[2]],
        "Interpretation": "Sensitivity of empirical one-sided lower-bound coverage to Monte Carlo truth uncertainty; not a confidence interval for coverage.",
    }


def nested_table(studies):
    output = []
    for name in NESTED:
        study = studies[name]
        row = study["summary"][0]
        entry = {"Experiment": name, **{key: row.get(key) for key in
                 ("N", "Jp", "Lambda", "DGP", "Effect", "Beta", "Reference", "Classification",
                  "TruthMean", "TargetMCSE", "Attempted", "Completed", "DatasetFailures", "OuterFailedDatasets")}}
        entry.update({key: study["config"][key] for key in ("Boot", "Outer", "InnerBoot", "Alpha")})
        entry["Methods"] = {}
        for suffix, label, _ in METHODS:
            values = [rep.get("p" + suffix) for rep in study["replicates"]]
            values = [value for value in values if finite(value)]
            stats = wilson(sum(value < entry["Alpha"] for value in values), len(values))
            reported = row.get("Reject" + suffix)
            require((stats["Rate"] is None and reported is None) or
                    (finite(stats["Rate"]) and finite(reported) and
                     math.isclose(stats["Rate"], reported, rel_tol=1e-12, abs_tol=1e-12)),
                    f"{name}: replicate-derived {suffix} rejection rate disagrees with summary")
            sensitivity = coverage_sensitivity(
                [rep.get("Lower" + suffix) for rep in study["replicates"]],
                row.get("TruthMean"), row.get("TargetMCSE"))
            point_coverage = sensitivity["CoverageAtTruthMean"]
            reported_coverage = row.get("LowerCoverage" + suffix)
            require((point_coverage is None and reported_coverage is None) or
                    (finite(point_coverage) and finite(reported_coverage) and
                     math.isclose(point_coverage, reported_coverage, rel_tol=1e-12, abs_tol=1e-12)),
                    f"{name}: replicate-derived {suffix} lower-bound coverage disagrees with summary")
            entry["Methods"][suffix] = {"Label": label, **stats,
                                        "LowerCoverage": point_coverage,
                                        "LowerCoverageMCTargetSensitivity": sensitivity}
        output.append(entry)
    return output


def paired_precision(studies):
    baseline = studies["nested_n80_sample_boundary"]
    refined = studies["nested_n80_sample_boundary_inner199"]
    for key in ("Ns", "Jps", "Lambdas", "DGPs", "Effects", "ReferenceModes", "Boot", "Outer", "Alpha",
                "Seed", "CalReplications", "CalIterations", "TruthReplications", "Tol", "MinRcond", "MaxIter"):
        require(baseline["config"].get(key) == refined["config"].get(key), f"Inner-count pairing changes config field {key}")
    require(baseline["config"]["InnerBoot"] == 99 and refined["config"]["InnerBoot"] == 199,
            "Expected inner bootstrap counts 99 and 199")
    first = sorted(baseline["replicates"], key=lambda row: row["Rep"])[:100]
    second = sorted(refined["replicates"], key=lambda row: row["Rep"])
    require(len(first) == len(second) == 100, "Inner-count comparison requires exactly 100 paired datasets")
    key = lambda row: (row.get("DataSeed"), row.get("BootstrapSeed"))
    lookup = {key(row): row for row in second}
    require(len(lookup) == 100 and set(lookup) == {key(row) for row in first},
            "First 100 inner99 rows do not match inner199 DataSeed/BootstrapSeed pairs")
    pairs = [(row, lookup[key(row)]) for row in first]
    original_fields = ("Raw", "Corrected", "Correction", "pRaw", "pShifted")
    for a, b in pairs:
        require(all(finite(value) for value in (*key(a), a.get("Beta"), b.get("Beta"))),
                "Pairing has missing seeds or beta")
        require(a["Beta"] == b["Beta"], "Paired precision comparison has unequal beta values")
        for field in ("Rep", "N", "Jp", "Lambda", "DGP", "Effect", "Reference"):
            require(a.get(field) == b.get(field), f"Paired datasets differ in {field}")
        for field in original_fields:
            require(same_value(a.get(field), b.get(field)),
                    f"Changing inner precision changed original-Boot statistic {field} for seeds {key(a)}")
    alpha = baseline["config"]["Alpha"]
    result = {"MatchedDatasets": 100, "PairingFields": ["DataSeed", "BootstrapSeed"],
              "BetaEqual": True, "Beta": pairs[0][0]["Beta"],
              "OriginalPointEstimatesIdentical": True,
              "OriginalBootstrapFieldsCompared": list(original_fields),
              "OriginalBootstrapFinitePairs": sum(all(finite(a.get(field)) and finite(b.get(field))
                                                       for field in original_fields) for a, b in pairs),
              "SerializationComparisonTolerance": {"Relative": 1e-12, "Absolute": 1e-14},
              "OriginalBootstrapComparisonNote": "All five original-Boot outputs match within serialization tolerance; missing values must also match. Finite-pair count is reported separately.",
              "Methods": {}}
    for suffix, label, _ in METHODS:
        valid = [(a, b) for a, b in pairs if finite(a.get("p" + suffix)) and finite(b.get("p" + suffix))]
        decisions = [(a["p" + suffix] < alpha, b["p" + suffix] < alpha) for a, b in valid]
        a_reject = sum(a for a, _ in decisions)
        b_reject = sum(b for _, b in decisions)
        discordant = sum(a != b for a, b in decisions)
        result["Methods"][suffix] = {
            "Label": label, "PairsWithFiniteInference": len(valid), "ExcludedPairs": 100 - len(valid),
            "Inner99": wilson(a_reject, len(valid)), "Inner199": wilson(b_reject, len(valid)),
            "BothReject": sum(a and b for a, b in decisions),
            "NeitherReject": sum(not a and not b for a, b in decisions),
            "Inner99OnlyReject": sum(a and not b for a, b in decisions),
            "Inner199OnlyReject": sum(not a and b for a, b in decisions),
            "DiscordantPairs": discordant,
            "DiscordanceRate": discordant / len(valid) if valid else None,
            "PairedRejectionRateDifference199Minus99": (b_reject - a_reject) / len(valid) if valid else None,
        }
    return result


def summarize(studies, summaries, truths, ignored, duplicates):
    primary = [row for row in summaries if row["Experiment"] in PRIMARY]
    learned = [row for row in primary if row["Jp"] == 3]
    require(len(primary) == 96 and len(learned) == 48, "Unexpected primary or Jp=3 row count")
    budget = []
    for name in EXPECTED:
        study = studies[name]
        cfg = study["config"]
        datasets = len(study["summary"]) * cfg["Replications"]
        budget.append({
            "Experiment": name, "SummaryRows": len(study["summary"]),
            "DatasetAnalysesRequested": datasets,
            "OriginalOOBFitsRequested": datasets * cfg["Boot"],
            "OuterEstimatorReplicatesRequested": datasets * cfg["Outer"],
            "OuterInnerOOBFitsRequested": datasets * cfg["Outer"] * cfg["InnerBoot"],
            "TotalOOBFitsRequested": datasets * (cfg["Boot"] + cfg["Outer"] * cfg["InnerBoot"]),
            "FullDatasetFitsRequested": datasets * (1 + cfg["Outer"]),
            "IndependentTruthFitsRequested": sum(row["TruthRequested"] for row in study["truth"]) * len(study.get("shard_configs", [cfg])),
        })
    budget_keys = [key for key in budget[0] if key not in ("Experiment", "SummaryRows")]
    budget_totals = {key: sum(row[key] for row in budget) for key in budget_keys}
    require(budget_totals["DatasetAnalysesRequested"] == 10300,
            "Requested dataset-analysis count differs from the prespecified 10,300")
    require(budget_totals["TotalOOBFitsRequested"] == 17830400,
            "Requested OOB-fit count differs from the prespecified 17,830,400")
    failures = {key: sum(row.get(key, 0) or 0 for row in summaries) for key in
                ("DatasetFailures", "BootstrapFitFailures", "OuterFailedDatasets", "OuterFailedFits")}
    failures.update({"IndependentTruthFitFailures": sum(row.get("TruthFailed", 0) or 0 for row in truths),
                     "InvalidTruthRows": sum(not row.get("TruthValid", False) for row in truths),
                     "InvalidCalibrationTruthRows": sum(not row.get("CalibrationValid", False) for row in truths)})
    return {
        "Estimand": "Delta EPE_p: fixed equal-weight target-indicator average prediction-risk difference.",
        "PrimaryObservedImprovements": improvement_counts(primary),
        "PrimaryJp3ObservedImprovements": improvement_counts(learned),
        "Failures": failures,
        "RequestedComputation": {"ByExperiment": budget,
                                 "Totals": budget_totals,
                                 "PrespecifiedDatasetAnalysisCount": 10300,
                                 "PrespecifiedOOBFitCount": 17830400,
                                 "PrespecifiedCountsVerifiedFromConfigs": True},
        "NestedComparison": nested_table(studies),
        "PairedInnerPrecisionComparison": paired_precision(studies),
        "ArtifactSources": {name: studies[name]["source"] for name in EXPECTED},
        "MergedShardCounts": {name: len(studies[name].get("shard_configs", [])) for name in EXPECTED},
        "IgnoredNonMainDirectories": ignored,
        "IgnoredIdenticalDuplicates": duplicates,
        "Notes": [
            "Counts of observed improvements are descriptive, not statistical evidence of uniform superiority.",
            "Oracle/sample reference rows reuse datasets and are paired analyses, not independent replications.",
            "Dataset analyses include repetitions across oracle/sample reference and the inner-count comparison.",
            "Requested OOB-fit budgets exclude full-dataset, boundary-calibration, and independent-truth fits; early failures may reduce actual execution.",
            "Independent truth-fit totals count repeated validation performed separately in each experiment/shard; merged truth rows retain the identical target once.",
            "Nested shards merge only when all four parts, all global replicate indices, matching metadata/truth, and the original seed formulas are verified. Aggregate summaries are recomputed from every saved replicate, including failure counts.",
            "BootstrapFitFailures includes failed original OOB draws/maps/losses, not solely optimizer failures.",
            "OuterFailedFits is the MATLAB summary field: it counts failed outer estimator replicates, not individual failed inner model fits.",
            "Boundary labels outside the analytic OLS bridge are Monte Carlo calibrated near-boundary candidates.",
            "Rejection intervals are Wilson 95% Monte Carlo intervals conditional on finite inference; examine failures and denominators.",
            "Target uncertainty affects reported bias/RMSE. Paired inner precision comparisons change simulation precision, not the inferential target.",
            "Nested lower-bound coverage sensitivity evaluates the same saved bounds at TruthMean and TruthMean plus/minus 1.96 TargetMCSE. Its range is not a confidence interval for coverage.",
        ],
    }


def plot_bias_rmse(summaries, output, plt):
    rows = [row for row in summaries if row["Experiment"] in PRIMARY and row["Jp"] == 3
            and row["Reference"] == "sample" and row["Effect"] == "boundary"]
    rows.sort(key=lambda row: (row["N"], row["Lambda"], row["DGP"] != "normal"))
    require(len(rows) == 8, "Expected eight Jp=3/sample near-boundary rows for the bias/RMSE plot")
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.9), sharey=True)
    y = list(range(8))
    labels = [f"N={row['N']}, ��={row['Lambda']:g}, " +
              ("normal" if row["DGP"] == "normal" else "heteroscedastic") for row in rows]
    colors = ("#0072b2", "#d55e00")
    for ax, fields, title in zip(axes, (("BiasRaw", "BiasCorrected"), ("RMSERaw", "RMSECorrected")),
                                 ("Bias", "RMSE")):
        for pos, row in zip(y, rows):
            values = [row.get(field) for field in fields]
            if all(finite(value) for value in values):
                ax.plot(values, [pos, pos], color="#b7b7b7", linewidth=1.4, zorder=1)
        for field, label, color, marker in zip(fields, ("Raw", "Variance corrected"), colors, ("o", "s")):
            xvals = [row.get(field) if finite(row.get(field)) else math.nan for row in rows]
            ax.scatter(xvals, y, label=label, color=color, marker=marker, s=40, zorder=3)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel(r"$\Delta\mathrm{EPE}_p$ units")
        ax.grid(axis="x", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].axvline(0, color="#555555", linewidth=0.9, linestyle="--", zorder=0)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[1].set_xlim(left=0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.61, 0.94), frameon=False)
    fig.suptitle("Fixed-indicator GSCA: point-estimate performance near the predictive boundary", x=0.03,
                 ha="left", fontsize=13, fontweight="bold")
    fig.text(0.03, 0.025, "Three indicators per block; sample-reference projection. Observed simulation bias/RMSE; independent Monte Carlo risk targets.",
             fontsize=9, color="#444444")
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.87))
    fig.savefig(output / "indicator_boundary_bias_rmse.png", dpi=220, facecolor="white")
    plt.close(fig)


def plot_nested(nested, output, plt):
    lookup = {row["Experiment"]: row for row in nested}
    order = ("nested_n80_sample_boundary", "nested_n80_oracle_boundary",
             "nested_n200_hetero_sample_boundary", "nested_n80_sample_boundary_inner199")
    boundary = [lookup[name] for name in order]
    positive = lookup["nested_n80_sample_positive"]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.6), gridspec_kw={"width_ratios": [2.5, 1]})
    offsets = (-0.24, -0.08, 0.08, 0.24)
    for suffix_info, offset in zip(METHODS, offsets):
        suffix, label, color = suffix_info
        for ax, rows in ((axes[0], boundary), (axes[1], [positive])):
            for i, row in enumerate(rows):
                stat = row["Methods"][suffix]
                if stat["Rate"] is None:
                    continue
                lo, hi = stat["Wilson95"]
                ax.errorbar(i + offset, stat["Rate"], yerr=[[stat["Rate"] - lo], [hi - stat["Rate"]]],
                            fmt="o", markersize=5, color=color, capsize=3, linewidth=1.1,
                            label=label if i == 0 else None)
    labels = []
    for row in boundary:
        dgp = "normal" if row["DGP"] == "normal" else "heteroscedastic"
        labels.append(f"N={row['N']}, {dgp}\n{row['Reference']}, inner={row['InnerBoot']}\nR={row['Attempted']}")
    axes[0].set_xticks(range(4), labels)
    axes[0].axhline(0.05, color="#555555", linestyle="--", linewidth=1, label="Nominal .05")
    axes[0].set_title("Near-boundary rejection", loc="left", fontweight="bold")
    axes[1].set_xticks([0], [f"N={positive['N']}, normal\nsample, inner={positive['InnerBoot']}\nR={positive['Attempted']}"])
    axes[1].set_title("Positive-contribution rejection", loc="left", fontweight="bold", fontsize=10)
    for ax, rows in ((axes[0], boundary), (axes[1], [positive])):
        highest = max((row["Methods"][suffix]["Wilson95"][1] or 0
                       for row in rows for suffix, _, _ in METHODS), default=0)
        ax.set_ylim(0, min(1, max(0.15, highest * 1.15)))
        ax.set_ylabel("Rejection probability")
        ax.grid(axis="y", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_xlim(-0.5, 3.5)
    axes[1].set_xlim(-0.5, 0.5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.935), frameon=False)
    fig.suptitle("Fixed-indicator GSCA: one-sided bootstrap rejection comparisons", x=0.035, ha="left",
                 fontsize=14, fontweight="bold")
    fig.text(0.035, 0.065, "Error bars: Wilson 95% Monte Carlo intervals. Near-boundary risks are estimated, not exact null values.",
             fontsize=9, color="#444444")
    fig.text(0.035, 0.037, "The inner=199 condition shares its 100 datasets with the first 100 of the inner=99 sample-reference condition.",
             fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, 0.105, 1, 0.84))
    fig.savefig(output / "indicator_nested_rejection.png", dpi=220, facecolor="white")
    plt.close(fig)


def write_json(path, value):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-shard-fixture", type=Path,
                        help="Validate four R=1 shards against a MATLAB unsplit R=4 fixture; no matplotlib required")
    args = parser.parse_args()
    if args.check_shard_fixture is not None:
        require(args.artifacts is None and args.output is None, "Use fixture validation separately from artifact collection")
        check_shard_fixture(args.check_shard_fixture)
        return
    require(args.artifacts is not None and args.output is not None, "Collection requires --artifacts and --output")
    require(args.artifacts.is_dir(), f"Artifact directory does not exist: {args.artifacts}")
    studies, ignored, duplicates = discover(args.artifacts)
    summaries = [{"Experiment": name, **row} for name in EXPECTED for row in studies[name]["summary"]]
    truths = [{"Experiment": name, **row} for name in EXPECTED for row in studies[name]["truth"]]
    analysis = summarize(studies, summaries, truths, ignored, duplicates)
    # Import only after data validation, so missing/incomplete artifacts fail clearly.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9})
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "summary.json", summaries)
    write_json(args.output / "truth.json", truths)
    write_json(args.output / "configs.json", {name: studies[name]["config"] for name in EXPECTED})
    write_json(args.output / "shard-configs.json", {name: studies[name]["shard_configs"]
                                                    for name in EXPECTED if "shard_configs" in studies[name]})
    write_json(args.output / "summary-analysis.json", analysis)
    fields = list(dict.fromkeys(key for row in summaries for key in row))
    with (args.output / "summary.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)
    plot_bias_rmse(summaries, args.output, plt)
    plot_nested(analysis["NestedComparison"], args.output, plt)
    print(json.dumps({"Output": str(args.output), "SummaryRows": len(summaries),
                      "PrimaryObservedImprovements": analysis["PrimaryObservedImprovements"],
                      "PrimaryJp3ObservedImprovements": analysis["PrimaryJp3ObservedImprovements"],
                      "Failures": analysis["Failures"],
                      "RequestedComputationTotals": analysis["RequestedComputation"]["Totals"]},
                     ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"collect_indicator_results: ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
