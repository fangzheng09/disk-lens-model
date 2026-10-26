"""Weighted posterior summaries for completed production fits."""

import csv
from io import TextIOWrapper
from pathlib import Path

from ._io import write_json
from ._result_files import open_result_file, result_checksum

import numpy as np


def weighted_quantile(values, quantiles, weights=None):
    """Return linearly interpolated quantiles of weighted samples."""
    values = np.asarray(values, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)
    weights = (
        np.ones_like(values) if weights is None else np.asarray(weights, dtype=float)
    )
    if values.ndim != 1 or weights.shape != values.shape or values.size == 0:
        raise ValueError("values and weights must be matching non-empty vectors")
    if (
        not np.all(np.isfinite(values))
        or not np.all(np.isfinite(weights))
        or np.any(weights < 0)
        or weights.sum() <= 0
        or not np.isfinite(weights.sum())
        or not np.all(np.isfinite(quantiles))
        or np.any((quantiles < 0) | (quantiles > 1))
    ):
        raise ValueError("invalid values, weights or quantiles")
    positive = weights > 0
    values, weights = values[positive], weights[positive]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= weights.sum()
    return np.interp(quantiles, cumulative, values, left=values[0], right=values[-1])


def summarize_samples_csv(samples_file, output_file, parameters=None):
    """Write q16/median/q84 using the exact model paths in an AutoFit CSV.

    By default, all model parameter columns are included. An explicit list can
    select a subset; padded statistical headers are stripped before use.
    """
    samples_path = Path(samples_file).resolve()
    with (
        open_result_file(samples_path) as binary,
        TextIOWrapper(binary, encoding="utf-8", newline="") as handle,
    ):
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("samples CSV has no header")
        names = [name.strip() for name in reader.fieldnames]
        if not all(names) or len(names) != len(set(names)):
            raise ValueError("samples CSV has empty or duplicate columns")
        reader.fieldnames = names
        rows = list(reader)
    if not rows:
        raise ValueError("samples CSV is empty")
    columns = set(names)
    if parameters is None:
        statistics = {
            "weight",
            "weights",
            "log_likelihood",
            "log_prior",
            "log_posterior",
        }
        parameters = [name for name in names if name not in statistics]
    if not parameters:
        raise ValueError("samples CSV has no model parameter columns")
    missing = set(parameters) - columns
    if missing:
        raise KeyError(f"missing requested parameter columns: {sorted(missing)}")
    weight_name = next(
        (name for name in ("weight", "weights") if name in columns), None
    )
    weights = (
        None
        if weight_name is None
        else np.asarray([float(row[weight_name]) for row in rows])
    )
    summary = {
        "schema_version": 1,
        "samples_file": str(samples_path),
        "samples_sha256": result_checksum(samples_path),
        "rows": len(rows),
        "weight_column": weight_name,
        "parameters": {},
    }
    for name in parameters:
        values = np.asarray([float(row[name]) for row in rows])
        q16, median, q84 = weighted_quantile(values, [0.15865, 0.5, 0.84135], weights)
        summary["parameters"][name] = {
            "q16": float(q16),
            "median": float(median),
            "q84": float(q84),
        }
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, summary)
    return summary
