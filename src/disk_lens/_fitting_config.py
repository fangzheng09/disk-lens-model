"""Shared five-stage configuration checks without the fitting stack."""

import json
from pathlib import Path

from ._paths import resource_path


STAGE_LIVE_POINTS = {
    "source_lp_n_live": 800,
    "source_pix_1_n_live": 800,
    "source_pix_2_n_live": 1000,
    "light_n_live": 1000,
    "mass_n_live": 1500,
}


def live_points(value):
    if value is None:
        return dict(STAGE_LIVE_POINTS)
    if isinstance(value, bool):
        raise ValueError("n_live must be an integer >= 2 or a stage mapping")
    if isinstance(value, int):
        if value < 2:
            raise ValueError("n_live must be >= 2")
        return {stage: value for stage in STAGE_LIVE_POINTS}
    if isinstance(value, dict) and set(value) == set(STAGE_LIVE_POINTS):
        if all(
            isinstance(n, int) and not isinstance(n, bool) and n >= 2
            for n in value.values()
        ):
            return dict(value)
    raise ValueError("n_live must be an integer >= 2 or a complete stage mapping")


def validate_fitting_config(config):
    """Reject changes to the historical model; allow validated live points."""
    required = {
        "schema_version",
        "experiment",
        "lens_redshift",
        "source_redshift",
        "lens_mge",
        "source_mge",
        "source_pixelization",
        "positions",
        "oversampling",
        "searches",
    }
    if (
        not isinstance(config, dict)
        or set(config) != required
        or config["schema_version"] != 1
        or config["experiment"] != "main"
    ):
        raise ValueError("fitting config is not the main experiment schema 1")
    template = json.loads(resource_path("configs", "fitting.template.json").read_text())
    # The original fit below spells out these choices; reject a config whose
    # model settings differ instead of silently ignoring it.
    for key in required - {"searches"}:
        if config[key] != template[key]:
            raise ValueError(f"original-style fit fixes {key}; configuration differs")
    live_points(config["searches"])
    return config


def load_fitting_config(path):
    path = Path(path).resolve()
    return path, validate_fitting_config(json.loads(path.read_text()))
