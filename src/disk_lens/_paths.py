"""Locate the same scientific resources in a checkout and an installed wheel."""

from pathlib import Path


def resource_path(directory, *parts):
    """Return a source or installed data path, independent of working directory."""
    directories = {"config": "config", "configs": "templates", "catalog": "catalog"}
    if directory not in directories:
        raise ValueError(f"unknown resource directory: {directory}")
    package = Path(__file__).resolve().parent
    checkout = package.parent.parent
    if (checkout / "pyproject.toml").is_file() and (
        checkout / directories[directory]
    ).is_dir():
        base = checkout / directories[directory]
    else:
        base = package / "resources" / directory
    target = base.joinpath(*parts)
    if not target.exists():
        raise FileNotFoundError(f"missing packaged scientific resource: {target}")
    return target
