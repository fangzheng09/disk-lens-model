"""Validated sample-set catalogues for reproducible batch selection."""

from dataclasses import dataclass
import json
from pathlib import Path
import re


_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _text(value, location):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{location} must be a non-empty string")
    return value


def _identifier(value, location):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{location} must be a non-negative integer")
    return value


def _unique(values, location):
    if len(values) != len(set(values)):
        raise ValueError(f"{location} contains duplicates")


@dataclass(frozen=True)
class SampleSet:
    """One named set of subhalo identifiers and its provenance."""

    name: str
    ids: tuple[int, ...]
    source: str
    status: str
    source_sha256: str | None = None


@dataclass(frozen=True)
class SampleCatalog:
    """Top-level catalogue metadata and named sample sets."""

    simulation: str
    snapshot: int
    status: str
    sets: dict[str, SampleSet]
    requires_review: tuple[dict, ...]
    unresolved_sets: tuple[str, ...]
    notes: tuple[str, ...]

    def ids(self, set_name):
        """Return immutable IDs for one explicitly named set."""
        try:
            return self.sets[set_name].ids
        except KeyError as exc:
            choices = ", ".join(sorted(self.sets))
            raise KeyError(
                f"unknown sample set {set_name!r}; available: {choices}"
            ) from exc


def load_sample_catalog(path):
    """Load schema version 1 and reject ambiguous or malformed membership."""
    path = Path(path)
    try:
        raw = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read sample catalogue {path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError("catalogue root must be an object")
    if raw.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    simulation = _text(raw.get("simulation"), "simulation")
    snapshot = _identifier(raw.get("snapshot"), "snapshot")
    status = _text(raw.get("status"), "status")
    raw_sets = raw.get("sets")
    if not isinstance(raw_sets, dict) or not raw_sets:
        raise ValueError("sets must be a non-empty object")

    sets = {}
    for name, item in raw_sets.items():
        _text(name, "set name")
        if not isinstance(item, dict):
            raise ValueError(f"sets.{name} must be an object")
        raw_ids = item.get("ids")
        if not isinstance(raw_ids, list):
            raise ValueError(f"sets.{name}.ids must be a list")
        ids = tuple(_identifier(value, f"sets.{name}.ids") for value in raw_ids)
        _unique(ids, f"sets.{name}.ids")
        digest = item.get("source_sha256")
        if digest is not None and (
            not isinstance(digest, str) or not _SHA256.fullmatch(digest)
        ):
            raise ValueError(f"sets.{name}.source_sha256 must be a lowercase SHA256")
        sets[name] = SampleSet(
            name=name,
            ids=ids,
            source=_text(item.get("source"), f"sets.{name}.source"),
            status=_text(item.get("status"), f"sets.{name}.status"),
            source_sha256=digest,
        )

    raw_review = raw.get("requires_review", [])
    if not isinstance(raw_review, list):
        raise ValueError("requires_review must be a list")
    review = []
    for index, item in enumerate(raw_review):
        if not isinstance(item, dict):
            raise ValueError(f"requires_review[{index}] must be an object")
        review.append(
            {
                "id": _identifier(item.get("id"), f"requires_review[{index}].id"),
                "status": _text(item.get("status"), f"requires_review[{index}].status"),
            }
        )
    _unique([item["id"] for item in review], "requires_review IDs")

    unresolved = raw.get("unresolved_sets", [])
    notes = raw.get("notes", [])
    if not isinstance(unresolved, list):
        raise ValueError("unresolved_sets must be a list")
    if not isinstance(notes, list):
        raise ValueError("notes must be a list")
    unresolved = tuple(_text(value, "unresolved_sets entry") for value in unresolved)
    notes = tuple(_text(value, "notes entry") for value in notes)
    _unique(unresolved, "unresolved_sets")
    return SampleCatalog(
        simulation, snapshot, status, sets, tuple(review), unresolved, notes
    )
