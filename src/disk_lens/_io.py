"""File checksums, atomic JSON records and completed imaging-bundle checks."""

import hashlib
import json
from pathlib import Path
import tempfile


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, data):
    """Replace a JSON record only after its complete contents have been written."""
    path = Path(path)
    content = json.dumps(data, indent=2, allow_nan=False) + "\n"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def validate_bundle(directory, stage, *, config_file=None, input_sha256=None):
    """Verify required products and their recorded provenance before reuse.

    A damaged or incomplete directory is left intact and reported to the caller.
    The checks apply to SPH foreground bundles and their production mocks.
    """
    required_files = {
        "lightmap": {"light.fits"},
        "mock": {
            "data.fits",
            "noise_map.fits",
            "psf.fits",
            "source_convolved.fits",
            "foreground_convolved.fits",
            "model_image.fits",
            "source_native.fits",
            "stellar_deflections_yx.npy",
        },
        "mock_truth": {"truth.json", "unit_stellar_deflections_unmasked_yx.npy"},
    }
    if stage not in required_files:
        raise ValueError(f"unknown bundle stage: {stage}")
    directory = Path(directory)
    try:
        manifest = json.loads((directory / "manifest.json").read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"missing or invalid bundle manifest: {directory}") from exc
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != 1
        or manifest.get("stage") != stage
        or manifest.get("status") != "complete"
    ):
        raise ValueError(f"input must be a complete {stage} bundle: {directory}")
    files = manifest.get("files_sha256")
    if not isinstance(files, dict) or not required_files[stage].issubset(files):
        raise ValueError(f"bundle is missing required file checksums: {directory}")
    for name, digest in files.items():
        if not name or Path(name).name != name:
            raise ValueError(f"invalid bundle filename: {name!r}")
        path = directory / name
        if not path.is_file() or checksum(path) != digest:
            raise ValueError(f"bundle file missing or checksum mismatch: {path}")
    if config_file is not None:
        provenance = manifest.get("provenance", {})
        recorded_config = manifest.get("config_sha256")
        if recorded_config is None and isinstance(provenance, dict):
            recorded_config = provenance.get("config_sha256")
        if recorded_config != checksum(config_file):
            raise ValueError(f"bundle configuration checksum mismatch: {directory}")
        try:
            saved = json.loads((directory / "config.json").read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"missing or invalid saved configuration: {directory}"
            ) from exc
        if saved != json.loads(Path(config_file).read_text()):
            raise ValueError(f"bundle saved configuration differs: {directory}")
    for name, digest in (input_sha256 or {}).items():
        if manifest.get(name) != digest:
            raise ValueError(f"bundle input checksum mismatch ({name}): {directory}")
    return manifest
