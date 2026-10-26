"""Version requirements and recorded environments for the production workflow."""

from importlib.metadata import PackageNotFoundError, version


REQUIRED_PACKAGES = {"autolens": "2025.1.18.7", "lenstronomy": "1.13.2"}
LEGACY_MOCK_PACKAGES = {"autolens": "2024.11.13.2", "lenstronomy": "1.13.2"}
MOCK_OVER_SAMPLE_SIZE = 1


def package_versions(names):
    """Read installed distribution metadata rather than copying requirement pins."""
    actual = {}
    for name in names:
        try:
            actual[name] = version(name)
        except PackageNotFoundError:
            actual[name] = None
    return actual


def require_environment():
    actual = package_versions(REQUIRED_PACKAGES)
    if actual != REQUIRED_PACKAGES:
        raise RuntimeError(
            f"wrong workflow environment: {actual}; expected {REQUIRED_PACKAGES}"
        )
    return actual


def validate_recorded_packages(manifest, required, *, stage):
    """Reject bundles from another environment even when their files are intact."""
    packages = manifest.get("packages", {})
    if not isinstance(packages, dict):
        packages = {}
    if stage == "mock_truth":
        actual = {
            name: packages.get(name, manifest.get(f"{name}_version"))
            for name in required
        }
    else:
        actual = {name: packages.get(name) for name in required}
    if actual != required:
        raise ValueError(
            f"{stage} bundle environment {actual} differs from plan requirements "
            f"{required}; prepare a new output directory"
        )
