"""Read AutoFit scientific results before or after automatic ZIP cleanup."""

from contextlib import contextmanager
import hashlib
from pathlib import Path
from zipfile import ZipFile


@contextmanager
def open_result_file(path):
    """Open a logical result path as binary data without extracting its ZIP."""
    path = Path(path).resolve()
    if path.is_file():
        with path.open("rb") as handle:
            yield handle
        return
    for directory in path.parents:
        if not directory.name:
            continue
        archive = directory.with_name(directory.name + ".zip")
        if not archive.is_file():
            continue
        member = path.relative_to(directory).as_posix()
        with ZipFile(archive) as zipped:
            matches = [info for info in zipped.infolist() if info.filename == member]
            if not matches:
                raise FileNotFoundError(f"missing result {member} in {archive}")
            if len(matches) != 1:
                raise ValueError(f"ambiguous duplicate result {member} in {archive}")
            with zipped.open(matches[0]) as handle:
                yield handle
        return
    raise FileNotFoundError(path)


def result_checksum(path):
    """Hash the file contents, independent of AutoFit's storage format."""
    digest = hashlib.sha256()
    with open_result_file(path) as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
