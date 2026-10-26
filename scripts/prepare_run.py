"""Prepare a per-system SPH workflow using the installed package."""

import sys
from disk_lens.cli import main

if __name__ == "__main__":
    main(["prepare-main", *sys.argv[1:]])
