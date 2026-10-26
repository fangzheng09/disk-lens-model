"""Measure one stellar SPH aperture M/L with the installed package."""

import sys
from disk_lens.cli import main

if __name__ == "__main__":
    main(["stellar-ml", *sys.argv[1:]])
