"""Run a prepared SPH workflow; install the package in both environments first."""

import sys
from disk_lens.cli import main

if __name__ == "__main__":
    main(["run-main", *sys.argv[1:]])
