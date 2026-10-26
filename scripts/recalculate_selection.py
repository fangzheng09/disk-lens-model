"""Calculate one system's particle-mass selection radii with the installed package."""

import sys
from disk_lens.cli import main

if __name__ == "__main__":
    main(["selection-radii", *sys.argv[1:]])
