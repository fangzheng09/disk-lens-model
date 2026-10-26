#!/usr/bin/env python3
"""Report bright-arc candidates and check stored coordinates in a mock."""

import sys

from disk_lens.cli import main

if __name__ == "__main__":
    main(["check-positions", *sys.argv[1:]])
