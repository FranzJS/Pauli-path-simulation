#!/usr/bin/env python3
"""Compatibility entry point for the generalized PPS convergence runner."""

from __future__ import annotations

import sys

from run_pps_convergence import main


if __name__ == "__main__":
    # Preserve the former positional support-only interface:
    #   BUILD MAX_SUPPORT MIN_MAGNITUDE PASSES OUTPUT [WORKERS]
    # New options, including --batches, use the generalized interface directly.
    if len(sys.argv) >= 4 and sys.argv[2] not in {"support", "l1"}:
        legacy = sys.argv[1:]
        translated = [
            sys.argv[0], legacy[0], "support", legacy[1], legacy[2],
        ]
        if len(legacy) >= 4:
            translated.extend(["--passes", legacy[3]])
        if len(legacy) >= 5:
            translated.extend(["--output", legacy[4]])
        if len(legacy) >= 6:
            translated.extend(["--workers", legacy[5]])
        if len(legacy) >= 7:
            translated.extend(legacy[6:])
        sys.argv = translated
    raise SystemExit(main())
