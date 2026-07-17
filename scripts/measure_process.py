#!/usr/bin/env python3
"""Run one command and report its child-process peak RSS in KiB."""

from __future__ import annotations

import resource
import subprocess
import sys


def main() -> int:
    if len(sys.argv) < 2:
        raise ValueError("a command is required")
    completed = subprocess.run(
        sys.argv[1:],
        capture_output=True,
        text=True,
    )
    sys.stdout.write(completed.stdout)
    sys.stderr.write(completed.stderr)
    peak_rss_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    print(f"__MAX_RSS_KB__={peak_rss_kb}", file=sys.stderr)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
