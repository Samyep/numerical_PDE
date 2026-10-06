"""Convenience entry point for the d=20 target-grid pilot and signal gate."""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from run_vb_high_budget import cli_main


if __name__ == "__main__":
    cli_main(forced_stage="pilot")
