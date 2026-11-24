#!/usr/bin/env python3
"""Compatibility entry point for all arrangement-specific diversity outputs."""
import subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
for script in ['infer_arrangements.py','estimate_tracks.py','plot_tracks.py','compare_inside_outside.py','summarize.py']:
    subprocess.run([sys.executable,str(HERE/script)],check=True)
