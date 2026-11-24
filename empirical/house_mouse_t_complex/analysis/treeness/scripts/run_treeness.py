#!/usr/bin/env python3
import runpy, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
sys.argv=['treeness_pipeline.py','--dataset','house_mouse_t_complex']
runpy.run_path(str(ROOT/'empirical/scripts/treeness_pipeline.py'),run_name='__main__')
