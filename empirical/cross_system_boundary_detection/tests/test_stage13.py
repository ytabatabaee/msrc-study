from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[3]
def test_shared_detector_outputs():
 p=ROOT/'empirical/cross_system_boundary_detection/results/fixed_model_predictions.tsv'; assert p.exists()
 rows=list(csv.DictReader(p.open(),delimiter='\t')); assert len(rows)>=20
 assert {r['system'] for r in rows}=={'Atlantic cod','Fire ants','Anopheles','House mouse'}
def test_multi_interval_output_and_figures():
 assert (ROOT/'empirical/cross_system_boundary_detection/results/multi_interval_predictions.tsv').exists()
 assert (ROOT/'empirical/cross_system_boundary_detection/figures/stage13_boundary_accuracy.png').exists()
