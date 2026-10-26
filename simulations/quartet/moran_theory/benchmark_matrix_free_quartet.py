"""Small reproducible sparse-versus-matrix-free benchmark."""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

try:
    from .matrix_free_quartet import matrix_free_shared_population_quartet_probabilities
    from .moran_quartet_operator import (
        integrated_quartet_probabilities, joint_operator,
        moran_generator, genealogy_generator,
    )
except ImportError:  # pragma: no cover
    from matrix_free_quartet import matrix_free_shared_population_quartet_probabilities
    from moran_quartet_operator import (
        integrated_quartet_probabilities, joint_operator,
        moran_generator, genealogy_generator,
    )


def run(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for M in (10, 20, 40):
        pi = np.zeros(M + 1)
        pi[M // 2] = 1.0
        for method, function in (
            ("sparse_assembled", integrated_quartet_probabilities),
            ("matrix_free", matrix_free_shared_population_quartet_probabilities),
        ):
            start = time.perf_counter()
            result = function(M, 0.17, 0.5, pi)
            elapsed = time.perf_counter() - start
            assembled_nnz = int(joint_operator(M, 0.17).nnz)
            block_nnz = int(moran_generator(M).nnz)
            block_nnz += sum(int(genealogy_generator(4, k, M, 0.17).nnz) for k in range(M + 1))
            rows.append({
                "M": M,
                "method": method,
                "state_dimension": (M + 1) * 19,
                "elapsed_seconds": f"{elapsed:.9f}",
                "assembled_operator_nnz": assembled_nnz,
                "matrix_free_stored_nnz": block_nnz,
                "probability_sum_error": f"{abs(float(result.quartet_probabilities.sum()) - 1.0):.3e}",
            })
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    destination = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results/exact_state_reduction_benchmark.tsv"
    run(destination)
