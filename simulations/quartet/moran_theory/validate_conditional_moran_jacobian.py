"""Validate and summarize the conditional Moran Jacobian certificate."""

from __future__ import annotations

import csv
from pathlib import Path

from .certify_conditional_moran_jacobian import run


def validate(output: Path) -> None:
    run(output)
    rows = list(csv.DictReader((output / "jacobian_certificate.tsv").open(), delimiter="\t"))
    if len(rows) != 6:
        raise RuntimeError("expected three finite-difference checks for each rooted shape")
    for row in rows:
        if abs(float(row["analytic_det"])) <= 1e-4:
            raise RuntimeError("determinant is not separated from zero numerically")
        if float(row["max_derivative_error"]) >= 1e-7:
            raise RuntimeError("analytic derivative check failed")


if __name__ == "__main__":
    validate(Path(__file__).parent / "results" / "conditional_moran_rigorous_coverage")
