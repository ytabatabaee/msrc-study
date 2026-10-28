"""Analytic duration derivatives for the exact conditional Moran quartet law.

The certificate is for the fixed finite-M family used by the local-image
analysis.  It differentiates the existing recursion using
``d exp(A t) v / d t = A exp(A t) v`` and the quotient rule.  The script also
compares the analytic derivatives with centered finite differences.  The
result is a strong numerical certificate; it does not claim interval-arithmetic
certification of the matrix exponential.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .full_species_tree import inheritance_kernel
from .moran_quartet_operator import (
    BACKGROUND_INDEX,
    BACKGROUND_STATES,
    PAIR_BACKGROUND_STATES,
    TRIPLE_BACKGROUND_INDEX,
    dual_lift,
    joint_operator_n,
    moran_generator,
    residual_vectors,
)


def _lift(M, rho, duration, value, n, selection=0.0, derivative=None,
          differentiate_duration=False):
    """Propagate a value and input derivative through one branch.

    ``differentiate_duration`` is separate from ``derivative``: when a
    terminal derivative is being propagated through an older branch, that
    branch's duration is held fixed and the ``A exp(A t)v`` term must not be
    added.
    """
    # ``joint_operator_n`` already contains the forward Moran block and the
    # transposed conditional genealogy blocks used by ``dual_lift``.
    A = joint_operator_n(M, rho, n, selection).tocsr()
    flat = np.asarray(value, dtype=float).reshape(-1)
    out = expm_multiply(A * duration, flat).reshape(np.asarray(value).shape)
    if derivative is None:
        din = np.zeros_like(flat)
    else:
        din = np.asarray(derivative, dtype=float).reshape(-1)
    dout = expm_multiply(A * duration, din)
    if differentiate_duration:
        dout = dout + A @ out.reshape(-1)
    dout = dout.reshape(out.shape)
    return out, dout


def _one(M, bit, duration, rho, selection=0.0, derivative_duration=False):
    value = np.zeros((M + 1, 2), dtype=float)
    for k in range(M + 1):
        value[k, bit] = (k / M) if bit else (1.0 - k / M)
    out, dout = _lift(M, rho, duration, value, 1, selection,
                      differentiate_duration=True)
    return out, dout if derivative_duration else np.zeros_like(out)


def _pair_with_terminal1_derivative(M_node, M_sister, M1, M2, lengths, sister,
                                    rho, bits, selection=0.0):
    u1, du1 = _one(M1, bits[0], lengths[0], rho, selection, True)
    u2, _ = _one(M2, bits[1], lengths[1], rho, selection, False)
    H1, H2 = inheritance_kernel(M_sister, M1), inheritance_kernel(M_sister, M2)
    young = np.zeros((M_sister + 1, 6), dtype=float)
    dy = np.zeros_like(young)
    for k in range(M_sister + 1):
        for b1, b2 in PAIR_BACKGROUND_STATES:
            a = H1[k] @ u1[:, b1]
            da = H1[k] @ du1[:, b1]
            b = H2[k] @ u2[:, b2]
            young[k, 2 * b1 + b2] = a * b
            dy[k, 2 * b1 + b2] = da * b
    old, dold = _lift(M_sister, rho, sister, young, 2, selection, dy)
    Hnode = inheritance_kernel(M_node, M_sister)
    return Hnode @ old, Hnode @ dold


def _pair_no_derivative(M_node, M_sister, M1, M2, lengths, sister,
                        rho, bits, selection=0.0):
    value, _ = _pair_with_terminal1_derivative(
        M_node, M_sister, M1, M2, lengths, sister, rho, bits, selection)
    # The helper above differentiates terminal 1, but the value is the same.
    return value


def _root_initial(pi, h):
    return np.kron(pi, h)


def _conditional_from_joint(joint, d_joint):
    den = float(np.sum(joint))
    dden = float(np.sum(d_joint))
    q = np.asarray(joint) / den
    dq = (np.asarray(d_joint) * den - np.asarray(joint) * dden) / den**2
    return q, dq, den


def balanced_conditional_derivatives(*, z="0101", M=10, rho=0.1,
                                     root_length=1.0, terminal1_length=0.2,
                                     selection=0.0):
    """Return q and exact derivatives wrt root and terminal-1 durations."""
    if z != "0101":
        raise ValueError("the certificate currently fixes Z=0101")
    pi = np.eye(M + 1)[M // 2]
    # Root count distribution is differentiated independently of genealogy.
    Q = moran_generator(M, selection).transpose().tocsr()
    pi_root = expm_multiply(Q * root_length, pi)
    dpi_root = Q @ pi_root

    left, dleft = _pair_with_terminal1_derivative(
        M, M, M, M, (terminal1_length, .1), .2, rho, (0, 1), selection)
    right = _pair_no_derivative(M, M, M, M, (.2, .3), .3, rho, (0, 1), selection)
    left_total, dleft_total = left.sum(axis=1), dleft.sum(axis=1)
    right_total = right.sum(axis=1)
    left_surv, dleft_surv = left[:, :4].sum(axis=1), dleft[:, :4].sum(axis=1)
    right_surv = right[:, :4].sum(axis=1)
    early = np.array([np.dot(pi_root, left_total * right_total - left_surv * right_surv), 0.0, 0.0])
    dearly_root = np.array([np.dot(dpi_root, left_total * right_total - left_surv * right_surv), 0.0, 0.0])
    dearly_t1 = np.array([np.dot(pi_root, dleft_total * right_total - dleft_surv * right_surv), 0.0, 0.0])

    F = np.zeros((M + 1, 19), dtype=float)
    dF = np.zeros_like(F)
    for il, zl in enumerate(PAIR_BACKGROUND_STATES):
        for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
            bg = BACKGROUND_STATES.index(zl + zr)
            F[:, bg] = left[:, il] * right[:, ir]
            dF[:, bg] = dleft[:, il] * right[:, ir]
    A4 = joint_operator_n(M, rho, 4, selection)
    propagated = expm_multiply(A4 * root_length, F.reshape(-1)).reshape(F.shape)
    dprop_t1 = expm_multiply(A4 * root_length, dF.reshape(-1)).reshape(F.shape)
    dprop_root = (A4 @ propagated.reshape(-1)).reshape(F.shape)
    initial = np.column_stack([_root_initial(pi, h) for h in residual_vectors()])
    root = initial.T @ propagated.reshape(-1)
    droot_t1 = initial.T @ dprop_t1.reshape(-1)
    droot_root = initial.T @ dprop_root.reshape(-1)
    J = early + root
    dJ_root = dearly_root + droot_root
    dJ_t1 = dearly_t1 + droot_t1
    q, dq_root, den = _conditional_from_joint(J, dJ_root)
    _, dq_t1, _ = _conditional_from_joint(J, dJ_t1)
    return q, np.column_stack((dq_root[:2], dq_t1[:2])), J, den


def _three_with_derivative(Mroot, Mint, pair, dpair, line3, intermediate, rho, selection):
    young = np.zeros((Mint + 1, 11), dtype=float)
    dy = np.zeros_like(young)
    line3_mass = line3.sum(axis=1)
    for k in range(Mint + 1):
        young[k, 8] = (pair[k, 4] + pair[k, 5]) * line3_mass[k]
        dy[k, 8] = (dpair[k, 4] + dpair[k, 5]) * line3_mass[k]
        for pair_index, pair_bits in enumerate(PAIR_BACKGROUND_STATES[:4]):
            for b3 in (0, 1):
                idx = TRIPLE_BACKGROUND_INDEX[pair_bits + (b3,)]
                young[k, idx] = pair[k, pair_index] * line3[k, b3]
                dy[k, idx] = dpair[k, pair_index] * line3[k, b3]
    old, dold = _lift(Mint, rho, intermediate, young, 3, selection, dy)
    H = inheritance_kernel(Mroot, Mint)
    return H @ old, H @ dold


def pectinate_conditional_derivatives(*, z="0101", M=10, rho=0.1,
                                      root_length=1.0, terminal1_length=0.2,
                                      selection=0.0):
    """Return q and exact derivatives for the pectinate recursion."""
    if z != "0101":
        raise ValueError("the certificate currently fixes Z=0101")
    pi = np.eye(M + 1)[M // 2]
    Q = moran_generator(M, selection).transpose().tocsr()
    pi_root = expm_multiply(Q * root_length, pi)
    dpi_root = Q @ pi_root
    pair, dpair = _pair_with_terminal1_derivative(
        M, M, M, M, (terminal1_length, .1), .3, rho, (0, 1), selection)
    line3_terminal, _ = _one(M, 0, .2, rho, selection, False)
    line3 = inheritance_kernel(M, M) @ line3_terminal
    triple, dtriple = _three_with_derivative(M, M, pair, dpair, line3, .2, rho, selection)
    line4_terminal, _ = _one(M, 1, .3, rho, selection, False)
    line4 = inheritance_kernel(M, M) @ line4_terminal
    mass4 = line4.sum(axis=1)
    early = np.array([
        np.dot(pi_root, triple[:, 8 + i] * mass4) for i in range(3)
    ])
    dearly_root = np.array([
        np.dot(dpi_root, triple[:, 8 + i] * mass4) for i in range(3)
    ])
    dearly_t1 = np.array([
        np.dot(pi_root, dtriple[:, 8 + i] * mass4) for i in range(3)
    ])
    F = np.zeros((M + 1, 19), dtype=float)
    dF = np.zeros_like(F)
    for z3, idx3 in TRIPLE_BACKGROUND_INDEX.items():
        if idx3 >= 8:
            continue
        for b4 in (0, 1):
            bg = BACKGROUND_INDEX[z3 + (b4,)]
            F[:, bg] = triple[:, idx3] * line4[:, b4]
            dF[:, bg] = dtriple[:, idx3] * line4[:, b4]
    A4 = joint_operator_n(M, rho, 4, selection)
    prop = expm_multiply(A4 * root_length, F.reshape(-1)).reshape(F.shape)
    dprop_t1 = expm_multiply(A4 * root_length, dF.reshape(-1)).reshape(F.shape)
    dprop_root = (A4 @ prop.reshape(-1)).reshape(F.shape)
    initial = np.column_stack([_root_initial(pi, h) for h in residual_vectors()])
    root = initial.T @ prop.reshape(-1)
    droot_t1 = initial.T @ dprop_t1.reshape(-1)
    droot_root = initial.T @ dprop_root.reshape(-1)
    J = early + root
    dJ_root = dearly_root + droot_root
    dJ_t1 = dearly_t1 + droot_t1
    q, dq_root, den = _conditional_from_joint(J, dJ_root)
    _, dq_t1, _ = _conditional_from_joint(J, dJ_t1)
    return q, np.column_stack((dq_root[:2], dq_t1[:2])), J, den


def _finite_difference(fun, x, h=1e-5):
    x = np.asarray(x, dtype=float)
    cols = []
    for j in range(2):
        e = np.zeros(2); e[j] = h
        cols.append((fun(x + e)[:2] - fun(x - e)[:2]) / (2 * h))
    return np.column_stack(cols)


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for shape, analytic, numeric in (
        ("balanced", balanced_conditional_derivatives, lambda x: __import__(
            "simulations.quartet.moran_theory.analyze_conditional_moran_image", fromlist=["balanced_q"]
        ).balanced_q("0101", x)),
        ("pectinate", pectinate_conditional_derivatives, lambda x: __import__(
            "simulations.quartet.moran_theory.analyze_conditional_moran_image", fromlist=["pectinate_q"]
        ).pectinate_q("0101", x, M=10)),
    ):
        q, jac, _, pz = analytic()
        fd = _finite_difference(numeric, (1.0, .2), 1e-5)
        det = float(np.linalg.det(jac))
        for h in (2e-4, 1e-4, 5e-5):
            fd_h = _finite_difference(numeric, (1.0, .2), h)
            rows.append((shape, h, det, pz, q[0], q[1], q[2], *jac.reshape(-1), *fd_h.reshape(-1),
                         float(np.max(np.abs(jac - fd_h))), float(np.linalg.det(fd_h))))
    with (output / "jacobian_certificate.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(("shape", "finite_difference_h", "analytic_det", "P_Z", "q1", "q2", "q3",
                    "dq1_droot", "dq1_dterminal1", "dq2_droot", "dq2_dterminal1",
                    "fd_dq1_droot", "fd_dq1_dterminal1", "fd_dq2_droot", "fd_dq2_dterminal1",
                    "max_derivative_error", "fd_det"))
        w.writerows(rows)


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "conditional_moran_coverage")
