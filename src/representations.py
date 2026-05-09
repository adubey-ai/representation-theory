"""Matrix representations of finite groups, characters, and decomposition.

Defines:
  - :class:`Representation` validating rho(g)rho(h) = rho(gh).
  - Common reps: trivial, regular, sign (S_n), permutation, standard (S_3).
  - Character machinery: characters from matrix traces and inner product.
  - Decomposition into irreducibles via inner-product multiplicities.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .groups import FiniteGroup


@dataclass
class Representation:
    group: FiniteGroup
    dim: int
    matrices: list[np.ndarray]  # matrices[g] = rho(e_g)

    def __post_init__(self) -> None:
        for i, m in enumerate(self.matrices):
            if m.shape != (self.dim, self.dim):
                raise ValueError(f"matrix {i} has shape {m.shape}, expected ({self.dim},{self.dim})")
        self.validate()

    def validate(self, tol: float = 1e-8) -> None:
        n = self.group.order
        # rho(e) = I
        I = np.eye(self.dim)
        if not np.allclose(self.matrices[self.group.identity], I, atol=tol):
            raise ValueError("rho(e) != I")
        for i in range(n):
            for j in range(n):
                lhs = self.matrices[i] @ self.matrices[j]
                rhs = self.matrices[self.group.mul(i, j)]
                if not np.allclose(lhs, rhs, atol=tol):
                    raise ValueError(
                        f"homomorphism violated: rho({i}) rho({j}) != rho({i}*{j})"
                    )

    def character(self) -> np.ndarray:
        return np.array([np.trace(m) for m in self.matrices], dtype=complex)


def trivial(group: FiniteGroup) -> Representation:
    n = group.order
    mats = [np.array([[1.0]]) for _ in range(n)]
    return Representation(group, dim=1, matrices=mats)


def regular(group: FiniteGroup) -> Representation:
    """Left regular representation: each g permutes group elements by left multiplication."""
    n = group.order
    mats: list[np.ndarray] = []
    for g in range(n):
        M = np.zeros((n, n))
        for h in range(n):
            gh = group.mul(g, h)
            M[gh, h] = 1.0
        mats.append(M)
    return Representation(group, dim=n, matrices=mats)


def sign_s3() -> Representation:
    """S_3 sign representation: +1 on e, r, r^2; -1 on transpositions."""
    from .groups import symmetric_s3
    G = symmetric_s3()
    # Name order: e, r, r2, s, rs, r2s  -> signs: +1 +1 +1 -1 -1 -1
    signs = [+1.0, +1.0, +1.0, -1.0, -1.0, -1.0]
    mats = [np.array([[s]]) for s in signs]
    return Representation(G, dim=1, matrices=mats)


def standard_s3() -> Representation:
    """The 2-dim irreducible of S_3 via its action on the plane x+y+z=0."""
    from .groups import symmetric_s3
    G = symmetric_s3()
    # Acting on R^3 by permutation, restricted to x+y+z=0 basis.
    # Basis of the 2D subspace: v1 = (1,-1,0), v2 = (0,1,-1).
    perms = [
        (0, 1, 2), (1, 2, 0), (2, 0, 1),
        (1, 0, 2), (2, 1, 0), (0, 2, 1),
    ]
    V = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])
    Vplus = np.linalg.pinv(V)  # 3x2
    mats: list[np.ndarray] = []
    for perm in perms:
        P = np.zeros((3, 3))
        for i, j in enumerate(perm):
            P[j, i] = 1.0
        rho = V @ P @ Vplus
        mats.append(rho)
    return Representation(G, dim=2, matrices=mats)


def character_inner_product(chi1: np.ndarray, chi2: np.ndarray, group: FiniteGroup) -> complex:
    n = group.order
    total = sum(chi1[g] * np.conjugate(chi2[g]) for g in range(n))
    return total / n


def decompose(rho: Representation, irreps: list[Representation]) -> list[int]:
    """Return multiplicities of each irrep in rho via character inner products.

    All irreps must be of the same group as rho.
    """
    chi = rho.character()
    mults: list[int] = []
    for irr in irreps:
        chi_i = irr.character()
        m = character_inner_product(chi, chi_i, rho.group)
        m_real = m.real
        if abs(m.imag) > 1e-6 or abs(m_real - round(m_real)) > 1e-6:
            raise ValueError(f"non-integer multiplicity {m}; irreps list incomplete or wrong?")
        mults.append(int(round(m_real)))
    return mults
