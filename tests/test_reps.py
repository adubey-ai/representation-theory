"""Tests verifying:
  - Cayley tables are valid groups (associativity, identity, inverses).
  - Representations are homomorphisms.
  - Orthogonality: sum of squared irrep dimensions equals |G|.
  - Regular rep decomposes as sum of d_i copies of each irrep V_i.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from src import (
    cyclic_group, symmetric_s3, klein_four,
    trivial, regular, sign_s3, standard_s3,
    decompose, character_inner_product,
)


def test_group_axioms():
    for G in (cyclic_group(4), klein_four(), symmetric_s3()):
        n = G.order
        for i in range(n):
            assert G.mul(G.identity, i) == i
            assert G.mul(i, G.identity) == i
            assert G.mul(i, G.inverse(i)) == G.identity
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    assert G.mul(G.mul(i, j), k) == G.mul(i, G.mul(j, k))
        print(f"  {G.name}: identity / inverses / associativity ✓")


def test_regular_rep_is_homomorphism():
    G = symmetric_s3()
    rho = regular(G)  # validation runs in __post_init__
    assert rho.dim == G.order
    print(f"  S3 regular rep validates as homomorphism (dim {rho.dim}) ✓")


def test_sum_of_squared_dims_equals_group_order():
    # S_3 irreps: trivial (d=1), sign (d=1), standard (d=2). Sum of squares = 6 = |S_3|.
    irreps = [trivial(symmetric_s3()), sign_s3(), standard_s3()]
    total = sum(r.dim ** 2 for r in irreps)
    assert total == symmetric_s3().order == 6
    print(f"  sum d_i^2 = {total} = |S_3| ✓")


def test_characters_orthonormal():
    irreps = [trivial(symmetric_s3()), sign_s3(), standard_s3()]
    G = symmetric_s3()
    for i, r1 in enumerate(irreps):
        for j, r2 in enumerate(irreps):
            inner = character_inner_product(r1.character(), r2.character(), G)
            expected = 1.0 if i == j else 0.0
            assert abs(inner - expected) < 1e-8, f"<chi_{i}, chi_{j}> = {inner}"
    print("  irrep characters are orthonormal (Schur) ✓")


def test_regular_rep_decomposition():
    """Regular rep has each irrep with multiplicity equal to its dimension."""
    G = symmetric_s3()
    irreps = [trivial(G), sign_s3(), standard_s3()]
    reg = regular(G)
    mults = decompose(reg, irreps)
    expected = [r.dim for r in irreps]
    assert mults == expected, f"got {mults}, expected {expected}"
    print(f"  regular rep of S3 decomposes as  {mults[0]}·triv + {mults[1]}·sign + {mults[2]}·std ✓")


def test_sign_rep_is_irreducible():
    G = symmetric_s3()
    rho = sign_s3()
    inner = character_inner_product(rho.character(), rho.character(), G)
    assert abs(inner - 1.0) < 1e-8
    print("  sign rep is irreducible (<chi, chi> = 1) ✓")


if __name__ == "__main__":
    print("Running representation-theory tests…")
    test_group_axioms()
    test_regular_rep_is_homomorphism()
    test_sum_of_squared_dims_equals_group_order()
    test_characters_orthonormal()
    test_sign_rep_is_irreducible()
    test_regular_rep_decomposition()
    print("All tests passed.")
