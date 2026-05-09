# Representation Theory of Finite Groups

A compact library for studying **finite-dimensional matrix representations**
of finite groups — with tests that verify classical results from the theory
(Schur orthogonality, Burnside's lemma on squared dimensions, decomposition
of the regular representation).

## What's here

- `src/groups.py` — `FiniteGroup` via explicit Cayley table, plus
  constructors for `C_n`, the Klein four-group, and `S_3`.
- `src/representations.py`:
    - `Representation` class that *validates* the homomorphism property
      `ρ(g)ρ(h) = ρ(gh)` at construction time.
    - Built-in reps: trivial, left-regular, sign (on `S_3`), standard
      2-dim (on `S_3`).
    - Characters via matrix traces; inner product; decomposition of a
      representation into irreducibles by computing character multiplicities.
- `tests/test_reps.py` — verifies every theorem numerically.
- `scripts/character_tables.py` — prints character tables.

## Run

```bash
python -m tests.test_reps              # all tests
python scripts/character_tables.py     # print S_3, V_4, C_4 tables
```

Deps: `numpy` only.

## What the tests actually prove (numerically)

| Test | Theorem |
| --- | --- |
| `test_group_axioms` | Identity, inverses, associativity hold in each Cayley table |
| `test_regular_rep_is_homomorphism` | `ρ(g)ρ(h) = ρ(gh)` for every pair in `S_3` |
| `test_sum_of_squared_dims_equals_group_order` | `∑ d_i² = |G|` (Burnside) |
| `test_characters_orthonormal` | `⟨χ_i, χ_j⟩ = δ_{ij}` (Schur orthogonality) |
| `test_regular_rep_decomposition` | Regular rep = ⊕ (d_i copies of V_i) |
| `test_sign_rep_is_irreducible` | `⟨χ, χ⟩ = 1` iff irreducible |

## Example — where this kind of thing shows up

Representation theory models **symmetries** — in crystallography, molecular
orbitals (Hückel, MO theory), selection rules in spectroscopy, and quantum
mechanics generally. Decomposing a reducible representation into
irreducibles is the algebraic version of "block-diagonalizing" a physical
system along its symmetry axes.
