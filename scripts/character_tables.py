"""Print character tables for S_3, Klein four, and cyclic groups."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from src import (
    cyclic_group, klein_four, symmetric_s3,
    trivial, regular, sign_s3, standard_s3,
)


def print_table(group, irreps: list) -> None:
    classes = group.conjugacy_classes()
    print(f"\n{group.name}   order = {group.order}   conjugacy classes = {len(classes)}")
    reps = [cls[0] for cls in classes]
    sizes = [len(cls) for cls in classes]

    header = ["chi"] + [f"[{group.elements[r]}]^{sz}" for r, sz in zip(reps, sizes)]
    print("  " + "  ".join(f"{h:>8}" for h in header))

    for i, rho in enumerate(irreps):
        chi = rho.character()
        vals = [f"{chi[r].real:>+.2f}" for r in reps]
        print("  " + "  ".join([f"chi_{i}(d={rho.dim})"] + [f"{v:>8}" for v in vals]))


def main():
    # S_3
    G = symmetric_s3()
    print_table(G, [trivial(G), sign_s3(), standard_s3()])

    # Klein four — Dim 1 only, 4 characters; here we show trivial + regular decomposition.
    V = klein_four()
    print(f"\nKlein four V_4: |G| = {V.order}, conjugacy classes = {len(V.conjugacy_classes())}")
    print("  (Abelian: all irreps are 1-dim; 4 characters exist)")

    # Cyclic Z/4
    C4 = cyclic_group(4)
    print(f"\nCyclic C_4: |G| = {C4.order}, regular rep is {C4.order}-dim and decomposes into 4 characters.")


if __name__ == "__main__":
    main()
