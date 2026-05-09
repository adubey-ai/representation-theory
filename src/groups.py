"""Small finite groups defined by their Cayley tables."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FiniteGroup:
    """A finite group with explicit Cayley table (rows/cols indexed by element id)."""
    name: str
    elements: tuple[str, ...]
    cayley: tuple[tuple[int, ...], ...]  # cayley[i][j] = index of e_i * e_j
    identity: int = 0

    @property
    def order(self) -> int:
        return len(self.elements)

    def mul(self, i: int, j: int) -> int:
        return self.cayley[i][j]

    def inverse(self, i: int) -> int:
        for j in range(self.order):
            if self.mul(i, j) == self.identity:
                return j
        raise ValueError(f"no inverse for element {i}")

    def conjugacy_classes(self) -> list[list[int]]:
        classes: list[list[int]] = []
        seen: set[int] = set()
        for x in range(self.order):
            if x in seen:
                continue
            cls = set()
            for g in range(self.order):
                ginv = self.inverse(g)
                cls.add(self.mul(self.mul(g, x), ginv))
            seen |= cls
            classes.append(sorted(cls))
        return classes


def cyclic_group(n: int) -> FiniteGroup:
    """Z/nZ under addition mod n."""
    cayley = tuple(tuple((i + j) % n for j in range(n)) for i in range(n))
    names = tuple(f"e{i}" if i else "e" for i in range(n))
    return FiniteGroup(name=f"C{n}", elements=names, cayley=cayley, identity=0)


def symmetric_s3() -> FiniteGroup:
    """S_3 with elements e, r, r^2, s, rs, r^2s where r^3 = e, s^2 = e, srs = r^{-1}."""
    # Permutations of {0,1,2}:
    perms = [
        (0, 1, 2),  # e
        (1, 2, 0),  # r
        (2, 0, 1),  # r^2
        (1, 0, 2),  # s
        (2, 1, 0),  # rs
        (0, 2, 1),  # r^2 s
    ]
    names = ("e", "r", "r2", "s", "rs", "r2s")
    n = len(perms)
    table = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            composed = tuple(perms[i][perms[j][k]] for k in range(3))
            table[i][j] = perms.index(composed)
    cayley = tuple(tuple(row) for row in table)
    return FiniteGroup(name="S3", elements=names, cayley=cayley, identity=0)


def klein_four() -> FiniteGroup:
    """Klein four-group V = Z/2 x Z/2."""
    # Elements: 00, 01, 10, 11 under componentwise XOR.
    table = [[(i ^ j) for j in range(4)] for i in range(4)]
    cayley = tuple(tuple(row) for row in table)
    return FiniteGroup(name="V4", elements=("e", "a", "b", "ab"), cayley=cayley, identity=0)
