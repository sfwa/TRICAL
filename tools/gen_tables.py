#!/usr/bin/env python3
"""Generate the constant index tables used by src/fit.c.

TRICAL accumulates the sums of all monomials x^a y^b z^c of total degree <= 4
(35 numbers). Every entry of the 10x10 quadric moment matrix is one of those
sums times a small integer, so the tables below are all that is needed to
work with it. Run this script and paste its output into src/fit.c if the
monomial ordering ever changes; test/test_fit.cpp checks the tables against
a brute-force computation.
"""
from itertools import combinations_with_replacement

# Monomials in graded lexicographic order: 1, x, y, z, xx, xy, xz, yy, ...
MONO = []
for degree in range(5):
    for combo in combinations_with_replacement(range(3), degree):
        e = [0, 0, 0]
        for axis in combo:
            e[axis] += 1
        MONO.append(tuple(e))
INDEX = {e: i for i, e in enumerate(MONO)}
assert len(MONO) == 35

# Quadric basis xi(x) = (xx, yy, zz, 2yz, 2xz, 2xy, 2x, 2y, 2z, 1)
XI_EXP = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (0, 1, 1), (1, 0, 1), (1, 1, 0),
          (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)]
XI_COEF = [1, 1, 1, 2, 2, 2, 2, 2, 2, 1]


def c_array(rows, per_line=10, indent="    "):
    out = []
    for row in rows:
        out.append(indent + "{ " + ", ".join("%2d" % v for v in row) + " },")
    return "\n".join(out)


print("/* exponents (a, b, c) of monomial k = x^a y^b z^c */")
print("static const unsigned char TRICAL_MONO_EXP[TRICAL_NUM_MOMENTS][3] = {")
print(c_array(MONO))
print("};\n")

print("/* index of the monomial that is the product of quadric basis terms i, j */")
print("static const unsigned char TRICAL_M_INDEX[10][10] = {")
print(c_array([[INDEX[tuple(a + b for a, b in zip(XI_EXP[i], XI_EXP[j]))]
                for j in range(10)] for i in range(10)]))
print("};\n")

print("/* integer coefficient of quadric basis term i */")
print("static const unsigned char TRICAL_XI_COEF[10] = { "
      + ", ".join(str(c) for c in XI_COEF) + " };\n")

print("/* index of monomial x^a y^b z^c for a + b + c <= 4 (255 otherwise) */")
print("static const unsigned char TRICAL_MONO_INDEX[5][5][5] = {")
for a in range(5):
    print("    {")
    for b in range(5):
        row = [INDEX.get((a, b, c), 255) for c in range(5)]
        print("        { " + ", ".join("%3d" % v for v in row) + " },")
    print("    },")
print("};")
