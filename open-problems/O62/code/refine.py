#!/usr/bin/env python3
"""Newton refinement of one exported fit at increasing working precision.

    python refine.py solutions/p05_H08_012.json

Reports res_2 as a function of the mpmath working precision.  The point of the
experiment:

  * if res_2 keeps falling roughly in step with the working precision, the
    stored point is a good approximation to a true solution of the system;
  * if it stalls at a floor independent of precision, the candidate is a
    near-miss and the floor is the distance to the nearest actual solution.

The Jacobian is rank-deficient (98 of 125 rows at p=5, H=8; the solution set is
22-dimensional), so a square solve is not available.  Steps are Tikhonov-damped
least squares,

    (J^T J + mu I) delta = -J^T F,       mu = 10^(-dps/2),

which is the minimum-norm Gauss-Newton step: mu sits far below the smallest
nonzero singular value (~0.29 at p=5) and far above the arithmetic noise at the
working precision, so it regularises only the exact 22-dimensional null space
and does not perturb the genuine directions.

THIS IS EVIDENCE, NOT PROOF.  See README.md: no amount of precision here
establishes that an exact fit of width 8 exists at p=5.  Approximate algorithms
of a given length can exist with error tending to zero when no exact algorithm
of that length does -- border rank exceeding rank -- and that possibility is
precisely what these numbers cannot exclude.
"""
import argparse
import json

from mpmath import mp, mpf, matrix, lu_solve, norm


def load(path):
    with open(path) as f:
        o = json.load(f)
    p, H = o["p"], o["H"]
    X = [[mpf(float.fromhex(s)) for s in r] for r in o["x_hex"]]
    Y = [[mpf(float.fromhex(s)) for s in r] for r in o["y_hex"]]
    V = [[mpf(float.fromhex(s)) for s in r] for r in o["V_hex"]]
    return o, p, H, X, Y, V


def pack(X, Y, V, p, H):
    w = []
    for j in range(H):
        w += list(X[j])
    for j in range(H):
        w += list(Y[j])
    for c in range(p):
        w += list(V[c])
    return w


def unpack(w, p, H):
    nx = H * p
    X = [[w[j * p + a] for a in range(p)] for j in range(H)]
    Y = [[w[nx + j * p + b] for b in range(p)] for j in range(H)]
    V = [[w[2 * nx + c * H + j] for j in range(H)] for c in range(p)]
    return X, Y, V


def resid(w, p, H):
    X, Y, V = unpack(w, p, H)
    out = []
    for c in range(p):
        for a in range(p):
            for b in range(p):
                acc = mpf(0)
                for j in range(H):
                    s = X[j][a] + Y[j][b]
                    acc += V[c][j] * s * s
                out.append(acc - (mpf(1) if (a + b) % p == c else mpf(0)))
    return out


def jac(w, p, H):
    X, Y, V = unpack(w, p, H)
    nx = H * p
    n = 2 * nx + p * H
    rows = []
    for c in range(p):
        for a in range(p):
            for b in range(p):
                r = [mpf(0)] * n
                for j in range(H):
                    s = X[j][a] + Y[j][b]
                    t = 2 * V[c][j] * s
                    r[j * p + a] += t
                    r[nx + j * p + b] += t
                    r[2 * nx + c * H + j] = s * s
                rows.append(r)
    return rows


def newton(w, p, H, dps, iters=12):
    mp.dps = dps
    w = [mpf(v) for v in w]
    mu = mpf(10) ** (-dps // 2)
    n = len(w)
    best = norm(matrix(resid(w, p, H)))
    for _ in range(iters):
        F = resid(w, p, H)
        J = jac(w, p, H)
        # normal equations with Tikhonov damping
        A = [[mpf(0)] * n for _ in range(n)]
        bb = [mpf(0)] * n
        for row, f in zip(J, F):
            nz = [(k, row[k]) for k in range(n) if row[k] != 0]
            for k, rk in nz:
                bb[k] -= rk * f
                Ak = A[k]
                for l, rl in nz:
                    Ak[l] += rk * rl
        for k in range(n):
            A[k][k] += mu
        d = lu_solve(matrix(A), matrix(bb))
        w = [w[k] + d[k] for k in range(n)]
        r = norm(matrix(resid(w, p, H)))
        if r >= best:
            break
        best = r
    return w, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--dps", type=int, nargs="+", default=[30, 60, 120])
    args = ap.parse_args()

    o, p, H, X, Y, V = load(args.path)
    w = pack(X, Y, V, p, H)
    mp.dps = 30
    r0 = norm(matrix(resid(w, p, H)))
    print(f"=== {o['id']}   p={p}  H={H}   d = 3pH = {3*p*H} ===")
    print(f"  as stored, evaluated at 30 digits:  res_2 = {mp.nstr(r0, 8)}\n")
    print(f"  {'dps':>6}  {'res_2 after Newton':>26}  {'log10(res_2)':>14}")
    for dps in args.dps:
        w, r = newton(w, p, H, dps)
        mp.dps = 20
        lg = mp.log(r, 10) if r > 0 else mpf('-inf')
        print(f"  {dps:>6}  {mp.nstr(r, 10):>26}  {mp.nstr(lg, 6):>14}")

    print("\n  Reading: res_2 falling in step with dps is consistent with a")
    print("  true nearby solution; a precision-independent floor would")
    print("  indicate a near-miss.  Neither outcome is a proof -- see README.md.")


if __name__ == "__main__":
    main()
