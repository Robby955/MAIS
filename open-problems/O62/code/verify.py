#!/usr/bin/env python3
"""Standalone verifier for the exported fits.  Requires numpy only.

    python verify.py solutions/

Reads the IEEE-754 hex fields (not the decimal ones), recomputes the residual
from scratch, and reports the diagnostics described in README.md.  Imports
nothing from the search code.

Exit status is 1 if any solution's res_inf exceeds --tol, so this is usable as
a CI check.

Index conventions (see README.md):
    x, y   H x p    row j = unit, column a = input residue
    V      p x H    row c = OUTPUT coordinate, column j = unit
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

# ----------------------------------------------------------------- defaults
TOL_INF = 1e-12          # res_inf above this fails the run
RANK_CUT = 1e-10         # singular values below cut * sigma_max are zero
DEAD_TOL = 1e-8          # ||A_j|| below this (normalised gauge) is dead
SUPP_TOL = 1e-8          # |s_j| below this is outside the support of s


def load(path):
    with open(path) as f:
        txt = f.read()
    o = json.loads(txt)
    X = np.array([[float.fromhex(s) for s in r] for r in o["x_hex"]])
    Y = np.array([[float.fromhex(s) for s in r] for r in o["y_hex"]])
    V = np.array([[float.fromhex(s) for s in r] for r in o["V_hex"]])
    p, H = o["p"], o["H"]
    if X.shape != (H, p) or Y.shape != (H, p) or V.shape != (p, H):
        raise SystemExit(f"{path}: index convention violated -- got x{X.shape} "
                         f"y{Y.shape} V{V.shape}, expected x({H},{p}) "
                         f"y({H},{p}) V({p},{H})")
    return o, txt, X, Y, V


def dump(obj):
    return json.dumps(obj, indent=2, sort_keys=False) + "\n"


def roundtrip_ok(o, txt):
    """Re-serialising what we read reproduces the file byte for byte."""
    return dump(o) == txt


def target(p):
    T = np.zeros((p, p, p))
    for a in range(p):
        for b in range(p):
            T[(a + b) % p, a, b] = 1.0
    return T


def residual(X, Y, V, p):
    n = (X[:, :, None] + Y[:, None, :]) ** 2
    R = np.einsum('cj,jab->cab', V, n) - target(p)
    return float(np.linalg.norm(R)), float(np.abs(R).max())


# ------------------------------------------------------------------- gauge

def normalise(X, Y, V):
    """Move each unit to the norm-minimising point of its gauge orbit.

    Two per-unit operations leave f_w unchanged:
        scale  (x,y,V) -> (t x, t y, t^-2 V)
        shift  (x,y)   -> (x + s, y - s)
    ||x+s||^2 + ||y-s||^2 is minimised at s = (ybar - xbar)/2, independently of
    t, so shift and scale do not interact.  Then t^2 u^2 + t^-4 v^2 (with
    u^2 = ||(x,y)||^2, v = ||V_j||) is minimised at t^6 = 2 v^2 / u^2, where
    ||V_j||^2 = (1/2)||(x_j,y_j)||^2.
    """
    X, Y, V = X.copy(), Y.copy(), V.copy()
    H = X.shape[0]
    for j in range(H):
        s = (Y[j].mean() - X[j].mean()) / 2.0
        X[j] += s
        Y[j] -= s
        u2 = (X[j] ** 2).sum() + (Y[j] ** 2).sum()
        v2 = (V[:, j] ** 2).sum()
        if u2 > 0 and v2 > 0:
            t = (2.0 * v2 / u2) ** (1.0 / 6.0)
            X[j] *= t
            Y[j] *= t
            V[:, j] *= t ** -2
    return X, Y, V


def pnorm(X, Y, V):
    return float(np.sqrt((X ** 2).sum() + (Y ** 2).sum() + (V ** 2).sum()))


# ------------------------------------------------------------- diagnostics

def units_table(X, Y, V, p, dead_tol):
    """Per-unit Fourier data.  ghat(k) = sum_a g(a) w^{-ka}, w = e^{2 pi i/p};
    this is numpy's fft convention."""
    rows = []
    for j in range(X.shape[0]):
        xh = np.fft.fft(X[j])
        yh = np.fft.fft(Y[j])
        nA = float(np.sqrt((np.abs(xh[1:]) ** 2).sum()))
        nB = float(np.sqrt((np.abs(yh[1:]) ** 2).sum()))
        g = float(xh[0].real + yh[0].real)
        kappa = float((((X[j][:, None] + Y[j][None, :]) ** 2)).sum())
        ident = g * g + nA ** 2 + nB ** 2
        dead = (nA < dead_tol) or (nB < dead_tol)
        rows.append(dict(j=j, nA=nA, nB=nB, g=g, kappa=kappa,
                         mass_lhs=kappa, mass_rhs=ident,
                         mass_err=abs(kappa - ident), dead=dead,
                         nV=float(np.linalg.norm(V[:, j]))))
    return rows


def jacobian(X, Y, V, p):
    """d/dw of F: w -> (f_w(a,b)_c).  Columns ordered [x, y, V] row-major."""
    H = X.shape[0]
    S = X[:, :, None] + Y[:, None, :]                 # H x p x p
    nx = H * p
    J = np.zeros((p * p * p, 2 * nx + p * H))
    idx = np.arange(p)
    for j in range(H):
        for c in range(p):
            base = c * p * p
            # d/dx[j,a]: rows (c,a,b) for all b
            blk = (2.0 * V[c, j] * S[j])              # p x p, [a,b]
            for a in range(p):
                J[base + a * p:base + (a + 1) * p, j * p + a] += blk[a]
            # d/dy[j,b]
            for b in range(p):
                J[base + idx * p + b, nx + j * p + b] += blk[:, b]
            # d/dV[c,j]
            J[base:base + p * p, 2 * nx + c * H + j] = (S[j] ** 2).ravel()
    return J


def rank_report(J, cut):
    sv = np.linalg.svd(J, compute_uv=False)
    thr = sv[0] * cut
    r = int((sv > thr).sum())
    lo = sv[max(0, r - 3):r]
    hi = sv[r:r + 3]
    return r, sv, thr, lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default="solutions")
    ap.add_argument("--tol", type=float, default=TOL_INF)
    ap.add_argument("--rank-cutoff", type=float, default=RANK_CUT)
    ap.add_argument("--dead-tol", type=float, default=DEAD_TOL)
    ap.add_argument("--supp-tol", type=float, default=SUPP_TOL)
    ap.add_argument("--full-spectrum", action="store_true")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.path, "*.json")))
    if not files:
        raise SystemExit(f"no .json files under {args.path}")

    print("Verifying exported fits for f_w(a,b)_c = sum_j V[c,j] "
          "(x[j,a] + y[j,b])^2  vs  T = [a+b == c]")
    print(f"arithmetic: IEEE double throughout; residual summed over all p^3 "
          f"equations (no subsampling)")
    print(f"thresholds: res_inf tol {args.tol:g}, rank cutoff {args.rank_cutoff:g}"
          f" (relative), dead {args.dead_tol:g}, support {args.supp_tol:g}\n")

    hdr = (f"{'id':<14}{'p':>3}{'H':>4}{'res_2':>11}{'res_inf':>11}"
           f"{'|w| found':>11}{'|w| gauge':>11}{'rank':>7}{'of':>6}"
           f"{'dead':>6}{'RT':>4}")
    print(hdr)
    print("-" * len(hdr))

    failures = []
    details = []
    for f in files:
        o, txt, X, Y, V = load(f)
        p, H = o["p"], o["H"]
        r2, ri = residual(X, Y, V, p)
        Xn, Yn, Vn = normalise(X, Y, V)
        r2n, rin = residual(Xn, Yn, Vn, p)
        rows = units_table(Xn, Yn, Vn, p, args.dead_tol)
        ndead = sum(r["dead"] for r in rows)
        J = jacobian(X, Y, V, p)
        r, sv, thr, lo, hi = rank_report(J, args.rank_cutoff)
        rt = roundtrip_ok(o, txt)
        print(f"{o['id']:<14}{p:>3}{H:>4}{r2:>11.2e}{ri:>11.2e}"
              f"{pnorm(X,Y,V):>11.4f}{pnorm(Xn,Yn,Vn):>11.4f}"
              f"{r:>7}{J.shape[0]:>6}{ndead:>6}{'ok' if rt else 'FAIL':>4}")
        if ri > args.tol:
            failures.append((o["id"], "res_inf", ri))
        if not rt:
            failures.append((o["id"], "roundtrip", 0.0))
        if abs(r2n - r2) > 1e-9 * max(1.0, r2):
            failures.append((o["id"], "gauge changed residual", r2n))
        details.append((o, X, Y, V, Xn, Yn, Vn, rows, r, sv, thr, lo, hi,
                        r2, r2n))

    for (o, X, Y, V, Xn, Yn, Vn, rows, r, sv, thr, lo, hi,
         r2, r2n) in details:
        p, H = o["p"], o["H"]
        print(f"\n{'='*74}\n{o['id']}   p={p}  H={H}   d = 3pH = {3*p*H}")
        print(f"  residual as-found      res_2 {r2:.4e}   res_inf "
              f"{residual(X,Y,V,p)[1]:.4e}")
        print(f"  residual after gauge   res_2 {r2n:.4e}   "
              f"(unchanged to {abs(r2n-r2):.2e} -- gauge self-test)")
        print(f"  ||w||  as-found {pnorm(X,Y,V):.6f}   "
              f"gauge-normalised {pnorm(Xn,Yn,Vn):.6f}")
        print(f"\n  per-unit (in the normalised gauge)")
        print(f"    {'j':>3}{'||A_j||':>10}{'||B_j||':>10}{'g_j':>11}"
              f"{'kappa_j':>11}{'mass err':>11}{'||V_j||':>10}  state")
        for q in rows:
            print(f"    {q['j']:>3}{q['nA']:>10.5f}{q['nB']:>10.5f}"
                  f"{q['g']:>11.5f}{q['kappa']:>11.5f}{q['mass_err']:>11.2e}"
                  f"{q['nV']:>10.5f}  {'DEAD' if q['dead'] else 'live'}")
        merr = max(q["mass_err"] for q in rows)
        print(f"    max discrepancy in kappa_j = g_j^2 + ||A_j||^2 + ||B_j||^2:"
              f"  {merr:.3e}")

        s = V.sum(axis=0)
        supp = [j for j in range(H) if abs(s[j]) > args.supp_tol]
        deadset = {q["j"] for q in rows if q["dead"]}
        print(f"\n  annihilator s_j = sum_c V[c,j]")
        print("    " + "  ".join(f"{v:+.3e}" for v in s))
        print(f"    support (|s_j| > {args.supp_tol:g}): {supp}"
              f"   ({len(supp)} of {H})")
        if supp:
            allsupp_dead = all(j in deadset for j in supp)
            print(f"    supporting units dead? {allsupp_dead}"
                  f"   dead units: {sorted(deadset)}")
            holds = "HOLDS" if (len(supp) == 1 and allsupp_dead) \
                else "DOES NOT HOLD as stated"
            print(f"    note claim (s supported on ONE unit, and that unit "
                  f"dead): {holds}")

        print(f"\n  Jacobian of F: R^{3*p*H} -> R^{p**3}")
        print(f"    rank {r} of {p**3} rows / {3*p*H} columns   "
              f"(cutoff {thr:.3e} = {args.rank_cutoff:g} * sigma_max)")
        print(f"    corank vs #equations {p**3 - r}   nullity {3*p*H - r}")
        print(f"    sigma just above cutoff: "
              + "  ".join(f"{v:.4e}" for v in lo))
        print(f"    sigma just below cutoff: "
              + ("  ".join(f"{v:.4e}" for v in hi) if len(hi) else "(none)"))
        gap = (lo[-1] / hi[0]) if len(hi) and hi[0] > 0 else float("inf")
        print(f"    ratio across the cutoff: {gap:.3e}"
              f"   (large => the cutoff is not doing delicate work)")
        if args.full_spectrum:
            print("    full spectrum:")
            for i in range(0, len(sv), 6):
                print("      " + "  ".join(f"{v:.4e}" for v in sv[i:i + 6]))

    print(f"\n{'='*74}")
    if failures:
        print("FAILURES:")
        for sid, what, val in failures:
            print(f"  {sid}: {what} = {val:g}")
        sys.exit(1)
    print(f"all {len(files)} solutions pass "
          f"(res_inf <= {args.tol:g}, hex round-trip exact, "
          f"gauge preserves residual)")
    print("\nNOTE: these are numerical candidates at the stated residual, not "
          "certified\nexact fits.  See README.md, section 'What this does not "
          "show'.")


if __name__ == "__main__":
    main()
