#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pure-Python variance decomposition of the SH (CGO+SPO) alpha constraint.

Reimplements phase6_phasor.mc_phasor's math with stdlib only (no numpy), then
toggles each uncertainty group on/off to attribute Var(alpha_13C_OH) for the
SH-only phasor-corrected constraint. This is the REAL decomposition, using the
actual harmonic fits / wetland phasors / sink params, unlike phase13's
hardcoded guess sigmas.
"""
import math
import random

# ---- constants (from phase6_phasor.py) ----
D13C_ATM, DD_ATM = -47.3, -86.0
D13C_WET, D13C_WET_SIG = -62.0, 5.0
Q_TOT, Q_TOT_SIG = 580.0, 50.0
WET_BC_FRAC_SIG = 0.20
ALPHA_D_OH, SIG_ALPHA_D_OH = 1.294, 0.01
A13C_CL, SIG_A13C_CL = 1.066, 0.005
AD_CL, SIG_AD_CL = 1.508, 0.05
A13C_SOIL, AD_SOIL = 1.022, 1.066
A13C_STRAT, AD_STRAT = 1.013, 1.16
F_OH, SIG_F_OH = 0.84, 0.04
F_CL, SIG_F_CL = 0.035, 0.01
F_SOIL, SIG_F_SOIL = 0.06, 0.02
F_STRAT = 0.065

def ci_sig(lo, hi): return (hi - lo) / (2 * 1.96)

def wet_to_mid(bq, cq):
    d = 2 * math.pi * 0.5 / 12.0
    return bq * math.cos(d) + cq * math.sin(d), -bq * math.sin(d) + cq * math.cos(d)

def ratio_to_alpha(R, adoh, foh, fcl, fsoil, fstrat, a13cl, adcl):
    epsD = (foh * (adoh - 1) + fcl * (adcl - 1) + fsoil * (AD_SOIL - 1) + fstrat * (AD_STRAT - 1)) * 1000.0
    eps13_need = R * epsD
    eps13_nonoh = (fcl * (a13cl - 1) + fsoil * (A13C_SOIL - 1) + fstrat * (A13C_STRAT - 1)) * 1000.0
    return 1.0 + (eps13_need - eps13_nonoh) / (foh * 1000.0)

# ---- SH site inputs (phase2 harmonic fits, wetland, dD source) ----
SITES = {
    "CGO": dict(B13=0.0503, C13=0.0498, amp13=[0.0293, 0.1305], pk13=[0.93, 4.56],
                BD=1.5005, CD=2.0113, ampD=[1.8929, 3.682], pkD=[1.0, 3.4],
                BQ=0.0044, CQ=0.0776, dDw=-301.0, dDsig=15.0),
    "SPO": dict(B13=0.0187, C13=0.0418, amp13=[0.0166, 0.0877], pk13=[0.7, 12.17],
                BD=1.4211, CD=1.3061, ampD=[0.8967, 3.3588], pkD=[1.03, 4.02],
                BQ=0.0044, CQ=0.0776, dDw=-301.0, dDsig=15.0),
}
GROUPS = ["obs_amp_13c", "obs_phase_13c", "obs_amp_dD", "obs_phase_dD",
          "wet_flux", "wet_d13C_sig", "wet_dD_sig", "sink"]

def draw_R(site, active, rng):
    """One R_corrected sample for a site, with only `active` groups varying."""
    s = SITES[site]
    on = lambda g: g in active
    A13 = math.hypot(s["B13"], s["C13"]); ph13 = math.atan2(s["C13"], s["B13"])
    AD = math.hypot(s["BD"], s["CD"]); phD = math.atan2(s["CD"], s["BD"])
    sa13, saD = ci_sig(*s["amp13"]), ci_sig(*s["ampD"])
    sp13 = ci_sig(*s["pk13"]) * (2 * math.pi / 12)
    spD = ci_sig(*s["pkD"]) * (2 * math.pi / 12)

    a13 = max(A13 + (rng.gauss(0, sa13) if on("obs_amp_13c") else 0), 1e-6)
    p13 = ph13 + (rng.gauss(0, sp13) if on("obs_phase_13c") else 0)
    b13, c13 = a13 * math.cos(p13), a13 * math.sin(p13)
    aD = max(AD + (rng.gauss(0, saD) if on("obs_amp_dD") else 0), 1e-6)
    pD = phD + (rng.gauss(0, spD) if on("obs_phase_dD") else 0)
    bD, cD = aD * math.cos(pD), aD * math.sin(pD)

    d13cw = D13C_WET + (rng.gauss(0, D13C_WET_SIG) if on("wet_d13C_sig") else 0)
    dDw = s["dDw"] + (rng.gauss(0, s["dDsig"]) if on("wet_dD_sig") else 0)
    if on("wet_flux"):
        qtot = max(rng.gauss(Q_TOT, Q_TOT_SIG), 300) / 12.0
        bq = s["BQ"] * (1 + rng.gauss(0, WET_BC_FRAC_SIG))
        cq = s["CQ"] * (1 + rng.gauss(0, WET_BC_FRAC_SIG))
    else:
        qtot = Q_TOT / 12.0; bq, cq = s["BQ"], s["CQ"]
    bqm, cqm = wet_to_mid(bq, cq)
    zfr_re, zfr_im = bqm / qtot, cqm / qtot
    zs13_re, zs13_im = (d13cw - D13C_ATM) * zfr_re, (d13cw - D13C_ATM) * zfr_im
    zsD_re, zsD_im = (dDw - DD_ATM) * zfr_re, (dDw - DD_ATM) * zfr_im
    zk13_re, zk13_im = b13 - zs13_re, c13 - zs13_im
    zkD_re, zkD_im = bD - zsD_re, cD - zsD_im
    A_k13 = math.hypot(zk13_re, zk13_im); A_kD = math.hypot(zkD_re, zkD_im)
    if A_kD <= 1e-6: return None
    R = A_k13 / A_kD
    return R if (0 < R < 1) else None

def base_sigma(site, N=60000, seed=42):
    rng = random.Random(seed)
    rs = [r for r in (draw_R(site, GROUPS, rng) for _ in range(N)) if r is not None]
    rs.sort()
    lo, hi = rs[int(0.025 * len(rs))], rs[int(0.975 * len(rs))]
    return ci_sig(lo, hi)

BASE_SIG = {c: base_sigma(c) for c in SITES}
w_raw = {c: 1 / BASE_SIG[c] ** 2 for c in SITES}
wsum = sum(w_raw.values())
W = {c: w_raw[c] / wsum for c in SITES}

def alpha_samples(active, N=120000, seed=7):
    rng = random.Random(seed)
    out = []
    for _ in range(N):
        Rs, ok = {}, True
        for c in SITES:
            r = draw_R(c, active, rng)
            if r is None: ok = False; break
            Rs[c] = r
        if not ok: continue
        Rcomb = sum(W[c] * Rs[c] for c in SITES)
        if "sink" in active:
            foh = min(max(rng.gauss(F_OH, SIG_F_OH), 0.5), 0.99)
            fcl = min(max(rng.gauss(F_CL, SIG_F_CL), 0), 0.1)
            fsoil = min(max(rng.gauss(F_SOIL, SIG_F_SOIL), 0), 0.15)
            fstrat = 1 - foh - fcl - fsoil
            adoh = rng.gauss(ALPHA_D_OH, SIG_ALPHA_D_OH)
            a13cl = rng.gauss(A13C_CL, SIG_A13C_CL)
            adcl = rng.gauss(AD_CL, SIG_AD_CL)
            a = ratio_to_alpha(Rcomb, adoh, foh, fcl, fsoil, fstrat, a13cl, adcl)
        else:
            a = ratio_to_alpha(Rcomb, ALPHA_D_OH, F_OH, F_CL, F_SOIL, F_STRAT, A13C_CL, AD_CL)
        out.append(a)
    return out

def stats(xs):
    xs = sorted(xs); n = len(xs)
    mean = sum(xs) / n
    var = sum((x - mean) ** 2 for x in xs) / n
    return dict(median=xs[n // 2], lo=xs[int(0.025 * n)], hi=xs[int(0.975 * n)],
                sd=math.sqrt(var), var=var, n=n)

print(f"Per-site baseline R sigma: " + ", ".join(f"{c}={BASE_SIG[c]:.4f}" for c in SITES))
print(f"Inverse-variance weights: " + ", ".join(f"{c}={W[c]:.3f}" for c in SITES))
print()
allst = stats(alpha_samples(GROUPS))
print(f"SH all-on: median={allst['median']:.5f}  95%CI=[{allst['lo']:.5f},{allst['hi']:.5f}]"
      f"  sd={allst['sd']:.5f}  var={allst['var']:.3e}  (n={allst['n']})")
print(f"reported : median=1.00460 95%CI=[0.99694,1.01578]\n")

rows = []
for g in GROUPS:
    st = stats(alpha_samples([g]))
    rows.append((g, st["sd"], st["var"]))
sv = sum(v for _, _, v in rows)
print(f"{'group':<16}{'sd(alpha)':>11}{'var':>12}{'%oat':>8}")
for g, s, v in sorted(rows, key=lambda r: -r[2]):
    print(f"{g:<16}{s:>11.5f}{v:>12.3e}{100 * v / sv:>7.1f}%")
print(f"{'SUM parts':<16}{math.sqrt(sv):>11.5f}{sv:>12.3e}{100:>7.1f}%")
print(f"{'ALL-ON':<16}{allst['sd']:>11.5f}{allst['var']:>12.3e}")
print(f"total_var/sum_parts = {allst['var']/sv:.3f}  (>1 => amplified by interaction)\n")

coarse = {
    "OBSERVATION total": ["obs_amp_13c", "obs_phase_13c", "obs_amp_dD", "obs_phase_dD"],
    "  d13C amp+phase": ["obs_amp_13c", "obs_phase_13c"],
    "  dD amp+phase": ["obs_amp_dD", "obs_phase_dD"],
    "WETLAND total": ["wet_flux", "wet_d13C_sig", "wet_dD_sig"],
    "SINK conversion": ["sink"],
}
print("Coarse groups (isolated):")
for name, gs in coarse.items():
    st = stats(alpha_samples(gs))
    print(f"  {name:<20} sd={st['sd']:.5f}  var={st['var']:.3e}")
