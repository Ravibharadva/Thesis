#!/usr/bin/env python3
"""
extract.py
----------
Build Raissi-format .mat datasets (x_star, y_star, t_star, U_star, V_star, P_star)
from Nektar++ .vtu checkpoint output, with a SINGLE consistent convective
non-dimensionalization and a built-in physics + sampling validator.

Confirmed from session_final.xml:
    TimeStep = 0.005 , IO_CheckSteps = 20  -> 0.1 s between checkpoints
    Kinvis (nu) = 0.01 = 1/Re  -> expected lambda2 = 0.01
    fields written as scalars u, v, p

U_INF is read from <BOUNDARYCONDITIONS> (inlet). If you set U_INF = None below,
the script MEASURES it from the inlet nodes and reports the implied Reynolds
number as a sanity check.
"""

import glob, os, re
import numpy as np
import scipy.io as sio
import meshio   # pip install meshio

# ----------------------------------------------------------------------
# 1. CONFIG
# ----------------------------------------------------------------------
HERE         = os.path.dirname(os.path.abspath(__file__))
VTU_GLOB     = os.path.join(HERE, "cyl_chk_*.vtu")
OUT_DIR      = HERE

D            = 1.0          # mesh is non-dimensional: U_inf=1, D=1, nu=0.01 -> Re=100
U_INF        = None         # inlet speed; None -> auto-measure from inlet nodes
RHO          = 1.0          # Nektar++ incompressible uses 1
RE           = 100.0        # expected lambda2 = 1/RE
NU           = 0.01         # Kinvis from session_final.xml (cross-check: nu==1/Re)
X_REF, Y_REF = 0.15, 0.16   # cylinder centre [mesh units]

TIMESTEP     = 0.005        # <PARAMETERS> TimeStep
IO_CHECKSTEPS= 20           # <PARAMETERS> IO_CheckSteps
DT_CKPT      = IO_CHECKSTEPS * TIMESTEP          # 0.1 s between checkpoint files

# which checkpoints to use (skip startup transient; don't read all 1975)
CKPT_START   = 300          # first checkpoint index to use (past transient)
CKPT_STRIDE  = 1            # take every Nth checkpoint
CKPT_COUNT   = 500          # how many snapshots to keep

# windows in NON-DIMENSIONAL diameters from the cylinder centre
WINDOWS = {
    "near": dict(xlim=(1.0,  8.0),  ylim=(-2.5, -0.5)),
    "mid" : dict(xlim=(8.0, 15.0),  ylim=(-2.0, 2.0)),
    "far" : dict(xlim=(15.0,25.0),  ylim=(-2.0, 2.0)),
}
TARGET_N  = None            # keep all raw mesh nodes in each window (no subsampling)
RNG_SEED  = 0

# ----------------------------------------------------------------------
# 2. READ THE CHECKPOINTS
# ----------------------------------------------------------------------
def _numkey(fn):
    m = re.findall(r"\d+", os.path.basename(fn))
    return int(m[-1]) if m else 0

def _get(pd, *names):
    for n in names:
        if n in pd:
            return np.asarray(pd[n]).reshape(-1)
    raise KeyError(f"none of {names} in point_data keys {list(pd.keys())}")

def read_checkpoints(vtu_glob):
    files = sorted(glob.glob(vtu_glob), key=_numkey)
    if not files:
        raise FileNotFoundError(f"no files match {vtu_glob}")
    files = files[CKPT_START : CKPT_START + CKPT_COUNT*CKPT_STRIDE : CKPT_STRIDE]
    print(f"  using {len(files)} checkpoints "
          f"({os.path.basename(files[0])} .. {os.path.basename(files[-1])}), "
          f"stride {CKPT_STRIDE}")
    m0 = meshio.read(files[0])
    xy = m0.points[:, :2].copy()
    N, T = xy.shape[0], len(files)
    U = np.empty((N, T)); V = np.empty((N, T)); P = np.empty((N, T))
    for k, f in enumerate(files):
        pd = meshio.read(f).point_data
        if "u" in pd and "v" in pd:
            U[:, k] = _get(pd, "u"); V[:, k] = _get(pd, "v")
        else:
            vel = np.asarray(pd.get("velocity", pd.get("U")))
            U[:, k] = vel[:, 0]; V[:, k] = vel[:, 1]
        P[:, k] = _get(pd, "p", "pressure")
        if (k+1) % 50 == 0:
            print(f"    read {k+1}/{T}")
    dt_phys = DT_CKPT * CKPT_STRIDE
    t_phys  = np.arange(T) * dt_phys
    return xy, t_phys, U, V, P, dt_phys

# ----------------------------------------------------------------------
# 3. MEASURE U_INF FROM THE INLET (if not given)
# ----------------------------------------------------------------------
def measure_uinf(xy, U):
    x = xy[:, 0]
    inlet = x < (x.min() + 0.02*(x.max()-x.min()))   # leftmost ~2% = inflow strip
    uinf  = float(np.median(U[inlet].mean(axis=1)))  # time-mean, then median over strip
    return uinf

# ----------------------------------------------------------------------
# 4. NON-DIMENSIONALIZE (one consistent convective scaling) + CROP
# ----------------------------------------------------------------------
def build_window(xy, t_phys, U, V, P, uinf, xlim, ylim):
    x_star = (xy[:, 0] - X_REF) / D
    y_star = (xy[:, 1] - Y_REF) / D
    u_star = U / uinf                  # full field, mean kept
    v_star = V / uinf
    p_star = P / (RHO * uinf**2)
    t_star = t_phys * uinf / D         # the scale that must not be skipped

    m = (x_star>=xlim[0])&(x_star<=xlim[1])&(y_star>=ylim[0])&(y_star<=ylim[1])
    idx = np.where(m)[0]
    if idx.size == 0:
        raise ValueError(f"window {xlim}x{ylim} empty; "
                         f"x_star[{x_star.min():.1f},{x_star.max():.1f}] "
                         f"y_star[{y_star.min():.1f},{y_star.max():.1f}]")
    raw_n = idx.size
    if TARGET_N and idx.size > TARGET_N:
        idx = np.sort(np.random.default_rng(RNG_SEED).choice(idx, TARGET_N, False))
        print(f"  raw mesh nodes in window: {raw_n}  -> subsampled to {TARGET_N}")
    elif TARGET_N:
        print(f"  raw mesh nodes in window: {raw_n}  -> BELOW TARGET_N={TARGET_N}, "
              f"keeping all {raw_n} (re-export hi-res VTUs to get more)")
    else:
        print(f"  raw mesh nodes in window: {raw_n}  -> keeping all (TARGET_N=None)")
    return dict(x_star=x_star[idx][:,None], y_star=y_star[idx][:,None],
                t_star=t_star[:,None],
                U_star=u_star[idx], V_star=v_star[idx], P_star=p_star[idx])

# ----------------------------------------------------------------------
# 5. VALIDATOR
# ----------------------------------------------------------------------
def validate(ds, dt_phys, uinf):
    x=ds["x_star"].ravel(); y=ds["y_star"].ravel(); t=ds["t_star"].ravel()
    U=ds["U_star"]; V=ds["V_star"]; dt=float(np.mean(np.diff(t))); ok=True

    umean=U.mean(); a=umean>0.3
    print(f"  [A] mean u = {umean:+.3f}  ({'PASS' if a else 'FAIL: mean-subtracted?'})"); ok&=a

    # short-baseline probe pairs (<~half a vortex wavelength) keep the
    # cross-correlation lag unambiguous in dense, multi-wavelength windows
    n0=np.where(np.abs(y)<0.3)[0]
    if n0.size>=2:
        n0=n0[np.argsort(x[n0])]; xs=x[n0]; Ucs=[]
        for i in range(len(n0)):
            for j in range(i+1,len(n0)):
                dx=xs[j]-xs[i]
                if dx<0.5 or dx>2.0: continue
                aa=V[n0[i]]-V[n0[i]].mean(); bb=V[n0[j]]-V[n0[j]].mean()
                lag=(np.argmax(np.correlate(bb,aa,"full"))-(len(aa)-1))*dt
                if lag>0: Ucs.append(dx/lag)
        if Ucs:
            Uc=float(np.median(Ucs)); r=Uc/max(abs(umean),1e-9); b=0.3<r<3.0
            print(f"  [B] convection {Uc:+.3f} vs mean u {umean:+.3f} "
                  f"(x{r:.1f}; {'PASS' if b else 'FAIL: time scale != space/velocity'})"); ok&=b

    # [D] temporal resolution: samples per shedding period (physical)
    vm=V.mean(0)-V.mean(); fr=np.fft.rfftfreq(len(vm),dt_phys); A=np.abs(np.fft.rfft(vm))
    if len(A)>2:
        fphys=fr[1:][np.argmax(A[1:])]; per=1/fphys if fphys else np.inf
        spp=per/dt_phys
        d=spp>=8
        print(f"  [D] shedding period {per:.3f}s, {spp:.1f} samples/period "
              f"({'PASS' if d else 'WARN: undersampled, increase IO output rate'})"); ok&=d

    re_imp=uinf*D/NU; lam2=NU/(uinf*D)
    print(f"  [C] U_inf={uinf:.3f} -> implied Re={re_imp:.1f} (target {RE:.0f}); "
          f"expected lambda2 = nu/(U*D) = {lam2:.4f}")
    print(f"  => {'READY TO TRAIN' if ok else 'CHECK FLAGS ABOVE'}")
    return ok

# ----------------------------------------------------------------------
# 6. MAIN
# ----------------------------------------------------------------------
if __name__ == "__main__":
    xy, t_phys, U, V, P, dt_phys = read_checkpoints(VTU_GLOB)
    uinf = U_INF if U_INF is not None else measure_uinf(xy, U)
    print(f"  U_inf = {uinf:.4f} ({'config' if U_INF else 'measured from inlet'})\n")
    for name in ("mid", "far"):
        ds = build_window(xy, t_phys, U, V, P, uinf, **WINDOWS[name])
        out = os.path.join(OUT_DIR, f"cylinder_nektar_{name}.mat")
        sio.savemat(out, ds)
        print(f"{name}: {ds['U_star'].shape[0]} pts x {ds['U_star'].shape[1]} snaps -> {out}")
        validate(ds, dt_phys, uinf); print()