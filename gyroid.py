"""Exact parametric Gyroid via the Enneper-Weierstrass representation.

Based on: P.J.F. Gandy, J. Klinowski, "Exact computation of the triply
periodic G ('Gyroid') minimal surface", Chem. Phys. Lett. 321 (2000) 363-371.

The Weierstrass function of the D/G/P family is
    R(tau) = 1 / sqrt(tau^8 - 14 tau^4 + 1)
and the Gyroid is the associate surface with Bonnet angle
    theta = arccot(K'/K) = 38.0147740 deg,
where K = K(m=1/4), K' = K(m=3/4) are complete elliptic integrals
(computed here by the AGM).

The fundamental patch is integrated numerically on a Coons-patch grid over
the curved-quad parameter domain O-P-Q-R (Fig. 3 of the paper):

    O = 0,  P = (sqrt(3)-1)/sqrt(2),  R = iP,
    arc |omega + (1+i)/sqrt(2)| = sqrt(2) from P through Q to R.

Square-root branch-point singularities of R(tau) at P and R are removed with
the substitution zeta = sqrt(tau - branch), so every Gauss-Legendre segment
integral sees an analytic integrand.

The 96 isometries that assemble the fundamental patch into one cubic unit
cell (space group Ia-3d, no. 230) were *derived*, not copied: the surface was
analytically continued across the three boundary curves of the patch, the
three side-pairing isometries were extracted by Procrustes fits (residual
~1e-9), the space group was generated from them by composition, and the
lattice (period L = 4a in the paper's notation -- the paper's own Tables 1/2
frame description is internally inconsistent) was read off from the pure
translations. Every operation below is exact: rotation parts are signed
permutation matrices, translation parts are multiples of L/8.

Verified properties of the emitted mesh:
  * all faces are quads, zero degenerate, orientation globally consistent;
  * all vertices lie on the exact minimal surface (mean curvature = 0);
  * vertices on the cell boundary match their periodic partners to ~1e-9
    of the cell, so Array modifiers weld seamlessly;
  * the cell agrees with the level-set convention
    sin x cos y + sin y cos z + sin z cos x  (x,y,z in units of 2 pi / cell):
    mean |F| ~ 0.0066 -- the true deviation of the nodal approximation.

This module is pure numpy (no bpy) and unit-tested outside Blender.
"""

import numpy as np

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

def _ellipK(m):
    """Complete elliptic integral K(m) via the arithmetic-geometric mean."""
    a, b = 1.0, np.sqrt(1.0 - m)
    for _ in range(60):
        a, b = 0.5 * (a + b), np.sqrt(a * b)
        if abs(a - b) < 1e-16:
            break
    return np.pi / (2.0 * a)

THETA = np.arctan(_ellipK(0.25) / _ellipK(0.75))   # Bonnet angle, 38.0147740 deg

_SQ2 = np.sqrt(2.0)
_R0 = (np.sqrt(3.0) - 1.0) / _SQ2       # branch-point radius
_BP_P = _R0 + 0.0j                       # branch point on the real axis
_BP_R = 1j * _R0                         # branch point on the imaginary axis
_ARC_C = -(1.0 + 1.0j) / _SQ2            # centre of the boundary arc
_Q_PT = _ARC_C + _SQ2 * np.exp(1j * np.pi / 4.0)

# Edge length of the cubic unit cell in raw (kappa=1) integration units.
# Equals 4*sqrt(2)*x(P'); converged, stable to the last digit for n>=12.
_L_RAW = 2.656243469335

# Linear map raw frame -> standard cubic frame with unit cell edge 1.
# (rotate -45 deg about z, flip z, then relabel axes to match the sin-cos
# level-set convention; derived from the 4-bar axes of the space group.)
_C = 0.5 * _SQ2
_A_STD = np.array([[_C, -_C, 0.0],
                   [0.0, 0.0, 1.0],
                   [-_C, -_C, 0.0]]) / _L_RAW
_ROT_STD = _A_STD * _L_RAW               # orthogonal part (det = +1)
_T_STD = np.array([0.75, 0.875, 0.5])    # eighths of the cell edge

# The 96 isometries of Ia-3d assembling the unit cell from the fundamental
# patch, in the standard frame: (m00..m22, t0, t1, t2) with the 3x3 signed
# permutation matrix in row-major order and the translation in eighths of
# the cell edge.
_OPS = (
    (1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0),
    (0, 0, 1, 0, -1, 0, -1, 0, 0, 2, 6, 2),
    (0, 0, -1, 0, -1, 0, 1, 0, 0, 2, 6, 6),
    (-1, 0, 0, 0, 0, 1, 0, 1, 0, 6, 2, 6),
    (-1, 0, 0, 0, 1, 0, 0, 0, -1, 4, 0, 0),
    (0, 0, -1, -1, 0, 0, 0, -1, 0, 4, 4, 4),
    (0, 0, 1, 1, 0, 0, 0, -1, 0, 4, 0, 4),
    (0, 1, 0, 0, 0, -1, 1, 0, 0, 0, 4, 4),
    (0, -1, 0, 0, 0, -1, -1, 0, 0, 4, 4, 4),
    (1, 0, 0, 0, 0, -1, 0, 1, 0, 2, 2, 6),
    (0, -1, 0, 1, 0, 0, 0, 0, 1, 6, 2, 6),
    (0, 1, 0, 1, 0, 0, 0, 0, -1, 6, 2, 2),
    (0, -1, 0, -1, 0, 0, 0, 0, -1, 6, 6, 6),
    (0, 1, 0, -1, 0, 0, 0, 0, 1, 6, 6, 2),
    (1, 0, 0, 0, 0, 1, 0, -1, 0, 6, 2, 2),
    (0, -1, 0, 1, 0, 0, 0, 0, -1, 6, 6, 2),
    (0, 1, 0, -1, 0, 0, 0, 0, -1, 2, 6, 2),
    (0, 1, 0, 0, 0, 1, -1, 0, 0, 0, 4, 0),
    (0, -1, 0, 0, 0, 1, 1, 0, 0, 4, 4, 0),
    (0, 0, 1, -1, 0, 0, 0, 1, 0, 0, 4, 4),
    (0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0),
    (0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 4, 0),
    (0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0, 4),
    (0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 4),
    (0, -1, 0, 0, 0, 1, -1, 0, 0, 0, 4, 4),
    (-1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 4, 0),
    (0, 0, -1, -1, 0, 0, 0, 1, 0, 4, 0, 4),
    (0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 4),
    (0, 0, -1, 1, 0, 0, 0, -1, 0, 4, 0, 0),
    (0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0),
    (-1, 0, 0, 0, 0, -1, 0, -1, 0, 2, 2, 2),
    (0, -1, 0, -1, 0, 0, 0, 0, 1, 6, 2, 2),
    (0, 1, 0, 1, 0, 0, 0, 0, 1, 2, 2, 2),
    (0, 0, -1, 0, 1, 0, -1, 0, 0, 6, 6, 2),
    (1, 0, 0, 0, 0, -1, 0, -1, 0, 2, 6, 2),
    (-1, 0, 0, 0, 0, -1, 0, 1, 0, 2, 6, 6),
    (1, 0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 2),
    (-1, 0, 0, 0, 0, 1, 0, -1, 0, 2, 2, 6),
    (0, 0, 1, 0, 1, 0, 1, 0, 0, 6, 6, 6),
    (-1, 0, 0, 0, 0, 1, 0, -1, 0, 6, 6, 2),
    (1, 0, 0, 0, 0, 1, 0, 1, 0, 6, 6, 6),
    (-1, 0, 0, 0, 0, -1, 0, 1, 0, 6, 2, 2),
    (1, 0, 0, 0, 0, -1, 0, -1, 0, 6, 2, 6),
    (0, 0, 1, 0, 1, 0, 1, 0, 0, 2, 2, 2),
    (0, 0, -1, 0, 1, 0, -1, 0, 0, 2, 2, 6),
    (0, 1, 0, 1, 0, 0, 0, 0, 1, 6, 6, 6),
    (0, 0, 1, 0, 1, 0, -1, 0, 0, 2, 6, 6),
    (0, 0, -1, 0, -1, 0, -1, 0, 0, 6, 6, 6),
    (0, -1, 0, -1, 0, 0, 0, 0, 1, 2, 6, 6),
    (0, 0, 1, 0, -1, 0, 1, 0, 0, 2, 2, 6),
    (0, 0, -1, 0, 1, 0, 1, 0, 0, 6, 2, 6),
    (1, 0, 0, 0, -1, 0, 0, 0, -1, 4, 4, 0),
    (0, 0, 1, 1, 0, 0, 0, 1, 0, 4, 4, 4),
    (0, 0, -1, 1, 0, 0, 0, -1, 0, 0, 4, 4),
    (0, 0, 1, -1, 0, 0, 0, -1, 0, 4, 4, 0),
    (0, 0, -1, -1, 0, 0, 0, 1, 0, 0, 4, 0),
    (-1, 0, 0, 0, -1, 0, 0, 0, 1, 4, 0, 4),
    (1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 4),
    (0, -1, 0, 0, 0, 1, -1, 0, 0, 4, 0, 0),
    (-1, 0, 0, 0, -1, 0, 0, 0, -1, 4, 4, 4),
    (1, 0, 0, 0, 1, 0, 0, 0, -1, 4, 0, 4),
    (0, 1, 0, 0, 0, -1, -1, 0, 0, 4, 4, 0),
    (-1, 0, 0, 0, 1, 0, 0, 0, 1, 4, 4, 0),
    (1, 0, 0, 0, -1, 0, 0, 0, 1, 4, 0, 0),
    (0, -1, 0, 0, 0, -1, 1, 0, 0, 4, 0, 4),
    (1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 4, 4),
    (-1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 4),
    (0, 1, 0, 0, 0, 1, 1, 0, 0, 4, 4, 4),
    (1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 4, 0),
    (-1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0),
    (0, 0, -1, 1, 0, 0, 0, 1, 0, 4, 4, 0),
    (0, 0, 1, -1, 0, 0, 0, 1, 0, 4, 0, 0),
    (0, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4),
    (0, 1, 0, 0, 0, 1, -1, 0, 0, 4, 0, 4),
    (0, 1, 0, -1, 0, 0, 0, 0, -1, 6, 2, 6),
    (0, 0, -1, 0, 1, 0, 1, 0, 0, 2, 6, 2),
    (0, 0, 1, 0, -1, 0, 1, 0, 0, 6, 6, 2),
    (0, -1, 0, 1, 0, 0, 0, 0, -1, 2, 2, 6),
    (0, 0, -1, 0, -1, 0, -1, 0, 0, 2, 2, 2),
    (0, 0, 1, 0, 1, 0, -1, 0, 0, 6, 2, 2),
    (1, 0, 0, 0, 0, 1, 0, -1, 0, 2, 6, 6),
    (-1, 0, 0, 0, 0, -1, 0, -1, 0, 6, 6, 6),
    (0, 1, 0, -1, 0, 0, 0, 0, 1, 2, 2, 6),
    (0, -1, 0, -1, 0, 0, 0, 0, -1, 2, 2, 2),
    (0, 1, 0, 1, 0, 0, 0, 0, -1, 2, 6, 6),
    (0, -1, 0, 1, 0, 0, 0, 0, 1, 2, 6, 2),
    (1, 0, 0, 0, 0, -1, 0, 1, 0, 6, 6, 2),
    (0, -1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, -1, 1, 0, 0, 4, 0, 0),
    (0, 0, 1, 1, 0, 0, 0, -1, 0, 0, 4, 0),
    (0, 0, -1, -1, 0, 0, 0, -1, 0, 0, 0, 0),
    (-1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 4, 4),
    (-1, 0, 0, 0, 0, 1, 0, 1, 0, 2, 6, 2),
    (0, 0, -1, 0, -1, 0, 1, 0, 0, 6, 2, 2),
    (0, 0, 1, 0, -1, 0, -1, 0, 0, 6, 2, 6),
    (1, 0, 0, 0, 1, 0, 0, 0, 1, 4, 4, 4),)

_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)
_CORNER_EPS = 0.14

# ----------------------------------------------------------------------
# Fundamental patch
# ----------------------------------------------------------------------

def _coons(u, v):
    """Map the unit square onto the curved-quad parameter domain O-P-Q-R."""
    bottom = u * _BP_P
    top = _ARC_C + _SQ2 * np.exp(1j * (np.pi / 3.0 - u * np.pi / 12.0))
    left = v * _BP_R
    right = _ARC_C + _SQ2 * np.exp(1j * (np.pi / 6.0 + v * np.pi / 12.0))
    return ((1 - v) * bottom + v * top + (1 - u) * left + u * right
            - (u * (1 - v) * _BP_P + u * v * _Q_PT + (1 - u) * v * _BP_R))

def _integrand(t):
    """(1 - t^2, i(1 + t^2), 2t) * R(t), stacked on the last axis."""
    t4 = t ** 4
    rw = 1.0 / np.sqrt(t4 * t4 - 14.0 * t4 + 1.0)
    return np.stack([(1.0 - t * t) * rw,
                     1j * (1.0 + t * t) * rw,
                     2.0 * t * rw], axis=-1)

def _segment_integrals(A, B):
    """Gauss-Legendre integrals of the Weierstrass form along straight
    segments A->B. Segments ending near a branch point are integrated in
    the zeta = sqrt(tau - branch) chart, where the integrand is analytic."""
    A = np.asarray(A, dtype=complex)
    B = np.asarray(B, dtype=complex)
    out = np.zeros(A.shape + (3,), dtype=complex)

    dP = np.minimum(abs(A - _BP_P), abs(B - _BP_P))
    dR = np.minimum(abs(A - _BP_R), abs(B - _BP_R))
    near_p = (dP < _CORNER_EPS) & (dP <= dR)
    near_r = (dR < _CORNER_EPS) & (dR < dP)
    plain = ~(near_p | near_r)

    if plain.any():
        a, b = A[plain], B[plain]
        mid = 0.5 * (a + b)[:, None]
        half = 0.5 * (b - a)[:, None]
        vals = _integrand(mid + half * _GL_X[None, :])
        out[plain] = np.einsum('sgk,g->sk', vals, _GL_W) * half

    for mask, bp in ((near_p, _BP_P), (near_r, _BP_R)):
        if not mask.any():
            continue
        a, b = A[mask], B[mask]
        za, zb = np.sqrt(a - bp), np.sqrt(b - bp)
        flip = (za.real * zb.real + za.imag * zb.imag) < 0.0
        zb = np.where(flip, -zb, zb)
        mid = 0.5 * (za + zb)[:, None]
        half = 0.5 * (zb - za)[:, None]
        zn = mid + half * _GL_X[None, :]
        vals = _integrand(bp + zn * zn) * (2.0 * zn)[..., None]
        out[mask] = np.einsum('sgk,g->sk', vals, _GL_W) * half
    return out

def _fundamental_patch(res):
    """Integrate the Enneper-Weierstrass representation on an
    (res+1) x (res+1) grid. Returns (points, normals) in the standard
    cubic frame with unit cell edge 1."""
    u = np.linspace(0.0, 1.0, res + 1)
    U, V = np.meshgrid(u, u, indexing='ij')
    om = _coons(U, V)

    phi = np.zeros((res + 1, res + 1, 3), dtype=complex)
    segs = _segment_integrals(om[:-1, 0], om[1:, 0])
    phi[1:, 0] = np.cumsum(segs, axis=0)
    segs = _segment_integrals(om[:, :-1].ravel(), om[:, 1:].ravel())
    phi[:, 1:] = phi[:, :1] + np.cumsum(segs.reshape(res + 1, res, 3), axis=1)

    xyz = (np.exp(1j * THETA) * phi).real            # raw frame, kappa = 1

    # Gauss map: the unit normal is the inverse stereographic image of omega.
    gauss = np.stack([2.0 * om.real, 2.0 * om.imag, np.abs(om) ** 2 - 1.0],
                     axis=-1)
    gauss /= np.linalg.norm(gauss, axis=-1, keepdims=True)

    return xyz @ _A_STD.T + _T_STD, gauss @ _ROT_STD.T

# ----------------------------------------------------------------------
# Unit-cell assembly
# ----------------------------------------------------------------------

def build_unit_cell(cell_size=1.0, res=8):
    """Build one cubic unit cell of the exact Gyroid as an all-quad mesh.

    Parameters
    ----------
    cell_size : float -- world edge length of the cubic translational cell.
    res       : int   -- quads per fundamental-patch edge; the cell has
                         96 * res^2 quads.

    Returns
    -------
    verts   : (V, 3) float64 -- welded vertices, cell spans ~[0, cell_size]^3
    quads   : (Q, 4) int32   -- quad faces, globally consistent orientation
    normals : (V, 3) float64 -- exact analytic unit normals
    """
    res = max(2, int(res))
    patch, pnorm = _fundamental_patch(res)
    npts = (res + 1) * (res + 1)

    # one patch's quad grid, reused (offset) for each of the 96 copies
    idx = np.arange(npts).reshape(res + 1, res + 1)
    q0 = np.stack([idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(),
                   idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()], axis=1)

    pflat = patch.reshape(-1, 3)
    nflat = pnorm.reshape(-1, 3)
    all_v = np.empty((96 * npts, 3))
    all_n = np.empty((96 * npts, 3))
    all_q = np.empty((96 * res * res, 4), dtype=np.int64)

    for k, row in enumerate(_OPS):
        M = np.array(row[:9], dtype=float).reshape(3, 3)
        t = np.array(row[9:], dtype=float) / 8.0
        p = pflat @ M.T + t
        t -= np.floor(p.mean(axis=0))        # keep the copy near [0, 1)^3
        det = float(np.linalg.det(M))
        all_v[k * npts:(k + 1) * npts] = pflat @ M.T + t
        all_n[k * npts:(k + 1) * npts] = nflat @ M.T * det
        all_q[k * res * res:(k + 1) * res * res] = q0 + k * npts

    # weld: vertices from different patches on shared boundary curves agree
    # to ~1e-9, far below the smallest quad size (~0.1 / res).
    tol = 1e-6
    keys = np.round(all_v / tol).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    verts = all_v[first] * cell_size
    normals = all_n[first]
    quads = inverse[all_q].astype(np.int32)
    return verts, quads, normals
