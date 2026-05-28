"""Scalar fields for the five supported TPMS families.

All fields are written as f(x,y,z) where coordinates are in *lattice space*
(one period = 2π). The operator converts world coordinates to lattice
coordinates by dividing by (cell_size / 2π).

References (level-set form):
  Gyroid:           sin x cos y + sin y cos z + sin z cos x
  Schwarz P:        cos x + cos y + cos z
  Schwarz D:        sin x sin y sin z + sin x cos y cos z
                    + cos x sin y cos z + cos x cos y sin z
  Schoen IWP:       2(cos x cos y + cos y cos z + cos z cos x)
                    − (cos 2x + cos 2y + cos 2z)
  Fischer-Koch S:   cos 2x sin y cos z + cos x cos 2y sin z
                    + sin x cos y cos 2z
"""

import numpy as np

TWO_PI = 2.0 * np.pi


def _gyroid(x, y, z):
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def _schwarz_p(x, y, z):
    return np.cos(x) + np.cos(y) + np.cos(z)


def _schwarz_d(x, y, z):
    sx, cx = np.sin(x), np.cos(x)
    sy, cy = np.sin(y), np.cos(y)
    sz, cz = np.sin(z), np.cos(z)
    return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz


def _schoen_iwp(x, y, z):
    cx, cy, cz = np.cos(x), np.cos(y), np.cos(z)
    return (
        2.0 * (cx * cy + cy * cz + cz * cx)
        - (np.cos(2.0 * x) + np.cos(2.0 * y) + np.cos(2.0 * z))
    )


def _fischer_koch_s(x, y, z):
    return (
        np.cos(2.0 * x) * np.sin(y) * np.cos(z)
        + np.cos(x) * np.cos(2.0 * y) * np.sin(z)
        + np.sin(x) * np.cos(y) * np.cos(2.0 * z)
    )


FIELDS = {
    'GYROID':         _gyroid,
    'SCHWARZ_P':      _schwarz_p,
    'SCHWARZ_D':      _schwarz_d,
    'SCHOEN_IWP':     _schoen_iwp,
    'FISCHER_KOCH_S': _fischer_koch_s,
}


def _euler_to_matrix(rx, ry, rz):
    """ZYX intrinsic Euler -> rotation matrix (numpy 3x3)."""
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    return np.array([
        [cy * cz,                 -cy * sz,               sy],
        [sx * sy * cz + cx * sz, -sx * sy * sz + cx * cz, -sx * cy],
        [-cx * sy * cz + sx * sz, cx * sy * sz + sx * cz,  cx * cy],
    ])


def sample_field(tpms_type, X, Y, Z, cell_size, origin=(0.0, 0.0, 0.0),
                 euler=(0.0, 0.0, 0.0)):
    """Evaluate the TPMS field on a numpy meshgrid.

    Parameters
    ----------
    tpms_type : str        — key into FIELDS.
    X, Y, Z   : ndarray    — world-space coordinate grids (same shape).
    cell_size : float      — world-unit length of one lattice period.
    origin    : (3,) float — world-space phase shift.
    euler     : (3,) float — lattice rotation (radians).
    """
    if tpms_type not in FIELDS:
        raise ValueError(f"Unknown TPMS type: {tpms_type!r}")

    ox, oy, oz = origin
    Xs = X - ox
    Ys = Y - oy
    Zs = Z - oz

    if any(abs(a) > 1e-9 for a in euler):
        R = _euler_to_matrix(*euler)
        # Transform world -> lattice frame by R^T
        flat = np.stack([Xs.ravel(), Ys.ravel(), Zs.ravel()], axis=0)
        rot = R.T @ flat
        Xs = rot[0].reshape(X.shape)
        Ys = rot[1].reshape(X.shape)
        Zs = rot[2].reshape(X.shape)

    k = TWO_PI / cell_size
    return FIELDS[tpms_type](Xs * k, Ys * k, Zs * k)
