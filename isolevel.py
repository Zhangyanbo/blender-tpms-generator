"""Level-set reshaping of the exact TPMS quad meshes.

The exact Enneper-Weierstrass meshes sample the minimal surfaces, which hug
the standard trigonometric TPMS fields F at level ~0 (coordinates in units
of one cell, angles in units of 2 pi / cell):

    Gyroid:    sin x cos y + sin y cos z + sin z cos x
    Schwarz P: cos x + cos y + cos z
    Schwarz D: sin x sin y sin z + sin x cos y cos z
               + cos x sin y cos z + cos x cos y sin z

A level surface F = t with t != 0 is not a minimal surface, so no
Weierstrass parametrization of it exists. Instead every vertex of the exact
quad mesh is transported along the field gradient onto F = t by a damped
Newton iteration. The transport map is smooth and commutes with the cell
lattice (F is 1-periodic in each coordinate), so mesh connectivity, the
genuine macro-quad patch structure, periodic seam welding and the all-quad
guarantee are all untouched; only vertex positions and normals change.

Valid levels are bounded by the first critical value of each field, where
the level set changes topology: |t| < sqrt(2) for the Gyroid, |t| < 1 for
Schwarz P and Schwarz D. Slightly before those limits the necks pinch so
hard that quads spanning them bend past the folded-quad acceptance check
(raising subdivisions does not help -- the bending is intrinsic), so the
practical ranges are |t| <~ 1.30 (Gyroid) and |t| <~ 0.99 (P/D); beyond
them generation fails with a clear error instead of emitting a dirty mesh.

This module is pure numpy (no bpy) and unit-tested outside Blender.
"""

import numpy as np

_TAU = 2.0 * np.pi

# First critical value of each field: level-set topology changes there.
LEVEL_LIMITS = {
    'GYROID': np.sqrt(2.0),
    'SCHWARZ_P': 1.0,
    'SCHWARZ_D': 1.0,
}

_LEVEL_MARGIN = 1.0e-3


def field(tpms_type, points):
    """Evaluate the standard trigonometric field in unit-cell coordinates."""
    a = _TAU * np.asarray(points, dtype=float)
    x, y, z = a[..., 0], a[..., 1], a[..., 2]
    if tpms_type == 'GYROID':
        return (np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z)
                + np.sin(z) * np.cos(x))
    if tpms_type == 'SCHWARZ_P':
        return np.cos(x) + np.cos(y) + np.cos(z)
    if tpms_type == 'SCHWARZ_D':
        sx, sy, sz = np.sin(x), np.sin(y), np.sin(z)
        cx, cy, cz = np.cos(x), np.cos(y), np.cos(z)
        return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz
    raise KeyError(f"Unknown TPMS type: {tpms_type}")


def gradient(tpms_type, points):
    """Gradient of `field` with respect to unit-cell coordinates."""
    a = _TAU * np.asarray(points, dtype=float)
    x, y, z = a[..., 0], a[..., 1], a[..., 2]
    sx, sy, sz = np.sin(x), np.sin(y), np.sin(z)
    cx, cy, cz = np.cos(x), np.cos(y), np.cos(z)
    if tpms_type == 'GYROID':
        g = (cx * cy - sz * sx, cy * cz - sx * sy, cz * cx - sy * sz)
    elif tpms_type == 'SCHWARZ_P':
        g = (-sx, -sy, -sz)
    elif tpms_type == 'SCHWARZ_D':
        g = (cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
             sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
             sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz)
    else:
        raise KeyError(f"Unknown TPMS type: {tpms_type}")
    return _TAU * np.stack(g, axis=-1)


def check_level(tpms_type, level):
    """Validate an isolevel against the field's first critical value."""
    limit = LEVEL_LIMITS[tpms_type]
    level = float(level)
    if not np.isfinite(level) or abs(level) > limit - _LEVEL_MARGIN:
        raise ValueError(
            f"{tpms_type} isolevel must satisfy |t| < {limit:.4g}; the level "
            "set pinches and changes topology at that critical value"
        )
    return level


def project_to_level(tpms_type, points, level, tol=1.0e-10, max_iter=200,
                     max_step=0.05):
    """Transport points along the field gradient onto F = level."""
    level = check_level(tpms_type, level)
    x = np.array(points, dtype=float)
    for _ in range(max_iter):
        residual = field(tpms_type, x) - level
        if np.max(np.abs(residual)) < tol:
            return x
        g = gradient(tpms_type, x)
        squared = np.einsum('ij,ij->i', g, g)
        step = -(residual / np.maximum(squared, 1.0e-12))[:, None] * g
        length = np.linalg.norm(step, axis=1)
        step *= np.minimum(1.0, max_step / np.maximum(length, 1.0e-30))[:, None]
        x += step
    raise RuntimeError(
        f"{tpms_type} isolevel projection did not converge at t = {level}; "
        "reduce |Iso Level|"
    )


def _unit_field_normals(tpms_type, points, reference_normals):
    """Normalized field gradient, sign-aligned with the mesh orientation."""
    g = gradient(tpms_type, points)
    g /= np.maximum(np.linalg.norm(g, axis=1, keepdims=True), 1.0e-30)
    if np.einsum('ij,ij->', g, reference_normals) < 0.0:
        g = -g
    return g


def _check_quads(tpms_type, verts, quads, context):
    p = verts[quads]
    n1 = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    n2 = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    if (np.min(np.linalg.norm(n1, axis=1)) <= 1.0e-12 or
            np.min(np.linalg.norm(n2, axis=1)) <= 1.0e-12 or
            np.min(np.einsum('ij,ij->i', n1, n2)) <= 0.0):
        raise ValueError(
            f"{tpms_type} {context} produced a degenerate or folded quad "
            "near a pinching neck; reduce |Iso Level|"
        )


def reshape_unit_cell(tpms_type, verts, quads, normals, iso_level=0.0):
    """Apply the isolevel parameter to an exact unit-cell quad mesh.

    All geometry is in unit-cell coordinates (cell edge 1). `iso_level = 0`
    keeps the exact minimal surface untouched; a nonzero level transports
    the mesh onto the standard trigonometric level set.
    """
    verts = np.asarray(verts, dtype=float)
    quads = np.asarray(quads)
    normals = np.asarray(normals, dtype=float)

    if iso_level != 0.0:
        verts = project_to_level(tpms_type, verts, iso_level)
        normals = _unit_field_normals(tpms_type, verts, normals)
        _check_quads(tpms_type, verts, quads, f"isolevel t = {iso_level:g}")

    return verts, quads, normals
