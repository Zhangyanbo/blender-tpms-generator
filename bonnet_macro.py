"""Genuine macro-quad meshes for the Schwarz P and Schwarz D associates."""

from collections import defaultdict

import numpy as np

try:
    from . import gyroid_macro as core
except ImportError:
    import gyroid_macro as core


_PAIR_CACHE = {}


def _surface_points(omega, surface):
    """Evaluate one exact fundamental triangle in standard cell coordinates."""
    shape = np.shape(omega)
    om = np.asarray(omega, dtype=complex).reshape(-1)
    t, weights = core._GL_T, core._GL_WEIGHTS
    tau = om[:, None] * t[None, :]
    polynomial = tau ** 8 - 14.0 * tau ** 4 + 1.0
    phase = np.unwrap(np.angle(polynomial), axis=1)
    root = np.sqrt(np.abs(polynomial)) * np.exp(0.5j * phase)
    common = weights[None, :] / root

    i0 = om * np.sum(common, axis=1)
    i1 = om ** 2 * np.sum(common * t[None, :], axis=1)
    i2 = om ** 3 * np.sum(common * t[None, :] ** 2, axis=1)
    rotation = np.exp(1j * surface.theta)
    raw = np.stack(((rotation * (i0 - i2)).real,
                    (rotation * 1j * (i0 + i2)).real,
                    (rotation * 2.0 * i1).real), axis=1)
    standard = raw @ surface.frame.T / surface.L_raw + surface.t_std
    return standard.reshape(shape + (3,))


def _operation_transforms(surface, reference_patch):
    """Resolve the periodic representative chosen by each space-group row."""
    reference = np.asarray(reference_patch).reshape(-1, 3)
    transforms = []
    for row in surface.ops:
        matrix = np.asarray(row[:9], dtype=float).reshape(3, 3)
        translation = np.asarray(row[9:12], dtype=float) / 16.0
        transformed = reference @ matrix.T + translation
        translation -= np.floor(transformed.mean(axis=0))
        transforms.append((matrix, translation))
    return tuple(transforms)


def _apply(points, transform):
    matrix, translation = transform
    return np.asarray(points) @ matrix.T + translation


def _macro_pairs(surface_name, surface, transforms):
    cached = _PAIR_CACHE.get(surface_name)
    if cached is not None:
        return cached

    phi = np.linspace(0.0, 0.5 * np.pi, 17)
    omega = core._rho_max(phi) * np.exp(1j * phi)
    arc = _surface_points(omega, surface)
    groups = defaultdict(list)
    curves = []
    for index, transform in enumerate(transforms):
        curve = _apply(arc, transform)
        curves.append(curve)
        groups[core._periodic_curve_key(curve)].append(index)

    expected = len(transforms) // 2
    if len(groups) != expected or any(len(group) != 2 for group in groups.values()):
        sizes = sorted(len(group) for group in groups.values())
        raise RuntimeError(
            f"{surface_name} circular-edge pairing failed: expected "
            f"{expected} groups of 2, got {len(groups)} with sizes {sizes}"
        )

    pairs = []
    for group in groups.values():
        lower, upper = sorted(group)
        lower_curve = curves[lower]
        upper_curve = curves[upper]
        periodic_error = np.mod(
            lower_curve - upper_curve[::-1] + 0.5, 1.0) - 0.5
        if np.max(np.linalg.norm(periodic_error, axis=1)) > 2.0e-5:
            raise RuntimeError(f"{surface_name} paired arcs do not run oppositely")
        shift = np.rint(np.mean(lower_curve - upper_curve[::-1], axis=0))
        seam_error = np.max(np.linalg.norm(
            lower_curve - (upper_curve[::-1] + shift), axis=1))
        if seam_error > 1.0e-9:
            raise RuntimeError(f"{surface_name} hidden-edge mismatch: {seam_error:.3g}")
        pairs.append((lower, upper, shift))

    pairs.sort(key=lambda item: item[:2])
    _PAIR_CACHE[surface_name] = tuple(pairs)
    return _PAIR_CACHE[surface_name]


def _macro_patch_points(u, v, pair, surface, transforms):
    """Evaluate one globally smooth analytic P/D macro patch."""
    omega, upper_mask = core._macro_square_omega(u, v)
    base = _surface_points(omega, surface).reshape(-1, 3)
    upper_mask = upper_mask.ravel()
    lower, upper, shift = pair
    points = np.empty_like(base)
    points[~upper_mask] = _apply(base[~upper_mask], transforms[lower])
    points[upper_mask] = _apply(base[upper_mask], transforms[upper]) + shift
    return points.reshape(np.shape(omega) + (3,))


def build_unit_cell(surface_name, surface, reference_patch, cell_size=1.0,
                    quad_subdivisions=2):
    """Build Schwarz P/D directly from genuine four-sided parameter patches."""
    cell_size = float(cell_size)
    if not np.isfinite(cell_size) or cell_size <= 0.0:
        raise ValueError("cell_size must be a positive finite number")
    n = max(1, int(quad_subdivisions))
    transforms = _operation_transforms(surface, reference_patch)
    pairs = _macro_pairs(surface_name, surface, transforms)

    axis = np.linspace(0.0, 1.0, n + 1)
    u, v = np.meshgrid(axis, axis, indexing='ij')
    index = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    local_faces = np.stack((index[:-1, :-1].ravel(),
                            index[1:, :-1].ravel(),
                            index[1:, 1:].ravel(),
                            index[:-1, 1:].ravel()), axis=1)

    vertices = []
    faces = []
    offset = 0
    for pair in pairs:
        patch = _macro_patch_points(
            u, v, pair, surface, transforms).reshape(-1, 3)
        vertices.append(patch)
        faces.append(local_faces + offset)
        offset += len(patch)
    vertices, faces = core._weld(np.vstack(vertices), np.vstack(faces))
    faces = core._orient_faces(vertices, faces)

    expected_faces = len(pairs) * n * n
    if len(faces) != expected_faces:
        raise RuntimeError(
            f"{surface_name} expected {expected_faces} macro quads, got {len(faces)}")
    metrics = core.validate_unit_cell(vertices, faces)
    expected_euler = {'SCHWARZ_P': -4, 'SCHWARZ_D': -16}[surface_name]
    required = (
        metrics['non_quad_faces'] == 0,
        metrics['duplicate_faces'] == 0,
        metrics['degenerate_quads'] == 0,
        metrics['folded_quads'] == 0,
        metrics['nonmanifold_edges'] == 0,
        metrics['all_periodic_boundary_edges_paired'],
        metrics['quotient_all_edges_incident_to_two_quads'],
        metrics['quotient_euler_characteristic'] == expected_euler,
    )
    if not all(required):
        raise RuntimeError(f"{surface_name} topology acceptance failed: {metrics}")
    normals = core._vertex_normals(vertices, faces)
    return vertices * cell_size, faces, normals


def macro_patch_count(surface_name, surface, reference_patch):
    transforms = _operation_transforms(surface, reference_patch)
    return len(_macro_pairs(surface_name, surface, transforms))
