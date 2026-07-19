"""Exact Gyroid unit cell made from 48 genuine quadrilateral patches.

Two fundamental triangles are analytically unfolded across their shared circle
into one genuine four-sided complex domain. A Coons map addresses that domain
directly; no runtime PDE, inverse map, lookup table, or fitted surrogate exists.
"""

from collections import defaultdict, deque

import numpy as np


_SQ2 = np.sqrt(2.0)
_ARC_CENTER = -(1.0 + 1.0j) / _SQ2
_ARC_RADIUS = _SQ2
_BRANCH_RADIUS = (np.sqrt(3.0) - 1.0) / _SQ2
_REFLECTED_ORIGIN = (1.0 + 1.0j) / _SQ2
_PAIR_CACHE = None


def _ellipk(m):
    a, b = 1.0, np.sqrt(1.0 - m)
    for _ in range(60):
        a, b = 0.5 * (a + b), np.sqrt(a * b)
        if abs(a - b) < 1.0e-16:
            break
    return np.pi / (2.0 * a)


_K = _ellipk(0.25)
_KP = _ellipk(0.75)
_THETA = np.arctan(_K / _KP)
_KAPPA = np.hypot(_K, _KP) / (_K * _KP)

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(32)
_GL_S = 0.5 * (_GL_NODES + 1.0)
_GL_T = 2.0 * _GL_S - _GL_S * _GL_S
_GL_WEIGHTS = 0.5 * _GL_WEIGHTS * 2.0 * (1.0 - _GL_S)


_PATCH_TRANSFORMS = (
    (([-1/_SQ2, -1/_SQ2, 0], [-1/_SQ2, 1/_SQ2, 0], [0, 0, 1]), [0, .5, .75]),
    (([ 1/_SQ2, -1/_SQ2, 0], [-1/_SQ2,-1/_SQ2, 0], [0, 0,-1]), [0, .5, .75]),
    (([-1/_SQ2,  1/_SQ2, 0], [0, 0, 1], [ 1/_SQ2, 1/_SQ2, 0]), [.5, .75, 1]),
    (([-1/_SQ2, -1/_SQ2, 0], [0, 0,-1], [-1/_SQ2, 1/_SQ2, 0]), [.5, .75, 1]),
    (((0, 0, 1), [ 1/_SQ2, 1/_SQ2, 0], [ 1/_SQ2,-1/_SQ2, 0]), [.75, 1, .5]),
    (((0, 0,-1), [-1/_SQ2, 1/_SQ2, 0], [ 1/_SQ2, 1/_SQ2, 0]), [.75, 1, .5]),
    (([ 1/_SQ2,  1/_SQ2, 0], [ 1/_SQ2,-1/_SQ2, 0], [0, 0,-1]), [1, .5, .25]),
    (([-1/_SQ2,  1/_SQ2, 0], [ 1/_SQ2, 1/_SQ2, 0], [0, 0, 1]), [1, .5, .25]),
    (([ 1/_SQ2, -1/_SQ2, 0], [0, 0,-1], [-1/_SQ2,-1/_SQ2, 0]), [.5, .25, 0]),
    (([ 1/_SQ2,  1/_SQ2, 0], [0, 0, 1], [ 1/_SQ2,-1/_SQ2, 0]), [.5, .25, 0]),
    (((0, 0,-1), [-1/_SQ2,-1/_SQ2, 0], [-1/_SQ2, 1/_SQ2, 0]), [.25, 0, .5]),
    (((0, 0, 1), [ 1/_SQ2,-1/_SQ2, 0], [-1/_SQ2,-1/_SQ2, 0]), [.25, 0, .5]),
)

_BLOCK_TRANSFORMS = (
    ((( 1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0)),
    (((-1, 0, 0), (0,-1, 0), (0, 0, 1)), (2, 1, 0)),
    (((-1, 0, 0), (0, 1, 0), (0, 0, 1)), (1, 1, 0)),
    ((( 1, 0, 0), (0,-1, 0), (0, 0, 1)), (1, 2, 0)),
    ((( 1, 0, 0), (0,-1, 0), (0, 0, 1)), (0, 1,-1)),
    (((-1, 0, 0), (0, 1, 0), (0, 0, 1)), (2, 0,-1)),
    (((-1, 0, 0), (0,-1, 0), (0, 0, 1)), (1, 2,-1)),
    ((( 1, 0, 0), (0, 1, 0), (0, 0, 1)), (1, 1,-1)),
)


def _triangle_transforms():
    """Return the 96 affine transforms in normalized lattice coordinates."""
    result = []
    for block_a, block_b in _BLOCK_TRANSFORMS:
        ba = np.asarray(block_a, dtype=float)
        bb = np.asarray(block_b, dtype=float)
        for patch_a, patch_b in _PATCH_TRANSFORMS:
            pa = np.asarray(patch_a, dtype=float)
            pb = np.asarray(patch_b, dtype=float)
            result.append((0.5 * (ba @ pa),
                           0.5 * (pb @ ba.T + bb) + (0.0, 0.0, 0.5)))
    return tuple(result)


_TRIANGLE_TRANSFORMS = _triangle_transforms()


def _rho_max(phi):
    cs = np.cos(phi) + np.sin(phi)
    return (-_SQ2 * cs + np.sqrt(2.0 * cs * cs + 4.0)) * 0.5


def _weierstrass(omega):
    """Evaluate the exact fundamental triangle along continuous radial roots."""
    shape = np.shape(omega)
    om = np.asarray(omega, dtype=complex).reshape(-1)
    t, weights = _GL_T, _GL_WEIGHTS
    tau = om[:, None] * t[None, :]
    z = tau ** 8 - 14.0 * tau ** 4 + 1.0
    phase = np.unwrap(np.angle(z), axis=1)
    root = np.sqrt(np.abs(z)) * np.exp(0.5j * phase)
    common = weights[None, :] / root
    i0 = _KAPPA * om * np.sum(common, axis=1)
    i1 = _KAPPA * om ** 2 * np.sum(common * t[None, :], axis=1)
    i2 = _KAPPA * om ** 3 * np.sum(common * t[None, :] ** 2, axis=1)
    rot = np.exp(1j * _THETA)
    xyz = np.stack([(rot * (i0 - i2)).real,
                    (rot * 1j * (i0 + i2)).real,
                    (rot * 2.0 * i1).real], axis=1)
    return xyz.reshape(shape + (3,))


def _apply_transform(points, transform):
    a, b = transform
    return np.asarray(points) @ a.T + b


def _periodic_curve_key(curve):
    c = np.mod(curve, 1.0)
    c[np.isclose(c, 1.0, atol=2.0e-5)] = 0.0
    forward = tuple(np.round(c, 5).ravel())
    backward = tuple(np.round(c[::-1], 5).ravel())
    return min(forward, backward)


def _macro_pairs():
    """Discover the 48 circular-edge pairs; pair numbers are never hard-coded."""
    global _PAIR_CACHE
    if _PAIR_CACHE is not None:
        return _PAIR_CACHE

    phi = np.linspace(0.0, 0.5 * np.pi, 17)
    arc = _weierstrass(_rho_max(phi) * np.exp(1j * phi))
    groups = defaultdict(list)
    curves = []
    for index, transform in enumerate(_TRIANGLE_TRANSFORMS):
        curve = _apply_transform(arc, transform)
        curves.append(curve)
        groups[_periodic_curve_key(curve)].append(index)
    if len(groups) != 48 or any(len(group) != 2 for group in groups.values()):
        sizes = sorted(len(group) for group in groups.values())
        raise RuntimeError(
            "Gyroid circular-edge pairing failed: expected 48 groups of 2, "
            f"got {len(groups)} groups with sizes {sizes}"
        )

    pairs = []
    for group in groups.values():
        lower, upper = sorted(group)
        c0, c1 = curves[lower], curves[upper]
        same = np.mod(c0 - c1 + 0.5, 1.0) - 0.5
        reverse = np.mod(c0 - c1[::-1] + 0.5, 1.0) - 0.5
        if np.max(np.linalg.norm(reverse, axis=1)) > 2.0e-5:
            raise RuntimeError("Paired Gyroid circular edges do not run oppositely")
        if np.mean(np.linalg.norm(same, axis=1)) <= np.mean(np.linalg.norm(reverse, axis=1)):
            raise RuntimeError("Paired Gyroid circular-edge orientation is ambiguous")
        # Put the upper copy next to the lower copy in Euclidean space.  This
        # is only non-zero for a macro patch crossing a periodic cell cut.
        shift = np.rint(np.mean(c0 - c1[::-1], axis=0))
        seam_error = np.max(np.linalg.norm(c0 - (c1[::-1] + shift), axis=1))
        if seam_error > 1.0e-9:
            raise RuntimeError(f"Gyroid hidden-edge mismatch: {seam_error:.3g}")
        pairs.append((lower, upper, shift))
    pairs.sort(key=lambda item: item[:2])
    _PAIR_CACHE = tuple(pairs)
    return _PAIR_CACHE


def _circle_reflect(omega):
    """Reflect across the circular edge; this is an anti-Mobius involution."""
    return (_ARC_CENTER + _ARC_RADIUS ** 2
            / (np.conj(omega) - np.conj(_ARC_CENTER)))


def _macro_square_domain(u, v):
    """Analytically map a square onto the unfolded genuine four-sided domain."""
    u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
    bottom = u * _BRANCH_RADIUS
    top = _circle_reflect(1j * (1.0 - u) * _BRANCH_RADIUS)
    left = 1j * v * _BRANCH_RADIUS
    right = _circle_reflect((1.0 - v) * _BRANCH_RADIUS)
    bilinear = (u * (1.0 - v) * _BRANCH_RADIUS
                + u * v * _REFLECTED_ORIGIN
                + (1.0 - u) * v * 1j * _BRANCH_RADIUS)
    return ((1.0 - v) * bottom + v * top
            + (1.0 - u) * left + u * right - bilinear)


def _macro_square_omega(u, v):
    """Return local triangle coordinates and the selected analytic branch."""
    unfolded = _macro_square_domain(u, v)
    upper = np.abs(unfolded - _ARC_CENTER) > _ARC_RADIUS + 1.0e-12
    omega = unfolded.copy()
    omega[upper] = 1j * np.conj(_circle_reflect(unfolded[upper]))
    return omega, upper


def _macro_patch_points(u, v, pair):
    """Evaluate one globally smooth analytic macro patch."""
    omega, upper_mask = _macro_square_omega(u, v)
    base = _weierstrass(omega).reshape(-1, 3)
    upper_mask = upper_mask.ravel()
    lower, upper, shift = pair
    points = np.empty_like(base)
    points[~upper_mask] = _apply_transform(
        base[~upper_mask], _TRIANGLE_TRANSFORMS[lower])
    points[upper_mask] = _apply_transform(
        base[upper_mask], _TRIANGLE_TRANSFORMS[upper]) + shift
    return points.reshape(np.shape(omega) + (3,))


def _weld(vertices, faces, tolerance=3.0e-6):
    vertices = np.asarray(vertices, dtype=float)
    keys = np.round(vertices / tolerance).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    return vertices[first], inverse[np.asarray(faces, dtype=np.int64)].astype(np.int32)


def _orient_faces(vertices, faces):
    """Orient every connected component consistently using edge adjacency."""
    faces = np.asarray(faces, dtype=np.int32).copy()
    edge_uses = defaultdict(list)
    for face_index, face in enumerate(faces):
        for k in range(4):
            a, b = int(face[k]), int(face[(k + 1) % 4])
            edge_uses[tuple(sorted((a, b)))].append((face_index, a, b))
    if any(len(uses) > 2 for uses in edge_uses.values()):
        raise RuntimeError("Gyroid mesh has a non-manifold edge")

    neighbors = defaultdict(list)
    for uses in edge_uses.values():
        if len(uses) == 2:
            (f0, a0, b0), (f1, a1, b1) = uses
            same_direction = a0 == a1 and b0 == b1
            neighbors[f0].append((f1, same_direction))
            neighbors[f1].append((f0, same_direction))
    flip = np.full(len(faces), -1, dtype=np.int8)
    for root in range(len(faces)):
        if flip[root] >= 0:
            continue
        flip[root] = 0
        queue = deque([root])
        while queue:
            current = queue.popleft()
            for other, toggle in neighbors[current]:
                wanted = flip[current] ^ int(toggle)
                if flip[other] < 0:
                    flip[other] = wanted
                    queue.append(other)
                elif flip[other] != wanted:
                    raise RuntimeError("Gyroid mesh is not consistently orientable")
    faces[flip == 1] = faces[flip == 1][:, [0, 3, 2, 1]]
    return faces


def _vertex_normals(vertices, faces):
    normals = np.zeros_like(vertices)
    for face in faces:
        p = vertices[face]
        normal = np.cross(p[1] - p[0], p[2] - p[0])
        normal += np.cross(p[2] - p[0], p[3] - p[0])
        normals[face] += normal
    lengths = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(lengths[:, None], 1.0e-30)
    return normals


def build_unit_cell(cell_size=1.0, quad_subdivisions=2):
    """Return a welded unit cell sampled from 48 genuine macro quads."""
    cell_size = float(cell_size)
    if not np.isfinite(cell_size) or cell_size <= 0.0:
        raise ValueError("cell_size must be a positive finite number")
    n = max(1, int(quad_subdivisions))
    axis = np.linspace(0.0, 1.0, n + 1)
    u, v = np.meshgrid(axis, axis, indexing='ij')
    pairs = _macro_pairs()

    local_idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    local_faces = np.stack((local_idx[:-1, :-1].ravel(),
                            local_idx[1:, :-1].ravel(),
                            local_idx[1:, 1:].ravel(),
                            local_idx[:-1, 1:].ravel()), axis=1)
    all_vertices, all_faces = [], []
    offset = 0
    for pair in pairs:
        patch = _macro_patch_points(u, v, pair).reshape(-1, 3)
        all_vertices.append(patch)
        all_faces.append(local_faces + offset)
        offset += len(patch)
    vertices, faces = _weld(np.vstack(all_vertices), np.vstack(all_faces))
    if len(faces) != 48 * n * n:
        raise RuntimeError("Incorrect Gyroid macro-quad face count")
    if np.any(np.apply_along_axis(lambda f: len(set(f)), 1, faces) != 4):
        raise RuntimeError("Degenerate Gyroid quad after welding")
    faces = _orient_faces(vertices, faces)

    p = vertices[faces]
    n1 = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    n2 = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    if (np.min(np.linalg.norm(n1, axis=1)) <= 1.0e-12 or
            np.min(np.linalg.norm(n2, axis=1)) <= 1.0e-12 or
            np.min(np.einsum('ij,ij->i', n1, n2)) <= 0.0):
        raise RuntimeError("Gyroid output contains a degenerate or folded quad")
    metrics = validate_unit_cell(vertices, faces)
    required = (
        metrics['non_quad_faces'] == 0,
        metrics['duplicate_faces'] == 0,
        metrics['degenerate_quads'] == 0,
        metrics['folded_quads'] == 0,
        metrics['nonmanifold_edges'] == 0,
        metrics['all_periodic_boundary_edges_paired'],
        metrics['quotient_all_edges_incident_to_two_quads'],
        metrics['quotient_euler_characteristic'] == -8,
        metrics['periodic_quotient_genus'] == 5,
    )
    if not all(required):
        raise RuntimeError(f"Gyroid topology acceptance failed: {metrics}")
    normals = _vertex_normals(vertices, faces)
    return vertices * cell_size, faces, normals


def macro_topology_counts():
    """Small public diagnostic used by tests and Blender-side reporting."""
    pairs = _macro_pairs()
    return {
        'fundamental_triangles': len(_TRIANGLE_TRANSFORMS),
        'circular_edge_groups': len(pairs),
        'triangles_per_group': tuple(len(pair[:2]) for pair in pairs),
        'genuine_quad_patches': len(pairs),
    }


def validate_unit_cell(vertices, faces, periodic_tolerance=2.0e-5):
    """Return the acceptance metrics from section 8/9 of the update spec."""
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    if faces.ndim != 2:
        raise ValueError("faces must be a rectangular index array")

    face_keys = [tuple(sorted(map(int, face))) for face in faces]
    spatial_edges = defaultdict(int)
    for face in faces:
        for k in range(len(face)):
            spatial_edges[tuple(sorted((int(face[k]),
                                        int(face[(k + 1) % len(face)]))))] += 1

    periodic = np.mod(vertices, 1.0)
    periodic[np.isclose(periodic, 1.0, atol=periodic_tolerance)] = 0.0
    _, quotient_vertex = np.unique(np.round(periodic, 5), axis=0,
                                   return_inverse=True)
    quotient_faces = quotient_vertex[faces]
    quotient_edges = defaultdict(list)
    for face in quotient_faces:
        for k in range(len(face)):
            directed = (int(face[k]), int(face[(k + 1) % len(face)]))
            quotient_edges[tuple(sorted(directed))].append(directed)
    quotient_vertices = len(np.unique(quotient_vertex))
    euler = quotient_vertices - len(quotient_edges) + len(faces)

    if faces.shape[1] == 4:
        p = vertices[faces]
        n1 = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        n2 = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
        area1 = np.linalg.norm(n1, axis=1)
        area2 = np.linalg.norm(n2, axis=1)
        folded = np.einsum('ij,ij->i', n1, n2) <= 0.0
        degenerate = (area1 <= 1.0e-12) | (area2 <= 1.0e-12)
    else:
        folded = np.ones(len(faces), dtype=bool)
        degenerate = np.ones(len(faces), dtype=bool)

    return {
        'vertices': len(vertices),
        'quad_faces': int(np.count_nonzero([len(face) == 4 for face in faces])),
        'non_quad_faces': int(np.count_nonzero([len(face) != 4 for face in faces])),
        'duplicate_faces': len(face_keys) - len(set(face_keys)),
        'degenerate_quads': int(np.count_nonzero(degenerate)),
        'folded_quads': int(np.count_nonzero(folded)),
        'nonmanifold_edges': int(sum(count > 2 for count in spatial_edges.values())),
        'all_periodic_boundary_edges_paired': all(
            len(uses) == 2 and uses[0] == uses[1][::-1]
            for uses in quotient_edges.values()),
        'quotient_all_edges_incident_to_two_quads': all(
            len(uses) == 2 for uses in quotient_edges.values()),
        'quotient_euler_characteristic': euler,
        'periodic_quotient_genus': 1 - euler // 2,
    }
