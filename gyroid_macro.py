"""Exact Gyroid unit cell made from 48 genuine quadrilateral patches.

The construction follows ``docs/Gyroid_Blender_Plugin_MacroQuad_Update_Spec.md``:
96 Weierstrass fundamental triangles are paired across their circular edges,
then a cached harmonic map removes the parameter kink at the hidden diagonal.
Only the final quadrilateral grid is returned to Blender.
"""

from collections import defaultdict, deque

import numpy as np


_SQ2 = np.sqrt(2.0)
_THETA = None
_KAPPA = None
_GL_CACHE = {}
_HARMONIC_CACHE = {}
_PAIR_CACHE = {}


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


def _quadrature(order):
    order = max(32, int(order))
    cached = _GL_CACHE.get(order)
    if cached is None:
        nodes, weights = np.polynomial.legendre.leggauss(order)
        s = 0.5 * (nodes + 1.0)
        weights = 0.5 * weights
        t = 2.0 * s - s * s
        weights = weights * 2.0 * (1.0 - s)
        cached = (t, weights)
        _GL_CACHE[order] = cached
    return cached


def _weierstrass(omega, order):
    """Evaluate the exact fundamental triangle along continuous radial roots."""
    shape = np.shape(omega)
    om = np.asarray(omega, dtype=complex).reshape(-1)
    t, weights = _quadrature(order)
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


def _omega_from_barycentric(lam_p, lam_r):
    r = lam_p + lam_r
    ratio = np.divide(lam_r, r, out=np.zeros_like(r, dtype=float), where=r > 0.0)
    phi = 0.5 * np.pi * ratio
    return r * _rho_max(phi) * np.exp(1j * phi)


def _apply_transform(points, transform):
    a, b = transform
    return np.asarray(points) @ a.T + b


def _periodic_curve_key(curve):
    c = np.mod(curve, 1.0)
    c[np.isclose(c, 1.0, atol=2.0e-5)] = 0.0
    forward = tuple(np.round(c, 5).ravel())
    backward = tuple(np.round(c[::-1], 5).ravel())
    return min(forward, backward)


def _macro_pairs(order):
    """Discover the 48 circular-edge pairs; pair numbers are never hard-coded."""
    cache_order = max(32, int(order))
    cached = _PAIR_CACHE.get(cache_order)
    if cached is not None:
        return cached

    phi = np.linspace(0.0, 0.5 * np.pi, 17)
    arc = _weierstrass(_rho_max(phi) * np.exp(1j * phi), cache_order)
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
    cached = tuple(pairs)
    _PAIR_CACHE[cache_order] = cached
    return cached


def _piecewise_points(s, t, pair, order):
    """Initial exact C0 map P0 from a square split along t=s."""
    s, t = np.broadcast_arrays(np.asarray(s, dtype=float), np.asarray(t, dtype=float))
    flat_s, flat_t = s.ravel(), t.ravel()
    lower_mask = flat_t <= flat_s
    points = np.empty((flat_s.size, 3), dtype=float)
    lower, upper, shift = pair

    lp = np.empty(flat_s.size)
    lr = np.empty(flat_s.size)
    lp[lower_mask] = 1.0 - flat_s[lower_mask]
    lr[lower_mask] = flat_t[lower_mask]
    lp[~lower_mask] = flat_s[~lower_mask]
    lr[~lower_mask] = 1.0 - flat_t[~lower_mask]
    base = _weierstrass(_omega_from_barycentric(lp, lr), order).reshape(-1, 3)
    points[lower_mask] = _apply_transform(
        base[lower_mask], _TRIANGLE_TRANSFORMS[lower])
    points[~lower_mask] = _apply_transform(
        base[~lower_mask], _TRIANGLE_TRANSFORMS[upper]) + shift
    return points.reshape(s.shape + (3,))


def _triangles_for_grid(resolution):
    idx = np.arange((resolution + 1) ** 2).reshape(resolution + 1,
                                                            resolution + 1)
    tris = []
    for i in range(resolution):
        for j in range(resolution):
            a, b = idx[i, j], idx[i + 1, j]
            c, d = idx[i, j + 1], idx[i + 1, j + 1]
            tris.append((a, b, d))
            tris.append((a, d, c))
    return np.asarray(tris, dtype=np.int32)


def _cotangent_adjacency(points, triangles):
    weights = defaultdict(float)
    for tri in triangles:
        ids = [int(x) for x in tri]
        p = points[ids]
        twice_area = np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        if twice_area <= 1.0e-15:
            raise RuntimeError("Degenerate auxiliary triangle in harmonic solver")
        for corner, edge in ((0, (1, 2)), (1, (2, 0)), (2, (0, 1))):
            v0 = p[edge[0]] - p[corner]
            v1 = p[edge[1]] - p[corner]
            cot = float(np.dot(v0, v1) / twice_area)
            a, b = sorted((ids[edge[0]], ids[edge[1]]))
            weights[(a, b)] += 0.5 * cot
    adjacency = [[] for _ in range(len(points))]
    for (a, b), weight in weights.items():
        adjacency[a].append((b, weight))
        adjacency[b].append((a, weight))
    return adjacency


def _cg(matvec, rhs, initial, tolerance=1.0e-11, max_iterations=10000):
    x = initial.copy()
    residual = rhs - matvec(x)
    direction = residual.copy()
    rr = float(np.dot(residual, residual))
    target = tolerance * max(1.0, float(np.linalg.norm(rhs)))
    if np.sqrt(rr) <= target:
        return x
    for _ in range(max_iterations):
        ad = matvec(direction)
        denom = float(np.dot(direction, ad))
        if denom <= 0.0:
            break
        alpha = rr / denom
        x += alpha * direction
        residual -= alpha * ad
        next_rr = float(np.dot(residual, residual))
        if np.sqrt(next_rr) <= target:
            return x
        direction = residual + (next_rr / rr) * direction
        rr = next_rr
    raise RuntimeError("Cotangent harmonic solve did not converge")


def _harmonic_template(solver_resolution, order):
    key = (max(8, int(solver_resolution)), max(32, int(order)))
    cached = _HARMONIC_CACHE.get(key)
    if cached is not None:
        return cached
    resolution = key[0]
    axis = np.linspace(0.0, 1.0, resolution + 1)
    ss, tt = np.meshgrid(axis, axis, indexing='ij')
    st = np.stack((ss, tt), axis=-1).reshape(-1, 2)
    pair = _macro_pairs(order)[0]
    points = _piecewise_points(ss, tt, pair, order).reshape(-1, 3)
    triangles = _triangles_for_grid(resolution)
    adjacency = _cotangent_adjacency(points, triangles)

    boundary = ((st[:, 0] == 0.0) | (st[:, 0] == 1.0) |
                (st[:, 1] == 0.0) | (st[:, 1] == 1.0))
    interior_ids = np.flatnonzero(~boundary)
    interior_lookup = np.full(len(st), -1, dtype=np.int32)
    interior_lookup[interior_ids] = np.arange(len(interior_ids))

    diagonal = np.empty(len(interior_ids))
    rhs = np.zeros((len(interior_ids), 2))
    neighbor_rows = []
    neighbor_weights = []
    for row, vertex in enumerate(interior_ids):
        diagonal[row] = sum(weight for _, weight in adjacency[vertex])
        rows, values = [], []
        for neighbor, weight in adjacency[vertex]:
            neighbor_row = interior_lookup[neighbor]
            if neighbor_row >= 0:
                rows.append(int(neighbor_row))
                values.append(float(weight))
            else:
                rhs[row] += weight * st[neighbor]
        neighbor_rows.append(np.asarray(rows, dtype=np.int32))
        neighbor_weights.append(np.asarray(values, dtype=float))

    def matvec(values):
        result = diagonal * values
        for row, (cols, weights) in enumerate(zip(neighbor_rows, neighbor_weights)):
            result[row] -= np.dot(weights, values[cols])
        return result

    uv = st.copy()
    for component in range(2):
        uv[interior_ids, component] = _cg(
            matvec, rhs[:, component], st[interior_ids, component])

    edge_1 = uv[triangles[:, 1]] - uv[triangles[:, 0]]
    edge_2 = uv[triangles[:, 2]] - uv[triangles[:, 0]]
    signed = edge_1[:, 0] * edge_2[:, 1] - edge_1[:, 1] * edge_2[:, 0]
    if np.min(signed) <= 1.0e-12:
        raise RuntimeError(
            f"Harmonic map contains a flipped triangle (min area {np.min(signed):.3g})")
    cached = (st, uv, triangles)
    _HARMONIC_CACHE[key] = cached
    return cached


def _inverse_harmonic_grid(subdivisions, solver_resolution, order):
    st, uv, triangles = _harmonic_template(solver_resolution, order)
    axis = np.linspace(0.0, 1.0, subdivisions + 1)
    uu, vv = np.meshgrid(axis, axis, indexing='ij')
    targets = np.stack((uu, vv), axis=-1).reshape(-1, 2)
    result = np.empty_like(targets)
    a = uv[triangles[:, 0]]
    ab = uv[triangles[:, 1]] - a
    ac = uv[triangles[:, 2]] - a
    denominator = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    lo = np.min(uv[triangles], axis=1) - 2.0e-12
    hi = np.max(uv[triangles], axis=1) + 2.0e-12
    for index, target in enumerate(targets):
        candidates = np.flatnonzero(np.all(target >= lo, axis=1) &
                                    np.all(target <= hi, axis=1))
        delta = target - a[candidates]
        b1 = (delta[:, 0] * ac[candidates, 1] -
              delta[:, 1] * ac[candidates, 0]) / denominator[candidates]
        b2 = (ab[candidates, 0] * delta[:, 1] -
              ab[candidates, 1] * delta[:, 0]) / denominator[candidates]
        bary = np.stack((1.0 - b1 - b2, b1, b2), axis=1)
        best = int(np.argmax(np.min(bary, axis=1)))
        if np.min(bary[best]) < -2.0e-8:
            raise RuntimeError(f"Could not invert harmonic map at {target}")
        tri = triangles[candidates[best]]
        result[index] = bary[best] @ st[tri]
    return result.reshape(subdivisions + 1, subdivisions + 1, 2)


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
    faces[flip == 1] = faces[flip == 1, ::-1]
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


def build_unit_cell(cell_size=1.0, quad_subdivisions=2,
                    solver_resolution=44, quadrature_order=200):
    """Return a welded unit cell sampled from 48 genuine macro quads."""
    cell_size = float(cell_size)
    if not np.isfinite(cell_size) or cell_size <= 0.0:
        raise ValueError("cell_size must be a positive finite number")
    n = max(1, int(quad_subdivisions))
    order = max(32, int(quadrature_order))
    inverse = _inverse_harmonic_grid(n, solver_resolution, order)
    s, t = inverse[..., 0], inverse[..., 1]
    pairs = _macro_pairs(order)

    local_idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    local_faces = np.stack((local_idx[:-1, :-1].ravel(),
                            local_idx[1:, :-1].ravel(),
                            local_idx[1:, 1:].ravel(),
                            local_idx[:-1, 1:].ravel()), axis=1)
    all_vertices, all_faces = [], []
    offset = 0
    for pair in pairs:
        patch = _piecewise_points(s, t, pair, order).reshape(-1, 3)
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


def macro_topology_counts(quadrature_order=200):
    """Small public diagnostic used by tests and Blender-side reporting."""
    pairs = _macro_pairs(quadrature_order)
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
