"""TPMS generation operator.

Pipeline:
    1. Compute the target's *world-space* axis-aligned bounding box and
       expand by `padding`.
    2. Build a regular voxel grid spanning that bbox; voxel pitch ≈
       `cell_size / resolution` along each axis (clamped so the grid never
       exceeds a sane max so we don't blow up memory).
    3. Evaluate the chosen TPMS scalar field f on the grid.
    4. Convert to a "solid-mode" field g:
            SHELL  : g = |f - iso| - thickness   →  iso-surface g=0 is closed
            VOLUME : g = f - iso                 →  iso-surface separates phases
    5. Run Naive Surface Nets to extract a closed quad mesh from g.
    6. Create a new Blender object with the mesh.
    7. If `clip_to_target`: add a Boolean (INTERSECT) modifier referencing
       the target. Optionally apply it.
"""

import time

import bmesh
import bpy
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree

from . import tpms_field, surface_nets


# Hard cap on total voxel count to avoid OOM. ~30 M voxels ≈ 240 MB float64.
_MAX_VOXELS = 30_000_000


def _world_bbox(obj):
    """Axis-aligned world-space bbox of `obj` from its 8 corner vertices."""
    corners = [obj.matrix_world @ Vector(c) for c in obj.bound_box]
    xs = [c.x for c in corners]
    ys = [c.y for c in corners]
    zs = [c.z for c in corners]
    return (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))


def _build_grid(bmin, bmax, cell_size, resolution):
    """Build a regular sampling grid for the TPMS field.

    The grid is extended by 2 voxels on each side beyond `bmin`/`bmax` so the
    outer boundary-cap layer is guaranteed to be *outside* the target. That
    way the Boolean INTERSECT cleanly removes the cap and we never see the
    flat AABB patches in the final result.
    """
    pitch = cell_size / max(int(resolution), 2)

    # Pad by 2 voxels beyond the user bbox so the cap layer is always outside.
    extra = 2.0 * pitch
    bmin = (bmin[0] - extra, bmin[1] - extra, bmin[2] - extra)
    bmax = (bmax[0] + extra, bmax[1] + extra, bmax[2] + extra)

    nx = max(int(np.ceil((bmax[0] - bmin[0]) / pitch)) + 1, 2)
    ny = max(int(np.ceil((bmax[1] - bmin[1]) / pitch)) + 1, 2)
    nz = max(int(np.ceil((bmax[2] - bmin[2]) / pitch)) + 1, 2)

    # Memory guard: rescale pitch if necessary.
    total = nx * ny * nz
    if total > _MAX_VOXELS:
        scale = (total / _MAX_VOXELS) ** (1.0 / 3.0)
        pitch *= scale
        nx = max(int(np.ceil((bmax[0] - bmin[0]) / pitch)) + 1, 2)
        ny = max(int(np.ceil((bmax[1] - bmin[1]) / pitch)) + 1, 2)
        nz = max(int(np.ceil((bmax[2] - bmin[2]) / pitch)) + 1, 2)

    xs = bmin[0] + np.arange(nx) * pitch
    ys = bmin[1] + np.arange(ny) * pitch
    zs = bmin[2] + np.arange(nz) * pitch
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    return X, Y, Z, (pitch, pitch, pitch), (xs[0], ys[0], zs[0]), (nx, ny, nz)


def _build_target_bvh(target, depsgraph):
    """Build a BVHTree of `target` in *world space* (modifiers applied)."""
    eval_obj = target.evaluated_get(depsgraph)
    mesh = eval_obj.to_mesh()
    try:
        M = target.matrix_world
        verts = [tuple((M @ v.co)[:]) for v in mesh.vertices]
        polys = [list(p.vertices) for p in mesh.polygons]
        bvh = BVHTree.FromPolygons(verts, polys, all_triangles=False)
    finally:
        eval_obj.to_mesh_clear()
    return bvh


# Five non-axis-aligned ray directions for majority-vote in/out test.
# More rays = more robust against small holes / open caps. Three rays were
# not enough for default Suzanne: her two earring rings + open neck caused
# correlated parity errors in roughly half the directions.
_RAY_DIRS = (
    Vector(( 1.000,  0.137,  0.043)).normalized(),
    Vector(( 0.073,  1.000,  0.310)).normalized(),
    Vector((-0.591,  0.183,  1.000)).normalized(),
    Vector(( 0.701, -0.503,  0.612)).normalized(),
    Vector(( 0.234,  0.787, -0.661)).normalized(),
)
_RAY_MAJORITY = (len(_RAY_DIRS) // 2) + 1  # ≥3 of 5


def _ray_inside(bvh, p, direction):
    """Parity test: cast ray in `direction`, count hits, odd => inside."""
    eps_offset = 1e-5
    hits = 0
    cur = p.copy()
    for _ in range(64):
        hloc, _, _, _ = bvh.ray_cast(cur, direction)
        if hloc is None:
            break
        hits += 1
        cur = hloc + direction * eps_offset
    return (hits & 1) == 1


def _sd_at_points(bvh, points, robust=False):
    """Signed distance from each row of `points` (N,3) to the target mesh.

    Negative inside, positive outside. Two strategies:

      robust=False : `find_nearest` + nearest-face outward-normal dot test.
                     Fast but unreliable in concave regions.
      robust=True  : 3-ray majority-vote parity test. Handles concavities
                     and small holes (Suzanne eye sockets, open neck) well.
    """
    n = points.shape[0]
    sd = np.empty(n, dtype=np.float64)
    if robust:
        for i in range(n):
            p = Vector((float(points[i, 0]),
                        float(points[i, 1]),
                        float(points[i, 2])))
            loc, _, _, dist = bvh.find_nearest(p)
            if loc is None:
                sd[i] = 1e9
                continue
            votes = 0
            for d in _RAY_DIRS:
                if _ray_inside(bvh, p, d):
                    votes += 1
            inside = votes >= _RAY_MAJORITY
            sd[i] = -float(dist) if inside else float(dist)
    else:
        for i in range(n):
            p = Vector((float(points[i, 0]),
                        float(points[i, 1]),
                        float(points[i, 2])))
            loc, nrm, _, dist = bvh.find_nearest(p)
            if loc is None:
                sd[i] = 1e9
                continue
            outward = (p - loc).dot(nrm)
            sd[i] = float(dist) if outward >= 0.0 else -float(dist)
    return sd


def _trilinear_upsample(field, target_shape):
    """Upsample a coarse 3D field to `target_shape` via trilinear interp."""
    snx, sny, snz = field.shape
    nx, ny, nz = target_shape
    xs = np.linspace(0.0, snx - 1.0, nx)
    ys = np.linspace(0.0, sny - 1.0, ny)
    zs = np.linspace(0.0, snz - 1.0, nz)
    ix = np.clip(xs.astype(np.int64), 0, snx - 2)
    iy = np.clip(ys.astype(np.int64), 0, sny - 2)
    iz = np.clip(zs.astype(np.int64), 0, snz - 2)
    fx = (xs - ix).astype(np.float64)
    fy = (ys - iy).astype(np.float64)
    fz = (zs - iz).astype(np.float64)
    IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing='ij')
    FX, FY, FZ = np.meshgrid(fx, fy, fz, indexing='ij')
    f000 = field[IX,     IY,     IZ    ]
    f100 = field[IX + 1, IY,     IZ    ]
    f010 = field[IX,     IY + 1, IZ    ]
    f001 = field[IX,     IY,     IZ + 1]
    f110 = field[IX + 1, IY + 1, IZ    ]
    f101 = field[IX + 1, IY,     IZ + 1]
    f011 = field[IX,     IY + 1, IZ + 1]
    f111 = field[IX + 1, IY + 1, IZ + 1]
    f00 = f000 * (1.0 - FX) + f100 * FX
    f10 = f010 * (1.0 - FX) + f110 * FX
    f01 = f001 * (1.0 - FX) + f101 * FX
    f11 = f011 * (1.0 - FX) + f111 * FX
    f0 = f00 * (1.0 - FY) + f10 * FY
    f1 = f01 * (1.0 - FY) + f11 * FY
    return f0 * (1.0 - FZ) + f1 * FZ


def _signed_distance_grid(bvh, X, Y, Z, pitch, robust=False, subsample=1,
                          refine_band=2.0):
    """Two-pass SDF on the voxel grid.

    Pass 1: SDF on a subsampled grid (every `subsample` voxels), then
            trilinear-upsample. Localisation error ≲ subsample·pitch.
    Pass 2: For voxels with |sd_upsampled| < refine_band·subsample·pitch
            (i.e. the band where the coarse SDF could still be on the wrong
            side of the boundary), recompute SDF exactly via BVH.

    The narrow refinement band has ~5–10× fewer voxels than the full grid,
    so the second pass is much cheaper than uniform full-resolution SDF.
    """
    s = max(int(subsample), 1)
    if s == 1:
        pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
        return _sd_at_points(bvh, pts, robust=robust).reshape(X.shape)

    # Pass 1: coarse SDF + trilinear upsample.
    sX = np.ascontiguousarray(X[::s, ::s, ::s])
    sY = np.ascontiguousarray(Y[::s, ::s, ::s])
    sZ = np.ascontiguousarray(Z[::s, ::s, ::s])
    pts_coarse = np.stack([sX.ravel(), sY.ravel(), sZ.ravel()], axis=1)
    sd_coarse = _sd_at_points(bvh, pts_coarse, robust=robust).reshape(sX.shape)
    sd_full = _trilinear_upsample(sd_coarse, X.shape)

    # Pass 2: refine the narrow band where the coarse SDF could be wrong.
    band = float(refine_band) * s * float(pitch)
    near_mask = np.abs(sd_full) < band
    if near_mask.any():
        near_pts = np.stack(
            [X[near_mask], Y[near_mask], Z[near_mask]], axis=-1
        )
        sd_near = _sd_at_points(bvh, near_pts, robust=robust)
        sd_full[near_mask] = sd_near

    return sd_full


def _project_to_isosurface(verts, tpms_type, solid_mode, iso, thickness,
                           cell_size, origin, euler, pitch, iters=3,
                           interior_mask=None):
    """Snap interior verts onto the actual TPMS iso-surface using Newton
    iterations against the analytic field.

    `interior_mask` (N,) bool: True for verts that lie on the TPMS surface
    (and should be projected). Verts on the target boundary (where the SDF
    constraint dominates the iso-surface) should be False so they stay put.
    """
    if iters <= 0 or verts.shape[0] == 0:
        return verts

    if interior_mask is None:
        interior = np.ones(verts.shape[0], dtype=bool)
    else:
        interior = interior_mask
    if not interior.any():
        return verts

    eps = max(cell_size * 1e-3, 1e-6)
    max_step = 0.4 * pitch

    pts = verts[interior].astype(np.float64, copy=True)
    for _ in range(iters):
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
        f   = tpms_field.sample_field(tpms_type, x,       y,       z,       cell_size, origin, euler)
        fxp = tpms_field.sample_field(tpms_type, x + eps, y,       z,       cell_size, origin, euler)
        fxm = tpms_field.sample_field(tpms_type, x - eps, y,       z,       cell_size, origin, euler)
        fyp = tpms_field.sample_field(tpms_type, x,       y + eps, z,       cell_size, origin, euler)
        fym = tpms_field.sample_field(tpms_type, x,       y - eps, z,       cell_size, origin, euler)
        fzp = tpms_field.sample_field(tpms_type, x,       y,       z + eps, cell_size, origin, euler)
        fzm = tpms_field.sample_field(tpms_type, x,       y,       z - eps, cell_size, origin, euler)
        inv = 1.0 / (2.0 * eps)
        gx = (fxp - fxm) * inv
        gy = (fyp - fym) * inv
        gz = (fzp - fzm) * inv

        if solid_mode == 'SHELL':
            target = np.where(f > iso, iso + thickness, iso - thickness)
        else:
            target = iso

        delta = f - target
        grad_sq = gx * gx + gy * gy + gz * gz
        grad_sq = np.where(grad_sq < 1e-12, 1.0, grad_sq)
        s = delta / grad_sq
        sx, sy, sz = s * gx, s * gy, s * gz

        step_len = np.sqrt(sx * sx + sy * sy + sz * sz)
        clip = np.minimum(1.0, max_step / np.maximum(step_len, 1e-30))
        pts[:, 0] -= sx * clip
        pts[:, 1] -= sy * clip
        pts[:, 2] -= sz * clip

    verts[interior] = pts
    return verts


def _remove_small_components(verts, quads, min_faces):
    """Drop connected components with fewer than `min_faces` quads.

    Used as a last-mile fix for boundary fragments that survive SDF clipping
    when the in/out test misclassifies a few voxels. The legitimate TPMS
    body is always one (or a few) huge component(s); spurious fragments are
    tiny isolated islands, so a face-count threshold cleanly separates them.

    Returns filtered (verts, quads) with vertex indices renumbered.
    """
    if min_faces <= 1 or quads.shape[0] == 0:
        return verts, quads

    n = int(verts.shape[0])
    parent = list(range(n))
    rank = [0] * n

    def find(a):
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    quads_list = quads.tolist()
    for q in quads_list:
        a, b, c, d = q
        union(a, b); union(b, c); union(c, d)

    quad_roots = np.fromiter((find(q[0]) for q in quads_list),
                             dtype=np.int64, count=len(quads_list))
    unique, counts = np.unique(quad_roots, return_counts=True)
    big_set = set(unique[counts >= int(min_faces)].tolist())
    if not big_set:
        return (np.zeros((0, 3), dtype=verts.dtype),
                np.zeros((0, 4), dtype=quads.dtype))

    keep = np.array([r in big_set for r in quad_roots], dtype=bool)
    kept_quads = quads[keep]

    used = np.zeros(n, dtype=bool)
    used[kept_quads.ravel()] = True
    new_idx = np.full(n, -1, dtype=np.int64)
    new_idx[used] = np.arange(int(used.sum()))

    return verts[used], new_idx[kept_quads]


def _make_mesh_object(name, verts, quads, collection, smooth):
    mesh = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bv = [bm.verts.new(v.tolist()) for v in verts]
    bm.verts.ensure_lookup_table()
    for q in quads:
        try:
            bm.faces.new((bv[int(q[0])], bv[int(q[1])],
                          bv[int(q[2])], bv[int(q[3])]))
        except ValueError:
            # Duplicate face from periodic topology — skip silently.
            pass
    bm.normal_update()
    if smooth:
        for f in bm.faces:
            f.smooth = True
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()

    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    return obj


class TPMS_OT_generate(bpy.types.Operator):
    """Generate a TPMS solid inside the selected target mesh."""
    bl_idname = "tpms.generate"
    bl_label = "Generate TPMS"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        props = context.scene.tpms_props
        return props.target is not None and props.target.type == 'MESH'

    def execute(self, context):
        props = context.scene.tpms_props
        target = props.target
        if target is None or target.type != 'MESH':
            self.report({'ERROR'}, "Pick a mesh as Target Mesh first.")
            return {'CANCELLED'}

        t0 = time.perf_counter()

        # 1. World bbox + padding
        bmin, bmax = _world_bbox(target)
        pad = props.padding
        bmin = (bmin[0] - pad, bmin[1] - pad, bmin[2] - pad)
        bmax = (bmax[0] + pad, bmax[1] + pad, bmax[2] + pad)

        # 2. Sampling grid
        X, Y, Z, spacing, origin, dims = _build_grid(
            bmin, bmax, props.cell_size, props.resolution
        )

        # 3. TPMS scalar field
        f = tpms_field.sample_field(
            props.tpms_type, X, Y, Z,
            cell_size=props.cell_size,
            origin=tuple(props.origin),
            euler=tuple(props.rotation),
        )

        # 4. Convert to solid-mode field
        if props.solid_mode == 'SHELL':
            g_tpms = np.abs(f - props.iso_value) - props.thickness
        else:  # 'VOLUME'
            g_tpms = f - props.iso_value

        # 4b. Clip in field-space using the target's signed distance, instead
        # of relying on a Boolean modifier (which is fragile when the target
        # mesh is non-closed or when the SN output has non-planar quads).
        pitch = spacing[0]
        sd = None
        bvh = None
        if props.clip_to_target:
            depsgraph = context.evaluated_depsgraph_get()
            try:
                bvh = _build_target_bvh(target, depsgraph)
            except Exception as ex:
                self.report({'ERROR'},
                            f"Could not build BVH from target: {ex}")
                return {'CANCELLED'}
            t_sd = time.perf_counter()
            sd = _signed_distance_grid(
                bvh, X, Y, Z, pitch,
                robust=props.robust_clip,
                subsample=props.sdf_subsample,
            )
            t_sd = time.perf_counter() - t_sd
            # Optional inset: shrink "inside" by N voxels to suppress any
            # residual boundary fuzz from imperfect in/out classification.
            if props.boundary_inset > 0.0:
                sd = sd + props.boundary_inset * pitch
            # g_combined = 0 is exactly: TPMS surface (inside target) ∪
            # target surface (inside shell). max() composes them as an SDF.
            g = np.maximum(g_tpms, sd)
        else:
            g = g_tpms
            # Without an SDF, still close the iso-surface at bbox so the
            # output is a closed manifold.
            cap = float(np.abs(g).max() + 1.0)
            g[ 0, :, :] = cap; g[-1, :, :] = cap
            g[:,  0, :] = cap; g[:, -1, :] = cap
            g[:, :,  0] = cap; g[:, :, -1] = cap

        # 5. Iso-surface (g = 0)
        verts, quads = surface_nets.surface_nets(
            g, iso=0.0, spacing=spacing, origin=origin
        )

        if verts.shape[0] == 0 or quads.shape[0] == 0:
            self.report({'WARNING'},
                        "No iso-surface found. Try adjusting Iso Value, "
                        "Thickness, or increasing Resolution, or check that "
                        "the target overlaps the TPMS pattern.")
            return {'CANCELLED'}

        # 5b. Project verts onto the analytic TPMS iso-surface. Skip verts
        # whose iso-surface is determined by the *target* SDF (those are on
        # the target's surface; they shouldn't be pulled toward the TPMS).
        interior_mask = None
        if bvh is not None:
            sd_v = _sd_at_points(bvh, verts, robust=False)
            # Vert is on the TPMS surface (not the target boundary) iff
            # the target SDF is well negative there.
            interior_mask = sd_v < -1.2 * pitch

        verts = _project_to_isosurface(
            verts.astype(np.float64, copy=True),
            tpms_type=props.tpms_type,
            solid_mode=props.solid_mode,
            iso=props.iso_value,
            thickness=props.thickness,
            cell_size=props.cell_size,
            origin=tuple(props.origin),
            euler=tuple(props.rotation),
            pitch=pitch,
            iters=props.project_iters,
            interior_mask=interior_mask,
        )

        # 5c. Drop tiny disconnected fragments (last-mile fix for any
        # voxels misclassified by the in/out test).
        n_before = quads.shape[0]
        verts, quads = _remove_small_components(
            verts, quads, props.min_component_faces
        )
        dropped = n_before - quads.shape[0]
        if quads.shape[0] == 0:
            self.report({'WARNING'},
                        "All components were smaller than the threshold. "
                        "Lower 'Min Component Faces' or check Iso/Thickness.")
            return {'CANCELLED'}

        # 6. Build object
        name = f"TPMS_{props.tpms_type.title()}"
        obj = _make_mesh_object(
            name, verts, quads, context.collection, props.smooth_shade
        )

        # Select the new object
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        context.view_layer.objects.active = obj

        dt = time.perf_counter() - t0
        self.report(
            {'INFO'},
            f"{name}: {dims[0]}×{dims[1]}×{dims[2]} voxels → "
            f"{len(verts)} verts, {len(quads)} quads "
            f"(dropped {dropped} fragment quads) in {dt:.2f}s"
        )
        return {'FINISHED'}


_classes = (TPMS_OT_generate,)


def register():
    for cls in _classes:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_classes):
        bpy.utils.unregister_class(cls)
