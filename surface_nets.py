"""Naive Surface Nets — dual iso-surface extraction.

Given a 3D scalar field, produces a closed quad mesh approximating the
iso-surface f = iso. Compared to Marching Cubes this:
  * Requires no 256-entry lookup table (≈100 lines instead of ≈1000).
  * Produces ~1 quad per cell with a sign change rather than ≤5 tris.
  * Handles the saddle topology of TPMS robustly.

The algorithm (Naive Surface Nets, S.F. Gibson 1998):

  1. For each *cell* (8 corners) whose corner signs are not uniform,
     place one *dual vertex* at the centroid of the field's
     zero-crossings on that cell's 12 edges.
  2. For each axis-aligned *grid edge* with a sign change, emit one
     quad joining the dual vertices of the (up to four) cells that
     share that edge. The quad's winding follows the sign-change
     direction so that all quads point outward (toward f > iso).

The output quads are all interior cells' quads, so the resulting mesh
is closed and orientable.
"""

import numpy as np


# Cube corner offsets, indexed by 3-bit (i + 2j + 4k).
# Kept as plain Python tuples so adding them to large array sizes (e.g.
# cnx >= 128) doesn't overflow an int8.
_CORNERS = (
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
    (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1),
)
_CORNERS_F = np.array(_CORNERS, dtype=np.float64)

# Cube edges, given as (corner_a, corner_b) pairs:
_EDGES = (
    # 4 edges along +X
    (0, 1), (2, 3), (4, 5), (6, 7),
    # 4 edges along +Y
    (0, 2), (1, 3), (4, 6), (5, 7),
    # 4 edges along +Z
    (0, 4), (1, 5), (2, 6), (3, 7),
)


def surface_nets(field, iso=0.0, spacing=(1.0, 1.0, 1.0),
                 origin=(0.0, 0.0, 0.0)):
    """Extract an iso-surface as quads.

    Parameters
    ----------
    field   : (nx,ny,nz) float ndarray — scalar samples on a regular grid.
    iso     : float                    — iso-value.
    spacing : (3,) float               — voxel size in world units.
    origin  : (3,) float               — world position of field[0,0,0].

    Returns
    -------
    verts : (V, 3) float ndarray
    quads : (Q, 4) int  ndarray  (CCW outward winding for f > iso)
    """
    f = np.asarray(field, dtype=np.float64) - iso
    nx, ny, nz = f.shape
    cnx, cny, cnz = nx - 1, ny - 1, nz - 1
    if cnx <= 0 or cny <= 0 or cnz <= 0:
        return np.zeros((0, 3)), np.zeros((0, 4), dtype=np.int64)

    sign = (f > 0.0).astype(np.uint8)              # corner signs, (nx,ny,nz)

    # --- 1. Sample 8 corner values per cell, build a "case index" ---------
    corner_vals = np.empty((8, cnx, cny, cnz), dtype=np.float64)
    corner_signs = np.empty((8, cnx, cny, cnz), dtype=np.uint8)
    for c, (di, dj, dk) in enumerate(_CORNERS):
        corner_vals[c]  = f   [di:di + cnx, dj:dj + cny, dk:dk + cnz]
        corner_signs[c] = sign[di:di + cnx, dj:dj + cny, dk:dk + cnz]

    case = np.zeros((cnx, cny, cnz), dtype=np.uint16)
    for c in range(8):
        case |= corner_signs[c].astype(np.uint16) << c
    active = (case != 0) & (case != 255)
    if not active.any():
        return np.zeros((0, 3)), np.zeros((0, 4), dtype=np.int64)

    # --- 2. Compute dual vertex per active cell --------------------------
    ai, aj, ak = np.nonzero(active)
    n_active = ai.size
    pos_sum = np.zeros((n_active, 3), dtype=np.float64)
    cnt = np.zeros(n_active, dtype=np.float64)

    for c0, c1 in _EDGES:
        v0 = corner_vals[c0, ai, aj, ak]
        v1 = corner_vals[c1, ai, aj, ak]
        crosses = (v0 > 0.0) != (v1 > 0.0)
        if not crosses.any():
            continue
        denom = v0 - v1
        # Where denom == 0, the edge is degenerate -> param midpoint:
        safe = np.where(np.abs(denom) > 1e-30, denom, 1.0)
        t = np.where(np.abs(denom) > 1e-30, v0 / safe, 0.5)
        t = np.clip(t, 0.0, 1.0)
        o0 = _CORNERS_F[c0]
        o1 = _CORNERS_F[c1]
        crossing = o0 + t[:, None] * (o1 - o0)        # (n_active, 3)
        m = crosses
        pos_sum[m] += crossing[m]
        cnt[m] += 1.0

    cnt_safe = np.where(cnt > 0.0, cnt, 1.0)
    local = pos_sum / cnt_safe[:, None]               # in cell-local [0,1]^3

    sx, sy, sz = spacing
    ox, oy, oz = origin
    verts = np.column_stack([
        (ai + local[:, 0]) * sx + ox,
        (aj + local[:, 1]) * sy + oy,
        (ak + local[:, 2]) * sz + oz,
    ])

    # Map (cell i,j,k) -> dual vertex index, -1 if cell is not active
    vert_idx = np.full((cnx, cny, cnz), -1, dtype=np.int64)
    vert_idx[ai, aj, ak] = np.arange(n_active)

    # --- 3. Build quads from interior grid edges with sign change --------
    quads_list = []

    # -- X-aligned grid edges:
    # An X-edge between corners (i,j,k) and (i+1,j,k) is shared by 4 cells:
    #   cells (i, j-1, k-1), (i, j, k-1), (i, j-1, k), (i, j, k).
    # For all 4 to exist we need 0<=i<=cnx-1, 1<=j<=cny, 1<=k<=cnz, i.e.
    # interior j in [1..ny-2], k in [1..nz-2].
    diff_x = sign[1:, :, :] != sign[:-1, :, :]           # (cnx, ny, nz)
    interior_x = diff_x[:, 1:-1, 1:-1]                   # (cnx, ny-2, nz-2)
    ix, jx, kx = np.nonzero(interior_x)
    if ix.size:
        j_corner = jx + 1
        k_corner = kx + 1
        c00 = vert_idx[ix, j_corner - 1, k_corner - 1]
        c10 = vert_idx[ix, j_corner,     k_corner - 1]
        c11 = vert_idx[ix, j_corner,     k_corner    ]
        c01 = vert_idx[ix, j_corner - 1, k_corner    ]
        # Orient: +X normal if sign goes from negative (at i) to positive (at i+1)
        flip = sign[ix, j_corner, k_corner] > sign[ix + 1, j_corner, k_corner]
        q = np.column_stack([c00, c10, c11, c01])
        # Flip winding where needed
        q_flip = np.column_stack([c00, c01, c11, c10])
        q = np.where(flip[:, None], q_flip, q)
        quads_list.append(q)

    # -- Y-aligned grid edges --
    diff_y = sign[:, 1:, :] != sign[:, :-1, :]           # (nx, cny, nz)
    interior_y = diff_y[1:-1, :, 1:-1]                   # (nx-2, cny, nz-2)
    iy, jy, ky = np.nonzero(interior_y)
    if iy.size:
        i_corner = iy + 1
        k_corner = ky + 1
        c00 = vert_idx[i_corner - 1, jy, k_corner - 1]
        c10 = vert_idx[i_corner,     jy, k_corner - 1]
        c11 = vert_idx[i_corner,     jy, k_corner    ]
        c01 = vert_idx[i_corner - 1, jy, k_corner    ]
        flip = sign[i_corner, jy, k_corner] > sign[i_corner, jy + 1, k_corner]
        # For Y edges we want opposite default winding so outward normal works:
        q       = np.column_stack([c00, c01, c11, c10])
        q_flip  = np.column_stack([c00, c10, c11, c01])
        q = np.where(flip[:, None], q_flip, q)
        quads_list.append(q)

    # -- Z-aligned grid edges --
    diff_z = sign[:, :, 1:] != sign[:, :, :-1]           # (nx, ny, cnz)
    interior_z = diff_z[1:-1, 1:-1, :]                   # (nx-2, ny-2, cnz)
    iz, jz, kz = np.nonzero(interior_z)
    if iz.size:
        i_corner = iz + 1
        j_corner = jz + 1
        c00 = vert_idx[i_corner - 1, j_corner - 1, kz]
        c10 = vert_idx[i_corner,     j_corner - 1, kz]
        c11 = vert_idx[i_corner,     j_corner,     kz]
        c01 = vert_idx[i_corner - 1, j_corner,     kz]
        flip = sign[i_corner, j_corner, kz] > sign[i_corner, j_corner, kz + 1]
        q       = np.column_stack([c00, c10, c11, c01])
        q_flip  = np.column_stack([c00, c01, c11, c10])
        q = np.where(flip[:, None], q_flip, q)
        quads_list.append(q)

    if not quads_list:
        return verts, np.zeros((0, 4), dtype=np.int64)

    quads = np.concatenate(quads_list, axis=0)
    # Drop any quad that touches a missing dual vertex (shouldn't happen for
    # interior edges, but guard against the edge case where the field is flat
    # in one of the four sharing cells):
    valid = (quads >= 0).all(axis=1)
    quads = quads[valid]
    return verts, quads
