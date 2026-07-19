"""TPMS generation operator.

Pipeline
--------
1. `weierstrass.build_unit_cell` evaluates the exact Enneper-Weierstrass
   parametrization of the chosen surface.  Gyroid uses 48 paired-triangle
   macro patches with cached harmonic square coordinates; P and D retain
   their existing exact patch grids.
   Every vertex lies on the exact minimal surface and boundary vertices
   match their periodic partners to ~1e-9 of the cell.
2. Three Array modifiers (X, Y, Z) with vertex merging tile the cell.
   Because the mesh is exactly periodic, the merge is seamless.
"""

import time

import bpy
import numpy as np
from mathutils import Vector

from . import weierstrass


def _make_quad_mesh_object(name, verts, quads, normals, collection, smooth):
    """Build a directly sampled quad mesh with smooth vertex normals."""
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(
        [(float(v[0]), float(v[1]), float(v[2])) for v in verts],
        [],
        [(int(q[0]), int(q[1]), int(q[2]), int(q[3])) for q in quads],
    )
    mesh.update()

    if smooth:
        for poly in mesh.polygons:
            poly.use_smooth = True

    # Exact per-vertex normals from the Gauss map of the Weierstrass data.
    # The custom split-normals API differs across Blender versions; fail
    # soft so the mesh still loads with computed normals.
    try:
        mesh.normals_split_custom_set_from_vertices(
            [(float(n[0]), float(n[1]), float(n[2])) for n in normals]
        )
    except Exception:
        pass

    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    return obj


def _add_array_modifier(obj, name, count, offset_vec, merge_threshold):
    m = obj.modifiers.new(name=name, type='ARRAY')
    m.fit_type = 'FIXED_COUNT'
    m.count = count
    m.use_relative_offset = False
    m.use_constant_offset = True
    m.constant_offset_displace = offset_vec
    m.use_object_offset = False
    m.use_merge_vertices = True
    m.use_merge_vertices_cap = False
    m.merge_threshold = merge_threshold
    return m


class TPMS_OT_generate(bpy.types.Operator):
    """Generate one exact TPMS unit cell + Array modifiers for tiling"""
    bl_idname = "tpms.generate"
    bl_label = "Generate TPMS"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        props = context.scene.tpms_props
        cs = float(props.cell_scale)

        t0 = time.perf_counter()
        verts, quads, normals = weierstrass.build_unit_cell(
            tpms_type=props.tpms_type, cell_size=cs,
            res=int(props.quad_subdivisions),
            solver_resolution=int(props.solver_resolution),
            quadrature_order=int(props.quadrature_order))
        dt_build = time.perf_counter() - t0

        name = f"TPMS_{props.tpms_type.title()}"
        obj = _make_quad_mesh_object(
            name, verts, quads, normals, context.collection,
            props.smooth_shade)

        merge_thr = max(cs * 1e-5, 1e-8)
        _add_array_modifier(obj, "Array X", int(props.cells_x),
                            Vector((cs, 0.0, 0.0)), merge_thr)
        _add_array_modifier(obj, "Array Y", int(props.cells_y),
                            Vector((0.0, cs, 0.0)), merge_thr)
        _add_array_modifier(obj, "Array Z", int(props.cells_z),
                            Vector((0.0, 0.0, cs)), merge_thr)

        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        context.view_layer.objects.active = obj

        dt = time.perf_counter() - t0
        self.report(
            {'INFO'},
            f"{name}: {len(verts)}v / {len(quads)} quads, all on the "
            f"exact minimal surface (build {dt_build:.2f}s, total {dt:.2f}s). "
            f"Tiled {props.cells_x}x{props.cells_y}x{props.cells_z}."
        )
        return {'FINISHED'}


_classes = (TPMS_OT_generate,)


def register():
    for cls in _classes:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_classes):
        bpy.utils.unregister_class(cls)
