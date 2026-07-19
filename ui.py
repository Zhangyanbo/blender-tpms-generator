"""UI panel for the TPMS Generator. Lives in the 3D View's N-panel."""

import bpy


class TPMS_PT_panel(bpy.types.Panel):
    bl_label = "TPMS Generator"
    bl_idname = "TPMS_PT_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TPMS"

    def draw(self, context):
        layout = self.layout
        p = context.scene.tpms_props

        layout.prop(p, "tpms_type")

        box = layout.box()
        box.label(text="Lattice")
        box.prop(p, "cell_scale")

        box = layout.box()
        box.label(text="Tiling")
        row = box.row(align=True)
        row.prop(p, "cells_x")
        row.prop(p, "cells_y")
        row.prop(p, "cells_z")
        box.label(text="(also editable on the Array modifiers)",
                  icon='INFO')

        box = layout.box()
        box.label(text="Quality")
        box.prop(p, "quad_subdivisions")
        patches = {'GYROID': 48, 'SCHWARZ_P': 48, 'SCHWARZ_D': 192}
        npatch = patches.get(p.tpms_type, 48)
        subdivisions = int(p.quad_subdivisions)
        if p.tpms_type != 'GYROID':
            subdivisions = max(2, subdivisions)
        box.label(text=f"{npatch * subdivisions ** 2} quads per cell",
                  icon='MESH_GRID')
        if p.tpms_type == 'GYROID':
            box.prop(p, "solver_resolution")
            box.prop(p, "quadrature_order")
        box.prop(p, "smooth_shade")

        layout.separator()
        col = layout.column(align=True)
        col.scale_y = 1.4
        col.operator("tpms.generate", icon='MESH_ICOSPHERE')


_classes = (TPMS_PT_panel,)


def register():
    for cls in _classes:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_classes):
        bpy.utils.unregister_class(cls)
