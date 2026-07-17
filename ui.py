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
        box.prop(p, "cell_size")

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
        box.prop(p, "resolution")
        patches = {'GYROID': 96, 'SCHWARZ_P': 48, 'SCHWARZ_D': 192}
        npatch = patches.get(p.tpms_type, 96)
        box.label(text=f"~{npatch * int(p.resolution) ** 2} quads per cell",
                  icon='MESH_GRID')
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
