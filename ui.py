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

        col = layout.column(align=True)
        col.prop(p, "target")
        col.prop(p, "tpms_type")

        box = layout.box()
        box.label(text="Lattice")
        box.prop(p, "cell_size")
        box.prop(p, "iso_value", slider=True)
        row = box.row(align=True)
        row.prop(p, "origin", text="")
        box.prop(p, "rotation", text="Rotation")

        box = layout.box()
        box.label(text="Solid")
        box.prop(p, "solid_mode", expand=True)
        sub = box.column()
        sub.enabled = (p.solid_mode == 'SHELL')
        sub.prop(p, "thickness", slider=True)

        box = layout.box()
        box.label(text="Sampling")
        box.prop(p, "resolution")
        box.prop(p, "project_iters")
        box.prop(p, "padding")

        box = layout.box()
        box.label(text="Output")
        box.prop(p, "clip_to_target")
        sub = box.column()
        sub.enabled = p.clip_to_target
        sub.prop(p, "robust_clip")
        sub.prop(p, "sdf_subsample")
        sub.prop(p, "boundary_inset")
        box.prop(p, "min_component_faces")
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
