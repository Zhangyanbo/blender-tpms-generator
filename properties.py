"""Scene-level properties for the Gyroid Generator add-on."""

import bpy
from bpy.props import (
    BoolProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
)


class TPMSProperties(bpy.types.PropertyGroup):
    cell_size: FloatProperty(
        name="Cell Size",
        description="World-unit edge length of one cubic Gyroid unit cell",
        default=1.0, min=0.001, soft_max=10.0,
        unit='LENGTH',
    )

    cells_x: IntProperty(
        name="Cells X",
        description="Number of unit cells along X (Array modifier count). "
                    "Can be edited live on the Array modifier after "
                    "generation",
        default=4, min=1, soft_max=64,
    )

    cells_y: IntProperty(
        name="Cells Y",
        description="Number of unit cells along Y",
        default=4, min=1, soft_max=64,
    )

    cells_z: IntProperty(
        name="Cells Z",
        description="Number of unit cells along Z",
        default=4, min=1, soft_max=64,
    )

    resolution: IntProperty(
        name="Resolution",
        description="Quads per fundamental-patch edge. The unit cell is "
                    "assembled from 96 patches, so it has 96 x res^2 quads "
                    "(res 8 = 6144 quads). Every vertex lies on the exact "
                    "minimal surface, so even low resolutions are accurate",
        default=8, min=2, soft_max=32,
    )

    smooth_shade: BoolProperty(
        name="Smooth Shading",
        description="Apply smooth shading with exact analytic normals",
        default=True,
    )


_classes = (TPMSProperties,)


def register():
    for cls in _classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.tpms_props = PointerProperty(type=TPMSProperties)


def unregister():
    del bpy.types.Scene.tpms_props
    for cls in reversed(_classes):
        bpy.utils.unregister_class(cls)
