"""Scene-level properties for the TPMS Generator add-on."""

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
)


TPMS_TYPES = [
    ('GYROID', "Gyroid",
     "Schoen G surface (Ia-3d): exact Enneper-Weierstrass parametrization, "
     "96 patches per cell"),
    ('SCHWARZ_P', "Schwarz Primitive",
     "Schwarz P surface (Im-3m): exact parametrization, 48 patches per cell"),
    ('SCHWARZ_D', "Schwarz Diamond",
     "Schwarz D surface (Fd-3m): exact parametrization, 192 patches per cell"),
]


class TPMSProperties(bpy.types.PropertyGroup):
    tpms_type: EnumProperty(
        name="Type",
        description="Which TPMS to generate (all from their exact "
                    "Enneper-Weierstrass parametrizations)",
        items=TPMS_TYPES,
        default='GYROID',
    )

    cell_size: FloatProperty(
        name="Cell Size",
        description="World-unit edge length of one cubic TPMS unit cell",
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
        description="Quads per fundamental-patch edge. The unit cell has "
                    "(patches x res^2) quads with 96/48/192 patches for "
                    "Gyroid/P/D. Every vertex lies on the exact minimal "
                    "surface, so even low resolutions are accurate",
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
