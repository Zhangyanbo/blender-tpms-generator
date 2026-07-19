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
     "48 genuine quadrilateral macro patches per cell"),
    ('SCHWARZ_P', "Schwarz Primitive",
     "Schwarz P surface (Im-3m): 24 genuine quadrilateral macro patches"),
    ('SCHWARZ_D', "Schwarz Diamond",
     "Schwarz D surface (Fd-3m): 96 genuine quadrilateral macro patches"),
]


class TPMSProperties(bpy.types.PropertyGroup):
    tpms_type: EnumProperty(
        name="Type",
        description="Which TPMS to generate (all from their exact "
                    "Enneper-Weierstrass parametrizations)",
        items=TPMS_TYPES,
        default='GYROID',
    )

    cell_scale: FloatProperty(
        name="Cell Scale",
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

    quad_subdivisions: IntProperty(
        name="Quad Subdivisions",
        description="Quads along each side of a genuine macro patch",
        default=2, min=1, soft_max=16,
    )

    solver_resolution: IntProperty(
        name="Harmonic Solver",
        description="Resolution of the cached harmonic reparameterization "
                    "shared by all three Bonnet associates",
        default=44, min=32, max=64,
    )

    quadrature_order: IntProperty(
        name="Quadrature Order",
        description="Gauss-Legendre order for exact Weierstrass integration",
        default=200, min=64, max=320,
    )

    smooth_shade: BoolProperty(
        name="Smooth Shading",
        description="Apply smooth shading with consistent vertex normals",
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
