"""Scene-level properties for the TPMS Generator add-on."""

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    PointerProperty,
)


TPMS_TYPES = [
    ('GYROID',         "Gyroid",            "Schoen Gyroid: sin·cos sum"),
    ('SCHWARZ_P',      "Schwarz Primitive", "Schwarz P: cos x + cos y + cos z"),
    ('SCHWARZ_D',      "Schwarz Diamond",   "Schwarz D (diamond) surface"),
    ('SCHOEN_IWP',     "Schoen IWP",        "Schoen I-WP surface"),
    ('FISCHER_KOCH_S', "Fischer-Koch S",    "Fischer-Koch S surface"),
]

SOLID_MODES = [
    ('SHELL',  "Shell",  "Thicken the iso-surface into a closed shell of uniform wall thickness"),
    ('VOLUME', "Volume", "Use one phase (f < c) of the TPMS as a solid volume"),
]


def _target_poll(self, obj):
    return obj is not None and obj.type == 'MESH'


class TPMSProperties(bpy.types.PropertyGroup):
    target: PointerProperty(
        name="Target Mesh",
        description="Mesh whose enclosed volume will be filled with TPMS",
        type=bpy.types.Object,
        poll=_target_poll,
    )

    tpms_type: EnumProperty(
        name="Type",
        description="Which TPMS family to generate",
        items=TPMS_TYPES,
        default='GYROID',
    )

    cell_size: FloatProperty(
        name="Cell Size",
        description="World-unit length of one TPMS period (2π in lattice space)",
        default=0.5, min=0.001, soft_max=10.0,
        unit='LENGTH',
    )

    iso_value: FloatProperty(
        name="Iso Value",
        description="Level-set constant c. 0 ≈ equal phase volumes. Shifts the volume fraction",
        default=0.0, soft_min=-1.2, soft_max=1.2,
    )

    solid_mode: EnumProperty(
        name="Solid Mode",
        items=SOLID_MODES,
        default='SHELL',
    )

    thickness: FloatProperty(
        name="Wall Thickness",
        description="Wall thickness for Shell mode (in scalar-field units, ≈ relative to range)",
        default=0.3, min=0.01, soft_max=2.0,
    )

    resolution: IntProperty(
        name="Resolution / Cell",
        description="Voxels per TPMS unit cell along each axis (quality vs. memory)",
        default=24, min=6, soft_max=80,
    )

    clip_to_target: BoolProperty(
        name="Clip to Target",
        description="Clip the TPMS to the target's interior using the target's "
                    "signed-distance field (robust to non-closed meshes; no "
                    "Boolean modifier needed)",
        default=True,
    )

    robust_clip: BoolProperty(
        name="Robust In/Out Test",
        description="Use 3-ray majority-vote parity for the inside/outside "
                    "test. Slower but mandatory for concave / non-closed "
                    "targets (Suzanne, anything with eye sockets / open neck)",
        default=True,
    )

    sdf_subsample: IntProperty(
        name="SDF Subsample",
        description="Compute the target SDF on every Nth voxel and trilinear-"
                    "interpolate the rest. 4 is a good balance; raise it if "
                    "generation is slow, drop to 1 for crisp boundary",
        default=4, min=1, soft_max=8,
    )

    boundary_inset: FloatProperty(
        name="Boundary Inset",
        description="Shrink the target by this much (in voxels) before "
                    "clipping. Use 0.5-1.5 to hide any residual fuzz at the "
                    "TPMS-target interface",
        default=0.0, min=0.0, soft_max=4.0,
    )

    min_component_faces: IntProperty(
        name="Min Component Faces",
        description="Drop isolated mesh fragments smaller than this many "
                    "quads. The main TPMS body has thousands of quads, so a "
                    "threshold of 50-500 cleanly removes boundary specks. "
                    "Set 0 to disable",
        default=200, min=0, soft_max=2000,
    )

    padding: FloatProperty(
        name="BBox Padding",
        description="Extra space added around the target's bounding box before sampling",
        default=0.0, min=0.0, soft_max=2.0,
        unit='LENGTH',
    )

    rotation: FloatVectorProperty(
        name="Lattice Rotation",
        description="Rotate the TPMS lattice relative to the target",
        subtype='EULER', size=3, default=(0.0, 0.0, 0.0),
    )

    origin: FloatVectorProperty(
        name="Lattice Origin",
        description="Translate (phase-shift) the TPMS lattice",
        subtype='TRANSLATION', size=3, default=(0.0, 0.0, 0.0),
        unit='LENGTH',
    )

    smooth_shade: BoolProperty(
        name="Smooth Shading",
        description="Apply smooth shading to the generated mesh",
        default=True,
    )

    project_iters: IntProperty(
        name="Surface Snap Iters",
        description="Number of Newton iterations that project each vertex "
                    "onto the analytic iso-surface. 0 = raw Surface Nets",
        default=3, min=0, soft_max=6,
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
