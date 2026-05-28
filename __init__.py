"""TPMS Generator — Blender 4.2+ extension.

Generates Triply Periodic Minimal Surfaces (Gyroid, Schwarz P/D, Schoen IWP,
Fischer-Koch S) as solid meshes inside a user-selected target mesh.
"""

from . import properties, operators, ui


def register():
    properties.register()
    operators.register()
    ui.register()


def unregister():
    ui.unregister()
    operators.unregister()
    properties.unregister()
