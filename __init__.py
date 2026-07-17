"""TPMS Generator — Blender 4.2+ extension.

Generates exact triply periodic minimal surfaces (Gyroid, Schwarz P,
Schwarz D) from their Enneper-Weierstrass parametrizations (Gandy,
Cvijovic, Mackay & Klinowski, Chem. Phys. Lett. 314 (1999) 543;
321 (2000) 363; 322 (2000) 579) as clean all-quad meshes: one cubic
unit cell each, tiled with Array modifiers.
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
