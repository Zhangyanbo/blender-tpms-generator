"""Gyroid Generator — Blender 4.2+ extension.

Generates the exact Gyroid triply periodic minimal surface from its
Enneper-Weierstrass parametrization (Gandy & Klinowski, Chem. Phys. Lett.
321 (2000) 363-371) as a clean all-quad mesh: one cubic unit cell, tiled
with Array modifiers.
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
