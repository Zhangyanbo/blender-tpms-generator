import unittest

import numpy as np

import gyroid_macro
import isolevel
import weierstrass


TYPES = ('GYROID', 'SCHWARZ_P', 'SCHWARZ_D')
PATCHES = {'GYROID': 48, 'SCHWARZ_P': 24, 'SCHWARZ_D': 96}


class IsolevelTests(unittest.TestCase):
    def test_exact_meshes_hug_their_fields(self):
        for name in TYPES:
            verts, _, normals = weierstrass.build_unit_cell(name, 1.0, 2)
            values = isolevel.field(name, verts)
            self.assertLess(np.abs(values).mean(), 0.05, name)
            g = isolevel.gradient(name, verts)
            g /= np.linalg.norm(g, axis=1, keepdims=True)
            alignment = np.einsum('ij,ij->i', g, normals)
            self.assertGreater(np.abs(alignment).min(), 0.9, name)

    def test_gradient_matches_finite_differences(self):
        rng = np.random.default_rng(7)
        points = rng.uniform(0.0, 1.0, size=(64, 3))
        eps = 1.0e-6
        for name in TYPES:
            numeric = np.stack([
                (isolevel.field(name, points + eps * e)
                 - isolevel.field(name, points - eps * e)) / (2.0 * eps)
                for e in np.eye(3)], axis=1)
            analytic = isolevel.gradient(name, points)
            self.assertLess(np.abs(numeric - analytic).max(), 1.0e-5, name)

    def test_zero_level_is_untouched(self):
        for name in TYPES:
            base = weierstrass.build_unit_cell(name, 1.0, 2)
            same = weierstrass.build_unit_cell(name, 1.0, 2, iso_level=0.0)
            for a, b in zip(base, same):
                self.assertTrue(np.array_equal(a, b), name)

    def test_isolevel_hits_the_level_and_keeps_topology(self):
        for name in TYPES:
            base_verts, base_quads, _ = weierstrass.build_unit_cell(
                name, 1.0, 3)
            reference = gyroid_macro.validate_unit_cell(
                base_verts, base_quads)
            for level in (-0.6, 0.35, 0.85):
                verts, quads, normals = weierstrass.build_unit_cell(
                    name, 1.0, 3, iso_level=level)
                self.assertLess(
                    np.abs(isolevel.field(name, verts) - level).max(),
                    1.0e-8, f"{name} t={level}")
                self.assertTrue(np.array_equal(quads, base_quads))
                self.assertEqual(len(verts), len(base_verts))
                metrics = gyroid_macro.validate_unit_cell(verts, quads)
                for key in ('non_quad_faces', 'duplicate_faces',
                            'degenerate_quads', 'folded_quads',
                            'nonmanifold_edges'):
                    self.assertEqual(metrics[key], 0, f"{name} t={level} {key}")
                self.assertTrue(
                    metrics['quotient_all_edges_incident_to_two_quads'])
                self.assertEqual(metrics['quotient_euler_characteristic'],
                                 reference['quotient_euler_characteristic'],
                                 f"{name} t={level}")
                lengths = np.linalg.norm(normals, axis=1)
                self.assertLess(np.abs(lengths - 1.0).max(), 1.0e-9)

    def test_isolevel_respects_cell_scale(self):
        verts1, _, _ = weierstrass.build_unit_cell(
            'GYROID', 1.0, 2, iso_level=0.5)
        verts3, _, _ = weierstrass.build_unit_cell(
            'GYROID', 3.0, 2, iso_level=0.5)
        self.assertLess(np.abs(verts3 - 3.0 * verts1).max(), 1.0e-9)

    def test_critical_levels_are_rejected(self):
        for name, bad in (('GYROID', 1.42), ('SCHWARZ_P', 1.0),
                          ('SCHWARZ_D', -1.0)):
            with self.assertRaises(ValueError):
                weierstrass.build_unit_cell(name, 1.0, 2, iso_level=bad)


if __name__ == '__main__':
    unittest.main()
