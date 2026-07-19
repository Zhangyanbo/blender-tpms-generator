import unittest

import numpy as np

import gyroid_macro


class GyroidMacroTests(unittest.TestCase):
    def test_macro_patch_pairing_is_discovered(self):
        counts = gyroid_macro.macro_topology_counts(quadrature_order=80)
        self.assertEqual(counts['fundamental_triangles'], 96)
        self.assertEqual(counts['circular_edge_groups'], 48)
        self.assertEqual(counts['triangles_per_group'], (2,) * 48)
        self.assertEqual(counts['genuine_quad_patches'], 48)

    def test_required_face_counts_are_all_quads(self):
        for subdivisions in (1, 2, 4, 8):
            vertices, faces, normals = gyroid_macro.build_unit_cell(
                quad_subdivisions=subdivisions,
                solver_resolution=32,
                quadrature_order=80,
            )
            self.assertEqual(faces.shape,
                             (48 * subdivisions * subdivisions, 4))
            self.assertEqual(len(vertices), len(normals))
            self.assertTrue(np.isfinite(vertices).all())

    def test_reference_n2_periodic_topology(self):
        vertices, faces, _ = gyroid_macro.build_unit_cell(
            quad_subdivisions=2,
            solver_resolution=44,
            quadrature_order=200,
        )
        metrics = gyroid_macro.validate_unit_cell(vertices, faces)
        self.assertEqual(metrics, {
            'vertices': 237,
            'quad_faces': 192,
            'non_quad_faces': 0,
            'duplicate_faces': 0,
            'degenerate_quads': 0,
            'folded_quads': 0,
            'nonmanifold_edges': 0,
            'all_periodic_boundary_edges_paired': True,
            'quotient_all_edges_incident_to_two_quads': True,
            'quotient_euler_characteristic': -8,
            'periodic_quotient_genus': 5,
        })

    def test_hidden_diagonal_matches_to_integration_accuracy(self):
        pair = gyroid_macro._macro_pairs(200)[0]
        parameter = np.linspace(0.05, 0.95, 11)
        epsilon = 1.0e-10
        below = gyroid_macro._piecewise_points(
            parameter, parameter - epsilon, pair, 200)
        above = gyroid_macro._piecewise_points(
            parameter, parameter + epsilon, pair, 200)
        error = np.max(np.linalg.norm(below - above, axis=1))
        self.assertLess(error, 1.0e-9)


if __name__ == '__main__':
    unittest.main()
