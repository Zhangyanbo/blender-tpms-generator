import unittest

import numpy as np

import gyroid_macro


class GyroidMacroTests(unittest.TestCase):
    def test_macro_patch_pairing_is_discovered(self):
        counts = gyroid_macro.macro_topology_counts()
        self.assertEqual(counts['fundamental_triangles'], 96)
        self.assertEqual(counts['circular_edge_groups'], 48)
        self.assertEqual(counts['triangles_per_group'], (2,) * 48)
        self.assertEqual(counts['genuine_quad_patches'], 48)

    def test_required_face_counts_are_all_quads(self):
        for subdivisions in (1, 2, 4, 8):
            vertices, faces, normals = gyroid_macro.build_unit_cell(
                quad_subdivisions=subdivisions,
            )
            self.assertEqual(faces.shape,
                             (48 * subdivisions * subdivisions, 4))
            self.assertEqual(len(vertices), len(normals))
            self.assertTrue(np.isfinite(vertices).all())

    def test_reference_n2_periodic_topology(self):
        vertices, faces, _ = gyroid_macro.build_unit_cell(
            quad_subdivisions=2,
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

    def test_analytic_square_map_has_positive_jacobian(self):
        axis = np.linspace(0.0, 1.0, 301)
        u, v = np.meshgrid(axis, axis, indexing='ij')
        unfolded = gyroid_macro._macro_square_domain(u, v)
        step = 2.0 / 300.0
        du = (unfolded[2:, 1:-1] - unfolded[:-2, 1:-1]) / step
        dv = (unfolded[1:-1, 2:] - unfolded[1:-1, :-2]) / step
        jacobian = du.real * dv.imag - du.imag * dv.real
        self.assertGreater(np.min(jacobian), 0.2)

    def test_analytic_branches_have_matching_tangents(self):
        pair = gyroid_macro._macro_pairs()[0]
        angles = []
        for v in np.linspace(0.05, 0.95, 31):
            low, high = 0.0, 1.0
            for _ in range(60):
                u = 0.5 * (low + high)
                zeta = gyroid_macro._macro_square_domain(u, v)
                if abs(zeta - gyroid_macro._ARC_CENTER) > gyroid_macro._ARC_RADIUS:
                    high = u
                else:
                    low = u
            u = 0.5 * (low + high)
            step = 1.0e-5
            points = gyroid_macro._macro_patch_points(
                np.array([u - step, u, u + step]),
                np.array([v, v, v]), pair)
            left = points[1] - points[0]
            right = points[2] - points[1]
            cosine = np.dot(left, right) / (
                np.linalg.norm(left) * np.linalg.norm(right))
            angles.append(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
        self.assertLess(max(angles), 0.01)


if __name__ == '__main__':
    unittest.main()
