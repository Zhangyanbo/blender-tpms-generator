import unittest

import numpy as np

import bonnet_macro
import gyroid_macro
import weierstrass


class BonnetMacroTests(unittest.TestCase):
    def _surface_data(self, surface_name):
        surface = weierstrass.SURFACES[surface_name]
        reference, _ = weierstrass._fundamental_patch(surface, 2)
        transforms = bonnet_macro._operation_transforms(surface, reference)
        return surface, reference, transforms

    def test_circular_edges_pair_without_hard_coded_indices(self):
        expected = {'SCHWARZ_P': 24, 'SCHWARZ_D': 96}
        for surface_name, patch_count in expected.items():
            surface, reference, transforms = self._surface_data(surface_name)
            pairs = bonnet_macro._macro_pairs(
                surface_name, surface, transforms)
            self.assertEqual(len(pairs), patch_count)
            self.assertEqual(
                sorted(index for pair in pairs for index in pair[:2]),
                list(range(2 * patch_count)),
            )

    def test_required_face_counts_are_all_quads(self):
        patch_counts = {'SCHWARZ_P': 24, 'SCHWARZ_D': 96}
        for surface_name, patch_count in patch_counts.items():
            for subdivisions in (1, 2, 4, 8):
                vertices, faces, normals = weierstrass.build_unit_cell(
                    surface_name,
                    res=subdivisions,
                )
                self.assertEqual(
                    faces.shape,
                    (patch_count * subdivisions * subdivisions, 4),
                )
                self.assertEqual(len(vertices), len(normals))
                self.assertTrue(np.isfinite(vertices).all())

    def test_reference_n2_periodic_topology(self):
        expected = {
            'SCHWARZ_P': (140, 96, -4, 3),
            'SCHWARZ_D': (459, 384, -16, 9),
        }
        for surface_name, values in expected.items():
            vertices, faces, _ = weierstrass.build_unit_cell(
                surface_name,
                res=2,
            )
            metrics = gyroid_macro.validate_unit_cell(vertices, faces)
            vertex_count, face_count, euler, genus = values
            self.assertEqual(metrics['vertices'], vertex_count)
            self.assertEqual(metrics['quad_faces'], face_count)
            self.assertEqual(metrics['non_quad_faces'], 0)
            self.assertEqual(metrics['duplicate_faces'], 0)
            self.assertEqual(metrics['degenerate_quads'], 0)
            self.assertEqual(metrics['folded_quads'], 0)
            self.assertEqual(metrics['nonmanifold_edges'], 0)
            self.assertTrue(metrics['all_periodic_boundary_edges_paired'])
            self.assertTrue(metrics['quotient_all_edges_incident_to_two_quads'])
            self.assertEqual(metrics['quotient_euler_characteristic'], euler)
            self.assertEqual(metrics['periodic_quotient_genus'], genus)

    def test_analytic_branches_match_on_internal_circle(self):
        phi = np.linspace(0.05, 0.5 * np.pi - 0.05, 17)
        omega = gyroid_macro._rho_max(phi) * np.exp(1j * phi)
        for surface_name in ('SCHWARZ_P', 'SCHWARZ_D'):
            surface, _, transforms = self._surface_data(surface_name)
            pair = bonnet_macro._macro_pairs(
                surface_name, surface, transforms)[0]
            lower, upper, shift = pair
            lower_points = bonnet_macro._apply(
                bonnet_macro._surface_points(omega, surface),
                transforms[lower])
            upper_omega = 1j * np.conj(omega)
            upper_points = bonnet_macro._apply(
                bonnet_macro._surface_points(upper_omega, surface),
                transforms[upper]) + shift
            error = np.max(np.linalg.norm(lower_points - upper_points, axis=1))
            self.assertLess(error, 1.0e-9)


if __name__ == '__main__':
    unittest.main()
