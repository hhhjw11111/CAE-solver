"""Protect frozen geometry while permitting separately controlled z refinement."""
import sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generate_hv_surface_mesh import generate,audit_geometry,write_mesh


class SurfaceMesh(unittest.TestCase):
    def test_original_level_is_byte_identical_to_frozen_mesh(self):
        root=Path(__file__).resolve().parents[1]
        coords,elements,_=generate(0,0)
        with tempfile.TemporaryDirectory() as d:
            mesh=Path(d)/'mesh';write_mesh(mesh,coords,elements)
            self.assertEqual(mesh.read_bytes(),(root/'cases/F1_20260906/H1_L0.mesh').read_bytes())

    def test_thickness_refinement_preserves_contacts_and_inplane_coordinates(self):
        coarse,_,_=generate(0,0); fine,elements,m=generate(0,2)
        self.assertEqual({tuple(p[:2]) for p in coarse},{tuple(p[:2]) for p in fine})
        self.assertEqual(m['first_z_m'],3.5e-9)
        a=audit_geometry(fine,elements)
        self.assertGreater(a['minimum_jacobian_3x3x3'],0)
        self.assertEqual(a['anchor_counts'],[1,1,1])


if __name__=='__main__':unittest.main()
