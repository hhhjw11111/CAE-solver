"""Protect frozen geometry while permitting separately controlled z refinement."""
import sys, tempfile, unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generate_hv_surface_mesh import generate,audit_geometry,write_mesh,refine_surface_intervals


class SurfaceMesh(unittest.TestCase):
    def test_original_level_is_byte_identical_to_frozen_mesh(self):
        root=Path(__file__).resolve().parents[1]
        for level in range(3):
            coords,elements,_=generate(level,level)
            with tempfile.TemporaryDirectory() as d:
                mesh=Path(d)/'mesh';write_mesh(mesh,coords,elements)
                self.assertEqual(mesh.read_bytes(),(root/f'cases/F1_20260906/H1_L{level}.mesh').read_bytes())

    def test_nested_refinement_splits_entire_near_surface_intervals(self):
        coarse,_,cm=generate(0,2)
        fine,elements,fm=generate(0,2,z_subdivisions=2)
        old=np.array(cm['z_corners_m']);new=np.array(fm['z_corners_m'])
        self.assertEqual({tuple(p[:2]) for p in coarse},{tuple(p[:2]) for p in fine})
        for a,b in zip(old,old[1:]):
            inside=new[(new>=a-1e-22)&(new<=b+1e-22)]
            distance=min(abs(a+5e-7),abs(a-5e-7),abs(b+5e-7),abs(b-5e-7))
            self.assertEqual(len(inside),3 if distance<1e-7 else 2)
            self.assertLessEqual(np.max(np.diff(inside)),(b-a)/(2 if distance<1e-7 else 1)+1e-21)
        self.assertGreater(audit_geometry(fine,elements)['minimum_jacobian_3x3x3'],0)
        for bad in [0,5]:
            with self.assertRaises(ValueError):refine_surface_intervals(list(old),bad,1e-7)

    def test_thickness_refinement_preserves_contacts_and_inplane_coordinates(self):
        coarse,_,_=generate(0,0); fine,elements,m=generate(0,2)
        self.assertEqual({tuple(p[:2]) for p in coarse},{tuple(p[:2]) for p in fine})
        self.assertEqual(m['first_z_m'],3.5e-9)
        a=audit_geometry(fine,elements)
        self.assertGreater(a['minimum_jacobian_3x3x3'],0)
        self.assertEqual(a['anchor_counts'],[1,1,1])

    def test_independent_angular_refinement_preserves_layers(self):
        coarse,ce,cm=generate(0,0);fine,fe,fm=generate(0,0,32)
        self.assertEqual(len(fe),2*len(ce))
        self.assertEqual(fm['z_corners_m'],cm['z_corners_m'])
        self.assertEqual(fm['first_xy_m'],cm['first_xy_m'])
        cset={tuple(np.round(p,18)) for p in coarse};fset={tuple(np.round(p,18)) for p in fine}
        self.assertTrue(cset<=fset)
        self.assertGreater(audit_geometry(fine,fe)['minimum_jacobian_3x3x3'],0)
        with self.assertRaises(ValueError):generate(0,0,30)


if __name__=='__main__':unittest.main()
