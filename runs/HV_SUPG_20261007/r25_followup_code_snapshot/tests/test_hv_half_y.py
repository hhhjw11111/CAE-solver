"""Half-domain reduction must preserve polynomials and remove only rigid motion."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from generate_hv_surface_mesh import generate
from generate_hv_half_y import restrict_mesh,rigid_gauge,verify_reflection


class HalfY(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        xyz,conn,_=generate(0,0)
        cls.xyz=np.array(xyz);cls.conn=np.array(conn)-1

    def test_restriction_retains_every_element_coordinate(self):
        ids,eids,xyz,conn=restrict_mesh(self.xyz,self.conn)
        self.assertEqual(len(conn)*2,len(self.conn))
        np.testing.assert_array_equal(xyz[conn],self.xyz[self.conn[eids]])
        self.assertTrue(np.all(xyz[:,1]>=-1e-15))

    def test_crossing_elements_are_rejected(self):
        xyz=self.xyz.copy();xyz[self.conn[0,0],1]=-1e-6
        with self.assertRaises(ValueError):restrict_mesh(xyz,self.conn)

    def test_reflection_requires_coordinates_and_materials(self):
        material=np.ones(len(self.conn),dtype=int)
        self.assertEqual(verify_reflection(self.xyz,self.conn,material)['matched_nodes'],len(self.xyz))
        material[0]=2
        with self.assertRaises(ValueError):verify_reflection(self.xyz,self.conn,material)
        xyz=self.xyz.copy();xyz[0,1]=1e-8
        with self.assertRaises(ValueError):verify_reflection(xyz,self.conn,np.ones(len(self.conn),dtype=int))

    def test_arbitrary_rigid_motion_recovers_symmetric_strain_field(self):
        xyz=self.xyz; strain=np.diag([2e-5,-3e-5,1e-5]);u=xyz@strain
        shift=np.array([3e-10,-2e-10,4e-10]);rotation=np.array([4e-5,-3e-5,2e-5])
        original,_=rigid_gauge(xyz,u)
        corrected,gauge=rigid_gauge(xyz,u+shift+np.cross(rotation,xyz))
        np.testing.assert_allclose(corrected,original,atol=1e-23,rtol=1e-13)
        self.assertLess(gauge['maximum_symmetry_uy_m'],1e-23)


if __name__=='__main__':unittest.main()
