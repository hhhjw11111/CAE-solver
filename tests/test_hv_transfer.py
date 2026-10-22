import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transfer_hv_initial_guess import shapes,transfer,NATURAL

class Transfer(unittest.TestCase):
    def test_hex20_nodal_identity(self):
        np.testing.assert_allclose(shapes(NATURAL),np.eye(20),atol=1e-15)
    def test_curved_linear_field_and_nodal_identity(self):
        c=(NATURAL+1)*.5e-6;c[:,0]+=NATURAL[:,1]*NATURAL[:,2]*.13e-6
        c[:,1]+=NATURAL[:,0]*NATURAL[:,2]*.08e-6
        vals=np.column_stack([np.ones(20),c])
        rng=np.random.default_rng(812);target=shapes(rng.uniform(-.95,.95,(30,3)))@c
        result,r=transfer(c,np.arange(20)[None,:],vals,target)
        np.testing.assert_allclose(result,np.column_stack([np.ones(30),target]),atol=1e-15)
        result,r=transfer(c,np.arange(20)[None,:],vals,c)
        np.testing.assert_allclose(result,vals,atol=1e-15)
    def test_large_extrapolation_rejected(self):
        c=(NATURAL+1)*.5e-6
        with self.assertRaises(ValueError):transfer(c,np.arange(20)[None,:],np.ones((20,1)),np.array([[2e-6,0,0]]))

if __name__=='__main__':unittest.main()
