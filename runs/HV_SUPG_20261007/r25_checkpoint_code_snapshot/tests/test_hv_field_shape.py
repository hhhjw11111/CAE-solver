import sys,unittest,json,tempfile
from unittest.mock import patch
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from audit_hv_field_shape import difference,audit
from run_hv_sweep import sha
class FieldShape(unittest.TestCase):
    def test_local_outlier_is_not_hidden_by_rms(self):
        a=np.ones(1000);b=a.copy();b[5]=1.03
        result=difference(a,b)
        self.assertLess(result['rms_relative'],.02)
        self.assertFalse(result['passed'])
    def test_zero_and_nonfinite(self):
        self.assertTrue(difference(np.zeros(3),np.zeros(3))['passed'])
        with self.assertRaises(ValueError):difference(np.array([float('nan')]),np.zeros(1))
    def test_small_perturbation_is_not_normalized_by_reference(self):
        a=np.expm1(np.array([1e-4]));b=np.expm1(np.array([1.1e-4]))
        self.assertFalse(difference(a,b)['passed'])
    def test_wrong_mapping_source_mesh_rejected_before_field_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';target=root/'target';mapping=root/'mapping'
            for directory in [source,target,mapping]:directory.mkdir()
            exe=root/'solver';exe.write_text('fixture')
            for directory in [source,target]:
                mesh=directory/'mesh';mesh.write_text(directory.name)
                inp=directory/'input.nml';inp.write_text('fixture')
                (directory/'nodes.csv').write_text('fixture')
                manifest=dict(mesh=str(mesh),mesh_sha256=sha(mesh),executable=str(exe),executable_sha256=sha(exe))
                (directory/'manifest.json').write_text(json.dumps(manifest))
                (directory/'execution.json').write_text(json.dumps(dict(exit_code=0,input_sha256=sha(inp),executable_sha256=sha(exe))))
                (directory/'field_audit.json').write_text(json.dumps(dict(passed=True)))
            report=dict(source_nodes_sha256=sha(source/'nodes.csv'),target_mesh_sha256=sha(target/'mesh'),source_mesh_sha256='wrong')
            (mapping/'transfer_report.json').write_text(json.dumps(report))
            metrics=dict(voltage_V=.5,model='NT',heat_source='original',source='fixture',transport_discretization='weighted',domain_symmetry='Y_HALF')
            with patch('audit_hv_field_shape.classify',return_value=('SOLVE_FLUX_PASS_MESH_PENDING',metrics)):
                with self.assertRaisesRegex(ValueError,'Mapping source mesh mismatch'):
                    audit(mapping,source,target)
if __name__=='__main__':unittest.main()
