"""Guard model identity and failed/incomplete prerequisites; no native solves."""
from pathlib import Path
import json
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_je_controlled_trace_batch_20261021 as batch
import run_je_trace_native_20261021 as native
import run_je_controlled_trace_batch_20261021_v2 as resume


class ControlledNumericsTests(unittest.TestCase):
    def test_local_comparison_mapping_is_voltage_specific_and_rejects_wrong_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            targets = [root/('A48_R22_Z18_S1_V'+v) for v in ('060','080')]
            source = root/'source'
            manifest = dict(mesh='synthetic',mesh_sha256='target')
            for target in targets:
                directory = root/'transfers'/('mesh_final_pair_'+target.name)
                directory.mkdir(parents=True)
                (directory/'transfer_report.json').write_text(json.dumps(dict(
                    boundary_extension_nodes=0,source_nodes_sha256='wrong-seed',
                    source_mesh_sha256='source',target_mesh_sha256='target')),encoding='utf8')
                with patch.object(resume,'OUT',root), \
                     patch.object(resume.b,'checked_half',return_value=(manifest,{})), \
                     patch.object(resume.r,'sha',return_value='correct-seed'), \
                     patch.object(resume,'field_audit') as compare:
                    with self.assertRaisesRegex(RuntimeError,'mapping identity'):
                        resume.compare(source,target,'mesh_final_pair')
                    compare.assert_not_called()
            self.assertEqual(len(list((root/'transfers').iterdir())),2)

    def test_resume_rejects_wrong_existing_supg_metric_before_reuse(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            root.mkdir(exist_ok=True)
            with patch.object(resume,'point',return_value=root), \
                 patch.object(resume.b,'checked_half',return_value=({},dict(supg_directional_tau='T'))), \
                 patch.object(resume.b.native,'run_native') as launch:
                with self.assertRaisesRegex(RuntimeError,'not this numerical variant'):
                    resume.run_point(0,6,1.,{})
                launch.assert_not_called()

    def test_missing_benchmark_blocks_launch_and_output_creation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with patch.object(batch,'ROOT',root), patch.object(batch,'OUT',root/'runs'), \
                 patch.object(batch.r,'validate_protocol',return_value={}), \
                 patch.object(batch.native,'run_native') as launch:
                with self.assertRaises(FileNotFoundError):
                    batch.execute()
                launch.assert_not_called()
                self.assertFalse((root/'runs').exists())

    def test_false_benchmark_pass_cannot_use_different_executable(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root/'native_trace_verification_attempt02'
            folder.mkdir()
            (folder/'acceptance.json').write_text(json.dumps(dict(passed=True,exit_code=0,
                original_core_unchanged=True,original_thresholds_unchanged=True,
                production_executable_sha256='different-binary')),encoding='utf8')
            with patch.object(batch,'ROOT',root), patch.object(batch.r,'validate_protocol',return_value={}):
                with self.assertRaisesRegex(RuntimeError,'prerequisite failed'):
                    batch.validate_prerequisites()

    def test_changed_benchmark_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root/'native_trace_verification_attempt02'
            folder.mkdir()
            source = root/'frozen.f90'
            source.write_text('changed',encoding='ascii')
            (folder/'acceptance.json').write_text(json.dumps(dict(passed=True,exit_code=0,
                original_core_unchanged=True,original_thresholds_unchanged=True,
                production_executable_sha256=batch.r.EXPECTED_EXE,
                verification_source_hashes={str(source):'wrong-hash'},copied_source_hashes={})),encoding='utf8')
            with patch.object(batch,'ROOT',root), patch.object(batch.r,'validate_protocol',return_value={}):
                with self.assertRaisesRegex(RuntimeError,'source changed'):
                    batch.validate_prerequisites()

    def test_trace_adapter_preserves_heat_load_and_tolerances_before_launch(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = root/'templates'
            base.mkdir()
            shutil.copy2(native.BASE/'input_V060.nml',base/'input_V060.nml')
            mesh, seed = root/'mesh', root/'seed'
            mesh.write_text('synthetic identity only',encoding='ascii')
            seed.write_text('synthetic identity only',encoding='ascii')
            out = root/'run'
            with patch.object(native,'BASE',base), patch.object(native,'preflight'), \
                 patch.object(native,'short',side_effect=lambda p:Path(p).as_posix()), \
                 patch.object(native,'runtime_env',return_value={}), \
                 patch.object(native.shutil,'disk_usage',return_value=shutil._ntuple_diskusage(100*1024**3,0,100*1024**3)), \
                 patch.object(native.subprocess,'Popen',side_effect=OSError('injected prelaunch failure')):
                with self.assertRaisesRegex(OSError,'injected'):
                    native.run_native(dict(name='test',voltage=6,supg_scale=2.,counts=[32,18,14]),
                                      mesh,seed,out,{})
            text = (out/'input.nml').read_text(encoding='ascii')
            for literal in ('heat_source_current_work=.true.','residual_tolerance=1e-9,update_tolerance=1e-9',
                            'target_joule_factor=1.0','restart_joule_factor=1.0',
                            'supg_directional_tau=.false.','supg_scale=2.0'):
                self.assertIn(literal,text)
            record = json.loads((out/'manifest.json').read_text(encoding='utf8'))
            self.assertEqual(record['runtime_env']['TPS_PARDISO_SPLIT_HOLES'],'1')
            self.assertEqual(record['executable_sha256'],batch.r.EXPECTED_EXE)
            self.assertFalse((out/'execution.json').exists())


if __name__ == '__main__':
    unittest.main()
