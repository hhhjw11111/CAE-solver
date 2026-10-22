"""Fail-closed orchestration and half-domain power accounting, no native solves."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import continue_je_half_mesh_supg_20261020 as batch


class HalfBatchTests(unittest.TestCase):
    def test_missing_equivalence_cannot_start_native_or_create_mesh_root(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with patch.object(batch.eq,'OUTPUT',root/'missing'), \
                 patch.object(batch,'OUTPUT',root/'mesh'), \
                 patch.object(batch.r,'validate_protocol',return_value={}), \
                 patch.object(batch.native,'run_native') as native:
                with self.assertRaises(FileNotFoundError):batch.execute()
                native.assert_not_called()
                self.assertFalse((root/'mesh').exists())

    def test_duplicate_equivalence_runs_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            record=dict(all_four_equivalence_passed=True,
                events=[dict(run=str(root/'A32_R18_Z14_V060_unperturbed'))]*4)
            (root/'acceptance.json').write_text(json.dumps(record))
            with patch.object(batch.eq,'OUTPUT',root):
                with self.assertRaisesRegex(RuntimeError,'Duplicate'):
                    batch.equivalence_gate()

    def test_stopped_equivalence_does_not_wait_or_advance(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'stopped.json').write_text('{}')
            with patch.object(batch.eq,'OUTPUT',root),patch.object(batch.time,'sleep') as sleep:
                with self.assertRaisesRegex(RuntimeError,'controller stopped'):
                    batch.wait_for_equivalence(30)
                sleep.assert_not_called()

    def power_fixture(self,temp,reference=2.8):
        folder=Path(temp);(folder/'power_audit').mkdir()
        b=dict(JE_W=4.,all_exterior_power_W=3.,
            boundary_sets={'y-':dict(power_W=.2)})
        met=dict(voltage_V=6.,physical_exterior_terminal_power_W=reference)
        return folder,b,met

    def test_cut_flux_is_recorded_and_full_integrals_are_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            folder,b,met=self.power_fixture(temp)
            with patch.object(batch,'checked_half',return_value=({},met)), \
                 patch.object(batch.r,'read',return_value=b), \
                 patch.object(batch.r,'sha',return_value='synthetic-unit-test-only'):
                result=batch.physical_power(folder,folder/'result.json')
            self.assertEqual(result['raw_half_symmetry_cut_power_W'],.2)
            self.assertAlmostEqual(result['raw_half_physical_exterior_power_W'],2.8)
            self.assertAlmostEqual(result['physical_exterior_gap_W'],1.2)
            self.assertEqual(result['reconstructed_full_JE_W'],8.)
            self.assertAlmostEqual(result['reconstructed_full_physical_exterior_power_W'],5.6)
            self.assertIsNone(result['power_gap_acceptance_threshold'])

    def test_wrong_native_boundary_definition_fails_before_report(self):
        with tempfile.TemporaryDirectory() as temp:
            folder,b,met=self.power_fixture(temp,reference=3.)
            with patch.object(batch,'checked_half',return_value=({},met)), \
                 patch.object(batch.r,'read',return_value=b):
                with self.assertRaisesRegex(RuntimeError,'does not match'):
                    batch.physical_power(folder,folder/'result.json')
            self.assertFalse((folder/'result.json').exists())


if __name__=='__main__':unittest.main()
