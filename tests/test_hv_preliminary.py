"""Preliminary delivery cannot hide failed flux checks or missing voltages."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from summarize_hv_preliminary import preliminary_status,summarize,REQUESTED

def record(v=6.,state='SOLVE_FLUX_PASS_MESH_PENDING'):
    return dict(voltage_V=v,state=state,evidence_errors=[],field_audit_recomputed={'passed':True})

class Preliminary(unittest.TestCase):
    def test_flux_failure_is_an_available_solution_but_never_a_pass(self):
        result=preliminary_status(record(state='SOLVED_FLUX_GATE_FAILED'))
        self.assertTrue(result['full_load_solution_available'])
        self.assertFalse(result['preliminary_flux_check_passed'])
        self.assertFalse(result['spatial_accuracy_certified'])

    def test_partial_or_invalid_field_is_not_complete(self):
        for r in [record(state='FAILED_NO_COMPLETE_STATE'),dict(record(),evidence_errors=['INPUT_HASH_MISMATCH']),
                  dict(record(),field_audit_recomputed={'passed':False})]:
            self.assertFalse(preliminary_status(r)['full_load_solution_available'])

    def test_missing_voltage_cannot_be_hidden_by_off_grid_points(self):
        result=summarize([record(6.),record(6.25)])
        self.assertFalse(result['all_requested_full_load_solutions_available'])
        self.assertEqual(len(result['points']),6)

    def test_complete_grid_and_flux_checks_are_separate(self):
        rows=[record(v) for v in REQUESTED];rows[-1]['state']='SOLVED_FLUX_GATE_FAILED'
        result=summarize(rows)
        self.assertTrue(result['all_requested_full_load_solutions_available'])
        self.assertFalse(result['all_requested_preliminary_checks_passed'])
        self.assertFalse(result['complete_quantitative_0_8V'])

    def test_later_failure_does_not_replace_a_valid_existing_solution(self):
        result=summarize([record(),record(state='FAILED_NO_COMPLETE_STATE')])
        self.assertTrue(result['points'][0]['preliminary_flux_check_passed'])

if __name__=='__main__':unittest.main()
