"""Protect separate runs from inherited OOC paths and implicit solver changes."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from hv_runtime import configure_storage


class RuntimeIsolation(unittest.TestCase):
    def test_runs_do_not_share_or_overwrite_inherited_scratch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            previous = root/'previous.cfg'
            previous.write_text('protected prior configuration')
            original = dict(MKL_PARDISO_OOC_FILE_NAME=str(root/'previous'),
                            MKL_PARDISO_OOC_CFG_FILE_NAME=str(previous))
            configurations = []
            for name in ['first', 'second']:
                out = root/name; out.mkdir()
                env = configure_storage(original, out, 'OOC', 2000)
                prefix = Path(env['MKL_PARDISO_OOC_FILE_NAME'])
                self.assertTrue(prefix.is_relative_to(out))
                self.assertTrue((prefix.parent/'run_ooc.cfg').is_file())
                configurations.append(prefix)
            self.assertNotEqual(*configurations)
            self.assertEqual(previous.read_text(), 'protected prior configuration')
            self.assertEqual(original['MKL_PARDISO_OOC_CFG_FILE_NAME'], str(previous))

    def test_default_does_not_inherit_experimental_block_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            env = configure_storage({'TPS_PARDISO_SPLIT_HOLES':'1',
                                     'MKL_PARDISO_OOC_FILE_NAME':'unrelated'}, out, 'IC', 2000)
            self.assertEqual(env['TPS_PARDISO_SPLIT_HOLES'], '0')
            self.assertEqual(env['TPS_PARDISO_STORAGE'], 'IC')
            self.assertFalse(list(out.iterdir()))
            self.assertNotIn('MKL_PARDISO_OOC_FILE_NAME', env)

    def test_invalid_budget_is_rejected_before_creating_files(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            for invalid in [0, -1, 0.5]:
                with self.assertRaises(ValueError):
                    configure_storage({}, out, 'OOC', invalid)
            self.assertFalse(list(out.iterdir()))


if __name__ == '__main__':
    unittest.main()
