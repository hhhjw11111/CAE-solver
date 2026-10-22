"""Run the unchanged JE half problem with the existing checked block solver.

Use a private instance of the frozen resource-monitored adapter.  Its sole
runtime difference is TPS_PARDISO_SPLIT_HOLES=1, which checks the original full
matrix to 1e-11 and falls back to the ordinary solve when that check fails.
Input equations, residual/Jacobian assembly and acceptance scripts are unchanged.
"""
from pathlib import Path
import importlib.util
from hv_runtime import configure_storage
import run_je_reconstruction_mesh_20261018 as r

ADAPTER = Path(__file__).with_name('run_je_half_native_20261020_v2.py')


def run_native(job, mesh, seed, directory, gate, ooc_mb=8000):
    spec = importlib.util.spec_from_file_location('_je_controlled_private_native', ADAPTER)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)

    def checked_storage(env, output, storage, core_mb, split_holes=False):
        if split_holes is not False:
            raise RuntimeError('Unexpected inherited block-solver control')
        return configure_storage(env, output, storage, core_mb, split_holes=True)

    adapter.configure_storage = checked_storage
    adapter.OOC_MB = ooc_mb
    adapter.run_native(job, mesh, seed, directory, gate)
    manifest = r.read(directory/'manifest.json')
    if manifest['runtime_env']['TPS_PARDISO_SPLIT_HOLES'] != '1':
        raise RuntimeError('Checked block mode not recorded')
    r.save(directory/'checked_block_provenance.json', dict(
        source_adapter=str(ADAPTER), adapter_sha256=r.sha(ADAPTER),
        wrapper_sha256=r.sha(Path(__file__)), executable_sha256=r.EXPECTED_EXE,
        original_full_matrix_backward_error_limit=1e-11,
        equations_unchanged=True, native_gates_unchanged=True,
        runtime_mode=manifest['runtime_env']['TPS_PARDISO_SPLIT_HOLES']))

