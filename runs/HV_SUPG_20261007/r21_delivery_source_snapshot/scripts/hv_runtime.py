"""Per-run PARDISO storage selection; does not change sparse equations."""
from pathlib import Path


def add_storage_arguments(parser):
    parser.add_argument('--storage', choices=['IC', 'OOC'], default='IC')
    parser.add_argument('--split-holes', action='store_true',
                        help='Experimental checked block solve; original full matrix is the fallback')
    parser.add_argument('--ooc-core-mb', type=int, default=12000,
                        help='PARDISO OOC RAM budget in MB, excluding assembly and other processes')


def configure_storage(env, output, storage, core_mb, split_holes=False):
    if storage not in ('IC', 'OOC') or not isinstance(core_mb, int) or core_mb <= 0:
        raise ValueError('Invalid PARDISO storage configuration')
    result = dict(env)
    # Never inherit another run's scratch prefix or configuration file.
    for key in list(result):
        if key.upper().startswith('MKL_PARDISO_OOC_'):
            del result[key]
    result['TPS_PARDISO_STORAGE'] = storage
    result['TPS_PARDISO_SPLIT_HOLES'] = '1' if split_holes else '0'
    if storage == 'OOC':
        scratch = Path(output).absolute() / 'pardiso_ooc'
        scratch.mkdir(exist_ok=False)
        result.update(MKL_PARDISO_OOC_FILE_NAME=str(scratch / 'factors'),
                      MKL_PARDISO_OOC_MAX_CORE_SIZE=str(core_mb),
                      MKL_PARDISO_OOC_MAX_SWAP_SIZE='0',
                      MKL_PARDISO_OOC_KEEP_FILE='1',
                      MKL_PARDISO_OOC_CFG_PATH=str(scratch),
                      MKL_PARDISO_OOC_CFG_FILE_NAME='run_ooc.cfg')
        # An explicit run-owned configuration prevents a working-directory
        # pardiso_ooc.cfg from silently overriding this run's memory budget.
        (scratch / 'run_ooc.cfg').write_text(
            f'MKL_PARDISO_OOC_PATH = {scratch}\n'
            f'MKL_PARDISO_OOC_MAX_CORE_SIZE = {core_mb}\n'
            'MKL_PARDISO_OOC_MAX_SWAP_SIZE = 0\n'
            'MKL_PARDISO_OOC_KEEP_FILE = 1\n', encoding='ascii')
    return result


def recorded_runtime(env):
    keys = ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'TPS_PARDISO_STORAGE', 'TPS_PARDISO_SPLIT_HOLES']
    keys.extend(sorted(k for k in env if k.startswith('MKL_PARDISO_OOC_')))
    return {k: env[k] for k in keys}
