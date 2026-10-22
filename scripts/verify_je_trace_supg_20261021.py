"""Verify the already implemented trace-metric SUPG using frozen native kernels.

Build an isolated verification executable from unchanged production sources.
Reuse the original SUPG high/low-Pe, current and monotonicity limits verbatim;
test both unit transverse lengths and 0.01 transverse lengths.  Only verification
copies select trace tau and account for the changed cross-sectional area.
"""
from pathlib import Path
import json
import argparse
import subprocess
import time
import run_je_reconstruction_mesh_20261018 as r

OUT = r.STUDY/'stage5_controlled_numerics_20261021/native_trace_verification_attempt02'
RUNTIME_CHECKS = True
SOURCES = ('core/kinds_m', 'core/constants_m', 'core/status_m', 'core/field_layout_m',
    'mesh/mesh_types_m', 'core/model_types_m', 'physics/equilibrium_m',
    'physics/reference_material_m', 'physics/transport_m', 'physics/recombination_m',
    'physics/thermo_piezo_m', 'physics/material_point_m', 'mesh/hex20_shape_m',
    'mesh/hex20_quadrature_m', 'mesh/dof_map_m', 'fem/transport_supg_m',
    'fem/element_kernel_m', 'verification/material_tangent_check_m',
    'verification/element_tangent_check_m')


def once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError('Expected one frozen verification fragment: '+old)
    return text.replace(old, new)


def seal_evidence():
    acceptance = r.read(OUT/'acceptance.json')
    if acceptance.get('passed') is not True or acceptance.get('exit_code') != 0:
        raise RuntimeError('Cannot seal a failed or incomplete trace verification')
    r.save(OUT/'evidence_identity.json', dict(
        acceptance_sha256=r.sha(OUT/'acceptance.json'),
        verification_log_sha256=r.sha(OUT/'verification.txt'),
        verification_executable_sha256=r.sha(OUT/'trace_check.exe')))


def execute():
    gate = r.validate_protocol('paper')
    OUT.mkdir(exist_ok=False)
    (OUT/'mod').mkdir()
    template = r.ROOT/'src/verification/supg_check_m.f90'
    text = template.read_text(encoding='utf-8')
    copied = []
    for name, aspect in [('iso', 1.), ('thin', .01)]:
        body = text.replace('supg_check_m', 'supg_trace_'+name+'_check_m')
        body = body.replace('verify_supg', 'verify_supg_trace_'+name)
        body = once(body, 'm%supg_directional_tau=.true.', 'm%supg_directional_tau=.false.')
        body = once(body, 'call verify_element_g03(status,.true.,.true.)',
                    'call verify_element_g03(status,.true.,.false.)')
        if aspect != 1.:
            body = once(body, 'coords=0.5_dp*(nat+1.0_dp);coords(1,:)=left+coords(1,:)*(right-left)',
                'coords=0.5_dp*(nat+1.0_dp);coords(2:3,:)=0.01_dp*coords(2:3,:)\n'
                '        coords(1,:)=left+coords(1,:)*(right-left)')
            body = once(body, 'scale=ELEMENTARY_CHARGE*m%reference_n*0.01_dp*thermal',
                'scale=ELEMENTARY_CHARGE*m%reference_n*0.01_dp*thermal*1.0e-4_dp')
        body = body.replace('verification_id=SUPG_HEX20_HIGH_PE',
                            'verification_id=TRACE_SUPG_'+name.upper()+'_HIGH_PE')
        path = OUT/('trace_'+name+'.f90')
        path.write_text(body, encoding='ascii')
        copied.append(path)
    driver = OUT/'driver.f90'
    driver.write_text('''program trace_controlled_verification
  use status_m, only: SolverStatusType
  use supg_trace_iso_check_m, only: verify_supg_trace_iso
  use supg_trace_thin_check_m, only: verify_supg_trace_thin
  implicit none
  type(SolverStatusType) :: status
  call verify_supg_trace_iso(status)
  if (.not.status%is_ok()) then
    write(*,'(a)') trim(status%message)
    stop 1
  end if
  call verify_supg_trace_thin(status)
  if (.not.status%is_ok()) then
    write(*,'(a)') trim(status%message)
    stop 2
  end if
  write(*,'(a)') 'trace_controlled_verification=PASS'
end program trace_controlled_verification
''', encoding='ascii')
    original = [r.ROOT/'src'/(p+'.f90') for p in SOURCES]+[template]
    hashes = {str(p):r.sha(p) for p in original}
    args = ['/nologo', '/O2', '/traceback']
    args += ['/check:all', '/check:noarg_temp_created'] if RUNTIME_CHECKS else ['/check:none']
    args += [
            '/warn:all', '/stand:f18', '/standard-semantics', '/module:mod', '/exe:trace_check.exe']
    # Shorten the Unicode workspace prefix only. Shortening each source file
    # makes TRANSPORT_M and TRANSPORT_SUPG_M both become TRANSP~1.obj at link time.
    args += ['"'+str(Path(r.short(r.ROOT))/p.relative_to(r.ROOT))+'"' for p in original[:-1]+copied+[driver]]
    (OUT/'ifx.rsp').write_text('\n'.join(args)+'\n', encoding='ascii')
    ifx = Path('C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe')
    vcvars = Path('C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build/vcvars64.bat')
    batch = '\n'.join(['@echo off', 'call "'+str(vcvars)+'" >nul',
        'if errorlevel 1 exit /b %errorlevel%',
        'set "PATH='+str(ifx.parent)+';%PATH%"',
        'set "LIB=C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\lib;%LIB%"',
        '"'+str(ifx)+'" @ifx.rsp', 'exit /b %errorlevel%'])+'\n'
    (OUT/'build.cmd').write_text(batch, encoding='ascii')
    start = time.monotonic()
    build = subprocess.run(['cmd', '/d', '/c', 'build.cmd'], cwd=r.short(OUT), capture_output=True)
    (OUT/'build.log').write_bytes(build.stdout+build.stderr)
    if build.returncode:
        r.save(OUT/'build_failed.json', dict(exit_code=build.returncode))
        raise RuntimeError('Isolated trace verification build failed; see '+str(OUT/'build.log'))
    result = subprocess.run([r.short(OUT/'trace_check.exe')], cwd=r.short(OUT),
                            env=r.runtime_env(), capture_output=True)
    (OUT/'verification.txt').write_bytes(result.stdout+result.stderr)
    unchanged = all(r.sha(Path(p)) == digest for p,digest in hashes.items())
    passed = result.returncode == 0 and b'trace_controlled_verification=PASS' in result.stdout and unchanged
    r.save(OUT/'acceptance.json', dict(passed=passed, exit_code=result.returncode,
        seconds=time.monotonic()-start, original_core_unchanged=unchanged,
        verification_source_hashes=hashes, copied_source_hashes={str(p):r.sha(p) for p in copied+[driver]},
        frozen_verification_gate=gate, verification_executable_sha256=r.sha(OUT/'trace_check.exe'),
        production_executable_sha256=r.sha(r.EXE), original_thresholds_unchanged=True,
        verification_log_sha256=r.sha(OUT/'verification.txt'),
        anisotropic_cross_section_scale=.01, anisotropic_residual_area_normalization=1e-4,
        high_Pe_values=[.1,1,10,50,200,1000], density_limit=1e-3, current_limit=1e-2,
        monotonicity_limit=1e-3, low_Pe_Galerkin_limit=1e-7,
        production_acceptance=False))
    print(json.dumps(dict(passed=passed, exit_code=result.returncode, output=str(OUT))), flush=True)
    if not passed:
        print(result.stdout.decode('utf8', errors='replace')[-3500:])
        raise RuntimeError('Trace numerical verification failed')
    seal_evidence()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized', action='store_true', help='Use the unchanged production O2/check:none flags')
    settings = parser.parse_args()
    if settings.optimized:
        OUT = r.STUDY/'stage5_controlled_numerics_20261021/native_trace_verification_optimized'
        RUNTIME_CHECKS = False
    execute()
