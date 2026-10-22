"""Frozen-matrix, equivalent-representation PARDISO experiment. Diagnostic only."""
import os
import sys
import json
import hashlib
import subprocess
import time
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = HERE.parent / 'a2p_matrix_diagnostic'
sys.path.insert(0, str(OLD / 'python_packages'))
sys.path.insert(0, str(OLD))
from analyze_matrices import load
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
sys.path.insert(0, str(ROOT / 'scripts'))
from run_hv_regression import runtime_env

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write(path, data):
    with path.open('x', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')

def maxabs(x):
    return float(np.max(np.abs(x), initial=0.0))

def backward(A, x, b, fixed):
    r = A @ x - b
    denominator = abs(A) @ np.abs(x) + np.abs(b)
    ratio = np.abs(r) / np.maximum(denominator, np.finfo(float).tiny)
    anorm = maxabs(np.asarray(abs(A).sum(axis=1)))
    return dict(normwise_inf=maxabs(r) / max(anorm*maxabs(x)+maxabs(b),np.finfo(float).tiny),
                relative_action_l2=float(np.linalg.norm(r)/max(np.linalg.norm(A@x)+np.linalg.norm(b),1e-300)),
                componentwise_all_rows=maxabs(ratio),componentwise_free_rows=maxabs(ratio[~fixed]),
                absolute_residual_max=maxabs(r),fixed_solution_max=maxabs(x[fixed]),
                electric_residual_max=maxabs(r[np.isin(np.arange(len(b))%7,[3,4,5])]))

def remove_exact_offdiag_zeros(A):
    row = np.repeat(np.arange(A.shape[0]), np.diff(A.indptr))
    keep = (A.data != 0.0) | (A.indices == row)
    counts = np.bincount(row[keep], minlength=A.shape[0])
    ip = np.r_[0,np.cumsum(counts)].astype(np.int32)
    return csr_matrix((A.data[keep].copy(),A.indices[keep].copy(),ip),shape=A.shape)

def representations(d):
    A = d['A']; free=np.flatnonzero(~d['fixed'])
    reduced=A[free,:][:,free].tocsr()
    return [('original',A,np.arange(d['n'])),
            ('drop_exact_zeros',remove_exact_offdiag_zeros(A),np.arange(d['n'])),
            ('eliminate_fixed',reduced,free),
            ('eliminate_fixed_drop_exact_zeros',remove_exact_offdiag_zeros(reduced),free)]

def write_input(path,A,b):
    assert A.has_sorted_indices and not A.has_canonical_format is False
    with path.open('xb') as f:
        np.asarray([A.shape[0],A.nnz],dtype='<i4').tofile(f)
        np.asarray(A.indptr+1,dtype='<i4').tofile(f)
        np.asarray(A.indices+1,dtype='<i4').tofile(f)
        np.asarray(A.data,dtype='<f8').tofile(f)
        np.asarray(b,dtype='<f8').tofile(f)

def read_output(path):
    with path.open('rb') as f:
        h=np.fromfile(f,dtype='<i4',count=10)
        version,n,nnz,error,release,status,wrapper_error,wrapper_pivots,wrapper_refinements,wrapper_nnz=map(int,h)
        assert version==1
        before=np.fromfile(f,dtype='<i4',count=64)
        after=np.fromfile(f,dtype='<i4',count=64)
        released=np.fromfile(f,dtype='<i4',count=64)
        x=np.fromfile(f,dtype='<f8',count=n)
        wx=np.fromfile(f,dtype='<f8',count=n)
        assert not f.read()
    return dict(n=n,nnz=nnz,error=error,release_error=release,wrapper_status=status,
                wrapper_error=wrapper_error,pivots=int(after[13]),wrapper_pivots=wrapper_pivots,
                refinements=int(after[6]),wrapper_refinements=wrapper_refinements,
                factor_nnz=int(after[17]),wrapper_factor_nnz=wrapper_nnz,
                iparm_before=before.tolist(),iparm_after=after.tolist(),iparm_after_release=released.tolist(),
                direct_wrapper_maxdiff=maxabs(x-wx)),x,wx

def main():
    t0=time.monotonic()
    oldbuild=ROOT/'build/JE_VERIFY_20261018_MATRIX'
    matrices=sorted((OLD/'native/A2PMAT_fields').glob('A2PMAT_gauge*_L*.bin'))
    assert len(matrices)==9
    build=HERE/'build';build.mkdir(exist_ok=False)
    inputs=HERE/'inputs';inputs.mkdir()
    results=HERE/'results';results.mkdir()
    env=runtime_env()
    env['TPS_PARDISO_SPLIT_HOLES']='0'
    # Use the same 4-thread native environment as the original A2PMAT exporter.
    assert env['TPS_PARDISO_STORAGE']=='IC'
    assert env['OMP_NUM_THREADS']=='4' and env['MKL_NUM_THREADS']=='4'
    original_source=ROOT/'src/linear/pardiso_solver_m.f90'
    driver=HERE/'pardiso_representation_driver.f90'
    pattern=r'iparm\((\d+)\)\s*=\s*(-?\d+)'
    production_text=original_source.read_text(encoding='utf-8')
    plain_text=production_text.split('subroutine solve_csr_pardiso_plain(')[1]
    plain_text=plain_text.split("! Optional bounded-memory")[0]
    expected={int(i):int(v) for i,v in re.findall(pattern,plain_text)}
    copied={int(i):int(v) for i,v in re.findall(pattern,driver.read_text(encoding='utf-8'))}
    assert copied==expected, (copied,expected)
    objects=[oldbuild/(s+'.obj') for s in ['kinds_m','field_layout_m','status_m','pardiso_solver_m']]
    protocol=dict(study='A2P actual matrix PARDISO equivalent representations',date='2026-10-18',
        scope='Independent numerical diagnosis only. Does not change A2P acceptance or establish an A2P PASS.',
        matrix_files={str(p):sha(p) for p in matrices},
        representations=['original CSR','remove exactly-zero off-diagonal entries; retain all diagonals',
                         'strict elimination of fixed zero increments; retain other explicit zeros',
                         'strict fixed elimination plus exact off-diagonal zero removal'],
        equations='All nonzero entries unchanged bitwise. Reduced solutions expanded with exact fixed increment 0.',
        ungauged_policy='Gauge case 0: three known-singular structural controls only; never call native factorization.',
        rhs=['original scaled Newton RHS','NumPy PCG64 seed 20261018+level, standard_normal(n); fixed entries set exactly 0'],
        rhs_pairing='Same pre-mask random vector on both gauge locations of each level; RHS differences only at their gauge-fixed rows.',
        settings={'mtype':11,'maxfct':1,'mnum':1,'nrhs':1,'msglvl':0,'phase':13,'release_phase':-1,
                  'one_based':True,'LP64':True,'iparm_explicit_1based':copied,
                  'other_iparm_initial':0,'thread_count':4,'storage':'IC'},
        metrics=['all 64 iparm before/after/release','PARDISO pivots, refinements, factor nonzeros, errors',
                 'full ORIGINAL scaled matrix normwise and componentwise backward errors',
                 'full ORIGINAL raw matrix backward errors with recovered physical increment',
                 'direct-vs-linked-original-wrapper difference',
                 'solution difference from original CSR and independent exact-elimination SuperLU reference'],
        safeguards=['No tolerance-based coefficient dropping','No new scaling or reordering by diagnostic code',
                    'No changed PARDISO options','No PDE, BC, test thresholds or core source edits',
                    'One native child at a time; 180-second hard timeout per solve pair'],
        linked_objects={str(p):sha(p) for p in objects},
        production_pardiso_source_sha256=sha(original_source),driver_sha256=sha(driver),runner_sha256=sha(Path(__file__)))
    write(HERE/'protocol.json',protocol)
    compiler=Path(r'C:\Program Files (x86)\Intel\oneAPI\compiler\latest')
    mkl=Path(r'C:\Program Files (x86)\Intel\oneAPI\mkl\latest')
    flags=['/nologo','/O2','/debug:full','/traceback','/check:none','/Qopenmp','/warn:all','/warn:interfaces',
           '/stand:f18','/standard-semantics',f'/I"{oldbuild / "mod"}"','/exe:"pardiso_representation_driver.exe"']
    libs=[mkl/'lib'/s for s in ['mkl_intel_lp64.lib','mkl_intel_thread.lib','mkl_core.lib']]
    libs.append(compiler/'lib/libiomp5md.lib')
    rsp='\n'.join(flags+[f'"{driver}"']+[f'"{p}"' for p in objects+libs])+'\n'
    (build/'ifx.rsp').write_text(rsp,encoding='utf-8')
    batch='\n'.join(['@echo off',
        'call "C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat" >nul',
        'if errorlevel 1 exit /b %errorlevel%',
        f'set "PATH={compiler / "bin"};%PATH%"',
        f'set "LIB={compiler / "lib"};{mkl / "lib"};%LIB%"',
        f'"{compiler / "bin/ifx.exe"}" @ifx.rsp','exit /b %errorlevel%',''])
    (build/'build.cmd').write_text(batch,encoding='utf-8')
    with (build/'build.log').open('xb') as log:
        command=[os.environ['ComSpec'],'/d','/c',str(build/'build.cmd')]
        proc=subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT,timeout=120)
    write(build/'build_result.json',dict(command=command,returncode=proc.returncode,
        compiler=str(compiler/'bin/ifx.exe'),compiler_sha256=sha(compiler/'bin/ifx.exe'),
        libraries={str(p):sha(p) for p in libs},response_file_sha256=sha(build/'ifx.rsp'),
        build_log_sha256=sha(build/'build.log')))
    if proc.returncode: raise SystemExit('Build failed; preserve logs and inspect.')
    exe=build/'pardiso_representation_driver.exe'
    write(HERE/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        environment={k:v for k,v in env.items() if k.startswith(('OMP_','MKL_','TPS_PARDISO_','KMP_'))},
        executable=str(exe),executable_sha256=sha(exe)))
    records=[];structural=[]
    expected_iparm=[0]*64
    for k,v in copied.items():expected_iparm[k-1]=v
    for path in matrices:
        d=load(path);A=d['A'];n=d['n'];free=np.flatnonzero(~d['fixed'])
        fidx=np.flatnonzero(d['fixed'])
        I=csr_matrix((np.ones(len(fidx)),(np.arange(len(fidx)),fidx)),shape=(len(fidx),n))
        assert maxabs((A[d['fixed'],:]-I).data)==0.0
        assert maxabs(d['rhs'][d['fixed']])==0.0
        reps=representations(d)
        structure=[]
        for label,B,selected in reps:
            target=A[selected,:][:,selected]
            assert maxabs((B-target).data)==0.0
            diag_count=np.count_nonzero(B.indices==np.repeat(np.arange(B.shape[0]),np.diff(B.indptr)))
            assert diag_count==B.shape[0]
            null=(selected%7==3).astype(float)
            structure.append(dict(representation=label,n=B.shape[0],stored_nnz=B.nnz,
                exact_zeros=int(np.count_nonzero(B.data==0)),retained_diagonals=int(diag_count),
                mathematical_operator_difference=0.0,constant_phi_action_max=maxabs(B@null)))
        structural.append(dict(matrix=path.stem,gauge=d['gauge'],level=d['level'],representations=structure))
        if d['gauge']==0: continue
        rng=np.random.Generator(np.random.PCG64(20261018+d['level']))
        random_rhs=rng.standard_normal(n);random_rhs[d['fixed']]=0.0
        lu=splu(A[free,:][:,free].tocsc())
        for rhslabel,b in [('original_rhs',d['rhs']),('random_rhs',random_rhs)]:
            original=None
            reference=np.zeros(n);reference[free]=lu.solve(b[free])
            np.save(results/(path.stem+'_'+rhslabel+'_rhs.npy'),b)
            np.save(results/(path.stem+'_'+rhslabel+'_reference.npy'),reference)
            for label,B,selected in reps:
                name=path.stem+'_'+rhslabel+'_'+label
                inp=inputs/(name+'.bin');out=results/(name+'.bin');log=results/(name+'.log')
                write_input(inp,B,b[selected])
                started=time.monotonic()
                with log.open('xb') as f:
                    p=subprocess.run([str(exe),str(inp),str(out)],cwd=build,env=env,stdout=f,
                                     stderr=subprocess.STDOUT,timeout=180)
                record,x,wx=read_output(out)
                assert record['iparm_before']==expected_iparm
                expanded=np.zeros(n);expanded[selected]=x
                expanded_wrapper=np.zeros(n);expanded_wrapper[selected]=wx
                if original is None: original=expanded.copy()
                electric=np.isin(np.arange(n)%7,[3,4,5])
                record.update(matrix=path.stem,gauge=d['gauge'],level=d['level'],rhs=rhslabel,
                    representation=label,seconds=time.monotonic()-started,returncode=p.returncode,
                    input_sha256=sha(inp),output_sha256=sha(out),log_sha256=sha(log),
                    full_scaled_backward=backward(A,expanded,b,d['fixed']),
                    full_raw_backward=backward(d['raw'],d['columns']*expanded,d['rows']*b,d['fixed']),
                    linked_wrapper_full_scaled_backward=backward(A,expanded_wrapper,b,d['fixed']),
                    relative_solution_difference_original=maxabs(expanded-original)/max(maxabs(original),1e-300),
                    relative_solution_difference_superlu=maxabs(expanded-reference)/max(maxabs(reference),1e-300),
                    electric_relative_solution_difference_superlu=maxabs((expanded-reference)[electric])/max(maxabs(reference[electric]),1e-300),
                    full_scaled_solution_max=maxabs(expanded),electric_scaled_solution_max=maxabs(expanded[electric]))
                if rhslabel=='original_rhs':
                    record['original_export_solution_maxdiff']=maxabs(expanded-d['native_solution'])
                records.append(record)
                write(results/(name+'.json'),record)
                print(json.dumps({k:record[k] for k in ['matrix','rhs','representation','pivots','wrapper_pivots',
                                                       'direct_wrapper_maxdiff','relative_solution_difference_superlu']}) ,flush=True)
                assert p.returncode==0 and record['error']==0 and record['release_error']==0
                assert record['wrapper_pivots']==record['pivots']
                assert record['direct_wrapper_maxdiff']==0.0
    write(HERE/'structural_controls.json',structural)
    write(HERE/'summary.json',dict(status='DIAGNOSTIC_COMPLETED_NOT_A2P_ACCEPTANCE',
        seconds=time.monotonic()-t0,count=len(records),records=records,
        original_pardiso_source_unchanged=sha(original_source)==protocol['production_pardiso_source_sha256'],
        all_matrix_hashes_unchanged=all(sha(Path(p))==h for p,h in protocol['matrix_files'].items()),
        all_linked_objects_unchanged=all(sha(Path(p))==h for p,h in protocol['linked_objects'].items())))
    with (HERE/'summary.csv').open('x',newline='',encoding='utf-8') as f:
        fields=['matrix','rhs','representation','n','nnz','pivots','wrapper_pivots','refinements',
                'factor_nnz','direct_wrapper_maxdiff','relative_solution_difference_original',
                'relative_solution_difference_superlu','electric_relative_solution_difference_superlu',
                'scaled_normwise_inf','scaled_componentwise_all','scaled_componentwise_free',
                'fixed_solution_max','seconds']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in records:
            row={k:r[k] for k in fields if k in r}
            row.update(scaled_normwise_inf=r['full_scaled_backward']['normwise_inf'],
                       scaled_componentwise_all=r['full_scaled_backward']['componentwise_all_rows'],
                       scaled_componentwise_free=r['full_scaled_backward']['componentwise_free_rows'],
                       fixed_solution_max=r['full_scaled_backward']['fixed_solution_max'])
            w.writerow(row)
    write(HERE/'artifact_hashes.json',{str(p.relative_to(HERE)):sha(p) for p in HERE.rglob('*') if p.is_file()})

if __name__=='__main__':main()
