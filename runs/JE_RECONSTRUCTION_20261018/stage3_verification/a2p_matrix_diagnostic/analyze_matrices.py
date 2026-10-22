"""Independent double-precision diagnostics of actual A2PMAT CSR snapshots."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '8'
import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / 'python_packages'))
import numpy as np
import scipy
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

FIELDS = ['u1', 'u2', 'u3', 'phi', 'eta_n', 'eta_p', 'theta']
def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
def write(p, data):
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
def load(p):
    with p.open('rb') as f:
        header = np.fromfile(f, dtype='<i4', count=11)
        version,n,nnz,nodes,anchor,gauge_dof,level,gauge,pivots,refinements,error = map(int, header)
        assert version == 1 and n == nodes*7
        ia = np.fromfile(f, dtype='<i4', count=n+1)-1
        ja = np.fromfile(f, dtype='<i4', count=nnz)-1
        data = dict(path=str(p), sha256=sha(p), n=n, nnz=nnz, nodes=nodes, anchor=anchor,
                    gauge_dof=gauge_dof, level=level, gauge=gauge, pivots=pivots,
                    refinements=refinements, native_error=error)
        raw = np.fromfile(f,dtype='<f8',count=nnz)
        scaled = np.fromfile(f,dtype='<f8',count=nnz)
        data['raw'] = csr_matrix((raw,ja,ia),shape=(n,n))
        data['A'] = csr_matrix((scaled,ja,ia),shape=(n,n))
        for name in ('rhs_raw','rhs','columns','rows','native_solution'):
            data[name] = np.fromfile(f,dtype='<f8',count=n)
        data['fixed'] = np.fromfile(f,dtype='<i4',count=n).astype(bool)
        data['coordinates'] = np.fromfile(f,dtype='<f8',count=3*nodes).reshape((3,nodes),order='F')
        assert len(f.read()) == 0
    return data
def maxabs(x):
    return float(np.max(np.abs(x),initial=0.))
def backward(A,x,b):
    residual = A@x-b
    denom = abs(A)@np.abs(x)+np.abs(b)
    return dict(normwise=float(np.linalg.norm(residual)/(np.linalg.norm(b)+np.linalg.norm(A@x)+1e-300)),
                componentwise=float(np.max(np.abs(residual)/np.maximum(denom,1e-300),initial=0.)),
                maximum_absolute=maxabs(residual))
def lu_probe(A,b,native=None):
    t0=time.monotonic()
    try:
        lu=splu(A.tocsc())
        x=lu.solve(b)
        diag=np.abs(lu.U.diagonal())
        report=dict(success=True,seconds=time.monotonic()-t0,
                    U_diagonal_min=float(diag.min()),U_diagonal_max=float(diag.max()),
                    residual=backward(A,x,b))
        if native is not None:
            report['native_difference_max']=maxabs(x-native)
        return report,x
    except Exception as exc:
        return dict(success=False,seconds=time.monotonic()-t0,error=str(exc)),None
def scales_by_field(d):
    A=d['A']; am=abs(A)
    rmax=np.asarray(am.max(axis=1).toarray()).ravel()
    cmax=np.asarray(am.max(axis=0).toarray()).ravel()
    result={}
    for field,name in enumerate(FIELDS):
        ind=np.arange(field,d['n'],7)
        result[name]=dict(unknowns=len(ind),fixed=int(d['fixed'][ind].sum()),
                         column_variable_scale=[float(d['columns'][ind].min()),float(d['columns'][ind].max())],
                         row_equilibration_scale=[float(d['rows'][ind].min()),float(d['rows'][ind].max())],
                         scaled_row_max=[float(rmax[ind].min()),float(rmax[ind].max())],
                         scaled_column_max=[float(cmax[ind].min()),float(cmax[ind].max())])
    return result
def analyze(d,out):
    name=Path(d['path']).stem
    A=d['A']; raw=d['raw']; n=d['n']; free=np.flatnonzero(~d['fixed'])
    electrical=free[np.isin(free%7,[3,4,5])]
    mechanical=free[free%7==2]
    report={k:d[k] for k in ('path','sha256','n','nnz','nodes','anchor','gauge_dof','level','gauge',
                             'pivots','refinements','native_error')}
    reconstruction = raw.data*d['columns'][raw.indices]/np.repeat(d['rows'],np.diff(raw.indptr))
    report['scaling_reconstruction_max_error']=maxabs(A.data-reconstruction)
    report['explicit_zero_entries']=int(np.count_nonzero(A.data==0))
    report['scales']=scales_by_field(d)
    fixed_rows=A[d['fixed'],:]
    fidx=np.flatnonzero(d['fixed'])
    report['fixed_rows_identity_error']=maxabs((fixed_rows-csr_matrix((np.ones(len(fidx)),
        (np.arange(len(fidx)),fidx)),shape=fixed_rows.shape)).data)
    report['electrical_free_dimension']=len(electrical)
    report['mechanical_free_dimension']=len(mechanical)
    report['mechanical_electrical_coupling_max']=max(maxabs(A[mechanical,:][:,electrical].data),
                                                    maxabs(A[electrical,:][:,mechanical].data))
    if d['gauge']:
        report['native_backward_error']=backward(A,d['native_solution'],d['rhs'])
        report['superlu_full_stored_zeros'],xfull=lu_probe(A,d['rhs'],d['native_solution'])
        compressed=A.copy();compressed.eliminate_zeros()
        report['superlu_full_exact_zeros_removed'],xcompressed=lu_probe(compressed,d['rhs'],d['native_solution'])
        Afree=A[free,:][:,free]
        report['superlu_strict_fixed_DOF_elimination'],xfree=lu_probe(Afree,d['rhs'][free],d['native_solution'][free])
        if xfree is not None:
            expanded=np.zeros(n);expanded[free]=xfree
            report['eliminated_solution_full_residual']=backward(A,expanded,d['rhs'])
            np.save(out/(name+'_eliminated_solution.npy'),expanded)
    B=A[electrical,:][:,electrical].toarray()
    print(json.dumps(dict(case=name,phase='dense_active_electrical_SVD',dimension=len(electrical))),flush=True)
    t0=time.monotonic()
    u,s,vh=linalg.svd(B,full_matrices=False,check_finite=False,lapack_driver='gesdd')
    tau=np.finfo(float).eps*len(electrical)*s[0]
    report['electrical_svd']=dict(seconds=time.monotonic()-t0,sigma_max=float(s[0]),
        smallest_singular_values=[float(x) for x in s[-8:]],rank_threshold=float(tau),
        numerical_rank=int(np.count_nonzero(s>tau)),dimension=len(s),
        condition2=float(s[0]/s[-1]) if s[-1]>0 else None,
        smallest_triplet_right_residual=maxabs(B@vh[-1,:]-s[-1]*u[:,-1]),
        smallest_triplet_left_residual=maxabs(B.T@u[:,-1]-s[-1]*vh[-1,:]))
    np.save(out/(name+'_electrical_singular_values.npy'),s)
    right=vh[-1,:];left=u[:,-1]
    expected=(electrical%7==3).astype(float);expected/=np.linalg.norm(expected)
    report['smallest_right_mode_constant_phi_alignment']=float(abs(right@expected))
    report['expected_constant_phi_residual']=maxabs(B@expected)
    report['right_mode_norm_by_field']={FIELDS[k]:float(np.linalg.norm(right[electrical%7==k])) for k in (3,4,5)}
    report['left_mode_norm_by_field']={FIELDS[k]:float(np.linalg.norm(left[electrical%7==k])) for k in (3,4,5)}
    with (out/(name+'_smallest_modes.csv')).open('x',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['dof_1based','node_1based','field','left_scaled','right_scaled',
                                             'left_raw_row_dual','right_physical'])
        writer.writerows((int(i+1),int(i//7+1),FIELDS[i%7],float(l),float(r),
                          float(l/d['rows'][i]),float(r*d['columns'][i]))
                         for i,l,r in zip(electrical,left,right))
    if not d['gauge']:
        pins=[([21,33,45][d['level']]-1)*7+3,([31,43,55][d['level']]-1)*7+3]
        report['left_null_at_candidate_gauge_rows']=[dict(dof_1based=pin+1,
            left_scaled=float(left[np.flatnonzero(electrical==pin)[0]]),
            left_raw=float(left[np.flatnonzero(electrical==pin)[0]]/d['rows'][pin])) for pin in pins]
    write(out/(name+'.json'),report)
    print(json.dumps(dict(case=name,rank=report['electrical_svd']['numerical_rank'],dimension=len(s),
                         native_pivots=d['pivots'],sigma_min=float(s[-1]))),flush=True)
    return report
def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    write(a.output/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
          executable=sys.executable,analysis_script_sha256=sha(Path(__file__)),protocol_sha256=sha(BASE/'protocol.json')))
    files=sorted(a.matrix_dir.glob('A2PMAT_*.bin'))
    if len(files)!=9: raise RuntimeError(f'Expected 9 matrices, found {len(files)}')
    reports=[]
    for file in files:
        reports.append(analyze(load(file),a.output))
        write(a.output/'progress.json',reports)
    compatibility=[]
    for level in range(3):
        ung=load(a.matrix_dir/f'A2PMAT_gauge0_L{level}.bin')
        for gauge in (1,2):
            d=load(a.matrix_dir/f'A2PMAT_gauge{gauge}_L{level}.bin')
            physical=d['native_solution']*d['columns']
            residual=ung['raw']@physical-ung['rhs_raw']
            pin=d['gauge_dof']-1
            compatibility.append(dict(level=level,gauge=gauge,
                replaced_gauss_row_raw_residual=float(residual[pin]),
                maximum_original_gauss_residual=maxabs(residual[3::7]),
                maximum_original_gauss_scaled_residual=maxabs(residual[3::7]/ung['rows'][3::7])))
    summary=dict(diagnostic_only=True,not_A2P_acceptance=True,reports=reports,gauss_compatibility=compatibility)
    write(a.output/'summary.json',summary)
    write(a.output/'artifact_hashes.json',{str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file()})
    print(json.dumps(dict(completed=True,matrix_count=len(reports),not_A2P_acceptance=True)),flush=True)
if __name__=='__main__': main()
