"""Continue unchanged preregistered solves after removing a diagnostic bitwise-equality assertion.

First solve's 7.9e-31 threaded difference is retained. No A2P acceptance gate changes.
"""
from run_diagnostic_retry01 import *

def main():
    t0=time.monotonic()
    protocol=json.loads((HERE/'protocol_build_retry01.json').read_text(encoding='utf-8'))
    exe=HERE/'build_retry01/pardiso_representation_driver.exe'
    environment=json.loads((HERE/'environment.json').read_text(encoding='utf-8'))
    assert sha(exe)==environment['executable_sha256']
    for path,digest in protocol['matrix_files'].items():assert sha(Path(path))==digest
    for path,digest in protocol['linked_objects'].items():assert sha(Path(path))==digest
    write(HERE/'continuation_record.json',dict(reason='Diagnostic runner stopped on exact floating-point direct/wrapper equality; the only observed difference was 7.888609052210118e-31. Report differences instead; no numerical setting or A2P gate changed.',
        runner_sha256=sha(Path(__file__)),protocol_sha256=sha(HERE/'protocol_build_retry01.json'),
        executable_sha256=sha(exe),existing_record='results/A2PMAT_gauge1_L0_original_rhs_original.json',
        continuation_policy='Reuse already recorded first solve. Execute all remaining preregistered representations and RHS sequentially.'))
    env=runtime_env();env['TPS_PARDISO_SPLIT_HOLES']='0'
    inputs=HERE/'inputs';results=HERE/'results';records=[];structural=[]
    expected=[0]*64
    for k,v in protocol['settings']['iparm_explicit_1based'].items():expected[int(k)-1]=v
    for pathstr in sorted(protocol['matrix_files']):
        path=Path(pathstr);d=load(path);A=d['A'];n=d['n'];free=np.flatnonzero(~d['fixed'])
        fidx=np.flatnonzero(d['fixed'])
        I=csr_matrix((np.ones(len(fidx)),(np.arange(len(fidx)),fidx)),shape=(len(fidx),n))
        assert maxabs((A[d['fixed'],:]-I).data)==0.0
        assert maxabs(d['rhs'][d['fixed']])==0.0
        reps=representations(d);structure=[]
        for label,B,selected in reps:
            assert maxabs((B-A[selected,:][:,selected]).data)==0.0
            diag_count=np.count_nonzero(B.indices==np.repeat(np.arange(B.shape[0]),np.diff(B.indptr)))
            assert diag_count==B.shape[0]
            structure.append(dict(representation=label,n=B.shape[0],stored_nnz=B.nnz,
                exact_zeros=int(np.count_nonzero(B.data==0)),retained_diagonals=int(diag_count),
                mathematical_operator_difference=0.0,
                constant_phi_action_max=maxabs(B@(selected%7==3).astype(float))))
        structural.append(dict(matrix=path.stem,gauge=d['gauge'],level=d['level'],representations=structure))
        if d['gauge']==0:continue
        rng=np.random.Generator(np.random.PCG64(20261018+d['level']))
        random_rhs=rng.standard_normal(n);random_rhs[d['fixed']]=0.0
        lu=splu(A[free,:][:,free].tocsc())
        for rhslabel,b in [('original_rhs',d['rhs']),('random_rhs',random_rhs)]:
            original=None;reference=np.zeros(n);reference[free]=lu.solve(b[free])
            for suffix,values in [('rhs',b),('reference',reference)]:
                savepath=results/(path.stem+'_'+rhslabel+'_'+suffix+'.npy')
                if savepath.exists():assert np.array_equal(np.load(savepath),values)
                else:np.save(savepath,values)
            for label,B,selected in reps:
                name=path.stem+'_'+rhslabel+'_'+label
                inp=inputs/(name+'.bin');out=results/(name+'.bin');log=results/(name+'.log')
                jsonpath=results/(name+'.json')
                if jsonpath.exists():
                    record=json.loads(jsonpath.read_text(encoding='utf-8'))
                    _,x,wx=read_output(out)
                    assert sha(inp)==record['input_sha256'] and sha(out)==record['output_sha256']
                else:
                    write_input(inp,B,b[selected]);started=time.monotonic()
                    with log.open('xb') as f:
                        p=subprocess.run([str(exe),str(inp),str(out)],cwd=exe.parent,env=env,stdout=f,
                                         stderr=subprocess.STDOUT,timeout=180)
                    record,x,wx=read_output(out)
                    record.update(matrix=path.stem,gauge=d['gauge'],level=d['level'],rhs=rhslabel,
                        representation=label,seconds=time.monotonic()-started,returncode=p.returncode,
                        input_sha256=sha(inp),output_sha256=sha(out),log_sha256=sha(log))
                assert record['iparm_before']==expected
                expanded=np.zeros(n);expanded[selected]=x
                ew=np.zeros(n);ew[selected]=wx
                if original is None:original=expanded.copy()
                electric=np.isin(np.arange(n)%7,[3,4,5])
                if not jsonpath.exists():
                    record.update(full_scaled_backward=backward(A,expanded,b,d['fixed']),
                        full_raw_backward=backward(d['raw'],d['columns']*expanded,d['rows']*b,d['fixed']),
                        linked_wrapper_full_scaled_backward=backward(A,ew,b,d['fixed']),
                        relative_solution_difference_original=maxabs(expanded-original)/max(maxabs(original),1e-300),
                        relative_solution_difference_superlu=maxabs(expanded-reference)/max(maxabs(reference),1e-300),
                        electric_relative_solution_difference_superlu=maxabs((expanded-reference)[electric])/max(maxabs(reference[electric]),1e-300),
                        full_scaled_solution_max=maxabs(expanded),electric_scaled_solution_max=maxabs(expanded[electric]))
                    if rhslabel=='original_rhs':record['original_export_solution_maxdiff']=maxabs(expanded-d['native_solution'])
                    write(jsonpath,record)
                records.append(record)
                print(json.dumps({k:record[k] for k in ['matrix','rhs','representation','pivots','wrapper_pivots',
                    'direct_wrapper_maxdiff','relative_solution_difference_superlu']}),flush=True)
                assert record['returncode']==0 and record['error']==0 and record['release_error']==0
    write(HERE/'structural_controls.json',structural)
    write(HERE/'summary.json',dict(status='DIAGNOSTIC_COMPLETED_NOT_A2P_ACCEPTANCE',
        continuation_seconds=time.monotonic()-t0,count=len(records),records=records,
        original_pardiso_source_unchanged=sha(ROOT/'src/linear/pardiso_solver_m.f90')==protocol['production_pardiso_source_sha256'],
        all_matrix_hashes_unchanged=all(sha(Path(p))==h for p,h in protocol['matrix_files'].items()),
        all_linked_objects_unchanged=all(sha(Path(p))==h for p,h in protocol['linked_objects'].items())))
    with (HERE/'summary.csv').open('x',newline='',encoding='utf-8') as f:
        fields=['matrix','rhs','representation','n','nnz','pivots','wrapper_pivots','refinements','factor_nnz',
                'direct_wrapper_maxdiff','relative_solution_difference_original','relative_solution_difference_superlu',
                'electric_relative_solution_difference_superlu','scaled_normwise_inf','scaled_componentwise_all',
                'scaled_componentwise_free','fixed_solution_max','seconds']
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
