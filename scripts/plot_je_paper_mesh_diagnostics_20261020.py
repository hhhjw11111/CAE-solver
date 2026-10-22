"""Read-only source audit and figures for completed independent JE mesh studies.

No native solver, source changes, or acceptance upgrades. Refuses incomplete
series and existing output directories. Original gate outputs remain explicit.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import matplotlib.pyplot as plt
import run_je_reconstruction_mesh_20261018 as r
from plot_je_paper_supplements_20261018 import configure, export


def collect(study):
    rows, audits, evidence = [], {}, {}
    for direction, counts in [('angular', [32, 40, 48]), ('radial', [18, 20, 22])]:
        for voltage in (6, 8):
            if direction == 'angular':
                levels = [r.baseline(voltage)] + [study/r.name((a,14,14),voltage) for a in counts[1:]]
            else:
                levels = [r.radial18(voltage)] + [study/r.name((32,n,14),voltage) for n in counts[1:]]
            audit = r.mesh_comparison(levels, [voltage])
            record = audit['voltages'][0]
            if audit['provenance_errors'] or not record['complete_solve_flux_all_levels']:
                raise ValueError(f'Incomplete or invalid {direction}/{voltage}: {audit}')
            audits[f'{direction}_{voltage}V'] = audit
            previous = None
            for count, run in zip(counts, levels):
                manifest = r.accepted(run, require_new_completion=run.parent==study)
                _, m = r.classify(run, r.read(run/'execution.json')['exit_code'])
                power = r.read(run/'power_audit/decomposition.json')
                prov = r.read(run/'power_audit/provenance.json')
                if prov['mesh_sha256'] != manifest['mesh_sha256'] or prov['nodes_sha256'] != r.sha(run/'nodes.csv'):
                    raise ValueError('Power-audit source identity changed')
                for path, digest in prov['scripts'].items():
                    if r.sha(Path(path)) != digest:
                        raise ValueError('Power-audit code changed: '+path)
                    evidence[path] = digest
                ring = m['ring_abs_p_relative_mean']
                change = None if previous is None else 100*abs(ring-previous)/max(abs(ring),abs(previous))
                if change is not None:
                    pair = 'L0_L1' if count==counts[1] else 'L1_L2'
                    expected = 100*record['changes']['ring_abs_p_relative_mean'][pair]
                    if not math.isclose(change, expected, abs_tol=1e-12):
                        raise ValueError('Figure value differs from original gate')
                row = dict(direction=direction, voltage_V=voltage, direction_count=count,
                    run=str(run), theta_max_K=m['theta_max_K'], JE_mW=m['joule_W']*1000,
                    power_gap_all_exterior_percent=100*power['gap_relative'],
                    power_gap_electrode_percent=100*power['electrode_gap_relative'],
                    ring_abs_p_relative_mean=ring, ring_p_adjacent_change_percent=change,
                    ring_p_limit_percent=100*record['changes']['ring_abs_p_relative_mean']['limit'],
                    complete_three_level_scalar_ring_pass=record['mesh_pass'],
                    local_n_change_percent=None, local_p_change_percent=None, local_all_fields_pass=None)
                if run.parent==study:
                    local = r.read(run/'mapped_field_comparison.json')
                    for path, digest in local['hashes'].items():
                        if r.sha(Path(path)) != digest:
                            raise ValueError('Field-comparison identity changed: '+path)
                        evidence[path] = digest
                    row.update(local_n_change_percent=100*local['comparisons']['eta_n_density_perturbation']['maximum_relative'],
                        local_p_change_percent=100*local['comparisons']['eta_p_density_perturbation']['maximum_relative'],
                        local_all_fields_pass=local['passed'])
                rows.append(row)
                previous = ring
                for path in [run/'metrics.txt',run/'execution.json',run/'manifest.json',run/'input.nml',
                             run/'field_audit.json',run/'power_audit/decomposition.json',run/'power_audit/provenance.json',
                             run/'nodes.csv',Path(manifest['mesh']),Path(manifest['executable'])]:
                    evidence[str(path)] = r.sha(path)
    return rows, audits, evidence


def plot(rows, output):
    configure()
    fig, axes = plt.subplots(2,2,figsize=(183/25.4,150/25.4))
    fig.subplots_adjust(left=.11,right=.975,bottom=.14,top=.88,wspace=.34,hspace=.46)
    fig.suptitle('Independent mesh studies | JE branch, r = 3 µm',y=.97,fontsize=10,fontweight='bold')
    for ax, letter in zip(axes.flat,'abcd'):
        ax.text(-.18,1.06,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
        ax.tick_params(direction='out',length=3,width=.7)
    for col,direction in enumerate(('angular','radial')):
        for voltage,color,marker in [(6,'#0072B2','o'),(8,'#D55E00','s')]:
            group=[x for x in rows if x['direction']==direction and x['voltage_V']==voltage]
            axes[0,col].plot([x['direction_count'] for x in group],
                [x['power_gap_all_exterior_percent'] for x in group],color=color,marker=marker,
                ms=4,mfc='white',label=f'{voltage} V')
            axes[1,col].plot([0,1],[x['ring_p_adjacent_change_percent'] for x in group[1:]],
                color=color,marker=marker,ms=4,mfc='white',label=f'{voltage} V')
        counts=[32,40,48] if direction=='angular' else [18,20,22]
        axes[0,col].set_xticks(counts)
        axes[0,col].set_xlabel('Circumferential elements' if direction=='angular' else 'Radial elements')
        axes[0,col].set_title('Angular refinement' if direction=='angular' else 'Radial refinement',pad=10)
        axes[0,col].set_ylim(bottom=0)
        axes[0,col].legend(loc='upper right')
        axes[1,col].set_xticks([0,1],[f'{counts[0]} → {counts[1]}',f'{counts[1]} → {counts[2]}'])
        axes[1,col].set_xlim(-.2,1.2)
        axes[1,col].set_ylim(bottom=0)
        axes[1,col].set_xlabel('Successive mesh comparison')
        axes[1,col].axhline(2,color='#666666',ls='--',label='Original 2% limit')
        axes[1,col].legend(loc='upper right')
        axes[0,col].set_ylabel('All-exterior power gap (%)')
        axes[1,col].set_ylabel('Ring hole-metric change (%)')
    # Equal row scales permit direct angular/radial comparison.
    for row in range(2):
        upper=max(ax.get_ylim()[1] for ax in axes[row])
        for ax in axes[row]: ax.set_ylim(0,upper)
    fig.text(.5,.045,'Power gap = (volume JE − physical exterior power) / volume JE.',ha='center',fontsize=7.5)
    fig.text(.5,.018,'Diagnostic evidence; failed mesh gates retained. Not a production or SUPG certificate.',ha='center',fontsize=7.5)
    export(fig,output/'mesh_power_and_hole_diagnostics')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=r.STUDY/'stage4_paper_mesh')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rows,audits,evidence=collect(a.study.resolve())
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'source_metrics.csv').open('x',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    r.save(a.output/'original_mesh_audits.json',audits)
    plot(rows,a.output)
    for path,digest in evidence.items():
        if r.sha(Path(path)) != digest: raise ValueError('Source changed while plotting: '+path)
    r.save(a.output/'provenance.json',dict(sources=evidence,script_sha256=r.sha(Path(__file__)),
        original_gates_unchanged=True,production_acceptance=False,
        interpretation='Deterministic independent-direction comparison. No extrapolated continuum error or statistical confidence intervals. '
                       'Ring means use absolute relative hole perturbation; local-field metrics use their separately frozen normalization.'))
    r.save(a.output/'artifact_hashes.json',{p.name:r.sha(p) for p in a.output.iterdir() if p.is_file()})
    print(json.dumps(dict(rows=len(rows),studies=len(audits),production_acceptance=False)))


if __name__=='__main__': main()
