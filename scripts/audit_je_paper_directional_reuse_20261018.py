"""Read-only replay of the existing JE Z14/16/18 directional evidence.

No native solver, mesh generator, transfer, or power quadrature is launched.
The original mesh/local-field gates are rerun; saved independent power integrals
are hash-verified and their arithmetic is reconstructed. Only a new output
directory below paper_scope/directional_reuse may be written.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from audit_hv_field_shape import audit as field_audit
from audit_hv_mesh import audit as mesh_audit
from audit_je_paper_reuse_20261018 import validate_acceptance
from run_hv_sweep import audit_nodes, classify, sha
from run_je_sensitivity_20261015 import replace

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'runs/JE_SENSITIVITY_20261015'
FINE = ROOT / 'runs/JE_FINE_CROSS_20261016'
DIRECTION = FINE / 'GLOBAL_DIRECTIONAL'
STUDY = ROOT / 'runs/JE_RECONSTRUCTION_20261018'
PAPER = STUDY / 'paper_scope'
OUTPUT = PAPER / 'directional_reuse'
EXE = ROOT / 'build/ARC_JE_FINAL_20261012/TPSolverV2.exe'
EXE_SHA = '405e953c8cb0b1eca2a61eccc011945892924b43fca6805e87fd8fec19006d92'
RUNNER = ROOT / 'scripts/run_je_reconstruction_mesh_20261018.py'
RUNNER_SHA = '087a9ff302a63f5e87c9d937657265702190eddf52117ae3e7e2fcb66fad0c9a'
MODEL = 'JE_CURRENT_WORK_DIAGNOSTIC_ZNO_NT'
HEAT_SOURCE = 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
ALLOWED_INPUT_CHANGES = ['case_id', 'mesh_file', 'restart_nodes_file',
                         'restart_voltage_factor', 'supg_scale']


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def require(value, message):
    if not value:
        raise ValueError(message)


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


class Evidence:
    def __init__(self):
        self.hashes = {}
        self.frozen = {}
        self.reference = {}

    def record(self, path, expected=None):
        path = Path(path).resolve()
        require(path.is_file(), f'Missing evidence: {path}')
        digest = sha(path)
        if expected is not None:
            require(digest == expected, f'Evidence hash mismatch: {path}')
        require(str(path) not in self.hashes or self.hashes[str(path)] == digest,
                f'Evidence changed during audit: {path}')
        self.hashes[str(path)] = digest
        return digest

    def inventory(self, path, base):
        self.record(path)
        for key, digest in read(path).items():
            target = (base / key).resolve()
            if target in self.frozen:
                require(self.frozen[target] == digest, f'Conflicting frozen inventories: {target}')
            self.frozen[target] = digest
            self.reference[str(target)] = str(path)

    def historical(self, path):
        path = Path(path).resolve()
        require(path in self.frozen, f'No historical inventory entry: {path}')
        return self.record(path, self.frozen[path])

    def subtree(self, path):
        path = Path(path).resolve()
        found = [p for p in self.frozen if p.is_relative_to(path)]
        require(found, f'No frozen evidence below: {path}')
        for target in found:
            self.historical(target)

    def verify_unchanged(self):
        for key, digest in self.hashes.items():
            require(sha(Path(key)) == digest, f'Changed during read-only audit: {key}')


def thickness_run(voltage, z):
    return (OLD / f'A32_S1_V{voltage * 10:03d}' if z == 14 else
            DIRECTION / f'A32_R14_Z{z}_S1_V{voltage * 10:03d}')


def replay_power(directory, metrics, manifest, evidence):
    """Reconstruct saved order-3/5 arithmetic without claiming new quadrature."""
    out = directory / 'power_audit'
    for name in ('order3.json', 'order5.json', 'decomposition.json', 'provenance.json'):
        evidence.historical(out / name)
    a, b, saved, provenance = [read(out / name) for name in
                              ('order3.json', 'order5.json', 'decomposition.json', 'provenance.json')]
    require(provenance['supg_scale'] == 1.0, 'Changed power-audit SUPG scale')
    require(provenance['mesh_sha256'] == manifest['mesh_sha256'], 'Power-audit mesh identity')
    evidence.record(directory / 'nodes.csv', provenance['nodes_sha256'])
    for path, digest in provenance['scripts'].items():
        evidence.record(Path(path), digest)
    parts = {'quadrature_W': b['JE_W'] - a['JE_W'],
             'SUPG_weighted_W': a['phi_weighted_SUPG_W'],
             'reaction_minus_physical_boundary_W': a['reaction_terminal_power_W'] - b['all_exterior_power_W'],
             'free_residual_correction_W': -a['phi_weighted_free_residual_W']}
    gap = b['JE_W'] - b['all_exterior_power_W']
    require(gap == b['gap_W'] and gap / b['JE_W'] == b['gap_relative'], 'Saved power-gap arithmetic differs')
    recomputed = {'gap_W': gap, 'gap_relative': gap / b['JE_W'],
                  'electrode_gap_relative': (b['JE_W'] - b['contact_only_power_W']) / b['JE_W'],
                  'parts': parts, 'decomposition_remainder_W': gap - sum(parts.values()),
                  'SUPG_fraction_of_gap': parts['SUPG_weighted_W'] / gap,
                  'free_residual_max_A': a['free_residual_max_A']}
    require(recomputed == saved, 'Independent power decomposition arithmetic differs')
    native_relative = {}
    for power_key, native_key in [('JE_W', 'joule_W'), ('all_exterior_power_W', 'terminal_power_W')]:
        native_relative[power_key] = abs(b[power_key] - metrics[native_key]) / abs(metrics[native_key])
        require(native_relative[power_key] < 1e-8, 'Original independent/native power agreement gate failed')
    return {'method': 'FROZEN_INTEGRALS_HASH_CHECK_AND_EXACT_ARITHMETIC_REPLAY',
            'quadrature_rerun': False, 'native_relative_agreement': native_relative,
            'original_native_relative_agreement_limit': 1e-8,
            'saved_affine_self_check': provenance['verification'],
            'JE_W': b['JE_W'], 'all_exterior_power_W': b['all_exterior_power_W'],
            'contact_only_power_W': b['contact_only_power_W'], **recomputed}


def replay_run(directory, voltage, z, evidence):
    evidence.subtree(directory)
    manifest, execution = read(directory / 'manifest.json'), read(directory / 'execution.json')
    require(manifest['executable_sha256'] == EXE_SHA == execution['executable_sha256'], 'Wrong executable identity')
    require(Path(manifest['executable']).resolve() == EXE.resolve(), 'Wrong executable path')
    evidence.record(EXE, EXE_SHA)
    mesh = Path(manifest['mesh'])
    evidence.historical(mesh)
    require(evidence.record(mesh) == manifest['mesh_sha256'], 'Manifest mesh hash differs')
    evidence.record(directory / 'input.nml', execution['input_sha256'])
    require(manifest['model'] == MODEL and manifest['supg_scale'] == 1.0 and
            manifest['domain_symmetry'] == 'FULL' and manifest['log_weighted_quadrature_order'] == 3 and
            manifest['voltages'] == [voltage] and manifest['heat_source_unchanged'] is True,
            'Changed model, domain, voltage, integration, or stabilization')
    require(manifest['changed_input_keys'] == ALLOWED_INPUT_CHANGES, 'Unexpected input change list')
    template = Path(manifest['input_template'])
    evidence.record(template)
    a, b = template.read_text(encoding='utf-8-sig'), (directory / 'input.nml').read_text(encoding='utf-8-sig')
    for key in ALLOWED_INPUT_CHANGES:
        a, b = replace(a, key, "'NORMALIZED'"), replace(b, key, "'NORMALIZED'")
    require(a == b, 'Input changed outside the historical declared keys')
    seed = Path(manifest['seed_nodes'])
    evidence.historical(seed)
    evidence.record(seed, manifest['seed_sha256'])
    state, metrics = classify(directory, execution['exit_code'])
    require(state == 'SOLVE_FLUX_PASS_MESH_PENDING', 'Original solve/flux gates failed')
    require(read(directory / 'summary.json')['state'] == state, 'Historical classification differs')
    # The frozen full-domain writer omits this optional metric; the existing
    # mesh audit uses FULL as its legacy default, and the input/manifest above
    # are independently identity checked.
    require(metrics['heat_source'] == HEAT_SOURCE and metrics.get('domain_symmetry', 'FULL') == 'FULL',
            'Different heat/domain')
    fields = audit_nodes(directory / 'nodes.csv', mesh, voltage)
    require(fields == read(directory / 'field_audit.json') and fields['passed'], 'Node audit differs or fails')
    power = replay_power(directory, metrics, manifest, evidence)
    row = {'run': directory.name, 'directory': str(directory), 'voltage_V': voltage,
           'angular': 32, 'radial': 14, 'thickness': z, 'supg_scale': 1,
           'state': state, 'node_count': fields['node_count'],
           'gap_percent': power['gap_relative'] * 100,
           'electrode_gap_percent': power['electrode_gap_relative'] * 100,
           'JE_mW': power['JE_W'] * 1000, 'all_exterior_mW': power['all_exterior_power_W'] * 1000,
           'contact_only_mW': power['contact_only_power_W'] * 1000,
           'theta_max_K': metrics['theta_max_K'], 'ring_n': metrics['ring_abs_n_relative_mean'],
           'ring_p': metrics['ring_abs_p_relative_mean'], 'original_seconds': execution['seconds'],
           'original_maximum_residual': metrics['maximum_accepted_residual'],
           'original_maximum_update': metrics['maximum_accepted_update']}
    return row, {'input_only_declared_changes': True, 'original_model_heat_and_thresholds_retained': True,
                 'fresh_node_audit': fields, 'power_replay': power, 'metrics': metrics}


def failure_index(old, fine, evidence):
    """Reference historical unresolved scopes; this does not recompute their FEM states."""
    rows = []
    for voltage in (6, 8):
        angular = old['angular_mesh'][str(voltage)]
        af = old['fields'][f'field_A24_A32_V{voltage:02d}']
        radial = fine['mesh'][f'radial_V{voltage}']
        rf = fine['fields'][f'A32_R18_Z14_S1_V{voltage * 10:03d}']
        sf = fine['fields'][f'A32_R14_Z14_S0.5_V{voltage * 10:03d}']
        for scope, levels, mesh, field, source in [
            ('angular', [16, 24, 32], angular, af, OLD / 'analysis/summary.json'),
            ('radial', [14, 16, 18], radial, rf, DIRECTION / 'analysis/summary.json')]:
            for path, digest in field['hashes'].items():
                evidence.historical(Path(path))
                evidence.record(Path(path), digest)
            rows.append({'scope': scope, 'voltage_V': voltage, 'levels': levels,
                         'last_pair': levels[-2:], 'source_report': str(source),
                         'source_pointer': (f'angular_mesh.{voltage}' if scope == 'angular' else f'mesh.radial_V{voltage}'),
                         'scalar_ring_pass': mesh['all_mesh_pass'], 'local_field_pass': field['passed'],
                         'local_electron': field['comparisons']['eta_n_density_perturbation'],
                         'local_hole': field['comparisons']['eta_p_density_perturbation'],
                         'ring_hole': mesh['voltages'][0]['changes']['ring_abs_p_relative_mean'],
                         'whole_study_certificate': False})
        scale2 = FINE / f'A32_R14_Z14_S2_V{voltage * 10:03d}'
        evidence.subtree(scale2)
        ex = read(scale2 / 'execution.json')
        failed_state, _ = classify(scale2, ex['exit_code'])
        require(failed_state == 'FAILED_NO_COMPLETE_STATE', 'SUPG=2 historical failure no longer matches')
        # The old same-mesh comparator records its two node hashes as a list,
        # unlike the mapped-field comparator's file-to-hash dictionary.
        same_mesh_runs = [thickness_run(voltage, 14),
                          FINE / f'A32_R14_Z14_S0.5_V{voltage * 10:03d}']
        require(sf['runs'] == [p.name for p in same_mesh_runs] and len(sf['nodes_hashes']) == 2,
                'Unexpected same-mesh comparison identity')
        for directory, digest in zip(same_mesh_runs, sf['nodes_hashes']):
            evidence.historical(directory / 'nodes.csv')
            evidence.record(directory / 'nodes.csv', digest)
        rows.append({'scope': 'SUPG', 'voltage_V': voltage, 'grid': [32, 14, 14],
                     'required_scales': [0.5, 1, 2], 'compared_completed_scales': [0.5, 1],
                     'source_report': str(DIRECTION / 'analysis/summary.json'),
                     'source_pointer': f'SUPG_pairs.V{voltage}_S0.5',
                     'scalar_ring_pass': fine['SUPG'][str(voltage)]['covered_metrics_passed'],
                     'local_field_pass': sf['passed'],
                     'local_electron': sf['comparisons']['eta_n_density_perturbation'],
                     'local_hole': sf['comparisons']['eta_p_density_perturbation'],
                     'ring_hole': fine['SUPG_pairs'][f'V{voltage}_S0.5']['ring_abs_p_relative_mean'],
                     'scale2': {'directory': str(scale2), 'state': failed_state,
                                'seconds': ex['seconds'], 'timed_out': ex['timed_out'],
                                'has_accepted_final_state': False},
                     'whole_study_certificate': False})
    return {'method': 'HISTORICAL_REPORT_REFERENCE_WITH_FILE_HASH_AND_FAILURE_STATUS_RECHECK',
            'recomputed_angular_radial_or_SUPG_fields': False, 'records': rows,
            'limitations': ['The 8 V radial last-pair pass covers the radial direction only.',
                            'SUPG=0.5 cannot be selected merely to reduce the power discrepancy.',
                            'No data point or boundary node was excluded.',
                            'Initial OOC -9 attempts and the retired RESOURCE_RETRY remain historical failures, not accepted evidence.']}


def report(result):
    lines = ['## Material Passport', '',
             '- Mode: reproducibility validation of existing numerical evidence.',
             '- Status: VERIFIED for exact mesh/local-field audit replay; original FEM was not rerun.',
             '- Model: FULL J·E current-work diagnostic, ZnO NT; A32/R14, SUPG=1, volume order=3.',
             '- Data: frozen JE_SENSITIVITY_20261015 and JE_FINE_CROSS_20261016/GLOBAL_DIRECTIONAL.',
             '', '## 厚度方向复用结果', '',
             '6 V、8 V 的 Z14/16/18 三层标量/环带审计及相邻两层局部场审计重新执行通过，数值与历史报告完全一致。',
             '仅复用并重审厚度方向：不构成周向、径向、联合网格、SUPG、应力精度或整个 6–8 V 区间证书。',
             '热源、验收门限和旧结果均未更改；未执行原生求解、网格生成、映射或新的功率积分。',
             '', '功率差额 = (∫J·E dΩ − 全物理外边界电功率)/∫J·E dΩ；电极差额仅使用接触边界。',
             '独立功率积分原记录和源代码 SHA 已核对，并重算其差额与分解算术；没有新增功率差额 PASS 门限。',
             '', '| 电压 | Z 层数 | 节点数 | J·E /mW | 全外边界 /mW | 全边界差额/% | 电极差额/% | 原求解/s |',
             '|---|---|---|---|---|---|---|---|']
    for row in result['runs']:
        lines.append(f"| {row['voltage_V']} V | {row['thickness']} | {row['node_count']} | {row['JE_mW']:.9f} | "
                     f"{row['all_exterior_mW']:.9f} | {row['gap_percent']:.9f} | {row['electrode_gap_percent']:.9f} | {row['original_seconds']:.3f} |")
    lines += ['', '局部载流子最大差使用两场 exp(eta)−1 的全场最大幅值归一化；不是逐点除以局部浓度，未删除边界节点。',
              '', '| 电压 | 厚度相邻层 | 局部电子最大差/% | 局部空穴最大差/% | 电子环带变化/% | 空穴环带变化/% | 门限 |',
              '|---|---|---|---|---|---|---|']
    for voltage in (6, 8):
        for low, high, pair in [(14, 16, 'L0_L1'), (16, 18, 'L1_L2')]:
            field = result['fields'][f'V{voltage}_Z{low}_Z{high}']['comparisons']
            changes = result['mesh'][str(voltage)]['voltages'][0]['changes']
            lines.append(f"| {voltage} V | {low}→{high} | {field['eta_n_density_perturbation']['maximum_relative']*100:.6f} | "
                         f"{field['eta_p_density_perturbation']['maximum_relative']*100:.6f} | "
                         f"{changes['ring_abs_n_relative_mean'][pair]*100:.6f} | {changes['ring_abs_p_relative_mean'][pair]*100:.6f} | 各 <2% |")
    lines += ['', '三层标量/环带完整门限：功率、最大温升、环带温升和电子端口电流 1%；环带电场、热源和两载流子 2%；热源/电场 P95 5%。',
              '三层验收使用最后两层；原始两个相邻层的所有指标均保存在 summary.json 的 mesh/fields 中。',
              '厚度最大角节点间距 279.5→244.5625→217.388889 nm；A32、R14 和环带范围不变。网格哈希与原几何/Jacobian证据一一对应，未重新进行昂贵的积分几何检查。',
              '', '## 保留的未解决范围索引', '',
              '下表引用原报告并校验其来源文件，不宣称重新求解或重新完成这些方向的场比较。',
              '', '| 范围 | 电压 | 末对局部电子/% | 末对局部空穴/% | 空穴环带变化/% | 结果 |',
              '|---|---|---|---|---|---|']
    for row in result['outstanding_index']['records']:
        ring = row['ring_hole']
        value = ring['relative_span'] if row['scope'] == 'SUPG' else ring['L1_L2']
        status = '方向通过，非联合证书' if row['scalar_ring_pass'] and row['local_field_pass'] else '未通过'
        lines.append(f"| {row['scope']} | {row['voltage_V']} V | {100*row['local_electron']['maximum_relative']:.6f} | "
                     f"{100*row['local_hole']['maximum_relative']:.6f} | {100*value:.6f}（门限 {100*ring['limit']:g}%） | {status} |")
    lines += ['', 'SUPG=2 在 6/8 V 均未得到完整状态；0.5 对 1 的空穴差异仍超限。系数 0.5 较小的功率差额不能代替稳定化敏感性通过。',
              '径向 6 V 和周向 6/8 V 失败保留；已有径向 8 V 通过也只覆盖该独立方向。',
              '', '## 身份和证据范围', '',
              f"- 旧 exe SHA-256: {EXE_SHA}",
              f"- 40 个物理/装配/求解核心源码、双构建身份以及论文范围 10 专项/16 旧 CLI 已调用只读验收重验证。",
              f"- 共记录并在结束前复查 {result['evidence_count']} 个输入/结果/源码文件哈希，其中历史清单匹配 {result['historically_matched_count']} 个。",
              '- 验证模块与当前构建的差别按 verification_identity.json 原样披露；未宣称当前全部验证源码均编译进旧 exe。',
              '- 原 A2P/原整体验证 FAIL 保留；本报告不会升级其状态。',
              '- 冻结 stage4 runner 在审计前后均校验为既定 SHA，运行中的 A40 probe 未被修改或控制。',
              '- uncertainty: this deterministic directional study does not supply extrapolated continuum error or confidence intervals.',
              '- statistics applicability: no stochastic sample, hypothesis tests, p-values, multiple comparisons, causal identification, population inference, power/sample-size claims, or confidence-interval claims are made.',
              '- 完整来源、原门限、所有相邻场指标和未解决范围路径见 summary.json、outstanding_failures.json、provenance.json。']
    return '\n'.join(lines) + '\n'


def run(output):
    output = Path(output).resolve()
    require(output == OUTPUT.resolve() or output.is_relative_to(OUTPUT.resolve()),
            'Output must be a new directory at or below paper_scope/directional_reuse')
    require(not output.exists(), 'Refusing to overwrite an existing audit')
    evidence = Evidence()
    evidence.record(Path(__file__))
    evidence.record(RUNNER, RUNNER_SHA)
    gate_path = PAPER / 'verification_reuse/acceptance.json'
    gate = validate_acceptance(gate_path, expected_exe=EXE)
    evidence.record(gate_path)
    for path, digest in gate['evidence_file_hashes'].items():
        evidence.record(Path(path), digest)
    evidence.inventory(OLD / 'analysis/evidence_hashes.json', OLD)
    evidence.inventory(DIRECTION / 'analysis/artifact_hashes.json', FINE)
    # Current audit implementation must equal the implementation used for old results.
    for name in ('audit_hv_mesh.py', 'audit_hv_field_shape.py', 'run_hv_sweep.py',
                 'run_je_sensitivity_20261015.py'):
        snapshot = OLD / 'analysis' / name
        expected = evidence.historical(snapshot)
        evidence.record(ROOT / 'scripts' / name, expected)
    evidence.historical(FINE / 'build_manifest.json')
    evidence.record(EXE.with_name('build_manifest.json'), evidence.record(FINE / 'build_manifest.json'))
    for source in (OLD / 'analysis/summary.json', DIRECTION / 'analysis/summary.json',
                   DIRECTION / 'geometry_validation.json', DIRECTION / 'merged_contract.json'):
        evidence.historical(source)
    old, previous = read(OLD / 'analysis/summary.json'), read(DIRECTION / 'analysis/summary.json')
    geometry = read(DIRECTION / 'geometry_validation.json')
    require(geometry['baseline_byte_identical'] is True, 'Historical baseline geometry differs')
    profiles = {}
    for z, directory in [(14, OLD / 'mesh_A32'), (16, DIRECTION / 'mesh_thickness_16'),
                         (18, DIRECTION / 'mesh_thickness_18')]:
        evidence.subtree(directory)
        profile = read(directory / 'manifest.json')
        require(profile['counts'] == {'circumferential': 32, 'radial': 14, 'thickness': z}, 'Wrong mesh profile')
        require(profile['geometry_audit']['minimum_jacobian_3x3x3'] > 0, 'Original geometry check failed')
        profiles[str(z)] = profile
        if z > 14:
            check = next(c for c in geometry['checks'] if c['direction'] == 'thickness' and c['counts']['thickness'] == z)
            require(check['other_direction_coordinates_identical'] is True, 'Changed independent directions')
            generator = ROOT / 'scripts/generate_je_global_direction_mesh_20261016.py'
            evidence.record(generator, profile['generator_sha256'])
    rows, runs, meshes, fields = [], {}, {}, {}
    for voltage in (6, 8):
        levels = [thickness_run(voltage, z) for z in (14, 16, 18)]
        for z, directory in zip((14, 16, 18), levels):
            row, check = replay_run(directory, voltage, z, evidence)
            rows.append(row)
            runs[directory.name] = check
        meshes[str(voltage)] = mesh_audit(levels, [voltage])
        require(meshes[str(voltage)] == previous['mesh'][f'thickness_V{voltage}'], 'Recomputed mesh audit differs from original')
        require(meshes[str(voltage)]['all_mesh_pass'], 'Thickness mesh gate failed')
        for low, high, source, target in [(14, 16, levels[0], levels[1]), (16, 18, levels[1], levels[2])]:
            mapping = DIRECTION / f'transfer_{target.name}'
            evidence.subtree(mapping)
            value = field_audit(mapping, source, target)
            expected = {k: v for k, v in previous['fields'][target.name].items() if k != 'carrier_locations'}
            require(value == expected, 'Recomputed mapped field audit differs from original')
            require(value['passed'], 'Thickness local-field gate failed')
            fields[f'V{voltage}_Z{low}_Z{high}'] = value
        print(f'{voltage} V: original mesh and both local-field comparisons replay exactly', flush=True)
    unresolved = failure_index(old, previous, evidence)
    evidence.verify_unchanged()
    result = {'scope': 'JE_PAPER_THICKNESS_DIRECTION_ONLY_20261018',
              'method': 'REUSED_FROZEN_FEM_STATES_WITH_EXACT_GATE_REPLAY',
              'native_solver_started': False, 'fully_validated': False, 'thickness_direction_passed': True,
              'voltage_points': [6, 8], 'fixed_angular_radial': [32, 14], 'supg_scale': 1,
              'runs': rows, 'run_rechecks': runs, 'mesh': meshes, 'fields': fields,
              'geometry': geometry, 'mesh_profiles': profiles, 'outstanding_index': unresolved,
              'verification_gate': str(gate_path), 'verification_gate_sha256': evidence.hashes[str(gate_path.resolve())],
              'evidence_count': len(evidence.hashes),
              'historically_matched_count': sum(Path(p) in evidence.frozen for p in evidence.hashes),
              'frozen_runner_unchanged': True, 'created_utc': datetime.now(timezone.utc).isoformat()}
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'summary.json', result)
    save(output / 'outstanding_failures.json', unresolved)
    save(output / 'verification_identity.json', gate)
    save(output / 'provenance.json', {'sha256': evidence.hashes,
                                     'historical_inventory_sources': {p: evidence.reference[p] for p in evidence.hashes if p in evidence.reference},
                                     'all_read_files_unchanged_at_end': True,
                                     'new_native_runs': 0, 'old_files_written': 0})
    (output / 'REPORT.md').write_text(report(result), encoding='utf-8')
    with (output / 'thickness_results.csv').open('x', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save(output / 'artifact_hashes.json', {p.name: sha(p) for p in output.iterdir() if p.is_file()})
    print(json.dumps({'output': str(output), 'thickness_direction_passed': True,
                      'evidence_count': result['evidence_count'], 'fully_validated': False,
                      'native_solver_started': False}, ensure_ascii=True), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == '__main__':
    main()
