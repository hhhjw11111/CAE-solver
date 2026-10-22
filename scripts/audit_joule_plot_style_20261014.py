"""Separate template drift from jumps in the saved, unaveraged FEM fields."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def pixel_identity(path):
    with Image.open(path) as image:
        return {'size':list(image.size),
                'rgb_sha256':hashlib.sha256(image.convert('RGB').tobytes()).hexdigest()}


def field_audit(path, stride=1):
    with np.load(path) as cache:
        data = {key:cache[key] for key in cache.files}
    xyz = data['xyz']
    if stride > 1:
        # Match the historical four subdivisions, retaining both sides of seams.
        selected = np.arange(len(xyz)).reshape(-1, 17, 17)[:, ::stride, ::stride].ravel()
    else:
        selected = np.arange(len(xyz))
    points, inverse, count = np.unique(np.round(xyz[selected]*1e6, 9),
                                     axis=0, return_inverse=True, return_counts=True)
    result = {'source':str(path.relative_to(ROOT)), 'sha256':digest(path),
              'source_samples':len(xyz), 'source_triangles':len(data['triangles']),
              'audited_samples':len(selected),
              'shared_coordinate_locations':int((count>1).sum()),
              'coordinate_rounding_um':1e-9, 'fields':{}}
    for key in ['Q', 'Q_drift_same_field', 'Q_diffusion_work', 'phi']:
        if key not in data:
            continue
        values = data[key][selected]
        lo = np.full(len(points), np.inf)
        hi = np.full(len(points), -np.inf)
        np.minimum.at(lo, inverse, values)
        np.maximum.at(hi, inverse, values)
        jumps = (hi-lo)[count>1]
        worst = int(np.argmax(hi-lo))
        peak = float(np.max(abs(values)))
        result['fields'][key] = {
            'min':float(values.min()), 'max':float(values.max()),
            'max_shared_coordinate_jump':float(jumps.max()),
            'p95_shared_coordinate_jump':float(np.quantile(jumps, .95)),
            'max_jump_over_absolute_peak_percent':float(jumps.max()/peak*100),
            'worst_coordinate_um':points[worst].tolist(),
            'worst_two_sided_range':[float(lo[worst]), float(hi[worst])],
        }
    if 'Q_diffusion_work' in data:
        result['max_diffusion_work_over_Q_peak_percent'] = float(
            np.max(abs(data['Q_diffusion_work']))/np.max(abs(data['Q']))*100)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-image', type=Path)
    parser.add_argument('--current-image', type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    report_path = args.output/'qa/style_and_field_audit.json'
    if report_path.exists():
        raise FileExistsError(report_path)
    reference = ROOT/'output/F3_焦耳热_下排再次放大_20260909/comparisons/compare_Q_midplane.png'
    original = ROOT/'output/JE_高电压_云图与放大曲面_20261014/figures/01_joule_mid.png'
    skill = Path.home()/'.codex/skills/analytic-solution-comparison'
    backup = Path.home()/'.codex/skills.backups/20261013-official-optimization/original/analytic-solution-comparison'
    old_skill = (backup/'SKILL.md').read_text(encoding='utf-8')
    new_skill = (skill/'SKILL.md').read_text(encoding='utf-8')
    report = {
        'scope':'Rendering repair and saved-field diagnosis only; no PDE or field changes.',
        'skill_comparison':{
            'backup':str(backup), 'current':str(skill),
            'body_identical':old_skill.split('---',2)[2]==new_skill.split('---',2)[2],
            'style_contract_identical':digest(backup/'references/style-contract.md')==digest(skill/'references/style-contract.md'),
            'old_skill_sha256':digest(backup/'SKILL.md'),
            'new_skill_sha256':digest(skill/'SKILL.md'),
            'style_contract_sha256':digest(skill/'references/style-contract.md'),
        },
        'reference_image':{'path':str(reference), **pixel_identity(reference)},
        'original_image':{'path':str(original), **pixel_identity(original)},
        'comparison':{
            'old':'3 V, radii 2/3/4 um, drift-conductivity heat, Q color limits 0..17 (1e12 W/m3)',
            'new':'6/7/8 V, radius 3 um, total-current work diagnostic, Q color limits 0..110 (1e12 W/m3)',
            'rendering_same':'Element-local triangulation, Gouraud 2D, unshaded trisurf 3D, same 256-level palette, common linear normalization.',
            'rendering_drift':'New script copied palette only; hardcoded axes, smaller ticks, unit font, tick padding and uncropped canvas replaced the saved layout.',
            'repair':'Read full frozen larger-surfaces layout, views, palette, fonts and export settings; preserve current data and common color limits.',
        },
        'field_audits':[],
    }
    for supplied, label in [(args.reference_image, 'reference_image'), (args.current_image, 'original_image')]:
        if supplied:
            report[label]['attachment_matches_pixels'] = pixel_identity(supplied)==pixel_identity(Path(report[label]['path']))
            assert report[label]['attachment_matches_pixels']
    report['field_audits'].append(field_audit(
        ROOT/'output/F3_孔洞半径扫描_20260909_v3/source_data/R03_NT_midplane.npz'))
    for v in [6,7,8]:
        path = args.output/'source_data'/f'V{v}_midplane.npz'
        report['field_audits'].append(field_audit(path))
        report['field_audits'].append(field_audit(path, stride=4))
    report['interpretation'] = (
        'Q jumps exist before color mapping, while potential is continuous to roundoff. '
        'The current field also has comparable jumps in drift-only heat evaluated on the same state; '
        'the small diffusion-work term alone does not explain these midplane seams. '
        'Increasing plotting subdivisions does not enforce inter-element derivative continuity. '
        'No averaging, clipping, interpolating across the hole or change to solver evidence is allowed by this repair. '
        'These checks diagnose visualization of saved fields, not mesh convergence or physical validity.')
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    shutil.copy2(__file__, args.output/Path(__file__).name)
    print(json.dumps(report['skill_comparison']))
    for item in report['field_audits']:
        print(Path(item['source']).name, item['audited_samples'],
              'Q jump/peak (%) =', item['fields']['Q']['max_jump_over_absolute_peak_percent'])


if __name__ == '__main__':
    main()
