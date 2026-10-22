"""Quantify a completed F1 byte-regression failure without changing its verdict.

Use the original input field scales, rather than dividing by nearly zero fields.
This diagnostic is not a replacement acceptance threshold or a mesh certificate.
"""
import argparse
import json
from pathlib import Path
import re
import numpy as np
from run_hv_sweep import read_metrics, sha


def compare(root):
    strict = json.loads((root / 'report.json').read_text(encoding='utf-8-sig'))
    case = Path(strict['case'])
    text = case.read_text(encoding='utf-8-sig')
    match = re.search(r'field_scales\s*=\s*([^\n/]+)', text, re.I)
    if match is None:
        raise ValueError('Explicit original field_scales are required')
    scales = [float(v.strip().replace('D','e').replace('d','e'))
              for v in match.group(1).strip().rstrip(',').split(',')]
    fields = ['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']
    if len(scales) != 7 or not all(np.isfinite(scales)) or min(scales) <= 0:
        raise ValueError('Exactly seven positive finite field scales are required')
    files = [case, root/'report.json', root/'executions.json']
    data, metrics = [], []
    for label in ['baseline','candidate']:
        directory = root / label
        nodes = directory / 'nodes.csv'
        table = np.genfromtxt(nodes, delimiter=',', names=True)
        if any(not np.isfinite(table[k]).all() for k in table.dtype.names):
            raise ValueError('Nonfinite production field')
        data.append(table)
        metrics.append(read_metrics(directory/'metrics.txt'))
        files.extend([nodes, directory/'metrics.txt', directory/'field_samples.csv'])
    same_topology = len(data[0]) == len(data[1]) and all(
        np.array_equal(data[0][k], data[1][k]) for k in ['node','x_m','y_m','z_m'])
    if not same_topology:
        raise ValueError('Different node order or coordinates cannot be compared pointwise')
    differences = {}
    for key, scale in zip(fields, scales):
        maximum = float(np.max(abs(data[0][key] - data[1][key])))
        differences[key] = dict(max_absolute=maximum, frozen_input_scale=scale,
                                scaled_by_frozen_input=maximum/scale)
    metric_differences = {}
    for key in ['theta_max_K','joule_W','electron_right_A','heat_balance_relative',
                'electron_balance_relative','hole_balance_diagnostic']:
        values = [m[key] for m in metrics]
        if not all(isinstance(v,(int,float)) and np.isfinite(v) for v in values):
            raise ValueError('Nonfinite or missing native metric')
        metric_differences[key] = dict(baseline=values[0], candidate=values[1],
                                       absolute_difference=abs(values[1]-values[0]))
    return dict(strict_byte_identical=all(v.get('byte_identical') is True for v in strict['files'].values()),
                original_strict_verdict=strict['passed'], same_topology=True,
                original_case_acceptance=[m.get('case_acceptance') for m in metrics],
                fields=differences, metrics=metric_differences,
                maximum_scaled_node_difference=max(v['scaled_by_frozen_input'] for v in differences.values()),
                file_hashes={str(p.absolute()):sha(p) for p in files},
                interpretation='Numeric differences and native acceptance are reported separately; '
                               'the original byte verdict is unchanged. This is one production case.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.production_run)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps({k:result[k] for k in ['strict_byte_identical','original_case_acceptance',
                                          'maximum_scaled_node_difference']}))


if __name__ == '__main__':
    main()
