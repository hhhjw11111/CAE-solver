"""Append each finished reduced-production case to 逐例归档.md (user request 2026-10-22: 每完成一个归档到md).

Read-only observer: it never starts, stops or edits a solve. A case is archived once its final state
is on disk: summary.json for failures, summary.json plus completion.json for accepted states.
Cases the runner skipped without a directory are archived from the batch progress events.
Exits after the batch finished record is archived, or when the controller is gone and nothing is pending.
"""
from pathlib import Path
import argparse
import json
import time
import run_je_production_reduced_20261022 as reduced

base, r = reduced.base, reduced.r
OUTPUT = reduced.OUTPUT
LOG = OUTPUT/'逐例归档.md'
LEDGER = OUTPUT/'逐例归档_ledger.json'
ACCEPTED = (base.PASS, base.ZERO_PASS)
HEADER = ('| # | 归档时间 | 算例 | 组 | 模型 | a (μm) | V (V) | 状态 | 耗时 (min) | P 全域 (mW) | θmax (K) | '
          'Rth,eff (K/W) | In 全域 (mA) | 物理外边界热平衡 | 外边界功率差额 | metrics.txt SHA-256 |\n'
          '|---:|---|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|\n')


def fmt(x, spec):
    return format(x, spec) if isinstance(x, float) else '—'


def controller_alive(pid):
    import ctypes
    handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)
    if not handle: return False
    code = ctypes.c_ulong()
    ctypes.windll.kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
    ctypes.windll.kernel32.CloseHandle(handle)
    return code.value == 259


def row(index, job, state, directory):
    d = directory
    m = base.sweep.read_metrics(d/'metrics.txt') if (d/'metrics.txt').exists() else {}
    summary = r.read(d/'summary.json') if (d/'summary.json').exists() else {}
    factor = 2 if job['domain'] == 'Y_HALF' else 1
    p = factor*m['joule_W']*1e3 if isinstance(m.get('joule_W'), float) else None
    theta = m.get('theta_max_K')
    rth = theta/(p*1e-3) if isinstance(theta, float) and p else None
    current = factor*m['electron_right_A']*1e3 if isinstance(m.get('electron_right_A'), float) else None
    balance = m.get('physical_exterior_heat_balance_relative', m.get('heat_balance_relative'))
    gap = None
    if (d/'physical_power.json').exists(): gap = r.read(d/'physical_power.json')['physical_exterior_gap_relative']
    minutes = summary.get('elapsed_seconds')/60 if isinstance(summary.get('elapsed_seconds'), (int, float)) else None
    digest = r.sha(d/'metrics.txt')[:16] if (d/'metrics.txt').exists() else '—'
    return (f"| {index} | {time.strftime('%m-%d %H:%M')} | {job['name']} | {job['group']} | {job['model']} | {job['radius']} | "
            f"{job['voltage']:g} | {state} | {fmt(minutes, '.1f')} | {fmt(p, '.6f')} | {fmt(theta, '.6f')} | {fmt(rth, '.2f')} | "
            f"{fmt(current, '.6f')} | {fmt(balance, '.3e')} | {fmt(gap, '.3e')} | {digest} |\n")


def batch_files():
    identities = sorted(OUTPUT.glob('batch_*_identity.json'))
    if not identities: return None, None, []
    stamp = identities[-1].name.split('_')[1]
    progress = sorted(OUTPUT.glob(f'batch_{stamp}_progress_*.json'))
    events = r.read(progress[-1])['events'] if progress else []
    return stamp, OUTPUT/f'batch_{stamp}_finished.json', events


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controller-pid', type=int, required=True)
    p.add_argument('--poll-seconds', type=int, default=60)
    a = p.parse_args()
    ledger = r.read(LEDGER) if LEDGER.exists() else dict(logged=[])
    if 'states' not in ledger:
        # Ledger v1 kept names only; recover each logged state from its archived table row.
        ledger['states'] = {}
        if LOG.exists():
            for line in LOG.read_text(encoding='utf-8').splitlines():
                cells = [c.strip() for c in line.strip('|').split('|')]
                if len(cells) == 16 and cells[2] in ledger['logged']: ledger['states'][cells[2]] = cells[7]
    if not LOG.exists():
        LOG.write_text('# 30 例精简生产：逐例归档\n\n'
                       f'矩阵 SHA-256 `{r.sha(reduced.INPUTS/"case_matrix.csv")}`；运行器 `run_je_production_reduced_20261022.py`。\n'
                       '每例到达最终状态即追加一行；P 与 In 为半域积分×2 的全域值，θmax 不倍增；'
                       '热平衡与功率差额为相对值（功率差额仅 NT，诊断量无门限）。空穴量按 2026-10-22 决定不列。\n'
                       '状态 SOLVE_FLUX_PASS_MESH_PENDING 表示求解与独立通量通过，网格验收另见审计；任何非 PASS 状态均原样记录。\n\n'
                       + HEADER, encoding='utf-8')
    jobs = {j['name']: j for j in base.jobs()}
    final_pass = False
    while True:
        stamp, finished, events = batch_files()
        skipped = {e['name']: e['state'] for e in events if not (OUTPUT/e['name']/'summary.json').exists()}
        for name, job in jobs.items():
            d = OUTPUT/name
            if (d/'summary.json').exists():
                state = r.read(d/'summary.json')['state']
                if state in ACCEPTED and not (d/'completion.json').exists():
                    if controller_alive(a.controller_pid): continue  # post-processing still running
                    state += '（completion.json 缺失，后处理未完成）'
            elif name in skipped:
                state = skipped[name]
            else:
                continue
            if ledger['states'].get(name) == state: continue
            shown = state if name not in ledger['states'] else f'{state}（更新，原记为 {ledger["states"][name]}）'
            with LOG.open('a', encoding='utf-8') as f: f.write(row(len(ledger['logged'])+1, job, shown, d))
            ledger['logged'].append(name); ledger['states'][name] = state
            LEDGER.write_text(json.dumps(ledger, ensure_ascii=False, indent=1), encoding='utf-8')
        done = finished is not None and finished.exists()
        if (done or not controller_alive(a.controller_pid)) and not final_pass:
            final_pass = True  # one more archiving pass after the controller has stopped
            continue
        if final_pass:
            stamp, finished, events = batch_files()
            final = r.read(finished) if finished and finished.exists() else None
            with LOG.open('a', encoding='utf-8') as f:
                f.write(f"\n批次结束于 {time.strftime('%Y-%m-%d %H:%M')}：已归档 {len(ledger['states'])}/{len(jobs)} 例（表中共 {len(ledger['logged'])} 行，含状态更新行）；"
                        + (f"失败或跳过 {len(final['failed'])} 例：{', '.join(final['failed']) or '无'}。\n" if final
                           else '未找到 finished 记录，控制进程已退出，需人工核查。\n'))
            return 0
        time.sleep(a.poll_seconds)


if __name__ == '__main__':
    raise SystemExit(main())
