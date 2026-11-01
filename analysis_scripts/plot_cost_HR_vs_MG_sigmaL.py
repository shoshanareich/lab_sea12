#!/usr/bin/env python3
"""
Cost reduction, HR-only vs multi-grid, 10-day SWOT runs with the sigmaL
(representation-error) obsfit file.  Poster figure.

y-axis: <z^2> = J_SWOT / N_SWOT, the swot_obsfit term from costfunction0*
divided by its own sample count (27960), as in plot_optim_results_HR_vs_MG.py.

Best-so-far, without marking accepted/rejected trials: at each driver
iteration the curve shows <z^2> of the evaluation with the lowest TOTAL cost
fc so far -- m1qn3's acceptance criterion.  Taking the running minimum of
<z^2> itself would credit rejected trials whose SWOT term happens to drop
while fc rises (MG it0011: J_SWOT 135226 < 137133 at it0010, fc 146327 >
146305; rejected).

Usage:  python3 plot_cost_HR_vs_MG_sigmaL.py
"""

import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
OUTDIR   = '/home/shoshi/MITgcm_c69j/lab_sea12/figs/'
OBSNAME  = 'swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL'

RUNS = [('assim_swot_MG_10days_sigmaL', 'Multi-grid', 'royalblue'),
        ('assim_swot_HR_10days_sigmaL', 'HR-only',    '#1baf7a')]

# Figure 1 only: MG continued past iteration 10 with the LR adjoint viscosity
# switched x1 -> x0.1 (submit_MG_10days_sigmaL_x1tox0.1.sh, cold m1qn3 restart
# from x1's accepted it0010).  Viscosity enters only the LR adjoint, so fc
# (HR forward) is comparable across the switch; the continuation's it0010 is
# bit-identical to x1's.  x1's own it0011 is a rejected trial and is dropped.
MG_X1, MG_X01 = 'assim_swot_MG_10days_sigmaL', 'assim_swot_MG_10days_sigmaL_x1tox0.1'
VISC_SWITCH = 10


def read_cost(run, itmin=0, itmax=None):
    """Per driver iteration: (iter, fc, J_SWOT, N_SWOT) from costfunction0*."""
    rows = []
    for fp in sorted(glob.glob(os.path.join(RUNSROOT, run,
                                            'run_adhi_it0*', 'costfunction0*'))):
        it = int(re.search(r'it(\d+)', fp).group(1))
        if it < itmin or (itmax is not None and it > itmax):
            continue
        fc = swot = num = np.nan
        for ln in open(fp):
            ln = ln.replace('D+', 'E+').replace('D-', 'E-')
            if ln.strip().startswith('fc ='):
                fc = float(ln.split('=')[1].split()[0])
            elif ln.split() and ln.split()[0] == OBSNAME and 'obs_val' in ln:
                swot, num = (float(x) for x in ln.split('=')[1].split()[:2])
        rows.append((it, fc, swot, num))
    rows = np.array(sorted(rows))
    assert np.all(np.diff(rows[:, 0]) == 1), f'{run}: gap in iterations'
    assert np.all(np.isfinite(rows)), f'{run}: unparsed costfunction'
    return rows


def read_mg_continued():
    """MG x1 for 0..VISC_SWITCH, then the x0.1 continuation."""
    a = read_cost(MG_X1, 0, VISC_SWITCH)
    b = read_cost(MG_X01, VISC_SWITCH)
    assert a[-1, 1] == b[0, 1], 'continuation does not start from the x1 iterate'
    rows = np.vstack([a, b[1:]])
    assert np.all(np.diff(rows[:, 0]) == 1)
    return rows


def best_so_far(rows):
    """<z^2> of the lowest-fc evaluation up to and including each iteration."""
    z2 = rows[:, 2] / rows[:, 3]
    best = np.minimum.accumulate(rows[:, 1])
    idx = [int(np.flatnonzero(rows[:k + 1, 1] == best[k])[0])
           for k in range(len(rows))]
    return rows[:, 0].astype(int), z2[idx]


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    plt.rcParams.update({'font.size': 13, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    plot_z2(visc_switch=False)
    plot_z2(visc_switch=True)
    plot_percent_decrease()


def plot_z2(visc_switch):
    """<z^2> best-so-far; with visc_switch, MG continues to 20 on x0.1."""
    fig, ax = plt.subplots(figsize=(8.5, 4.6) if visc_switch else (7.5, 4.6))

    itmax = 0
    for run, label, color in RUNS:
        if run == MG_X1 and visc_switch:
            rows = read_mg_continued()
        else:
            rows = read_cost(run)
        it, z2 = best_so_far(rows)
        itmax = max(itmax, int(it.max()))
        ax.plot(it, z2, '-o', color=color, lw=2.4, ms=5.5, label=label)
        print(f'{label:11s} ' + '  '.join(f'{i}:{v:.2f}' for i, v in zip(it, z2)))
        ax.annotate(f'{z2[-1]:.1f}', (it[-1], z2[-1]), xytext=(8, 0),
                    textcoords='offset points', va='center', color=color,
                    fontsize=12, fontweight='bold')

    if visc_switch:
        ax.axvline(VISC_SWITCH, color='0.35', ls='--', lw=1.2, zorder=0)
        ax.text(VISC_SWITCH + 0.2, 0.97,
                'MG adjoint viscosity\n' + r'$\times$1 $\rightarrow$ $\times$0.1',
                transform=ax.get_xaxis_transform(), va='top', ha='left',
                fontsize=11, color='0.25')

    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'SWOT misfit  $\langle z^2 \rangle = J_\mathrm{SWOT}/N_\mathrm{SWOT}$')
    ax.set_xticks(range(0, itmax + 1))
    ax.set_xlim(-0.3, itmax + 1.0)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()

    name = 'cost_reduction_HR_vs_MG_sigmaL' + ('_viscswitch' if visc_switch else '')
    out = os.path.join(OUTDIR, name + '.pdf')
    fig.savefig(out)
    print('wrote', out)


def plot_percent_decrease():
    """Same best-so-far iterates, as % decrease of J_SWOT from iteration 0."""
    fig, ax = plt.subplots(figsize=(7.5, 4.6))

    itmax = 0
    for run, label, color in RUNS:
        it, z2 = best_so_far(read_cost(run))
        pct = 100.0 * (1.0 - z2 / z2[0])
        itmax = max(itmax, int(it.max()))
        ax.plot(it, pct, '-o', color=color, lw=2.4, ms=5.5, label=label)
        print(f'{label:11s} ' + '  '.join(f'{i}:{v:.1f}%' for i, v in zip(it, pct)))
        ax.annotate(f'{pct[-1]:.0f}%', (it[-1], pct[-1]), xytext=(8, 0),
                    textcoords='offset points', va='center', color=color,
                    fontsize=12, fontweight='bold')

    ax.set_xlabel('Iteration')
    ax.set_ylabel('SWOT cost decrease from iteration 0 (%)')
    ax.set_xticks(range(0, itmax + 1))
    ax.set_xlim(-0.3, itmax + 1.0)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, loc='lower right')
    fig.tight_layout()

    out = os.path.join(OUTDIR, 'cost_percent_decrease_HR_vs_MG_sigmaL.pdf')
    fig.savefig(out)
    print('wrote', out)


if __name__ == '__main__':
    main()
