#!/usr/bin/env python3
"""
mg_19 -- control adjustment by field, HR-only vs multi-grid, ctrlWyearly runs.

Same figure as mg_07 / mg_08 / mg_15, but the two columns are the two
configurations being compared here:

    assim_swot_HR_10days_ctrlWyearly   single-resolution 1/12 deg
    assim_swot_MG_10days_ctrlWyearly   1/12 deg forward + 1/3 deg adjoint

and it runs to it0020 instead of it0010.

Top row     rms control  sqrt(sum v^2 / N)  per field, over HR wet points.  v is
            the raw dimensionless control; under WC01 the penalty is sum(v^2)
            with no weight, so v is in units of the prior sigma and the fields
            are directly comparable.
Bottom row  each field's share of sum(v^2).

Filled marker = accepted m1qn3 iterate, open marker = rejected line-search
trial.  Both runs warm-start (coldstart=.FALSE.) so OPTIM/m1qn3_output.txt holds
the full history and the labels are recoverable.

CAVEAT for the MG column.  xx_* in run_adhi_it####/ are the HR fields the
forward model actually consumed, i.e. the injected LR control, so the rms is
directly comparable to the HR column.  But the MG penalty the optimizer sees is
sum over LR points, not HR -- the bottom-row shares for MG are the HR-space
proxy for it.  Same convention as mg_08.

Usage:  python3 plot_ctrl_by_field_HR_vs_MG.py [outfile.png] [--itermax 20]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
from plot_ctrl_by_field import MAIN, COL, collect
from plot_mg_diagnostics import _style, _save, SURFACE, TEXT_PRIMARY, TEXT_SECOND

RUNS   = ['assim_swot_HR_10days_ctrlWyearly',
          'assim_swot_MG_10days_ctrlWyearly']
LABELS = ['HR-only  (1/12$^\\circ$ fwd + 1/12$^\\circ$ adj)',
          'Multi-grid  (1/12$^\\circ$ fwd + 1/3$^\\circ$ adj)']
ITERMAX = 20


def main(outfile, runs, labels, itermax):
    data = [collect(r, itermax=itermax) for r in runs]

    for r, (its, rms, pen, kinds, _) in zip(runs, data):
        tot = np.array([sum(pen[f][k] for f in MAIN + ['other'])
                        for k in range(len(its))])
        print('\n%s' % r)
        print('  %-5s %-11s %-9s %-9s %-9s %-9s %s'
              % ('iter', 'kind', 'theta', 'salt', 'uwind', 'vwind', 'sum v^2'))
        for k, i in enumerate(its):
            print('  %-5d %-11s %-9.5f %-9.5f %-9.5f %-9.5f %.4e'
                  % (i, kinds[k], rms['theta'][k], rms['salt'][k],
                     rms['uwind'][k], rms['vwind'][k], tot[k]))

    fig, ax = plt.subplots(2, len(runs), figsize=(7.6 * len(runs), 8.4),
                           facecolor=SURFACE, sharey='row', sharex=True)
    ymax = 1.08 * max(max(max(d[1][f]) for f in MAIN) for d in data)

    for j, ((its, rms, pen, kinds, classified), lab) in enumerate(zip(data, labels)):
        a = ax[0, j]
        for f in MAIN:
            a.plot(its, rms[f], color=COL[f], lw=2.0, zorder=3)
            acc = [k for k in range(len(its)) if kinds[k] != 'rejected']
            rej = [k for k in range(len(its)) if kinds[k] == 'rejected'] if classified else []
            a.plot([its[k] for k in acc], [rms[f][k] for k in acc], 'o',
                   color=COL[f], ms=6, mec=SURFACE, mew=1.0, zorder=4)
            a.plot([its[k] for k in rej], [rms[f][k] for k in rej], 'o',
                   mfc=SURFACE, mec=COL[f], mew=1.8, ms=6, zorder=4)
        a.set_ylim(0, ymax); a.set_xlim(-0.5, itermax + 0.5)
        a.set_xticks(range(0, itermax + 1, 2))
        _style(a, title=lab if classified else
               lab + '\n(m1qn3 log is a partial history — accepted/rejected not recoverable)',
               ylabel=r'rms control  $\sqrt{\Sigma v^2/N}$   [prior $\sigma$]'
               if j == 0 else None)
        # direct labels: <= 4 series, so identity is never colour-alone.
        # theta and salt end within 0.5% of each other, so nudge the labels
        # apart -- a collided pair reads as one word and defeats the point.
        gap = 0.030 * ymax
        order = sorted(MAIN, key=lambda f: rms[f][-1])
        ypos, prev = {}, -1e30
        for f in order:
            y = max(rms[f][-1], prev + gap)
            ypos[f], prev = y, y
        for f in MAIN:
            a.annotate(f, (its[-1], ypos[f]), xytext=(6, 0),
                       textcoords='offset points', color=COL[f], fontsize=9,
                       va='center', fontweight='medium', annotation_clip=False)
        if j == 0:
            h = [plt.Line2D([], [], color=COL[f], lw=2.0, marker='o', ms=6) for f in MAIN]
            h += [plt.Line2D([], [], color=TEXT_SECOND, ls='', marker='o', ms=6),
                  plt.Line2D([], [], mfc=SURFACE, mec=TEXT_SECOND, mew=1.8, ls='',
                             marker='o', ms=6)]
            a.legend(h, MAIN + ['accepted iterate', 'rejected trial'], fontsize=9,
                     ncol=2, frameon=False, loc='upper left')

        a = ax[1, j]
        tot = np.array([sum(pen[f][k] for f in MAIN + ['other'])
                        for k in range(len(its))])
        ok = tot > 0                       # iteration 0 has v = 0, share undefined
        share = np.array([100 * np.array(pen[f])[ok] / tot[ok] for f in MAIN + ['other']])
        a.stackplot(np.array(its)[ok], share,
                    colors=[COL[f] for f in MAIN + ['other']], alpha=0.8,
                    edgecolor=SURFACE, linewidth=0.8)
        a.set_ylim(0, 100); a.set_xlim(-0.5, itermax + 0.5)
        a.set_xticks(range(0, itermax + 1, 2))
        _style(a, xlabel='outer iteration',
               ylabel=r'share of total penalty $\Sigma v^2$   [%]' if j == 0 else None)
        if j == 0:
            hh = [plt.Rectangle((0, 0), 1, 1, fc=COL[f], alpha=0.8)
                  for f in MAIN + ['other']]
            a.legend(hh, MAIN + ['other (6 fields)'], fontsize=9, ncol=5,
                     frameon=False, loc='lower center')

    fig.suptitle('Control adjustment by field, 10-day SWOT window, yearly control priors: '
                 'HR-only vs multi-grid', color=TEXT_PRIMARY, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    _save(fig, outfile)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?',
                    default='mg_19_ctrl_by_field_HR_vs_MG_ctrlWyearly.png')
    ap.add_argument('--itermax', type=int, default=ITERMAX)
    ap.add_argument('--runs', nargs=2, default=RUNS)
    ap.add_argument('--labels', nargs=2, default=LABELS)
    a = ap.parse_args()
    main(a.outfile, a.runs, a.labels, a.itermax)
