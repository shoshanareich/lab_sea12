#!/usr/bin/env python3
"""
Control adjustment by field, per outer iteration -- the mg_07/mg_08 figure,
regenerated for an arbitrary pair of runs.

Top row     rms control  sqrt(sum v^2 / N)  per field.  v is the raw
            dimensionless control the optimizer works in; under WC01 the
            penalty is sum(v^2) with no weight, so v is already in units of the
            prior sigma and the fields are directly comparable.
Bottom row  each field's share of the total penalty sum(v^2).

Filled marker = accepted m1qn3 iterate, open marker = rejected trial.  m1qn3
evaluates the cost at trial points as well as at iterates, so not every it####
directory is an iterate -- read the marker before comparing two runs at the
same it number.

N counts wet points only (hFacC > 0), because ctrl_cost_gen.F:317 gates the
penalty on maskC.  2D time-varying controls use the surface mask across all
their records.

theta, salt, uwind and vwind carry 99.97% of the penalty; the other six
controls (atemp, precip, swdown, lwdown, aqh, diffkr) are plotted as one
'other' band and are invisible at this scale.

Usage:  python3 plot_ctrl_by_field.py [outfile.png] [--runs A B] [--labels A B]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds
from plot_mg_diagnostics import classify_iterates, read_m1qn3, _style, _save, SURFACE, \
                                TEXT_PRIMARY, TEXT_SECOND, RUNSROOT, HR_GRID

# ============================== configuration ==============================
MAIN   = ['theta', 'salt', 'uwind', 'vwind']            # 99.97% of the penalty
OTHER  = ['atemp', 'precip', 'swdown', 'lwdown', 'aqh', 'diffkr']
FLD3D  = {'theta', 'salt', 'diffkr'}
# validated all-pairs in light mode: worst CVD dE 9.1, worst normal dE 22.9
COL    = {'theta': '#dc143c', 'salt': '#4169e1',
          'uwind': '#1baf7a', 'vwind': '#eda100', 'other': '#8a8880'}
ITERMAX = 10
# ===========================================================================


def load_ctrl(rundir, field, it, hf):
    """Raw control v, masked to wet points.  Returns the 1-D vector of values."""
    a = np.asarray(rdmds(os.path.join(rundir, 'xx_' + field), it), dtype=np.float64)
    if field in FLD3D:
        return a[hf > 0]
    a = a.reshape(-1, hf.shape[1], hf.shape[2])
    return a[np.broadcast_to(hf[0] > 0, a.shape)]


def collect(run, subpat='run_adhi_it%04d', itermax=ITERMAX):
    """Per-iteration rms and penalty by field, plus accepted/rejected labels."""
    root = os.path.join(RUNSROOT, run)
    hf = np.asarray(rdmds(HR_GRID + 'hFacC'), dtype=np.float64)
    # optim.x rewrites OPTIM/m1qn3_output.txt on every outer iteration.  With
    # coldstart=.TRUE. m1qn3 does not replay its history from OPWARM, so the
    # surviving log covers only the LAST call and accepted/rejected cannot be
    # recovered.  Detect that instead of mislabelling every point 'rejected'.
    classified = True
    try:
        iterates, _ = read_m1qn3(os.path.join(root, 'OPTIM'))
        if len(iterates) < 0.5 * itermax:
            raise ValueError('m1qn3 log holds only %d iterate(s) for %d '
                             'directories -- partial history' % (len(iterates), itermax))
        kinds = {i: k for i, _, k in
                 classify_iterates(root, subpat, os.path.join(root, 'OPTIM'), itermax)}
    except Exception as e:
        print('  %s: iterates not classifiable (%s)' % (run, e))
        kinds = {i: 'unknown' for i in range(itermax + 1)}
        classified = False

    rms, pen = {f: [] for f in MAIN + ['other']}, {f: [] for f in MAIN + ['other']}
    its = []
    for i in range(itermax + 1):
        d = os.path.join(root, subpat % i)
        if not os.path.isdir(d):
            break
        its.append(i)
        oth = 0.0
        for f in MAIN + OTHER:
            v = load_ctrl(d, f, i, hf)
            s = float((v ** 2).sum())
            if f in MAIN:
                rms[f].append(np.sqrt(s / v.size)); pen[f].append(s)
            else:
                oth += s
        rms['other'].append(np.nan); pen['other'].append(oth)
    return its, rms, pen, [kinds.get(i, 'unknown') for i in its], classified


def main(outfile, runs, labels):
    data = [collect(r) for r in runs]
    for r, (its, rms, pen, kinds, _) in zip(runs, data):
        tot = np.array([sum(pen[f][k] for f in MAIN + ['other']) for k in range(len(its))])
        print('\n%s' % r)
        print('  %-5s %-11s %-9s %-9s %-9s %-9s %s'
              % ('iter', 'kind', 'theta', 'salt', 'uwind', 'vwind', 'sum v^2'))
        for k, i in enumerate(its):
            print('  %-5d %-11s %-9.5f %-9.5f %-9.5f %-9.5f %.4e'
                  % (i, kinds[k], rms['theta'][k], rms['salt'][k],
                     rms['uwind'][k], rms['vwind'][k], tot[k]))

    fig, ax = plt.subplots(2, len(runs), figsize=(7.0 * len(runs), 8.0),
                           facecolor=SURFACE, sharey='row')
    ax = np.atleast_2d(ax.T).T if len(runs) > 1 else ax.reshape(2, 1)

    ymax = 1.08 * max(max(max(d[1][f]) for f in MAIN) for d in data)
    for j, ((its, rms, pen, kinds, classified), lab) in enumerate(zip(data, labels)):
        a = ax[0, j]
        for f in MAIN:
            a.plot(its, rms[f], color=COL[f], lw=1.8, zorder=3,
                   label=f if j == 0 else None)
            acc = [k for k in range(len(its)) if kinds[k] != 'rejected']
            rej = [k for k in range(len(its)) if kinds[k] == 'rejected'] if classified else []
            a.plot([its[k] for k in acc], [rms[f][k] for k in acc], 'o',
                   color=COL[f], ms=6, zorder=4)
            a.plot([its[k] for k in rej], [rms[f][k] for k in rej], 'o',
                   mfc=SURFACE, mec=COL[f], mew=1.6, ms=6, zorder=4)
        a.set_ylim(0, ymax); a.set_xlim(-0.4, ITERMAX + 0.4)
        _style(a, title=lab if classified else lab + '\n(m1qn3 log is a partial history \u2014 accepted/rejected not recoverable)', ylabel=r'rms control  $\sqrt{\Sigma v^2/N}$   [prior $\sigma$]'
               if j == 0 else None)
        if j == 0:
            h = [plt.Line2D([], [], color=COL[f], lw=1.8, marker='o', ms=6) for f in MAIN]
            h += [plt.Line2D([], [], color=TEXT_SECOND, ls='', marker='o', ms=6),
                  plt.Line2D([], [], mfc=SURFACE, mec=TEXT_SECOND, mew=1.6, ls='',
                             marker='o', ms=6)]
            a.legend(h, MAIN + ['accepted', 'rejected trial'], fontsize=9, ncol=2,
                     frameon=False, loc='lower right')

        a = ax[1, j]
        tot = np.array([sum(pen[f][k] for f in MAIN + ['other']) for k in range(len(its))])
        ok = tot > 0                       # iteration 0 has v = 0, share undefined
        share = np.array([100 * np.array(pen[f])[ok] / tot[ok] for f in MAIN + ['other']])
        a.stackplot(np.array(its)[ok], share,
                    colors=[COL[f] for f in MAIN + ['other']], alpha=0.75)
        a.set_ylim(0, 100); a.set_xlim(-0.4, ITERMAX + 0.4)
        _style(a, xlabel='iteration',
               ylabel=r'share of total penalty $\Sigma v^2$   [%]' if j == 0 else None)

    fig.suptitle('Control adjustment by field, HR 1/12$^\\circ$, v2 prior: '
                 'warm start vs cold start', color=TEXT_PRIMARY, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    _save(fig, outfile)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?', default='ctrl_by_field.png')
    ap.add_argument('--runs', nargs=2,
                    default=['assim_swot_HR_10days_newCtrlW_v2',
                             'assim_swot_HR_10days_newCtrlW_v2_coldstart'])
    ap.add_argument('--labels', nargs=2,
                    default=['HR v2  (warm start, coldstart=.FALSE.)',
                             'HR v2  (cold start, coldstart=.TRUE.)'])
    a = ap.parse_args()
    main(a.outfile, a.runs, a.labels)
