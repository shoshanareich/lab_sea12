#!/usr/bin/env python3
"""
adxx_theta gradient, HR-only assimilation vs multi-grid, iteration by iteration.

    row 0   HR-only run's HR adjoint gradient, RESTRICTED to the LR grid
    row 1   MG run's LR adjoint gradient, native LR

Both rows are therefore on the same 120x96 grid and share a colour scale within
each column, so a column is a like-for-like comparison of the same quantity.

The restriction is the wet-area-weighted block mean

    (R g)_c = sum_block(hFacC*rA * g) / sum_block(hFacC*rA)

which is the exact adjoint of piecewise-constant injection under the wet-area
inner product -- the same operator make_cost_cp_v2.py's 'area' mode and
make_lr_weights_from_hr.py use.  A plain 4x4 mean is wrong at coastlines
(varying wet-cell count) and HR cell area varies 3.3x across this domain.

Colour scale: symmetric, clipped at a robust percentile of |g| over wet cells
pooled across both panels of the column (default 99th), with the colourbar
extended at both ends so saturated cells stay visible as saturated.  A raw
max-abs scale is dominated by a handful of coastal spikes and renders every
other cell white.

Usage:
    python3 plot_grad_hr_vs_mg.py [outfile.png] [-k LEVEL] [-p PCTL] [-n NITER]
                                  [--field adxx_theta] [--shared]

    --shared   one colour scale for the whole figure instead of per column
               (honest about the ~10x decay across iterations, but late
               iterations then read as blank)
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

# ------------------------------------------------------------------ config
RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
HR_RUN   = RUNSROOT + 'assim_swot_HR_10days_newCtrlW_v2/'
MG_RUN   = RUNSROOT + 'assim_swot_MG_10days_xxmult16_newCtrlW_v2/'
HR_GRID  = RUNSROOT + 'grid_hires_cleanbathy/'
LR_GRID  = RUNSROOT + 'grid_lores_cleanbathy_v2/'

NX, NY, NZ = 120, 96, 21     # LR
FAC        = 4               # HR cells per LR cell, per direction

# diverging pair: blue <-> red, neutral warm-gray midpoint (never a hue at the
# midpoint, never a rainbow).  Land is a mid gray, clearly darker than the
# midpoint so "no data" never reads as "near zero".
DIVERGING = ['#0d366b', '#1c5cab', '#2a78d6', '#86b6ef', '#cde2fb',
             '#f0efec',
             '#fbd6d3', '#f09a95', '#e34948', '#b02f2f', '#7a1f1f']
LAND      = '#9e9c97'
INK       = '#0b0b0b'
INK2      = '#52514e'

CMAP = LinearSegmentedColormap.from_list('blue_red_div', DIVERGING, N=256)
CMAP.set_bad(LAND)


def native(a, dtype='f8'):
    """MITgcm binaries are big-endian; keep everything native inside numpy."""
    return np.asarray(np.asarray(a), dtype=dtype)


def load_grids():
    hf_hr = native(rdmds(HR_GRID + 'hFacC'))
    rA_hr = native(rdmds(HR_GRID + 'RAC'))
    hf_lr = native(rdmds(LR_GRID + 'hFacC'))
    assert hf_hr.shape == (NZ, NY * FAC, NX * FAC), hf_hr.shape
    assert hf_lr.shape == (NZ, NY, NX), hf_lr.shape
    W_hr = hf_hr * rA_hr[None, :, :]                      # wet area per HR cell
    W_lr = W_hr.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    return W_hr, W_lr, hf_lr


def restrict(g_hr, W_hr, W_lr):
    """Wet-area-weighted block mean, HR -> LR.  Land carries W=0 and drops out."""
    num = (np.nan_to_num(g_hr) * W_hr).reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    return np.where(W_lr > 0, num / np.where(W_lr > 0, W_lr, 1.0), 0.0)


def read_iter(run, sub, field, it):
    path = run + (sub % it) + '/' + field
    if not os.path.exists(path + '.%010d.data' % 0) and not os.path.exists(path + '.meta'):
        import glob
        if not glob.glob(path + '.*.data'):
            return None
    return native(rdmds(path, itrs=np.nan))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?', default='grad_hr_vs_mg.png')
    ap.add_argument('-k', '--level', type=int, default=0, help='vertical level (0-based)')
    ap.add_argument('-p', '--pctl', type=float, default=99.0,
                    help='percentile of |g| used as the symmetric colour limit')
    ap.add_argument('-n', '--niter', type=int, default=10, help='iterations 0..n-1')
    ap.add_argument('--field', default='adxx_theta')
    ap.add_argument('--shared', action='store_true',
                    help='one colour scale for the whole figure')
    a = ap.parse_args()

    k = a.level
    W_hr, W_lr, hf_lr = load_grids()
    wet = hf_lr[k] > 0
    print('level %d: %d wet LR cells of %d' % (k, int(wet.sum()), wet.size))

    # ---- read everything first, so the colour limits can be computed up front
    cols = []
    for it in range(a.niter):
        g_hr = read_iter(HR_RUN, 'run_adhi_it%04d', a.field, it)
        g_mg = read_iter(MG_RUN, 'run_adlo_it%04d', a.field, it)
        if g_hr is None or g_mg is None:
            print('iter %d: missing %s -- stopping at %d columns' % (it, a.field, len(cols)))
            break
        hr_lr = restrict(g_hr, W_hr, W_lr)[k]
        mg    = g_mg[k]
        cols.append((it, hr_lr, mg))

    if not cols:
        sys.exit('no iterations found')

    def lim(fields):
        v = np.concatenate([np.abs(f[wet]) for f in fields])
        v = v[np.isfinite(v)]
        return max(np.percentile(v, a.pctl), np.finfo(float).tiny)

    if a.shared:
        vmax_all = lim([f for _, hr, mg in cols for f in (hr, mg)])
        limits = [vmax_all] * len(cols)
    else:
        limits = [lim([hr, mg]) for _, hr, mg in cols]

    print('\n%-5s %12s %12s %12s %12s %12s'
          % ('iter', 'HR|max|', 'HR rms', 'MG|max|', 'MG rms', 'clim'))
    for (it, hr, mg), vm in zip(cols, limits):
        print('%-5d %12.4e %12.4e %12.4e %12.4e %12.4e'
              % (it, np.abs(hr[wet]).max(), np.sqrt((hr[wet] ** 2).mean()),
                 np.abs(mg[wet]).max(), np.sqrt((mg[wet] ** 2).mean()), vm))

    # ---- figure: 2 map rows + a thin row of per-column colourbars.
    # Laid out in inches so the equal-aspect panels tile without gutters and the
    # colourbar labels cannot collide with the neighbouring column.
    n   = len(cols)
    pw  = 1.55                       # panel width, inches
    ph  = pw * NY / NX               # equal aspect
    gx, gy = 0.07, 0.10              # gaps between panels
    ml, mr = 1.12, 0.08              # margins: row labels / right
    mt, mb = 0.78, 0.40              # title block / colourbar labels
    cbh = 0.10                       # colourbar height
    cbgap = 0.20

    fw = ml + n * pw + (n - 1) * gx + mr
    fh = mt + 2 * ph + gy + cbgap + cbh + mb
    fig = plt.figure(figsize=(fw, fh), facecolor='#fcfcfb')

    def box(x, y, w, h):             # inches from bottom-left -> figure fraction
        return [x / fw, y / fh, w / fw, h / fh]

    y_row = [fh - mt - ph, fh - mt - 2 * ph - gy]
    y_cb  = mb

    for j, ((it, hr, mg), vm) in enumerate(zip(cols, limits)):
        x = ml + j * (pw + gx)
        for i, fld in enumerate((hr, mg)):
            ax = fig.add_axes(box(x, y_row[i], pw, ph))
            m = np.ma.masked_invalid(np.where(wet, fld, np.nan))
            mesh = ax.pcolormesh(m, cmap=CMAP, vmin=-vm, vmax=vm,
                                 rasterized=True, shading='auto')
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color('#d5d3ce'); s.set_linewidth(0.6)
            if i == 0:
                ax.set_title('it %d' % it, fontsize=9.5, color=INK, pad=3)

        cax = fig.add_axes(box(x + 0.10 * pw, y_cb + cbh, 0.80 * pw, cbh))
        cb = fig.colorbar(mesh, cax=cax, orientation='horizontal', extend='both')
        cb.set_ticks([]); cb.outline.set_visible(False)
        fig.text((x + 0.5 * pw) / fw, (y_cb + 0.30 * cbh) / fh,
                 u'\u00b1%.2g' % vm, ha='center', va='top',
                 fontsize=7.5, color=INK2)

    for i, lab in enumerate(('HR-only\n(coarsened)', 'multi-grid\n(LR adjoint)')):
        fig.text((ml - 0.10) / fw, (y_row[i] + 0.5 * ph) / fh, lab,
                 ha='right', va='center', fontsize=9.5, color=INK)

    fig.text(ml / fw, (fh - 0.30) / fh,
             '%s, level k = %d   \u2014   HR-only adjoint restricted to 1/3\u00b0 '
             'vs multi-grid LR adjoint' % (a.field, k),
             ha='left', va='center', fontsize=12, color=INK)
    fig.text(ml / fw, (fh - 0.52) / fh,
             'colour limits %s, symmetric at the %g-th percentile of |g| over wet cells; '
             'arrowheads mark cells beyond the limit; gray is land'
             % ('per iteration' if not a.shared else 'shared across the figure', a.pctl),
             ha='left', va='center', fontsize=8.5, color=INK2)

    fig.savefig(a.outfile, dpi=170, facecolor=fig.get_facecolor())
    print('\nwrote %s' % os.path.abspath(a.outfile))


if __name__ == '__main__':
    main()
