#!/usr/bin/env python3
"""
Diagnostics for the 10-day multi-grid SWOT optimizations.

Reproduces every claim made about assim_swot_MG_10days{,_newinterp,_xxmult16}
against assim_swot_HR_10days / assim_swot_LR_10days.  Each figure is generated
by its own function so you can run them one at a time from an interpreter:

    from plot_mg_diagnostics import *
    fig_cost_history()          # 1. cost per outer iteration, accepted vs rejected
    fig_step_lengths()          # 2. accepted line-search step length
    fig_salt_depth_profile()    # 3. applied salinity adjustment vs depth
    fig_isolated_point()        # 4. the single-grid-point adjustment at k=8
    fig_weight_support()        # 5. LR control cells with no HR weight support
    fig_penalty_concentration() # 6. how concentrated the LR salt penalty is

or all of them:

    python3 plot_mg_diagnostics.py [outdir]

Every function also prints the numbers it plots, so the figures have a table
view (needed: three light-mode palette slots sit below 3:1 contrast).

Iterations used, and why they matter -- see NOTE_ON_ITERATIONS below.
"""

import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
HR_GRID  = RUNSROOT + 'grid_hires_cleanbathy/'
LR_GRID  = RUNSROOT + 'grid_lores_cleanbathy_v2/'
LR_WEIGHTS = '/home/shoshi/MITgcm_c69j/lab_sea12/input_weights_lores/'
OUTDIR   = '.'

NX, NY, NZ   = 120, 96, 21          # LR
FAC          = 4                    # LR * FAC = HR
NXH, NYH     = NX * FAC, NY * FAC
ITERMAX      = 10

# (label, run dir, subdir pattern) -- order fixes the colour slots
RUNS = [
    ('MG old (bilinear)', 'assim_swot_MG_10days',            'run_adhi_it%04d'),
    ('MG newinterp',      'assim_swot_MG_10days_newinterp',  'run_adhi_it%04d'),
    ('MG xxmult16',       'assim_swot_MG_10days_xxmult16',   'run_adhi_it%04d'),
    ('HR only',           'assim_swot_HR_10days',            'run_adhi_it%04d'),
    ('LR only',           'assim_swot_LR_10days',            'run_adlo_it%04d'),
]

NOTE_ON_ITERATIONS = """
All field comparisons are at it0010 unless stated.  Check the accepted/rejected
label first: m1qn3 evaluates the cost at trial points as well as at iterates, so
not every it####/ directory is an iterate.

    MG old        it0010 is a REJECTED TRIAL.  Last accepted iterate: it0006.
    MG newinterp  it0010 is an accepted iterate (m1qn3 iter 8).
    MG xxmult16   it0010 is an accepted iterate (m1qn3 iter 8).
    HR only       it0010 is an accepted iterate (m1qn3 iter 11).

So MG-old field numbers at it0010 describe a point the optimizer rejected.  Use
it0006 for a like-for-like comparison against the other two.
"""

# --------------------------------------------------------------------------
# palette -- reference instance, fixed slot order, adjacent pairlist (lines).
# Not re-derived here; substitute your own and re-run the skill's validator.
# --------------------------------------------------------------------------
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
SURFACE      = '#fcfcfb'
TEXT_PRIMARY = '#0b0b0b'
TEXT_SECOND  = '#52514e'
TEXT_MUTED   = '#8a8880'
GRID         = '#e8e7e3'
DIVERGING    = ['#0d366b', '#256abf', '#86b6ef', '#f0efec',   # blue -> gray
                '#f2a09f', '#e34948', '#a11f1e']              # -> red
SEQUENTIAL   = ['#cde2fb', '#9ec5f4', '#5598e7', '#2a78d6', '#184f95', '#0d366b']

DIVCMAP = LinearSegmentedColormap.from_list('blue_gray_red', DIVERGING)
SEQCMAP = LinearSegmentedColormap.from_list('blue_seq', SEQUENTIAL)


def _style(ax, title=None, xlabel=None, ylabel=None):
    """Recessive grid and axes; thin marks carry the data."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=TEXT_SECOND, labelsize=9, length=3, width=0.8)
    if title:
        ax.set_title(title, color=TEXT_PRIMARY, fontsize=11, loc='left', pad=10)
    if xlabel:
        ax.set_xlabel(xlabel, color=TEXT_SECOND, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT_SECOND, fontsize=9)


def _save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, dpi=150, bbox_inches='tight', facecolor=SURFACE)
    plt.close(fig)
    print('  wrote %s' % path)
    return path


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def read_costfinal(rundir):
    with open(os.path.join(rundir, 'costfinal')) as fh:
        return float(fh.read().strip())


def read_m1qn3(optimdir):
    """Parse OPTIM/m1qn3_output.txt.

    Returns (iterates, trials):
      iterates  list of (m1qn3_iter, simul, f, hprime0)
      trials    list of (t, ffn, accepted)   in file order
    'accepted' is decided by column position: mlis3 prints the accepted trial
    flush after the routine name and pads a rejected one far to the right.
    """
    with open(os.path.join(optimdir, 'm1qn3_output.txt')) as fh:
        lines = fh.read().split('\n')

    def num(x):
        return float(x.rstrip(',').replace('D', 'e'))

    iterates, trials = [], []
    for ln in lines:
        m = re.match(r"\s*m1qn3: iter\s+(\d+), simul\s+(\d+), "
                     r"f=\s*([0-9.DE+-]+), h'\(0\)=\s*([0-9.DE+-]+)", ln)
        if m:
            iterates.append((int(m.group(1)), int(m.group(2)),
                             num(m.group(3)), num(m.group(4))))
            continue
        m = re.match(r'\s+mlis3\s{2,12}(\d\.\d+D[-+]\d+)\s+(\S+)\s+(\S+)', ln)
        if m:
            trials.append((num(m.group(1)), num(m.group(2)), True))
            continue
        m = re.match(r'\s+mlis3\s{20,}(\d\.\d+D[-+]\d+)\s+(\S+)\s+(\S+)', ln)
        if m:
            trials.append((num(m.group(1)), num(m.group(2)), False))
    return iterates, trials


def classify_iterates(rundir_root, subpat, optimdir, itermax=ITERMAX):
    """Label each it#### as initial / accepted iterate / rejected trial.

    it000i is the (i+1)-th cost evaluation, so it matches an m1qn3 header whose
    'simul' equals i+1 -- provided the f value agrees, which is the actual test.
    """
    iterates, _ = read_m1qn3(optimdir)
    by_simul = {s: f for _, s, f, _ in iterates}
    out = []
    for i in range(itermax + 1):
        c = read_costfinal(os.path.join(rundir_root, subpat % i))
        if i == 0:
            kind = 'initial'
        else:
            f = by_simul.get(i + 1)
            kind = ('accepted' if f is not None
                    and abs(f - c) < 1e-3 * max(1.0, abs(c)) else 'rejected')
        out.append((i, c, kind))
    return out


def read_xx_hr(rundir, field, it):
    """HR control: rdmds, because the .meta sidecar is written."""
    return np.asarray(rdmds(os.path.join(rundir, 'xx_' + field), it),
                      dtype=np.float64)


def read_xx_lr(rundir, field, it, nz=NZ):
    """LR control: make_cost_cp_v2.py writes .data only, so no rdmds."""
    f = os.path.join(rundir, 'xx_%s.%010d.data' % (field, it))
    a = np.fromfile(f, dtype='>f8')
    return a.reshape(-1, NY, NX)


def applied_adjustment(run, field='salt', it0=0, it1=ITERMAX):
    """The perturbation the model actually received, in physical units.

    xx_*.effective holds the resulting field, so differencing two iterations
    isolates the adjustment and cancels the background.  Do not read the
    absolute .effective values: their extremes sit on points that never change
    and lie outside the namelist bounds.
    """
    d0 = os.path.join(RUNSROOT, run, 'run_adhi_it%04d' % it0)
    d1 = os.path.join(RUNSROOT, run, 'run_adhi_it%04d' % it1)
    e0 = np.asarray(rdmds(os.path.join(d0, 'xx_%s.effective' % field), it0),
                    dtype=np.float64)
    e1 = np.asarray(rdmds(os.path.join(d1, 'xx_%s.effective' % field), it1),
                    dtype=np.float64)
    return e1 - e0


def hfac(hr=True):
    return rdmds((HR_GRID if hr else LR_GRID) + 'hFacC')


def salt_weights():
    """(HR weight, LR weight).  Both are raw binaries: no .meta, big-endian f8."""
    wh = np.fromfile(os.path.join(
        RUNSROOT, 'assim_swot_MG_10days_xxmult16', 'run_adhi_it0010',
        'S_wt_fromASTE_HR_removeboundary.bin'), dtype='>f8')
    wl = np.fromfile(LR_WEIGHTS + 'S_wt_fromASTE_removeboundary_v4.bin',
                     dtype='>f8')
    return wh.reshape(NZ, NYH, NXH), wl.reshape(NZ, NY, NX)


def hr_weight_support_per_block(mask_by_hfac=True):
    """Count of HR points per LR block that carry a usable salt control.

    'Usable' means nonzero weight AND wet.  ctrl_map_genarr.F:571-577 zeroes the
    perturbation wherever the HR weight is zero, so a block with a count of 0
    cannot receive any adjustment at all.
    """
    wh, _ = salt_weights()
    sup = wh > 0
    if mask_by_hfac:
        sup &= hfac(hr=True) > 0
    return sup.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))


# ==========================================================================
# 1. cost history
# ==========================================================================
def fig_cost_history(name='mg_01_cost_history.png'):
    print('\n[1] cost per outer iteration')
    print('    %-20s %-7s %12s  %s' % ('run', 'iter', 'costfinal', 'kind'))
    fig, ax = plt.subplots(figsize=(9, 5.2))
    _style(ax, 'Cost per outer iteration, 10-day SWOT-only',
           'outer iteration', 'costfinal  (misfit + HR-space control penalty)')

    for slot, (label, run, subpat) in enumerate(RUNS):
        root = os.path.join(RUNSROOT, run)
        rows = classify_iterates(root, subpat, os.path.join(root, 'OPTIM'))
        for i, c, kind in rows:
            print('    %-20s %-7d %12.1f  %s' % (label, i, c, kind))
        col = SERIES[slot]
        keep = [(i, c) for i, c, k in rows if k in ('initial', 'accepted')]
        rej  = [(i, c) for i, c, k in rows if k == 'rejected']
        ax.plot([i for i, _ in keep], [c for _, c in keep],
                color=col, linewidth=2, marker='o', markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=2,
                label=label, zorder=3 + slot)
        if rej:
            # open markers, deliberately unconnected: a rejected trial is not
            # on the optimizer's path
            ax.plot([i for i, _ in rej], [c for _, c in rej], linestyle='none',
                    marker='o', markersize=5, markerfacecolor='none',
                    markeredgecolor=col, markeredgewidth=1.4, zorder=2)
    # 5 series: legend only.  Direct labels are capped at 4, and the printed
    # table above is the table view the low-contrast slots require.
    ax.set_yscale('log')
    ax.set_xlim(-0.4, 10.4)
    ax.set_xticks(range(0, 11))
    ax.set_yticks([3e5, 4e5, 6e5, 1e6, 1.5e6, 2.5e6])
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: '%.2fM' % (v / 1e6)))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    leg = ax.legend(frameon=False, fontsize=8.5, loc='upper right',
                    labelcolor=TEXT_SECOND)
    leg.set_title('filled = iterate,  open = rejected trial',
                  prop={'size': 8})
    leg.get_title().set_color(TEXT_MUTED)
    return _save(fig, name)


# ==========================================================================
# 2. accepted step length
# ==========================================================================
def fig_step_lengths(name='mg_02_step_lengths.png'):
    print('\n[2] accepted line-search step length t')
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    _style(ax, 'Accepted line-search step length',
           'accepted step (m1qn3 iteration)', 'step length  t')

    for slot, (label, run, _) in enumerate(RUNS):
        _, trials = read_m1qn3(os.path.join(RUNSROOT, run, 'OPTIM'))
        acc = [t for t, _, ok in trials if ok]
        print('    %-20s %s' % (label, ', '.join('%.4g' % t for t in acc)))
        ax.plot(range(1, len(acc) + 1), acc, color=SERIES[slot], linewidth=2,
                marker='o', markersize=5, markeredgecolor=SURFACE,
                markeredgewidth=2, label=label)
    ax.axhline(1.0, color=TEXT_MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate('t = 1: full quasi-Newton step accepted', (0.6, 1.0),
                textcoords='offset points', xytext=(0, 7),
                color=TEXT_MUTED, fontsize=8.5)
    ax.set_yscale('log')
    ax.legend(frameon=False, fontsize=8.5, loc='lower left',
              labelcolor=TEXT_SECOND)
    return _save(fig, name)


# ==========================================================================
# 3. applied salinity adjustment vs depth
# ==========================================================================
def fig_salt_depth_profile(name='mg_03_salt_depth.png',
                           runs=(('MG xxmult16', 'assim_swot_MG_10days_xxmult16'),
                                 ('HR only',     'assim_swot_HR_10days')),
                           thresh=2.0):
    print('\n[3] applied salinity adjustment by level  (effective it%02d - it0000)'
          % ITERMAX)
    rc = rdmds(HR_GRID + 'RC').squeeze()
    nshow = 13                      # top ~1160 m; deeper levels are flat
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.4), sharey=True)
    _style(axes[0], 'RMS of applied adjustment', 'psu', 'depth  (m)')
    _style(axes[1], 'points exceeding %.0f psu' % thresh, 'count')

    print('    %-6s %-9s %-11s %-11s' % ('k', 'depth m', 'rms', 'n>%.0f' % thresh))
    for slot, (label, run) in enumerate(runs):
        dd = applied_adjustment(run)
        rms = np.array([dd[k][dd[k] != 0].std() if (dd[k] != 0).any() else 0.0
                        for k in range(NZ)])
        cnt = np.array([(np.abs(dd[k]) > thresh).sum() for k in range(NZ)])
        for k in range(12):
            print('    %-6d %-9.1f %-11.4f %-11d  <- %s'
                  % (k, rc[k], rms[k], cnt[k], label))
        col = SERIES[slot]
        for ax, series in ((axes[0], rms), (axes[1], cnt)):
            ax.plot(series[:nshow], rc[:nshow], color=col, linewidth=2,
                    marker='o', markersize=4.5, markeredgecolor=SURFACE,
                    markeredgewidth=1.5, label=label)
            # the two series converge at depth, so labels would collide:
            # identity comes from the legend plus the printed table
    for ax in axes:
        ax.set_ylim(rc[nshow - 1] * 1.05, 0)
    axes[0].legend(frameon=False, fontsize=9, labelcolor=TEXT_SECOND,
                   loc='lower right')
    fig.suptitle('Where in the water column each run puts its salinity '
                 'adjustment', color=TEXT_PRIMARY, fontsize=12, x=0.09,
                 ha='left', y=1.0)
    return _save(fig, name)


# ==========================================================================
# 4. the isolated single-grid-point adjustment
# ==========================================================================
def fig_isolated_point(name='mg_04_isolated_point.png',
                       run='assim_swot_MG_10days_xxmult16',
                       k=8, jlr=55, ilr=82, halo=3):
    """Zoom on one LR block where only 1 of 16 HR points carries a weight.

    SMOOTH_CORREL3D runs at ctrl_map_genarr.F:555, BEFORE the pointwise
    1/sqrt(weight) scaling at 571-577.  So the correlation operator smooths v
    and the scaling then re-imprints the weight field's zero pattern, which the
    operator no longer gets a chance to regularize.
    """
    print('\n[4] isolated adjustment, LR block (k=%d, j=%d, i=%d)' % (k, jlr, ilr))
    dd = applied_adjustment(run)
    cnt = hr_weight_support_per_block()
    j0, i0 = jlr * FAC, ilr * FAC
    js, je = j0 - halo * FAC, j0 + (halo + 1) * FAC
    is_, ie = i0 - halo * FAC, i0 + (halo + 1) * FAC
    sub = dd[k, js:je, is_:ie]

    print('    HR points with a usable weight in this block: %d of %d'
          % (cnt[k, jlr, ilr], FAC * FAC))
    print('    applied adjustment inside the block (psu):')
    print(np.array2string(dd[k, j0:j0 + FAC, i0:i0 + FAC], precision=3,
                          suppress_small=True, prefix='      '))

    lim = float(np.abs(sub).max()) or 1.0
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    _style(ax, 'Applied salinity adjustment, level k=%d' % k,
           'HR i index', 'HR j index')
    im = ax.pcolormesh(np.arange(is_, ie), np.arange(js, je), sub,
                       cmap=DIVCMAP, norm=TwoSlopeNorm(0.0, -lim, lim),
                       shading='nearest')
    # LR block boundaries, so the 4x4 structure is visible
    for jj in range(js, je + 1, FAC):
        ax.axhline(jj - 0.5, color=TEXT_MUTED, linewidth=0.6, alpha=0.55)
    for ii in range(is_, ie + 1, FAC):
        ax.axvline(ii - 0.5, color=TEXT_MUTED, linewidth=0.6, alpha=0.55)
    ax.add_patch(plt.Rectangle((i0 - 0.5, j0 - 0.5), FAC, FAC, fill=False,
                               edgecolor=TEXT_PRIMARY, linewidth=2))
    cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.03)
    cb.set_label('psu', color=TEXT_SECOND, fontsize=9)
    cb.ax.tick_params(colors=TEXT_SECOND, labelsize=8)
    cb.outline.set_visible(False)
    ax.annotate('LR cell j=%d i=%d\n%d/16 HR points weighted'
                % (jlr, ilr, cnt[k, jlr, ilr]), (i0 + FAC, j0 + FAC),
                textcoords='offset points', xytext=(8, 8),
                color=TEXT_PRIMARY, fontsize=9)
    return _save(fig, name)


# ==========================================================================
# 5. LR control cells with no HR weight support
# ==========================================================================
def fig_weight_support(name='mg_05_weight_support.png'):
    """Cells where the LR adjoint hands back a gradient the HR run discards."""
    print('\n[5] LR salt control cells with no HR weight support')
    wh, wl = salt_weights()
    hh, hl = hfac(True) > 0, hfac(False) > 0
    for tag, (SH, SL) in [('unmasked', (wh > 0, wl > 0)),
                          ('hFacC-masked', ((wh > 0) & hh, (wl > 0) & hl))]:
        c = SH.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
        print('    %-14s LRw>0 & no HR support = %5d    HR support & LRw=0 = %5d'
              % (tag, int(((c == 0) & SL).sum()), int(((c > 0) & ~SL).sum())))
    cnt = hr_weight_support_per_block()
    live = (wl > 0) & hl
    bad = (cnt == 0) & live
    print('    inert cells: %d of %d wet LR control cells (%.2f%%)'
          % (bad.sum(), live.sum(), 100 * bad.sum() / live.sum()))

    per_k = np.array([[live[k].sum(), bad[k].sum()] for k in range(NZ)])
    rc = rdmds(HR_GRID + 'RC').squeeze()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0),
                             gridspec_kw={'width_ratios': [1, 1.35]})
    _style(axes[0], 'Inert share of LR salt controls, by level',
           'percent of wet LR cells', 'depth  (m)')
    frac = 100 * per_k[:, 1] / np.maximum(per_k[:, 0], 1)
    axes[0].plot(frac, rc, color=SERIES[0], linewidth=2, marker='o',
                 markersize=4.5, markeredgecolor=SURFACE, markeredgewidth=1.5)
    for k in range(NZ):
        print('    k=%2d  depth %8.1f m   wet %5d   inert %4d  (%.2f%%)'
              % (k, rc[k], per_k[k, 0], per_k[k, 1], frac[k]))

    kmap = int(np.argmax(per_k[:, 1]))
    _style(axes[1], 'Where they are, level k=%d (%.0f m)' % (kmap, rc[kmap]),
           'LR i', 'LR j')
    axes[1].grid(False)
    # three states, drawn as explicit masked layers so nothing is ambiguous
    ok = live[kmap] & ~bad[kmap]
    axes[1].imshow(np.ma.masked_where(~ok, np.ones_like(ok, float)),
                   origin='lower', cmap=matplotlib.colors.ListedColormap(
                       [SEQUENTIAL[1]]), vmin=0, vmax=1, interpolation='nearest')
    axes[1].imshow(np.ma.masked_where(~bad[kmap], np.ones_like(ok, float)),
                   origin='lower', cmap=matplotlib.colors.ListedColormap(
                       [SERIES[1]]), vmin=0, vmax=1, interpolation='nearest')
    axes[1].legend(handles=[
        plt.Rectangle((0, 0), 1, 1, facecolor=SEQUENTIAL[1], edgecolor='none',
                      label='control acts on the HR model (%d)'
                            % int(ok.sum())),
        plt.Rectangle((0, 0), 1, 1, facecolor=SERIES[1], edgecolor='none',
                      label='LR gradient, HR perturbation discarded (%d)'
                            % int(bad[kmap].sum())),
        plt.Rectangle((0, 0), 1, 1, facecolor=SURFACE, edgecolor=GRID,
                      label='dry / not a control')],
        frameon=False, fontsize=8.5, loc='upper left', labelcolor=TEXT_SECOND)
    axes[1].set_aspect('equal')
    return _save(fig, name)


# ==========================================================================
# 6. concentration of the LR control penalty
# ==========================================================================
def fig_penalty_concentration(name='mg_06_penalty_concentration.png',
                              field='salt', topn=15):
    """How much of sum(v^2) comes from how few LR cells."""
    print('\n[6] concentration of the LR %s penalty' % field)
    cnt = hr_weight_support_per_block()
    fig, ax = plt.subplots(figsize=(9, 4.8))
    _style(ax, 'Share of the LR %s control penalty, largest cells first' % field,
           'LR control cell, ranked', 'percent of sum(v^2)')

    for slot, (label, run, subpat) in enumerate(RUNS[:3]):
        d = os.path.join(RUNSROOT, run, subpat % ITERMAX)
        v = read_xx_hr(d, field, ITERMAX)
        # injection makes the block uniform up to the hFacC factor, so the block
        # max recovers the LR value; the ::FAC corner can land on land
        vb = np.abs(v).reshape(NZ, NY, FAC, NX, FAC).max(axis=(2, 4))
        sq = (vb ** 2).ravel()
        tot = sq.sum()
        order = np.argsort(sq)[::-1]
        share = 100 * sq[order[:topn]] / tot
        ax.plot(range(1, topn + 1), share, color=SERIES[slot], linewidth=2,
                marker='o', markersize=5, markeredgecolor=SURFACE,
                markeredgewidth=2, label=label)
        print('    %-20s total sum(v^2) = %.1f' % (label, tot))
        cflat = cnt.ravel()
        for r in range(5):
            idx = order[r]
            k, j, i = np.unravel_index(idx, vb.shape)
            print('        rank %d: k=%2d j=%2d i=%3d  v=%7.2f  %.1f%% of total'
                  '  HR weighted pts=%2d %s'
                  % (r + 1, k, j, i, vb[k, j, i], share[r], cflat[idx],
                     '<-- INERT' if cflat[idx] == 0 else ''))
    ax.legend(frameon=False, fontsize=8.5, labelcolor=TEXT_SECOND)
    return _save(fig, name)


# ==========================================================================
# 7. mask provenance -- run this FIRST when a mask result looks wrong
# ==========================================================================
def fig_verify_masks(name='mg_07_verify_masks.png', klevels=(0, 8)):
    """Three separate masks get confused with each other.  Keep them apart:

      grid_hires_cleanbathy/hFacC       the HR bathymetry the model runs on
      grid_lores_cleanbathy_v2/hFacC    the LR bathymetry -- verified below to be
                                        the exact 4x4 coarsening of the HR one
      the xx_genarr3d_weight binary     the control PRIOR, a separate file whose
                                        support need not match either

    Two traps this figure exists to catch.  First, a level-k mask is not the
    surface coastline: at k=8 (222 m) the shelves are gone, so it looks nothing
    like hFacC[0].  Second, the LR weight file's support is grid_lores, NOT
    grid_lores_cleanbathy_v2 -- so it is nonzero on cells that are land in the
    grid the model actually uses.
    """
    print('\n[7] mask provenance')
    hh = hfac(hr=True) > 0
    hl = hfac(hr=False) > 0
    hlo = rdmds(RUNSROOT + 'grid_lores/hFacC') > 0
    _, wl = salt_weights()
    wl = wl > 0

    print('    LR bathymetry == 4x4 coarsening of HR bathymetry?')
    for k in range(NZ):
        c = hh[k].reshape(NY, FAC, NX, FAC).any(axis=(1, 3))
        flag = '' if (c == hl[k]).all() else '   <-- MISMATCH'
        print('      k=%2d  HR wet %7d   LR wet %6d   coarsened %6d   match %s%s'
              % (k, hh[k].sum(), hl[k].sum(), c.sum(),
                 bool((c == hl[k]).all()), flag))
    print('    LR weight support == grid_lores wet:            %s  (%d vs %d)'
          % (bool((wl == hlo).all()), wl.sum(), hlo.sum()))
    print('    LR weight support == grid_lores_cleanbathy_v2:  %s  (%d vs %d)'
          % (bool((wl == hl).all()), wl.sum(), hl.sum()))
    print('    weight nonzero on cells that are LAND in cleanbathy_v2: %d'
          % int((wl & ~hl).sum()))

    cmap = matplotlib.colors.ListedColormap([SURFACE, SEQUENTIAL[1]])
    fig, axes = plt.subplots(len(klevels), 4,
                             figsize=(15, 4.2 * len(klevels)))
    axes = np.atleast_2d(axes)
    for r, k in enumerate(klevels):
        coarse = hh[k].reshape(NY, FAC, NX, FAC).any(axis=(1, 3))
        panels = [
            (hh[k],  'HR hFacC>0  (%d wet)' % hh[k].sum()),
            (hl[k],  'LR cleanbathy_v2 hFacC>0  (%d)' % hl[k].sum()),
            (coarse, 'HR coarsened 4x4  (%d)  match=%s'
                     % (coarse.sum(), bool((coarse == hl[k]).all()))),
            (wl[k],  'LR WEIGHT support = grid_lores  (%d)' % wl[k].sum()),
        ]
        for c, (arr, title) in enumerate(panels):
            ax = axes[r, c]
            ax.imshow(arr.astype(float), origin='lower', cmap=cmap,
                      vmin=0, vmax=1, interpolation='nearest')
            ax.set_title('k=%d  %s' % (k, title), color=TEXT_PRIMARY,
                         fontsize=9, loc='left')
            ax.set_aspect('equal')
            ax.tick_params(colors=TEXT_SECOND, labelsize=8)
    fig.suptitle('Mask provenance: bathymetry agrees across resolutions; the '
                 'weight support does not match either grid',
                 color=TEXT_PRIMARY, fontsize=12, x=0.07, ha='left')
    return _save(fig, name)


def report_inert_breakdown():
    """Split the inert-control count by cause.  Prints, does not plot.

    An LR cell can look inert for two very different reasons, and they need
    different fixes:
      (i)  its HR block is entirely dry     -> the bathymetries disagree
      (ii) its HR block is wet but carries  -> the weight masks disagree
           no nonzero weight
    """
    print('\n[*] why are the inert cells inert?')
    hh, hl = hfac(True) > 0, hfac(False) > 0
    wh, wl = salt_weights()
    nwet = hh.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    nusable = (hh & (wh > 0)).reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    live = (wl > 0) & hl
    bad = live & (nusable == 0)
    print('    LR wet but HR block entirely dry : %5d'
          % int((hl & (nwet == 0)).sum()))
    print('    LR dry but HR block has wet cells: %5d'
          % int(((~hl) & (nwet > 0)).sum()))
    print('    inert total                                        %5d' % bad.sum())
    print('      (i)  bathymetry mismatch (HR block all dry)      %5d'
          % int((bad & (nwet == 0)).sum()))
    print('      (ii) weight-mask mismatch (wet, none weighted)   %5d'
          % int((bad & (nwet > 0)).sum()))


ALL = [fig_verify_masks, fig_cost_history, fig_step_lengths,
       fig_salt_depth_profile, fig_isolated_point, fig_weight_support,
       fig_penalty_concentration]


if __name__ == '__main__':
    if len(sys.argv) > 1:
        OUTDIR = sys.argv[1]
    print(NOTE_ON_ITERATIONS)
    for f in ALL:
        f()
    report_inert_breakdown()
