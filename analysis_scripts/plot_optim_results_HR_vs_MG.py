#!/usr/bin/env python3
"""
Optimization results panel for the ctrlWyearly 10-day SWOT runs -- the
optim_results.ipynb figure, restructured to show HR-only and multi-grid at once.

    row 1   Cost reduction, both runs, full 0-20 history
    row 2   SWOT SSH  |  MG model SSH  |  MG normalized misfit
    row 3             |  HR model SSH  |  HR normalized misfit

with row 3 column-aligned under row 2, and a summary block where the SWOT panel
would repeat.

What is plotted, and why
-----------------------
pkg/obsfit removes ONE global scalar offset before forming the cost
(obsfit_cost.F:229-271):

    offset = <mod_val> - <obs_val>        over samples where both are valid
    misfit = mod_val - obs_val - offset
    J      = sum  misfit^2 / obs_uncert^2        (obs_uncert > 0 only)

so the misfit panels show the DIMENSIONLESS

    z = misfit / obs_uncert

whose mean square is exactly the quantity on the cost axis.  |z| < 1 is
"within the stated observational uncertainty", which is why z is plotted
rather than the misfit in metres.  Verified against the model: recomputing J
this way reproduces the reported swot_obsfit cost to 12 significant digits
(2343397.486669 at it0).

The z colour range is +/-5.  |z| is large early on and clipping it too hard hides
the structure that is actually changing -- fraction of samples outside the range:

    iter        0     1     4     8    14    20
    MG  +/-3  78.2  64.2  33.1  31.1  29.5  26.1  %
    MG  +/-5  61.0  44.9  13.9  13.0  12.7  10.9  %
    MG +/-10  29.4  12.4   2.7   2.5   2.4   1.9  %
    HR  +/-3  78.2  76.8  46.2  33.4  23.6  19.2  %
    HR  +/-5  61.0  59.3  24.7  14.0   7.7   6.1  %
    HR +/-10  29.4  27.2   5.2   2.2   1.1   0.6  %

+/-5 keeps ~90% (MG) and ~94% (HR) of the converged iterate inside the range
while still resolving the |z| ~ 1-3 interior, which +/-10 compresses into the
near-white middle of the ramp.  Iteration 0 is heavily saturated at any of these
-- rms z there is 9.15, so that is the data, not the choice of scale.  Use
--vmax-z to retune.  The colourbar carries explicit ticks at +/-1 so the "within
uncertainty" band stays findable.

  !! obs_uncert here is 2.0-4.0 cm (median 2.003), which is KaRIn noise after
  !! the sigma=2px Gaussian and NOTHING ELSE -- no representativeness error.
  !! A 1/12 deg model cannot match individual 10-km samples to 2 cm no matter
  !! how good the controls are, so mean z^2 = 1 is NOT the attainable floor and
  !! no target line is drawn.  Treat the axis as relative progress.

The SSH panels are anomalies about their own mean over the plotted samples, so
SWOT and model are directly comparable and the difference of the two panels is
misfit = z * obs_uncert by construction.

Data notes
----------
* The per-run global <obsfitfile>.equi.nc has mod_mask identically 0 (it is
  never filled in the global file, only in the per-tile ones), so the sample
  mask used here is np.isfinite(mod_val).  That gives exactly 27960 samples,
  matching num_obsfit in costfunction*.
* That mask is bit-identical across all 42 (run, iteration) pairs -- checked --
  so the scatter offsets are built ONCE from the 27960 usable samples and only
  the colour array is updated per frame.  Do not pad back to 170776 with NaN:
  sub-pixel markers over a mostly-masked array render unreliably.
* <obsfitfile>.misfit.nc, which obsfit also writes, is all NaN in these runs --
  the NF_PUT_VARA_DOUBLE writes are never synced.  Do not read it.
* The cost axis is the swot_obsfit term divided by its own num_obsfit, i.e.
  EXACTLY <z^2> over the misfit panels -- not fc/swot_num.  fc also carries the
  control penalty (2.0% of fc at HR it0020, 2.1% at MG it0020), which is a
  penalty on N_ctrl control elements and has no business being divided by
  N_SWOT.  Read the axis as chi^2 per observation: the target is ~1, and going
  well below 1 would mean overfitting, not success.
* Filled marker = accepted m1qn3 iterate, open = rejected line-search trial.
  Without that, the bumps in the cost curve read as instability rather than as
  the line search probing.

Usage
-----
    python3 plot_optim_results_HR_vs_MG.py            # both PNGs and the GIF
    python3 plot_optim_results_HR_vs_MG.py --iter 8   # one PNG
    python3 plot_optim_results_HR_vs_MG.py --gif-only
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
from plot_mg_diagnostics import classify_iterates

# ============================== configuration ==============================
RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
HR_GRID  = RUNSROOT + 'grid_hires_cleanbathy/'
OUTDIR   = '/home/shoshi/MITgcm_c69j/lab_sea12/figures/'
OBSNAME  = 'swot_obsfit_cycles_9thru11_labsea_L3v3'
OBS_FILE = ('/home/shoshi/MITgcm_c69j/lab_sea12/input_binaries_hires/'
            + OBSNAME + '.nc')

RUNS = [('MG', 'assim_swot_MG_10days_ctrlWyearly', 'Multi-grid', 'royalblue'),
        ('HR', 'assim_swot_HR_10days_ctrlWyearly', 'HR-only',    '#1baf7a')]
ITERMAX = 20
BEST    = {'MG': 20, 'HR': 20}          # lowest costfinal in each run

VMAX_SSH = 0.6
VMAX_Z   = 5.0        # see the saturation table in the docstring; --vmax-z retunes
LOG_COST = False      # linear, as in the notebook; --log-cost switches back
CMAP  = 'RdBu_r'
MSIZE = 1.6          # pt^2; 0.1 as in the notebook is invisible in a saved PNG
# ===========================================================================


def read_cost(run):
    """fc/swot_num per iteration, from costfunction0*, plus the raw pieces."""
    out = {}
    for fp in sorted(glob.glob(os.path.join(RUNSROOT, run,
                                            'run_adhi_it0*', 'costfunction0*'))):
        it = int(re.search(r'it(\d+)', fp).group(1))
        fc = swot = num = np.nan
        for ln in open(fp):
            ln = ln.replace('D+', 'E+').replace('D-', 'E-')
            if ln.strip().startswith('fc ='):
                fc = float(ln.split('=')[1].split()[0])
            elif OBSNAME in ln and 'obs_val' in ln:
                swot = float(ln.split('=')[1].split()[0])
                num = float(ln.split('=')[1].split()[1])
        out[it] = (fc, swot, num)
    its = np.array(sorted(out))
    fc = np.array([out[i][0] for i in its])
    swot = np.array([out[i][1] for i in its])
    num = np.array([out[i][2] for i in its])
    return its, fc / num, swot / num, num


def load_obs():
    ds = xr.open_dataset(OBS_FILE)
    return (ds.obs_val.values, ds.obs_uncert.values,
            ds.sample_lon.values, ds.sample_lat.values)


def load_mod(run, it):
    """mod_val from the global equi file; mask = finite (mod_mask is unfilled)."""
    fn = os.path.join(RUNSROOT, run, 'run_adhi_it%04d' % it, OBSNAME + '.equi.nc')
    with xr.open_dataset(fn) as d:
        mv = d.mod_val.values.copy()
    return mv, np.isfinite(mv)


def fields(mv, use, obs, unc):
    """SSH anomalies about their own mean, and z -- compact, over `use` only.

    offset = <mod> - <obs> is the model's own, so z * unc == (mod - <mod>)
    - (obs - <obs>): the two SSH panels differ by exactly the misfit panel.
    """
    o, m, u = obs[use], mv[use], unc[use]
    om, mm = o.mean(), m.mean()
    return o - om, m - mm, (m - o - (mm - om)) / u


def land(ax, hf, xc, yc):
    ax.pcolormesh(xc, yc, hf, cmap='Greys_r', vmin=-3, vmax=1,
                  shading='nearest', rasterized=True)
    ax.set_xlim(xc.min(), xc.max()); ax.set_ylim(yc.min(), yc.max())


def build(iteration_of, title_suffix):
    """Build the whole figure.  Returns (fig, update, suptitle-setter).

    `update(i)` accepts an int iteration (both rows move together, which is what
    the GIF and both PNGs want) and returns the artists it touched.
    """
    import xmitgcm
    g = xmitgcm.open_mdsdataset(HR_GRID, iters=None)
    hf = g.hFacC[0].values
    xc = np.asarray(g.XC.values, dtype='f8')
    yc = np.asarray(g.YC.values, dtype='f8')

    obs, unc, lon, lat = load_obs()
    _, use = load_mod(RUNS[0][1], 0)            # constant across all run/iter
    slon, slat, sunc = lon[use], lat[use], unc[use]

    cost = {tag: read_cost(run) for tag, run, _, _ in RUNS}
    kinds = {}
    for tag, run, _, _ in RUNS:
        root = os.path.join(RUNSROOT, run)
        kinds[tag] = {i: k for i, _, k in
                      classify_iterates(root, 'run_adhi_it%04d',
                                        os.path.join(root, 'OPTIM'), ITERMAX)}

    fig = plt.figure(figsize=(16.5, 15.0))
    gs = gridspec.GridSpec(3, 6, height_ratios=[1.05, 2, 2],
                           hspace=0.26, wspace=0.50,
                           left=0.065, right=0.955, top=0.925, bottom=0.100)

    # ------------------------------------------------------------- cost panel
    axc = fig.add_subplot(gs[0, 1:5])
    marks = {}
    for tag, run, label, color in RUNS:
        its, _, fcn, _ = cost[tag]          # fcn = SWOT term / N_SWOT = <z^2>
        acc = [k for k in range(len(its)) if kinds[tag].get(its[k]) != 'rejected']
        rej = [k for k in range(len(its)) if kinds[tag].get(its[k]) == 'rejected']
        axc.plot(its, fcn, '-', color=color, lw=2.0, alpha=0.9, zorder=3,
                 label=label)
        axc.plot(its[acc], fcn[acc], 'o', color=color, ms=6, mec='white',
                 mew=0.9, zorder=4)
        axc.plot(its[rej], fcn[rej], 'o', mfc='white', mec=color, mew=1.7,
                 ms=6, zorder=4)
        marks[tag], = axc.plot([], [], 'o', mfc='none', mec='crimson', mew=2.4,
                               ms=14, zorder=6)
    if LOG_COST:
        axc.set_yscale('log')
    else:
        axc.set_ylim(0, None)
    axc.set_xticks(range(0, ITERMAX + 1))
    axc.set_xlim(-0.5, ITERMAX + 0.5)
    axc.set_xlabel('Iteration Number', fontsize=13)
    axc.set_ylabel(r'$J_{\rm SWOT}/N_{\rm SWOT} = \langle z^2\rangle$',
                   fontsize=13)
    axc.set_title('Cost Reduction', size=18, fontweight='bold', pad=8)
    axc.grid(True, which='both' if LOG_COST else 'major', alpha=0.35)
    h, l = axc.get_legend_handles_labels()
    h += [plt.Line2D([], [], color='0.35', ls='', marker='o', ms=6),
          plt.Line2D([], [], mfc='white', mec='0.35', mew=1.7, ls='', marker='o', ms=6),
          plt.Line2D([], [], mfc='none', mec='crimson', mew=2.0, ls='', marker='o', ms=11)]
    l += ['accepted iterate', 'rejected trial', 'shown below']
    axc.legend(h, l, fontsize=10, frameon=False, ncol=2, loc='upper right')

    # ---------------------------------------------------------- spatial panels
    ax_swot = fig.add_subplot(gs[1, 0:2])
    axm = {'MG': fig.add_subplot(gs[1, 2:4]), 'HR': fig.add_subplot(gs[2, 2:4])}
    axz = {'MG': fig.add_subplot(gs[1, 4:6]), 'HR': fig.add_subplot(gs[2, 4:6])}
    ax_note = fig.add_subplot(gs[2, 0:2]); ax_note.axis('off')

    for a in [ax_swot] + list(axm.values()) + list(axz.values()):
        land(a, hf, xc, yc)
        a.tick_params(labelsize=10)

    ssh_o, _, _ = fields(*(load_mod(RUNS[0][1], 0)[0],
                           use, obs, unc))
    ax_swot.scatter(slon, slat, c=ssh_o, cmap=CMAP, vmin=-VMAX_SSH,
                    vmax=VMAX_SSH, s=MSIZE, linewidths=0, rasterized=True)
    ax_swot.set_title('SWOT SSH', size=16)
    ax_swot.set_ylabel(r'Latitude $^\circ$N', fontsize=13)
    ax_swot.set_xlabel(r'Longitude $^\circ$E', fontsize=13)

    # seed every animated scatter with real data, so nothing depends on
    # set_array being able to revive an all-masked collection
    sc_m, sc_z, ttl_m, ttl_z = {}, {}, {}, {}
    for tag, run, label, color in RUNS:
        mv0, _ = load_mod(run, 0)
        _, m0, z0 = fields(mv0, use, obs, unc)
        sc_m[tag] = axm[tag].scatter(slon, slat, c=m0, cmap=CMAP,
                                     vmin=-VMAX_SSH, vmax=VMAX_SSH, s=MSIZE,
                                     linewidths=0, rasterized=True)
        sc_z[tag] = axz[tag].scatter(slon, slat, c=z0, cmap=CMAP,
                                     vmin=-VMAX_Z, vmax=VMAX_Z, s=MSIZE,
                                     linewidths=0, rasterized=True)
        ttl_m[tag] = axm[tag].set_title('', size=16)
        ttl_z[tag] = axz[tag].set_title('', size=16)
        axm[tag].set_ylabel(label, fontsize=14, color=color,
                            fontweight='bold', labelpad=8)
    for a in (axm['HR'], axz['HR']):
        a.set_xlabel(r'Longitude $^\circ$E', fontsize=13)

    note = ax_note.text(0.02, 0.95, '', transform=ax_note.transAxes, va='top',
                        ha='left', fontsize=11.5, linespacing=1.6)

    # shared horizontal colourbars, below the columns they describe
    cax1 = fig.add_axes([0.40, 0.048, 0.26, 0.012])
    cb1 = fig.colorbar(sc_m['HR'], cax=cax1, orientation='horizontal',
                       extend='both')
    cb1.set_label('SSH anomaly about the sample mean   [m]', fontsize=11)
    cax2 = fig.add_axes([0.675, 0.048, 0.28, 0.012])
    cb2 = fig.colorbar(sc_z['HR'], cax=cax2, orientation='horizontal',
                       extend='both')
    cb2.set_label(r'Normalized misfit  $z = (\rm mod-obs-offset)\,/\,'
                  r'\sigma_{obs}$', fontsize=11)
    cb2.set_ticks([-5, -3, -1, 0, 1, 3, 5])
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=10)

    # the footnote lives in the free lower-left panel, not at the figure edge,
    # where it would sit on top of the colourbar labels
    ax_note.text(0.02, 0.44,
                 '\n'.join([
                     r'$z = (\rm mod - obs - offset)/\sigma_{obs}$.  The cost',
                     r'axis above is $\langle z^2\rangle$ over exactly these',
                     r'samples.  $|z|<1$ = within the stated uncertainty;',
                     r'$\langle z^2\rangle \approx 1$ is the nominal target.',
                     '',
                     r'$\sigma_{\rm obs}$ = obs_uncert, median 2.0 cm: KaRIn',
                     r'noise after the $\sigma$=2px smoothing and nothing',
                     r'else — no representativeness error.  So',
                     r'$\langle z^2\rangle=1$ is NOT an attainable floor for a',
                     r'1/12$^\circ$ model, and no target line is drawn.',
                 ]),
                 transform=ax_note.transAxes, va='top', ha='left',
                 fontsize=10, color='0.38', linespacing=1.45)

    def update(frame):
        it = int(frame)
        lines = []
        for tag, run, label, color in RUNS:
            mv, _ = load_mod(run, it)
            _, ssh_m, z = fields(mv, use, obs, unc)
            sc_m[tag].set_array(ssh_m); sc_m[tag].set_clim(-VMAX_SSH, VMAX_SSH)
            sc_z[tag].set_array(z);     sc_z[tag].set_clim(-VMAX_Z, VMAX_Z)
            its, _, fcn, _ = cost[tag]
            k = int(np.where(its == it)[0][0])
            ttl_m[tag].set_text('%s model SSH — iter %d' % (label, it))
            ttl_z[tag].set_text('%s misfit — iter %d' % (label, it))
            marks[tag].set_data([it], [fcn[k]])
            lines.append('%s   iter %d   (%s)\n'
                         '    $\\langle z^2\\rangle$ = %.2f       rms $z$ = %.2f\n'
                         '    rms misfit = %.1f cm     $|z|<1$ on %.1f%% of '
                         'samples'
                         % (label, it, kinds[tag].get(it, 'unknown'), fcn[k],
                            np.sqrt(np.mean(z ** 2)),
                            100 * np.sqrt(np.mean((z * sunc) ** 2)),
                            100.0 * np.mean(np.abs(z) < 1)))
        note.set_text('\n\n'.join(lines))
        return (list(sc_m.values()) + list(sc_z.values())
                + list(marks.values()) + [note])

    def sup():
        fig.suptitle('SWOT assimilation, Labrador Sea, 10-day window, '
                     'yearly control priors' + title_suffix,
                     size=17, fontweight='bold', y=0.978)

    return fig, update, sup


def one_png(it, name):
    fig, update, sup = build({'MG': it, 'HR': it}, '')
    sup()
    update(it)
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, dpi=140, facecolor='white')
    plt.close(fig)
    print('wrote %s' % path)


def png_best(name):
    assert BEST['MG'] == BEST['HR'], 'layout assumes one frame index'
    one_png(BEST['MG'], name)


def make_gif(name, fps=1.2):
    fig, update, sup = build({'MG': 0, 'HR': 0}, '')
    sup()
    ani = FuncAnimation(fig, update, frames=range(ITERMAX + 1), blit=False)
    path = os.path.join(OUTDIR, name)
    ani.save(path, writer=PillowWriter(fps=fps), dpi=70,
             savefig_kwargs=dict(facecolor='white'))
    plt.close(fig)
    print('wrote %s' % path)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--iter', type=int, default=None,
                    help='render one PNG at this iteration and stop')
    ap.add_argument('--gif-only', action='store_true')
    ap.add_argument('--vmax-z', type=float, default=VMAX_Z)
    ap.add_argument('--log-cost', action='store_true')
    a = ap.parse_args()
    VMAX_Z = a.vmax_z
    LOG_COST = a.log_cost
    os.makedirs(OUTDIR, exist_ok=True)
    if a.iter is not None:
        one_png(a.iter, 'optim_results_HR_vs_MG_ctrlWyearly_it%04d.png' % a.iter)
    elif a.gif_only:
        make_gif('optim_results_HR_vs_MG_ctrlWyearly.gif')
    else:
        one_png(8, 'optim_results_HR_vs_MG_ctrlWyearly_it0008.png')
        png_best('optim_results_HR_vs_MG_ctrlWyearly_best_it0020.png')
        make_gif('optim_results_HR_vs_MG_ctrlWyearly.gif')
