import numpy as np
import os
import re
import sys
from machine_config import DIRROOT, add_mitgcmutils
add_mitgcmutils()
from MITgcmutils import rdmds
from file_utils import *
from read_write import *


nx=120
ny=96
nz=21

factor = 4 # lowres * factor = hires

nxh=nx*factor
nyh=ny*factor

dirroot = DIRROOT

#iter = '0000'
iter = sys.argv[1] # sys.argv[0] is name of python file
print(iter)
ext = sys.argv[2]
rundirs = sys.argv[3]
rads_cost = sys.argv[4]

# --- optional resolution overrides (defaults = 1/3 deg, as the 63-day drivers use) ---
#     argv 5=factor  6=nx  7=ny  8=LR grid dir (holding maskCtrlC)
factor     = int(sys.argv[5]) if len(sys.argv) > 5 else factor
nx         = int(sys.argv[6]) if len(sys.argv) > 6 else nx
ny         = int(sys.argv[7]) if len(sys.argv) > 7 else ny
lr_maskdir = sys.argv[8]      if len(sys.argv) > 8 else 'grid_lores_cleanbathy_v2'
# argv 9 = restriction operator applied when coarsening HR -> LR:
#   'area' : wet-area-weighted block mean,
#              (R w)_c = sum_block(hFacC*rA*w) / sum_block(hFacC*rA)
#            This is the EXACT adjoint of piecewise-constant injection in the
#            wet-area-weighted inner product, so with interp_xx run in
#            'injection' mode the LR gradient is the true gradient of F(v)=f(Pv).
#   'mean' : original unweighted mean over non-zero cells.  The divisor varies
#            per block, so it is not proportional to P^T at coastlines.
# DEFAULT is 'area' (changed Sep 2026), paired with interp_xx_lores_to_hires_itX_v2.py
# argv 7 = 'injection'.  Pass 'mean' explicitly to reproduce pre-Sep-2026 runs.
restrict_method = sys.argv[9] if len(sys.argv) > 9 else 'area'
# argv 10 = HR grid dir supplying RAC and hFacC for the 'area' weights
hr_griddir      = sys.argv[10] if len(sys.argv) > 10 else 'grid_hires_cleanbathy'
# argv 11 = which space the cost handed to the optimizer is expressed in:
#   'hr' : leave costfinal alone; the LR run picks up the HR run's costfinal,
#          i.e. HR misfit + HR-space control penalty ||P v||^2.  Default, and
#          the right choice if the LR data.ctrl carries mult_gen* factors that
#          rescale the LR penalty up to the HR one.
#   'lr' : write costfinal_lo = HR misfit + LR-space control penalty ||v||^2,
#          which is the function the LR adjoint actually differentiates.
cost_space      = (sys.argv[11] if len(sys.argv) > 11 else 'hr').lower()
if cost_space not in ('hr', 'lr'):
    raise SystemExit(f'make_cost_cp: argv 11 (cost space) must be hr|lr, got {cost_space!r}')
if cost_space == 'lr' and rads_cost == 'True':
    raise SystemExit('make_cost_cp: cost_space=lr and rads_cost=True both write '
                     'costfinal_lo from different misfits; pick one')
# argv 12 = multiplier applied to the LR-space control penalty in costfinal_lo.
# It MUST equal the mult_gentim2d / mult_genarr* values in the LR data.ctrl,
# because those scale the penalty term of the gradient (ecco_cost_final.F:262-278
# is inside the differentiated code) while this scales the reported value.  Set
# both from one variable in the driver so they cannot drift apart.
#   1.0 (default) : costfinal_lo = misfit + sum_LR(v^2)          -- mult_gen* = 1
#   16.           : costfinal_lo = misfit + 16*sum_LR(v^2)       -- mult_gen* = 16
pen_mult        = float(sys.argv[12]) if len(sys.argv) > 12 else 1.0
if pen_mult != 1.0 and cost_space != 'lr':
    raise SystemExit(f'make_cost_cp: argv 12 (penalty multiplier) = {pen_mult} has no '
                     'effect unless argv 11 (cost space) is lr')
if pen_mult <= 0:
    raise SystemExit(f'make_cost_cp: argv 12 (penalty multiplier) must be > 0, got {pen_mult}')
nxh = nx*factor
nyh = ny*factor
print(f'make_cost_cp: factor={factor} nx={nx} ny={ny} nxh={nxh} nyh={nyh} maskdir={lr_maskdir}')

dirrun_lr = rundirs + '/run_adlo_it' + iter + ext + '/' #'_5day_ASTEwts/'
dirrun_hr = rundirs + '/run_adhi_it' + iter + ext + '/' #'_5day_ASTEwts/'


## NEW BATHY
# hfacc_hr = rdmds(dirrun_hr + 'hFacC')
maskc_lr = rdmds(dirroot + lr_maskdir + '/maskCtrlC')

# Wet-area weights for the 'area' restriction.
if restrict_method == 'area':
    _rA_hr    = rdmds(dirroot + hr_griddir + '/RAC')
    _hfacc_hr = rdmds(dirroot + hr_griddir + '/hFacC')
    # W[k] = hFacC(k)*rA  -> wet horizontal area of each HR cell at level k
    _W_hr = _hfacc_hr * _rA_hr[None, :, :]
    _W_lr = _W_hr.reshape(nz, ny, factor, nx, factor).sum(axis=(2, 4))
    print(f'bin_avg: wet-area-weighted restriction, '
          f'grid={hr_griddir} rA_hr {_rA_hr.min():.4g}..{_rA_hr.max():.4g}')


def bin_avg(fld, etan=False):
    """Restriction HR -> LR.

    'area' mode is the exact adjoint of piecewise-constant injection under the
    wet-area-weighted inner product:  (R w)_c = sum_block(W*w) / sum_block(W)
    with W = hFacC*rA.  Verified numerically to give <Pv,w>_HR = <v,Rw>_LR to
    machine precision.  Land cells carry W=0 so they drop out without needing
    any zero-to-NaN surgery, and the divisor is the block's wet area rather than
    a varying non-zero count.
    """
    fld = np.asarray(fld, dtype=np.float64)
    nk = fld.shape[0]

    if restrict_method == 'area':
        W = _W_hr if nk == nz else _W_hr[:1]
        Wl = _W_lr if nk == nz else _W_lr[:1]
        if nk != nz and nk != 1:
            # 2D-in-time control: one record per time level, all use the surface weight
            W = np.repeat(_W_hr[:1], nk, axis=0)
            Wl = np.repeat(_W_lr[:1], nk, axis=0)
        num = (np.nan_to_num(fld) * W).reshape(nk, ny, factor, nx, factor).sum(axis=(2, 4))
        den = Wl
        fld_lr = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
        return fld_lr

    fld = fld.reshape(nk, ny, factor, nx, factor)
    fld[fld == 0] = np.nan
    fld_lr = np.nanmean(fld, axis=(2,4))
    fld_lr[np.isnan(fld_lr)] = 0
    return fld_lr


#def write_float64(fout,fld):
#    with open(fout, 'wb') as f:
#        np.array(fld, dtype=">f8").tofile(f)


# read in high-res xx, adxx, m_, adm_

xx_files = get_files(dirrun_hr, 'xx')
adxx_files = get_files(dirrun_hr, 'adxx')
m_files = get_files(dirrun_hr, 'm')
adm_files = get_files(dirrun_hr, 'adm')


# bin-avg to low-res
xx_cost = 0
xx_cost_hr = 0
for f in xx_files.union(xx_files, adxx_files, m_files, adm_files):
    print(f)
    xx = rdmds(dirrun_hr + f)
    if f in xx_files:
        # HR-space control penalty, exactly as ctrl_cost_gen.F computes it:
        # a plain sum of squares over every record of the raw xx_* field, no
        # area weight, no smoothing applied first.  Taken before the zero->NaN
        # surgery below so the value is unambiguous.
        xx_cost_hr += float(np.sum(np.asarray(xx, dtype=np.float64)**2))
    xx[xx == 0] = np.nan
    
    etan = ('etan' in f) or ('m_eta' in f)
    if xx.shape[0] == nz:
        xx_lr = bin_avg(xx, etan=etan) * maskc_lr
    else:
        xx_lr = bin_avg(xx, etan=etan) * maskc_lr[0]

    xx_lr[np.isnan(xx_lr)] = 0
    write_float64(dirrun_lr + f + '.data', xx_lr)

    if f in xx_files:
        print(np.nansum(xx_lr**2))
        xx_cost += np.nansum(xx_lr**2)


#
#####

#misfit_sst = rdmds(dirrun_hr + 'misfit_sst')
#
#misfit_sst[misfit_sst == 0] = np.nan
#
#misfit_sst_lr = bin_avg(misfit_sst) * maskc_lr[0]
#
#misfit_sst_lr[np.isnan(misfit_sst_lr)] = 0
#
#
#write_float64(dirrun_lr + 'misfit_sst.data', misfit_sst_lr)
#
#
#

# if running ecco pkg

if os.path.exists(dirrun_hr + 'm_eta.000000' + iter + '.data' ):

    m_eta = rdmds(dirrun_hr + 'm_eta.000000' + iter )
    cost_hr = 0
    cost_lr = 0

    for f in ['tp', 'ers', 'gfo']:
        
        misfit = rdmds(dirrun_hr + 'sladiff_' + f + '_raw')

        ### confirm high-res cost
        sigma = 0.03
        cost = np.nansum((misfit / sigma)**2)
        print(f'{f} cost {cost}')
        cost_hr += cost

        misfit_lr = bin_avg(misfit) * maskc_lr[0]
        cost_lr += np.nansum((misfit_lr * maskc_lr[0] / sigma) **2)

        slaobs = rdmds(dirrun_hr + 'slaobs_' + f + '_raw')
        slaobs_lr = bin_avg(slaobs) * maskc_lr[0]
       
        write_float64(dirrun_lr + f + '.data', xx_lr)
        write_float32(dirrun_lr + 'sladiff_' + f + '_raw' + '.data', misfit_lr)
        write_float32(dirrun_lr + 'slaobs_' + f + '_raw' '.data', slaobs_lr)

    if rads_cost == 'True':
        print(f'rads total HR cost {cost_hr}')
        print(f'rads total LR cost = {cost_lr}')
        print(f'xx total LR cost = {xx_cost}')

        with open(dirrun_hr + 'costfinal_lo', 'w') as file:
            file.write(str(cost_lr + xx_cost))


# ===========================================================================
# Cost-space bookkeeping.
#
#   J_HR(v) = misfit + || P v ||^2_HR      <- what the HR run writes to costfinal
#   J_LR(v) = misfit + ||   v ||^2_LR      <- what the LR adjoint differentiates
#
# The misfit term is identical in the two.  obsfit_cost.F reads the model
# equivalents back from the staged .equi.nc files (ACTIVE_READ_OBS_TILE) and
# never re-interpolates, so the LR adjoint is seeded with the HR misfit
# bit-for-bit -- verified: sum(mod_val*mod_mask) over the 27991 SWOT samples
# agrees to 12 digits between the HR and LR run dirs.  The two objectives
# therefore differ in the control penalty and nothing else, and both penalties
# are computable right here:
#
#   costfinal_lo = costfinal_HR - sum_HR(xx^2) + sum_LR(xx^2)
#
# Subtracting the HR penalty rather than re-deriving the misfit keeps this
# agnostic to whichever observation multipliers are in force (mult_obsfit,
# mult_profiles, mult_gencost).  The one assumption is that the HR
# mult_gentim2d / mult_genarr2d / mult_genarr3d are all 1, so the HR penalty
# enters the HR costfinal with coefficient 1.
#
# Why LR is the right space: under piecewise-constant injection the 16 HR cells
# of a block carry one number, so ||P v||^2 = sum_c n_c * v_c^2 counts each
# independent parameter n_c times (n_c = wet HR cells in the block, 16 in the
# interior).  P has rank N_LR, so that is not a valid prior on the space being
# controlled -- it asserts a prior sqrt(n_c) tighter than intended.
# ===========================================================================

if cost_space == 'lr':
    fcost_hr = dirrun_hr + 'costfinal'
    if not os.path.exists(fcost_hr):
        raise SystemExit(f'make_cost_cp: cost_space=lr but no {fcost_hr}')

    with open(fcost_hr) as fh:
        costfinal_hr = float(fh.read().strip())

    misfit_hr    = costfinal_hr - xx_cost_hr
    costfinal_lo = misfit_hr + pen_mult * xx_cost

    # Cross-check against the HR run's own printout.  Informational only: the
    # value written is derived from costfinal, which is authoritative.  A
    # mismatch means either a non-unit control multiplier at HR or another
    # active cost term, both worth knowing about.
    try:
        with open(dirrun_hr + 'STDOUT.0000') as fh:
            hits = re.findall(r'obsfit_cost\(\s*\d+\)\s*=\s*(\S+)', fh.read())
        printed = [float(h.replace('D', 'e').replace('E', 'e')) for h in hits]
    except (IOError, OSError):
        printed = []

    print('--- cost space ---')
    print(f'  costfinal (HR)            = {costfinal_hr:.10g}')
    print(f'  control penalty at HR     = {xx_cost_hr:.10g}')
    print(f'  control penalty at LR     = {xx_cost:.10g}')
    if pen_mult != 1.0:
        print(f'  penalty multiplier        = {pen_mult:g}  '
              f'-> {pen_mult * xx_cost:.10g}   (must match mult_gen* in the LR data.ctrl)')
    if xx_cost > 0:
        print(f'  penalty ratio HR/LR       = {xx_cost_hr / xx_cost:.4f}   '
              f'(exactly {factor**2} for interior blocks under injection)')
    print(f'  misfit (HR, = costfinal - penalty_HR) = {misfit_hr:.10g}')
    if printed:
        print(f'  obsfit_cost printed by HR run         = '
              + ', '.join(f'{v:.10g}' for v in printed))
        if abs(misfit_hr - printed[0]) > 1e-8 * max(1.0, abs(printed[0])):
            print('  WARNING: derived misfit does not match obsfit_cost(1).  Either a '
                  'control multiplier at HR is not 1, or another cost term is active.')
    print(f'  costfinal_lo written      = {costfinal_lo:.10g}')

    with open(dirrun_hr + 'costfinal_lo', 'w') as fh:
        fh.write(str(costfinal_lo))
