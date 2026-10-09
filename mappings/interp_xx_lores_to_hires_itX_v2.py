import numpy as np
from scipy.interpolate import griddata, Rbf 
import sys
from machine_config import DIRROOT, add_mitgcmutils
add_mitgcmutils()
from MITgcmutils import rdmds

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

# --- optional resolution overrides (defaults = 1/3 deg, as the 63-day drivers use) ---
#     argv 3=factor  4=nx  5=ny  6=LR grid dir
factor  = int(sys.argv[3]) if len(sys.argv) > 3 else factor
nx      = int(sys.argv[4]) if len(sys.argv) > 4 else nx
ny      = int(sys.argv[5]) if len(sys.argv) > 5 else ny
lr_grid = sys.argv[6]      if len(sys.argv) > 6 else 'grid_lores'
# argv 7 = prolongation operator:
#   'injection' : replicate each LR value into its factor x factor HR block.
#                 This makes the HR->LR block mean in make_cost_cp_v2.py the exact
#                 adjoint of this operator (P^T = factor^2 * bin_avg), so the LR
#                 gradient is the true gradient of the reduced cost F(v)=f(Pv).
#   'linear'    : original scipy griddata bilinear + biharmonic inpaint + clip.
#                 NOT adjoint to bin_avg -- kept only for reproducing old runs.
# DEFAULT is 'injection' (changed Sep 2026).  It must stay paired with
# make_cost_cp_v2.py argv 9 = 'area', whose default was flipped at the same time.
interp_method = sys.argv[7] if len(sys.argv) > 7 else 'injection'
# argv 8 = HR grid dir supplying XC/YC and hFacC
hr_grid = sys.argv[8] if len(sys.argv) > 8 else 'grid_hires_cleanbathy'
nxh = nx*factor
nyh = ny*factor
print(f'interp_xx: factor={factor} nx={nx} ny={ny} nxh={nxh} nyh={nyh} lr_grid={lr_grid} hr_grid={hr_grid} method={interp_method}')

griddir_lr = dirroot + lr_grid + '/'
griddir_hr = dirroot + hr_grid + '/'

#dirrun_lr = dirroot + 'run_adlo_it' + iter + '/'
#dirrun_hr = dirroot + 'run_adhi_it' + iter + '/'

dirrun_pup = sys.argv[2] 
dir_out = dirrun_pup + '/xx_hires/'

# read in high-res grid and low-res grid
xc_lr = rdmds(griddir_lr + 'XC')
yc_lr = rdmds(griddir_lr + 'YC')

xc_hr = rdmds(griddir_hr + 'XC')
yc_hr = rdmds(griddir_hr + 'YC')

## NEW BATHY
#hfacc_hr = rdmds(griddir_hr + 'hFacC')
hfacc_hr = rdmds(griddir_hr + 'hFacC')

# use linear interpolation and inpaint nans for edges

from skimage.restoration import inpaint

def linear_interp(xx_lr, etan=False):

    points = np.array([xc_lr.ravel(), yc_lr.ravel()]).T

    if etan:
        z = 1
    else:
        z = xx_lr.shape[0]

    xx_hr = np.zeros((z, nyh, nxh))
    for i in range(z):

        if etan:
            values = xx_lr.ravel()
        else:
            values = xx_lr[i].ravel()
        values[values == 0] = np.nan

        if np.count_nonzero(~np.isnan(values)) == 0:
            xx_hr[i] = np.nan
        else:

            if np.count_nonzero(~np.isnan(values)) < 4:
                method = 'nearest'
            else:
                method = 'linear'

            # Create a mask to filter out NaN values in the low-resolution data
            mask = ~np.isnan(values)

            # Interpolate onto the high-resolution grid using linear interpolation
            tmp = griddata(points[mask], values[mask], (xc_hr, yc_hr), method=method)

            # Identify where the high-resolution grid still has NaN values after linear interpolation
            nan_mask = np.isnan(tmp)

            if np.any(nan_mask):
                tmp = inpaint.inpaint_biharmonic(tmp, nan_mask, multichannel=False)

            # set floor and ceiling to that of low-res  
            if etan:
                lmin = np.min(xx_lr)
                lmax = np.max(xx_lr)
            else:
                lmin = np.min(xx_lr[i,:,:])
                lmax = np.max(xx_lr[i,:,:])
            tmp = np.clip(tmp, lmin, lmax)

            if xx_lr.shape[0] == nz:
                xx_hr[i] = tmp * hfacc_hr[i]
            else:
                xx_hr[i] = tmp * hfacc_hr[0]
    print('done')

    return xx_hr


def inject(xx_lr, etan=False):
    """Piecewise-constant prolongation P: copy each LR cell into its factor x factor
    HR block, then apply the HR wet mask.

    Adjointness: with w=0 on land (which adxx_* satisfies), the HR->LR block mean
    used by make_cost_cp_v2.py::bin_avg equals (1/factor**2) * P^T exactly.  The
    remaining constant is a uniform rescaling of the gradient and does not change
    the descent direction.

    Unlike the bilinear path there is no inpainting and no clipping -- both are
    nonlinear operations with no adjoint counterpart in bin_avg.
    """
    arr = xx_lr[None, :, :] if etan else xx_lr
    z = arr.shape[0]

    # replicate each LR cell into a factor x factor block
    xx_hr = np.repeat(np.repeat(arr.astype(np.float64), factor, axis=-2), factor, axis=-1)

    for i in range(z):
        k = i if z == nz else 0
        xx_hr[i] = xx_hr[i] * hfacc_hr[k]
    print('done (injection)')
    return xx_hr


def write_float64(fout,fld):
    with open(fout, 'wb') as f:
        np.array(fld, dtype=">f8").tofile(f)


# read in all xx files

import glob
import os

files = glob.glob(os.path.join(dirrun_pup, 'xx_*.000000' + iter + '*'))
xx_files = []
for f in files:
    # split the basename, not the full path: a '.' in a directory name
    # (e.g. .../assim_..._viscx0.25/...) would otherwise match nothing
    parts = os.path.basename(f).split('.')
    if len(parts) == 3:
        if parts[-1] == 'meta':
            xx_files.append(f[:-5])


for xx_name in xx_files:
    
    xx = rdmds(xx_name)
    prolong = inject if interp_method == 'injection' else linear_interp
    if 'etan' in xx_name:
        xx_hr = prolong(xx, etan=True)[0]
    else:
        xx_hr = prolong(xx)
    xx_hr[np.isnan(xx_hr)] = 0
    
    write_float64(dir_out + xx_name.split('/')[-1] + '.data', xx_hr)



