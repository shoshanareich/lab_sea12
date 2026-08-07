#!/bin/bash -x
#SBATCH -J labseaHRswot_adj
#SBATCH -o labseaHRswot_adj.%j.out
#SBATCH -e labseaHRswot_adj.%j.err
#SBATCH -t 12:00:00
#SBATCH -N 6 
#SBATCH -n 103
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


iter=0
optimext=0000
nprocs_hr=103

#--- 2.set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
datadir=/home/shoshi/MITgcm_obsfit/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_swot_HR/
#dirrun_pup=${scratchdir}/run_adlo_packunpack/

hr_dir_iter0=${scratchdir}/run_adhi_it0000

builddir=${rootdir}/build_adhi_2lev_seaice_update_mpi_32x32

#echo $(printf "%04d" $((iter)))  
#source $(conda info --base)/etc/profile.d/conda.sh
#conda activate py38
#python3 ${rootdir}/mappings/interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter))) $dirrun_pup
#conda deactivate
#

ext2=$(printf "%04d" $iter)


workdir_hi=${scratchdir}/run_adhi_it${ext2}


cd $workdir_hi
  

set -x

date > run.MITGCM.timing
mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout
date >> run.MITGCM.timing


cd ..


