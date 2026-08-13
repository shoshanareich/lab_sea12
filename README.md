based on MITgcm c69j 

to submit a job on any of sverdrup, stampede3, or pfe:

./submit.sh run_optimization.sh      # detects machine → sbatch or qsub
This calls submit_slurm_sverdrup.sh, submit_slurm_stampede3.sh, or submit_pbs_pfe.sh


OR:
# sverdrup / stampede3
sbatch --export=ALL,DRIVER=run_optimization.sh submit_slurm_sverdrup.sh

# pfe
qsub -v DRIVER=run_optimization.sh submit_pbs_pfe.sh

other helpful options:
./submit.sh -n run_optimization.sh   # dry run
./submit.sh run_swot_optimization.sh -t 12:00:00   # extra args passed to sbatch/qsub
LABSEA_MACHINE=pfe ./submit.sh -n run_optimization.sh   # test another machine's config
./run_optimization.sh                # directly, inside an existing allocation
./submit.sh                                 # no driver → lists available ones



to check the linking:
./check_staging.sh                    # full check against this machine
./check_staging.sh -s link_lores.sh   # one script
