#!/bin/bash
#SBATCH --job-name=stress20
#SBATCH -o stress20.o%j 
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --constraint=OPA
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH -t 72:00:00            # run time (hh:mm:ss)
#SBATCH --mail-user=jli436@buffalo.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0 intel/19.5 mkl/2019.5 intel-mpi/2020.1
module list

export gmx="/projects/academic/vmonje/software/gromacs-ls/gromacs-ls-2016.3/my_install/bin/gmx_LS"
export tensor="/projects/academic/vmonje/software/gromacs-ls/gromacs-4.5.5-ls-5.0/my_install/bin/tensortools"


# production
srun -n 1 $gmx mdrun -s ../traj20ns.tpr -rerun ../traj_20.trr -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsmindihang 0.0005
