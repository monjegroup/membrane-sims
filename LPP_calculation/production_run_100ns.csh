#!/bin/bash
#SBATCH --job-name=cbdlpp
#SBATCH -o cbdlpp.o%j 
#SBATCH -N 4
#SBATCH --ntasks-per-node=32
#SBATCH --constraint=OPA
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH -t 48:00:00            # run time (hh:mm:ss)
#SBATCH --mail-user=jli436@buffalo.edu
#SBATCH --mail-type=ALL

module load foss gromacs/2021.5
module list

# production
srun -n 1 $gmx grompp -f prod_mem.mdp -o traj.tpr -c prod.gro -r prod.gro -t prod.cpt -p topol.top -maxwarn -1
srun $gmx mdrun -v -deffnm traj -dlb yes -resethway
