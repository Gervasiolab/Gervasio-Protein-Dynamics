#!/bin/bash -l

#SBATCH --job-name=1jwp_r0
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=6
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=pr126

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

# Load your version of Gromacs
module load daint-gpu
source $HOME/software/gromacs-2021.3/gromacs/bin/GMXRC

i='0'
cd rep_$i

# EM
if [ ! -e em.tpr ]; then
  
  gmx grompp -f ../em.mdp -p ../1jwp.top -c ../1jwp_ions.gro -r ../1jwp_ions.gro -o em.tpr

fi
if [ ! -e em.gro ]; then
  
  srun gmx_mpi mdrun -deffnm em -pin on
  
fi
# NVT
if [ ! -e nvt.tpr ]; then
    
    gmx grompp -f ../nvt.mdp -p ../1jwp.top -c em.gro -r em.gro -o nvt.tpr -n ../1jwp.ndx 
    
fi
if [ ! -e nvt.gro ]; then
  
  srun gmx_mpi mdrun -deffnm nvt -notunepme -cpi -dlb yes -pin on
  
fi
# NPT
if [ ! -e npt.tpr ]; then
   
    gmx grompp -f ../npt.mdp -p ../1jwp.top -c nvt.gro -r nvt.gro -o npt.tpr -n ../1jwp.ndx 

fi
if [ ! -e npt.gro ]; then
  
  srun gmx_mpi mdrun -deffnm npt -notunepme -cpi -dlb yes -pin on

fi
# NPT2
if [ ! -e npt2.tpr ]; then
    
    gmx grompp -f ../npt2.mdp -p ../1jwp.top -c npt.gro -r npt.gro -o npt2.tpr -n ../1jwp.ndx   
    
fi
if [ ! -e npt2.gro ]; then

  srun gmx_mpi mdrun -deffnm npt2 -notunepme -cpi -dlb yes -pin on

fi
# PRODUCTION
if [ ! -e run.tpr ]; then
    
    gmx grompp -f ../prod.mdp -p ../1jwp.top -c npt2.gro -o run.tpr -n ../1jwp.ndx
  
fi
if [ ! -e run.gro ]; then
 
 srun gmx_mpi mdrun -deffnm run -cpi -dlb yes -maxh 11.9 -pin on

fi

# Resubmit script if unfinished
if [ ! -e run.gro ]; then
cd ../
sbatch unbiased_r0_1jwp.sh 
fi
