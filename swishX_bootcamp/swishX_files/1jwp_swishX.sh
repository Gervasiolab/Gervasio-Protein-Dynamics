#!/bin/bash -l

# Job settings
#SBATCH --job-name=tem1_swishX           # Job name
#SBATCH --time=12:00:00                   # Walltime limit (hh:mm:ss)
#SBATCH --nodes=6                         # Number of nodes
#SBATCH --ntasks-per-node=1               # Number of tasks per node
#SBATCH --cpus-per-task=12                # Number of CPU cores per task
#SBATCH --partition=normal                # Partition/queue name
#SBATCH --constraint=gpu                  # Specify GPU constraint

# Environment variables
export OMP_NUM_THREADS=6                  # Set number of OpenMP threads
export CRAY_CUDA_MPS=1                    # Enable CUDA Multi-Process Service (MPS)

# Load user-specific configurations and modules
source $HOME/.bashrc

# SWISHX
# Check if the production run files exist, if not, generate them
if [ ! -e rep_0/run.tpr ]; then
  for i in $(seq 0 5); do
    cd rep_$i
    gmx grompp -f ../prod.mdp -p 1jwp_swish$i.top -c npt2.gro -o run.tpr -n ../1jwp_benz.ndx
    cd ..
  done
fi

# Perform the production run if the output files are not present
if [ ! -e rep_0/run.gro ]; then
  srun gmx_mpi mdrun -deffnm run -plumed ../plumed.dat -hrex -multidir $(ls -d rep_*/) -replex 5000 -cpi -maxh 11.9 -pin on
fi

# Resubmit script if unfinished
if [ ! -e rep_0/run.gro ]; then
  # Check if the run files have been generated, modify plumed.dat accordingly
  if [ -e rep_0/DeltaFs.0.data ]; then
    sed -i "s/#RESTART/RESTART/g" plumed.dat
  fi
  # Resubmit the job script
  sbatch 1jwp_swishX.sh
fi
