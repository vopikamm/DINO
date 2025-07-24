#!/bin/bash
#
#SBATCH --job-name=DINO_demonstrator
#SBATCH --ntasks=13                  # Total MPI tasks (NEMO + XIOS)
#SBATCH --time=2:00:00               # Job time limit (HH:MM:SS or MM)
#SBATCH --output=slurm-%j.out        # Stdout log
#SBATCH --error=slurm-%j.err         # Stderr log

# Load necessary modules (edit as needed)
module purge

module load intel/2021.4.0
module load openmpi/4.0.7
module load hdf5/1.10.7-mpi
module load netcdf-fortran/4.5.3-mpi
module load netcdf-c/4.7.4-mpi

# Enable debug output

# Launch NEMO (12 ranks) and XIOS server (1 rank)
mpirun -np 12 ./nemo -np 1 /home/dkamm/XIOS3/bin/xios_server.exe