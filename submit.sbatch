#!/bin/bash
#SBATCH --export=ALL
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --time=00:30:00  
#SBATCH --ntasks-per-node=16 
#SBATCH --mem-per-cpu=512
#SBATCH --constraint=opath
#SBATCH --partition=scavenge
#SBATCH --mem-per-cpu=512
#SBATCH --export=ALL

# sample sbatch script for running the provided MPI program on 4 ranks.
# when customizing the file for your submissions, just vary nodes and tasks.  

# for high performance, request nodes that have Intel's ommipath interconnect
# by replacing the word OPTIONAL above with SBATCH
# I've left this constraint off by default because it limits the nodes in the cluster
# that can be used to execute the job, which will result in longer queueing times.

# run one for timing
lscpu
N=7560
time srun -n 64 ./matrixMul $N 64 1
time srun -n 49 ./matrixMul $N 49 1
time srun -n 36 ./matrixMul $N 36 1
time srun -n 25 ./matrixMul $N 25 1
time srun -n 16 ./matrixMul $N 16 1
time srun -n 9  ./matrixMul $N 9  1
time srun -n 4  ./matrixMul $N 4  1
time srun -n 72 ./matrixMul $N 72 2
time srun -n 50 ./matrixMul $N 50 2
time srun -n 32 ./matrixMul $N 32 2
time srun -n 18 ./matrixMul $N 18 2
time srun -n 8  ./matrixMul $N 8  2
time srun -n 48 ./matrixMul $N 48 3
time srun -n 27 ./matrixMul $N 27 3
time srun -n 64 ./matrixMul $N 64 4

# run one to collect a trace


