#!/bin/bash
#SBATCH --export=ALL
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=16 
#SBATCH --cpus-per-task=1 
#SBATCH --exclusive
#SBATCH --time=01:00:00 
#SBATCH --mem-per-cpu=512
#OPTIONAL --constraint=opath

# sample sbatch script for running the provided MPI program on 4 ranks.
# when customizing the file for your submissions, just vary nodes and tasks.  

# for high performance, request nodes that have Intel's ommipath interconnect
# by replacing the word OPTIONAL above with SBATCH
# I've left this constraint off by default because it limits the nodes in the cluster
# that can be used to execute the job, which will result in longer queueing times.

# run one for timing
srun ./reduce_avg 2000

# run one to collect a trace
srun ./reduce_avg-trace 2000
