#!/bin/bash 
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-gauss-mpi-intel
#@ job_type = MPICH
#@ output = par_io.out
#@ error = par_io.out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = gauss
#@ minimize_time_to_solution = yes
#@ notification = never
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
. $HOME/.bashrc

module unload mpi.ibm
module load mpi.intel

for i in `seq 1 64`; do  mpirun -n $i ./par_io ./ge_data/size8192x8192; done
