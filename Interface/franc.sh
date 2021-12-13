#!/bin/bash
#PBS -N franc
#PBS -q normal
#PBS -P CBBI0818
#PBS -l select=2:ncpus=24
#PBS -l walltime=10:00:00	
#PBS -o /mnt/lustre/users/egeza/LOCAL_DATA/FRANC/franc/franc.err
#PBS -e /mnt/lustre/users/egeza/LOCAL_DATA/FRANC/franc/franc.out
#PBS -V
#PBS -M ephie@aims.ac.za

###Loading modules required by FRANC providing their paths (these depend on the cluster one is using:PLEASE ADD THEM THE WAY THEY ARE ACCESSIBLE)
module add chpc/BIOMODULES
module add /apps/chpc/scripts/modules/bio/app/
module add /apps/chpc/scripts/modules/bio/app/gcc/7.2.0
module add chpc/python/2.7.12_gcc-6.1.0
module add chpc/perl/5.24.0
module add chpc/perl/perlbrew
module add chpc/R/3.2.3-gcc5.1.0
module add /apps/chpc/scripts/modules/bio/lib/png/1.6.21
module add /apps/chpc/scripts/modules/bio/lib/gnu/gsl_2.1

python %sfrancserver.py

