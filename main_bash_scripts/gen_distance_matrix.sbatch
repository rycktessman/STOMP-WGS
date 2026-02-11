#!/bin/bash

#SBATCH --output=%A.out
#SBATCH --error=%A.err 
#SBATCH --time=12:00:00 
#SBATCH -p shared
#SBATCH --mem=100000
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tryckma1@jhmi.edu

module purge
ml av
ml gcc/9.3.0
ml r/4.0.2

#Rscript code/find_low_readdepths.R
#Rscript code/calculate-pairs-dist.R
Rscript code/calculate-pairs-dist-tSNPs-alignments.R





