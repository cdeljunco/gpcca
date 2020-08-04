#!/bin/sh

#SBATCH --job-name=gpcca-step2klist-MATRIX
#SBATCH --partition=svaikunt
#SBATCH --constraint=ib
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=cdeljunco.midway@gmail.com
module load matlab/2016a

matlab -r "step_2_klist "../G-PCCA-Results/MATRIX/MATRIX" KLIST; exit"
