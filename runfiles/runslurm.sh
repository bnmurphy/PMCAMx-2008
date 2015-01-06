#!/bin/bash
#
#SBATCH -J PMCAMx-2008
#SBATCH -t 2-3:00:00
#SBATCH --mem=4000
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=benjamin.murphy@itm.su.se
#SBATCH --constraint=vtune
#SBATCH -v
#
        # Run a single task in the foreground
        ./run_days_EastUS --verbose >&PMCAMx_LOG.$$ log.$$
#
#end
