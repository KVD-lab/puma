#!/bin/bash

#SBATCH -J puma
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -A iPlant-Collabs

module load tacc-singularity

IMG="/work/05066/imicrobe/singularity/puma-1.0.0.img"

if [[ ! -e "$IMG" ]]; then
    echo "Missing Singularity image \"$IMG\""
    exit 1
fi

DATA_DIR="/work/05066/imicrobe/iplantc.org/data/puma/data_dir"
singularity exec $IMG run_puma.py "$@" -o puma_out -d ${DATA_DIR} -D info

echo "Comments to kyclark@email.arizona.edu"
