#!/bin/bash

#SBATCH -J puma
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -A iPlant-Collabs

module load tacc-singularity

IMG="/work/05066/imicrobe/singularity/puma-0.1.0.img"

if [[ ! -e "$IMG" ]]; then
    echo "Missing Singularity image \"$IMG\""
    exit 1
fi

singularity exec $IMG puma.py "$@" -o puma-out

echo "Comments to kyclark@email.arizona.edu"
