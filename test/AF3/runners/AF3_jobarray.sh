#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e %j.err
#SBATCH -o %j.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -n 32
#SBATCH --array=1-1%2

case $SLURM_ARRAY_TASK_ID in
1) ./AF3/runners/AF3_runner.sh carboxamide_short.json ;;
esac