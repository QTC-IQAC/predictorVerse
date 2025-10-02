#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e %j.err
#SBATCH -o %j.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -n 32
#SBATCH --array=1-4%2

case $SLURM_ARRAY_TASK_ID in
1) ./AF3/runners/AF3_runner.sh sys1.json ;;
2) ./AF3/runners/AF3_runner.sh sys2.json ;;
3) ./AF3/runners/AF3_runner.sh sys3.json ;;
4) ./AF3/runners/AF3_runner.sh sys4.json ;;
esac