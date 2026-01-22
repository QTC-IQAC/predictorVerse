#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e %j.err
#SBATCH -o %j.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -n 32
#SBATCH --array=1-11%2

case $SLURM_ARRAY_TASK_ID in
1) ./AF3/runners/AF3_runner.sh carboxamide.json ;;
2) ./AF3/runners/AF3_runner.sh Gly2NDI.json ;;
3) ./AF3/runners/AF3_runner.sh HCatNDI.json ;;
4) ./AF3/runners/AF3_runner.sh MeCatPDI.json ;;
5) ./AF3/runners/AF3_runner.sh Ala2PDI.json ;;
6) ./AF3/runners/AF3_runner.sh PEG72PDI.json ;;
7) ./AF3/runners/AF3_runner.sh Gly2NDI_short.json ;;
8) ./AF3/runners/AF3_runner.sh HCatNDI_short.json ;;
9) ./AF3/runners/AF3_runner.sh MeCatPDI_short.json ;;
10) ./AF3/runners/AF3_runner.sh Ala2PDI_short.json ;;
11) ./AF3/runners/AF3_runner.sh PEG72PDI_short.json ;;
esac