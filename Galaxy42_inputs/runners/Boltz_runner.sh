#!/bin/bash
#SBATCH -J boltz
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1

echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"

inputs_dir=$cwd/Galaxy42_inputs/Boltz
outputs_dir=$cwd/Galaxy42_outputs/Boltz

for input_file in $inputs_dir/*.fasta; do
    boltz predict $input_file --use_msa_server --diffusion_samples 5 --out_dir $outputs_dir
done

echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"

