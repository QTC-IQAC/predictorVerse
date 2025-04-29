#!/bin/bash
#SBATCH -J Galaxy42_RFAA
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

cwd=$(pwd)
echo $cwd
inputs_dir=$cwd/Galaxy42_inputs/RFAA

echo "Starting RFAA..."

cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed

for input_file in $inputs_dir/*.yaml; do 
    # Split input file into file name and folder
    filename=$(basename "$input_file")
    foldername=$(dirname "$input_file")

    python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername
done

cd $cwd

