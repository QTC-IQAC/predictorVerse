#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e %j.err
#SBATCH -o %j.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -n 32


cwd=$(pwd)


inputs_dir=$cwd/AF3/inputs
outputs_dir=$cwd/AF3/outputs


input_name=$1
input=$inputs_dir/$input_name

# Check if each argument is provided
if [ -z "$1" ] ; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <input.json> "
    exit 1
fi


module load alphafold/3

#### Main program ####

start_date=$(date)

singularity exec  \
      \
     --bind $inputs_dir:/root/$inputs_dir \
     --bind $outputs_dir:/root/$outputs_dir \
     --bind $cwd/models:/root/models \
     --bind /data/ddbb/alphafold3:/root/public_databases \
           /prod/container/alphafold3/alphafold3.sif \
     python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \
              --json_path=/root/$input \
              --output_dir=/root/$outputs_dir --num_recycles 13


final_date=$(date)

cd $cwd

echo "START DATE: $start_date"
echo "FINAL DATE: $final_date"

