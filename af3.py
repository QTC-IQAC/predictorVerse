af3_json_template = """{
    "name": "{0}",
    "modelSeeds": [10, 42],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "{1}"
        }
      },
      {
        "ligand": {
          "id": "Z",
          "smiles": "{2}"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 2

  }
"""
#.format(system.name, system.seq, system.smiles)


#.format(workpath.name, workpath.inputs_predictor_unix, workpath.outputs_predictor_unix)
af3_runner_template = """#!/bin/bash
#SBATCH --job-name={0}
#SBATCH -e AF3.err
#SBATCH -o AF3.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --reservation=cuda_12.6
#SBATCH --gres=gpu:1
#SBATCH -n 32

module load alphafold/3

# alphafold3.sh /data/ucsqab/aortega/inserted_mutants


#RD=<YOUR_WORKDIR_PATH>
input_folder={1}
input_name=$2
input=$1/$2
output={2}
RD=$(pwd)


# Check if each argument is provided
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <inputs_path> <input.json> <output_path>"
    exit 1
fi


singularity exec  \
      \
     --bind $RD/$input_folder:/root/$input_folder \
     --bind $RD/$output:/root/$output \
     --bind $RD/models:/root/models \
     --bind /data/ddbb/alphafold3:/root/public_databases \
           /prod/container/alphafold3/alphafold3.sif \
     python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \
              --json_path=/root/$input \
              --output_dir=/root/$output


"""



