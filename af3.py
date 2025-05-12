import os
from utils import System, Workspace


# Double braces to not confuse with {0}
af3_json_template = """{{
    "name": "{0}",
    "modelSeeds": [10, 42],
    "sequences": [
      {{
        "protein": {{
          "id": "A",
          "sequence": "{1}"
        }}
      }},
      {{
        "ligand": {{
          "id": "B",
          "smiles": "{2}"
        }}
      }}
    ],
    "dialect": "alphafold3",
    "version": 2

}}
"""
#.format(system.name, system.seq, system.smiles)


af3_runner_template = """#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e AF3.err
#SBATCH -o AF3.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --reservation=cuda_12.6
#SBATCH --gres=gpu:1
#SBATCH -n 32

module load alphafold/3


#RD=<YOUR_WORKDIR_PATH>
input_folder={0} # Relative path
input_name=$1
input=$input_folder/$input_name
output={1} # Relative path
RD=$(pwd)


# Check if each argument is provided
if [ -z "$1" ] ; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <input.json> "
    exit 1
fi


singularity exec  \\
      \\
     --bind $RD/$input_folder:/root/$input_folder \\
     --bind $RD/$output:/root/$output \\
     --bind $RD/models:/root/models \\
     --bind /data/ddbb/alphafold3:/root/public_databases \\
           /prod/container/alphafold3/alphafold3.sif \\
     python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \\
              --json_path=/root/$input \\
              --output_dir=/root/$output


"""
# Double / to avoid writing in the same line everything and make it more readable
#.format(workspace.inputs_predictor_unix, workspace.outputs_predictor_unix)

af3_jobarray_template="""#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH -e AF3.err
#SBATCH -o AF3.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --reservation=cuda_12.6
#SBATCH --gres=gpu:1
#SBATCH -n 32
#SBATCH --array=1-{0}%{1}



"""
"""
case $SLURM_ARRAY_TASK_ID in
1) ./jobscript0 ;;
2) ./jobscript1 ;;
3) ./jobscript2 ;;
4) ./jobscript3 ;;
5) ./jobscript4 ;;
6) ./jobscript5 ;;
7) ./jobscript6 ;;
8) ./jobscript7 ;;
9) ./jobscript8 ;;
10) ./jobscript9 ;;
11) ./jobscript10 ;;
12) ./jobscript11 ;;
esac



"""
#.format(num_of_jobs,max_num_of_jobs_running_at_the_same_time)


def gen_af3_input(system:System, workspace:Workspace)-> None:
    """
    Generate a fasta for all the proteins in system.
    OF is only a predictor for proteins. 
    """

    input_file = os.path.join(workspace.inputs_predictor,system.name+".json")
    af3_str = af3_json_template.format(system.name, system.seq, system.smiles)
    
    with open(input_file, "w") as aa:
          aa.write(af3_str)



def gen_af3_runner(workspace:Workspace) -> None:
    runner_file = os.path.join(workspace.runners,"AF3_runner.sh")

    runner_str = af3_runner_template.format(workspace.inputs_predictor_unix, workspace.outputs_predictor_unix)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)



def gen_af3_jobarray(system_list:list[System], workspace:Workspace, max_cap_jobs=None | int) -> None:
    """
    
    max_cap_jobs : max number of jobs that can run at the same time.
                  If None, set to number of jobs
    """
    # Set variables
    num_jobs = len(system_list)
    num_cap_jobs = num_jobs if max_cap_jobs is None else max_cap_jobs
    af3_jobarr_header = af3_jobarray_template.format(num_jobs, num_cap_jobs)

    # Runner file relative path
    runner_file_rel_path = os.path.join(".",workspace.runners,"AF3_runner.sh")

    # Open jobarr file
    jobarr_file = os.path.join(workspace.runners,"AF3_jobarray.sh")
    with open(jobarr_file, "w") as jobarr:
        # Write header with correct parameters
        jobarr.write(af3_jobarr_header)

        # Write case start
        jobarr.write("case $SLURM_ARRAY_TASK_ID in\n")

        # Start looping through systems
        for ii, system in enumerate(system_list):
            jobarr.write(f"{ii+1}) {runner_file_rel_path} {system.name}.json ;;\n")
        
        jobarr.write("esac")
            


def main(system_list:list(System), workspace:Workspace):
    # Change current predictor in in Workspace
    workspace.predictor = "AF3"

    # Create directories
    os.makedirs(workspace.inputs_predictor,exist_ok=True)
    os.makedirs(workspace.outputs_predictor,exist_ok=True)

    # Generate input
    for system in system_list:
      gen_af3_input(system, workspace)

    # Generate runner
    gen_af3_runner(workspace)

    # Generate jobarray
    gen_af3_jobarray(system_list,workspace,max_cap_jobs=None)
