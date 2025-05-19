import os
from utils import System, Predictor, RunnerParams


# Double braces to not confuse with {0}

af3_prot_json = """      {{
        "protein": {{
          "id": "A",
          "sequence": "{system.seq}"
        }}
      }}"""



af3_lig_json = """      {{
        "ligand": {{
          "id": "B",
          "smiles": "{system.smiles}"
        }}
      }}"""


af3_json_template = """{{
    "name": "{system.name}",
    "modelSeeds": [10, 42],
    "sequences": [
      {input}
    ],
    "dialect": "alphafold3",
    "version": 2

}}
"""

extra_cmds = "module load alphafold/3"

exec_command = """singularity exec  \\
      \\
     --bind $inputs_dir:/root/$inputs_dir \\
     --bind $outputs_dir:/root/$outputs_dir \\
     --bind $cwd/models:/root/models \\
     --bind /data/ddbb/alphafold3:/root/public_databases \\
           /prod/container/alphafold3/alphafold3.sif \\
     python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \\
              --json_path=/root/$input \\
              --output_dir=/root/$outputs_dir
"""

runner_params = RunnerParams(header="csuc",
                             extra_cmds=True,
                             extra_inputs=True,
                             looper=False,
                             jobarray=True)


af3_data = Predictor(name= "AF3",
            prot_temp= af3_prot_json,
            lig_temp= af3_lig_json,
            joiner=",\n",
            prot_lig_temp= af3_json_template,
            input_extension= ".json",
            extra_cmds=extra_cmds,
            main_cmds=exec_command,
            runner_params=runner_params
)


#.format(system.name, system.seq, system.smiles)


# af3_runner_template = """#!/bin/bash
# #SBATCH --job-name=AF3
# #SBATCH -e AF3.err
# #SBATCH -o AF3.log
# #SBATCH -t 00-01:00
# #SBATCH -p gpu
# #SBATCH --reservation=cuda_12.6
# #SBATCH --gres=gpu:1
# #SBATCH -n 32

# module load alphafold/3


# #RD=<YOUR_WORKDIR_PATH>
# input_folder={0} # Relative path
# input_name=$1
# input=$input_folder/$input_name
# output={1} # Relative path
# RD=$(pwd)


# # Check if each argument is provided
# if [ -z "$1" ] ; then
#     echo "Error: Missing arguments."
#     echo "Usage: $0 <input.json> "
#     exit 1
# fi


# singularity exec  \\
#       \\
#      --bind $RD/$input_folder:/root/$input_folder \\
#      --bind $RD/$output:/root/$output \\
#      --bind $RD/models:/root/models \\
#      --bind /data/ddbb/alphafold3:/root/public_databases \\
#            /prod/container/alphafold3/alphafold3.sif \\
#      python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \\
#               --json_path=/root/$input \\
#               --output_dir=/root/$output


# """
# Double / to avoid writing in the same line everything and make it more readable
#.format(predictor.inputs_predictor_unix, predictor.outputs_predictor_unix)

# af3_jobarray_template="""#!/bin/bash
# #SBATCH --job-name=AF3
# #SBATCH -e AF3.err
# #SBATCH -o AF3.log
# #SBATCH -t 00-01:00
# #SBATCH -p gpu
# #SBATCH --reservation=cuda_12.6
# #SBATCH --gres=gpu:1
# #SBATCH -n 32
# #SBATCH --array=1-{0}%{1}
# """
#.format(num_of_jobs,max_num_of_jobs_running_at_the_same_time)




# def gen_af3_input(system:System, predictor:Predictor)-> None:
#     """
#     Generate a fasta for all the proteins in system.
#     """

#     input_file = os.path.join(predictor.inputs_predictor,system.name+".json")
#     af3_str = af3_json_template.format(system.name, system.seq, system.smiles)
    
#     with open(input_file, "w") as aa:
#           aa.write(af3_str)



# def gen_af3_runner(predictor:Predictor) -> None:
#     runner_file = os.path.join(predictor.runners,"AF3_runner.sh")

#     runner_str = af3_runner_template.format(predictor.inputs_predictor_unix, predictor.outputs_predictor_unix)
    
#     with open(runner_file,"w") as rr:
#         rr.write(runner_str)



# def gen_af3_jobarray(system_list:list[System], predictor:Predictor, max_cap_jobs=None | int) -> None:
#     """
    
#     max_cap_jobs : max number of jobs that can run at the same time.
#                   If None, set to number of jobs
#     """
#     # Set variables
#     num_jobs = len(system_list)
#     num_cap_jobs = num_jobs if max_cap_jobs is None else max_cap_jobs
#     af3_jobarr_header = af3_jobarray_template.format(num_jobs, num_cap_jobs)

#     # Runner file relative path
#     runner_file_rel_path = os.path.join(".",predictor.runners,"AF3_runner.sh")

#     # Open jobarr file
#     jobarr_file = os.path.join(predictor.runners,"AF3_jobarray.sh")
#     with open(jobarr_file, "w") as jobarr:
#         # Write header with correct parameters
#         jobarr.write(af3_jobarr_header)

#         # Write case start
#         jobarr.write("case $SLURM_ARRAY_TASK_ID in\n")

#         # Start looping through systems
#         for ii, system in enumerate(system_list):
#             jobarr.write(f"{ii+1}) {runner_file_rel_path} {system.name}.json ;;\n")
        
#         jobarr.write("esac")
            


# def main(system_list:list[System], predictor:Predictor):
#     # Change current predictor in in Predictor
#     predictor.predictor = "AF3"

#     # Create directories
#     os.makedirs(predictor.inputs_predictor,exist_ok=True)
#     os.makedirs(predictor.outputs_predictor,exist_ok=True)

#     # Generate input
#     for system in system_list:
#       gen_af3_input(system, predictor)

#     # Generate runner
#     gen_af3_runner(predictor)

#     # Generate jobarray
#     gen_af3_jobarray(system_list,predictor,max_cap_jobs=None)


