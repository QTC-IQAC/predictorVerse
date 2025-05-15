"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
import os
from utils import System, Predictor, RunnerParams

chai_prot_fasta = """>protein|{system.name}_prot\n{system.seq}"""
chai_lig_fasta = """>ligand|{system.name}_lig\n{system.smiles}"""
chai_fasta_template = """{input}"""

exec_command = """# Get name of system

system=$(basename $file | sed "s/.fasta//g")
mkdir -p $outputs_dir/$system
chai-lab fold --seed 42 --num-diffn-samples 10 $file $outputs_dir/$system
"""

runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)

chai_data = Predictor(name = "Chai-1",
            prot_temp= chai_prot_fasta,
            lig_temp= chai_lig_fasta,
            joiner="\n",
            prot_lig_temp= chai_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)



# chai_runner_template = """ #!/bin/bash

# inputs_folder={0}
# outputs_folder={1}

# for file in $inputs_folder/*.fasta; 
#     do
#     # Get name of system
#     system=$(basename $file | sed "s/.fasta//g")
#     mkdir -p $outputs_folder/$system
#     chai-lab fold --seed 42 --num-diffn-samples 10 $file $outputs_folder/$system
# done

# """
##################### RUNNER PARAMS #######################



"""
header = False # depends on local or cluster. Cluster should have an option to choose cluster header
extra_cmds = False # If there are extra commands
extra_inputs = False # Only True for job arrays
looper = True # True if looper indicat, false if not (in certain cases and in jobarrays)
"""

# def gen_chai_runner(predictor:Predictor) -> None:
#     runner_file = os.path.join(predictor.runners,"Chai_runner.sh")

#     runner_str = chai_runner_template.format(predictor.inputs_predictor_unix, predictor.outputs_predictor_unix)
    
#     with open(runner_file,"w") as rr:
#         rr.write(runner_str)


# def main(system_list:System, predictor:Predictor)-> None:
#     # Change current predictor in in Predictor
#     predictor.predictor = "Chai"

#     # Create directories
#     os.makedirs(predictor.inputs_predictor,exist_ok=True)
#     os.makedirs(predictor.outputs_predictor,exist_ok=True)

#     # Generate input (the default fasta file)
#     for system in system_list:
#         gen_fasta(system, predictor.inputs_predictor,mode=None)

#     # Generate runner
#     gen_chai_runner(predictor)

