"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from utils import System, Predictor, RunnerParams
import os

boltz_prot_fasta = """>A|protein||{system.name}_prot\n{system.seq}"""
boltz_lig_fasta = """>B|smiles||{system.name}_lig\n{system.smiles}"""
boltz_fasta_template = """{input}"""

exec_command = "boltz predict $file --use_msa_server --diffusion_samples 10 --out_dir $outputs_dir"

runner_params = RunnerParams(header="clusteriqac",
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)


boltz_data = Predictor(name= "Boltz",
            prot_temp= boltz_prot_fasta,
            lig_temp= boltz_lig_fasta,
            joiner="\n",
            prot_lig_temp= boltz_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)
    
# def gen_boltz_input(system:System , predictor:Predictor)-> None:
#     """
#     Generate fasta for system Boltz-1 format.
#     out_path: path were to write the fastas

#     """

#     fasta_file = os.path.join(predictor.inputs_predictor,f"{system.name}.fasta")
#     boltz_txt = boltz_prot_fasta.format(system=system, predictor=predictor) + boltz_lig_fasta.format(system=system, predictor=predictor)

#     with open(fasta_file,"w") as bb:
#         bb.write(boltz_txt)



# boltz_runner_temp = """#!/bin/bash
# #SBATCH -J boltz
# #SBATCH -e %J.err
# #SBATCH -o %J.out
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --partition=long
# #SBATCH --gres=gpu:1

# echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"
# cwd=$(pwd)

# inputs_dir=$cwd/{predictor.inputs_predictor_unix}
# outputs_dir=$cwd/{predictor.outputs_predictor_unix}

# for input_file in $inputs_dir/*.fasta; do
#     boltz predict $input_file --use_msa_server --diffusion_samples 5 --out_dir $outputs_dir
# done

# echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"

# """
# #.format(predictor.inputs_predictor_unix, predictor.outputs_predictor_unix)



# def gen_boltz_runner(predictor:Predictor) -> None:
#     runner_file = os.path.join(predictor.runners,"Boltz_runner.sh")

#     runner_str = boltz_runner_temp.format(predictor=predictor)
    
#     with open(runner_file,"w") as rr:
#         rr.write(runner_str)


# def main(system_list: list[System], predictor: Predictor):
#     # Change current predictor in in Predictor
#     predictor.predictor = "Boltz"

#     # Create directories
#     os.makedirs(predictor.inputs_predictor,exist_ok=True)
#     os.makedirs(predictor.outputs_predictor,exist_ok=True)

#     # Generate input
#     for system in system_list:
#       gen_boltz_input(system, predictor)

#     # Generate runner
#     gen_boltz_runner(predictor)

