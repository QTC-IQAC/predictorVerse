import os
from utils import System


def gen_of_input(system_list:list[System], workspace:str)-> None:
    """
    Generate a fasta for all the proteins in system.
    OF is only a predictor for proteins. 
    """
    inputs_dir = os.path.join(workspace+"_inputs","OmegaFold")
    output_dir = os.path.join(workspace+"_outputs","OmegaFold")

    os.makedirs(inputs_dir,exist_ok=True)
    os.makedirs(output_dir,exist_ok=True)
    
    fasta_file = os.path.join(inputs_dir,"OF_input.fasta")

    with open(fasta_file, "w") as oo:
        for system in system_list:
            oo.write(f">{system.name}_prot\n")
            oo.write(f"{system.seq}\n")

of_runner_temp = """#!/bin/bash
#SBATCH -J OF
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

input={0}/OF_input.fasta
output_dir={1}

omegafold $input $output_dir

"""
#.format(inputs_dir,outputs_dir)


def gen_of_runner(workspace:str) -> None:
    inputs_dir = workspace+"_inputs/OmegaFold"
    outputs_dir = workspace+"_outputs/OmegaFold"
    runner_file = os.path.join(workspace+"_inputs","runners","OF_runner.sh")

    # TODO: generalize this in utils

    runner_str = of_runner_temp.format(inputs_dir,outputs_dir)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)
   