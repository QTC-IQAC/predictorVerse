import os
from utils import System, Workpath


def gen_of_input(system_list:list[System], workspace:Workpath)-> None:
    """
    Generate a fasta for all the proteins in system.
    OF is only a predictor for proteins. 
    """

    fasta_file = os.path.join(workspace.inputs_predictor,"OF_input.fasta")

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


def gen_of_runner(workspace:Workpath) -> None:
    
    runner_file = os.path.join(workspace.runners,"OF_runner.sh")


    runner_str = of_runner_temp.format(workspace.inputs_predictor_unix,
                                       workspace.outputs_predictor_unix)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)

def main(system_list:System, workpath:Workpath):
    # Change current predictor in in Workpath
    workpath.predictor = "OmegaFold"

    # Create directories
    os.makedirs(workpath.inputs_predictor,exist_ok=True)
    os.makedirs(workpath.outputs_predictor,exist_ok=True)

    # Generate input
    gen_of_input(system_list, workpath)

    # Generate runner
    gen_of_runner(workpath)
