"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from utils import System
import os

    
def gen_boltz_input(system:System , workspace:str)-> None:
    """
    Generate fasta for system Boltz-1 format.
    out_path: path were to write the fastas

    """
    inputs_dir = os.path.join(workspace+"_inputs","Boltz")
    output_dir = os.path.join(workspace+"_outputs","Boltz")

    os.makedirs(inputs_dir,exist_ok=True)
    os.makedirs(output_dir,exist_ok=True)

    fasta_file = os.path.join(inputs_dir,f"{system.name}.fasta")
    
    with open(fasta_file,"w") as bb:
        bb.write(f">A|protein||{system.name}_prot\n{system.seq}\n")
        bb.write(f">B|smiles||{system.name}_lig\n{system.smiles}\n")



boltz_runner_temp = """#!/bin/bash
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
cwd=$(pwd)

inputs_dir=$cwd/{0}
outputs_dir=$cwd/{1}

for input_file in $inputs_dir/*.fasta; do
    boltz predict $input_file --use_msa_server --diffusion_samples 5 --out_dir $outputs_dir
done

echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"

"""
#.format(inputs_dir, output_dir)



def gen_boltz_runner(workspace:str) -> None:
    inputs_dir = workspace+"_inputs/Boltz"
    outputs_dir = workspace+"_outputs/Boltz"
    runner_file = os.path.join(workspace+"_inputs","runners","Boltz_runner.sh")

    # TODO: generalize this in utils

    runner_str = boltz_runner_temp.format(inputs_dir,outputs_dir)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)