"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from utils import System, Workpath
import os

boltz_prot_fasta = """>A|protein||{system.name}_prot
{system.seq}
"""

boltz_lig_fasta = """>B|smiles||{system.name}_lig
{system.smiles}
"""
    
def gen_boltz_input(system:System , workpath:Workpath)-> None:
    """
    Generate fasta for system Boltz-1 format.
    out_path: path were to write the fastas

    """

    fasta_file = os.path.join(workpath.inputs_predictor,f"{system.name}.fasta")
    boltz_txt = boltz_prot_fasta.format(system=system, workpath=workpath) + boltz_lig_fasta.format(system=system, workpath=workpath)

    with open(fasta_file,"w") as bb:
        bb.write(boltz_txt)



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

inputs_dir=$cwd/{workpath.inputs_predictor_unix}
outputs_dir=$cwd/{workpath.outputs_predictor_unix}

for input_file in $inputs_dir/*.fasta; do
    boltz predict $input_file --use_msa_server --diffusion_samples 5 --out_dir $outputs_dir
done

echo "%%%%%%%%%%%% $(date) %%%%%%%%%%%%%"

"""
#.format(workspace.inputs_predictor_unix, workspace.outputs_predictor_unix)



def gen_boltz_runner(workpath:Workpath) -> None:
    runner_file = os.path.join(workpath.runners,"Boltz_runner.sh")

    runner_str = boltz_runner_temp.format(workpath=workpath)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)


def main(system_list: list[System], workpath: Workpath):
    # Change current predictor in in Workpath
    workpath.predictor = "Boltz"

    # Create directories
    os.makedirs(workpath.inputs_predictor,exist_ok=True)
    os.makedirs(workpath.outputs_predictor,exist_ok=True)

    # Generate input
    for system in system_list:
      gen_boltz_input(system, workpath)

    # Generate runner
    gen_boltz_runner(workpath)