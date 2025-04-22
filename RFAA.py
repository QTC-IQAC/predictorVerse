import sys
import os
import utils as uu
from utils import System


lig_prot_yaml =""" defaults:
 - base
 - _self_
 job_name: {0}
 output_path: {2}
 protein_inputs:
   A:
     prot_file: {1}/{0}_prot.fasta
 sm_inputs:
   B:
     input: {1}/{0}_lig.fasta
     input_type: 'smiles'
"""
#.format(system.name, fastas_dir, output_dir)



def gen_RFAA_input(system: System, workspace:str) -> None:
    """
    Generate RFAA yaml inputs


    # save input as workspace+"_inputs/RFAA/system.name"
    # save output in workspace+"_outputs/RFAA/

    Also needs a FUCKING fasta for each protein and ligand
    """
    
    # Handle input/output directories
    fastas_dir = os.path.join(workspace+"_inputs","RFAA","fastas")
    yamls_dir = os.path.join(workspace+"_inputs","RFAA")
    output_dir = os.path.join(workspace+"_outputs","RFAA")
    
    # TODO: generalize this in utils
    os.makedirs(fastas_dir,exist_ok=True)
    os.makedirs(yamls_dir,exist_ok=True)
    os.makedirs(output_dir,exist_ok=True)

    # Generate a fasta for protein and ligand (in fasta folder in RFAA)
    uu.gen_fasta(system, fastas_dir, mode="protein")
    uu.gen_fasta(system, fastas_dir, mode="ligand")

    # Generate yaml
    yaml_file = os.path.join(yamls_dir, system.name+".yaml")
    yaml_str = lig_prot_yaml.format(system.name, fastas_dir, output_dir)
    with open(yaml_file,"w") as yy:
        yy.write(yaml_str)


#TODO make a gestor for slurm options
runner_temp ="""#!/bin/bash
#SBATCH -J {0}_RFAA
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

cwd=$(pwd)
echo $cwd
inputs_dir=$cwd/{1}

echo "Starting RFAA..."

cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed

for input_file in $inputs_dir/*.yaml; do 
    # Split input file into file name and folder
    filename=$(basename "$input_file")
    foldername=$(dirname "$input_file")

    python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername
done

cd $cwd

"""
#.format(workspace)

def gen_RFAA_runner(workspace:str) -> None:
    runner_file = os.path.join(workspace+"_inputs","runners","RFAA_runner.sh")
    rfaa_dir = workspace+"_inputs/RFAA" #done manually because this is for cluster in Linux

    runner_str = runner_temp.format(workspace,rfaa_dir)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)