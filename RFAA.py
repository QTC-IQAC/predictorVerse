import sys
import os
import utils as uu
from utils import System, Workpath


lig_prot_yaml ="""defaults:
- base
- _self_
job_name: {0}
output_path: {2}
protein_inputs:
  A:
    fasta_file: {1}/{0}_prot.fasta
sm_inputs:
  B:
    input: {1}/{0}_lig.sdf
    input_type: 'sdf'
"""
#.format(system.name, fastas_dir, output_dir)



def gen_RFAA_input(system: System, workspace:Workpath, paths:list) -> None:
    """
    Generate RFAA yaml inputs


    # save input as workspace+"_inputs/RFAA/system.name"
    # save output in workspace+"_outputs/RFAA/

    Also needs a FUCKING fasta for each protein and ligand
    """
    
    # Handle input/output directories
    fastas_dir, fastas_dir_unix, output_dir_unix = paths

    # Generate a fasta for protein and and sdf for the ligand (in fasta folder in RFAA)
    uu.gen_fasta(system, fastas_dir, mode="protein")
    uu.lig_smiles_to_sdf(system, fastas_dir)

    # Generate yaml
    yaml_file = os.path.join(workspace.inputs_predictor, system.name+".yaml")
    yaml_str = lig_prot_yaml.format(system.name, fastas_dir_unix, output_dir_unix)
    with open(yaml_file,"w") as yy:
        yy.write(yaml_str)


# TODO make a gestor for slurm options
#TODO make it a job array
runner_temp ="""#!/bin/bash
#SBATCH -J RFAA
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

cwd=$(pwd)
echo $cwd
inputs_dir=$cwd/{0}

echo "Starting RFAA..."

cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed

for input_file in $inputs_dir/*.yaml; do 
    # Change the relative path for absolute ones (will need to run this from the folder before Galaxy42_inputs)
    sed -i "s|\./|$cwd/|g" $file


    # Split input file into file name and folder
    filename=$(basename "$input_file")
    foldername=$(dirname "$input_file")

    python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername
done

cd $cwd

"""
#.format(workspace)

def gen_RFAA_runner(workspace:Workpath) -> None:
    runner_file = os.path.join(workspace.runners,"RFAA_runner.sh")

    runner_str = runner_temp.format(workspace.inputs_predictor_unix)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)


def main(system_list:System, workpath:Workpath):
    # Change current predictor in in Workpath
    workpath.predictor = "RFAA"

    # Create directories
    os.makedirs(workpath.inputs_predictor,exist_ok=True)
    os.makedirs(workpath.outputs_predictor,exist_ok=True)

    fastas_dir = os.path.join(".",workpath.inputs_predictor,"fastas")
    fastas_dir_unix = "./"+workpath.inputs_predictor_unix+"/fastas"
    output_dir_unix = "./"+workpath.outputs_predictor_unix

    path_list = [fastas_dir,fastas_dir_unix,output_dir_unix]
    
    os.makedirs(fastas_dir,exist_ok=True)

    # Generate input
    for system in system_list:
      gen_RFAA_input(system, workpath, path_list)

    # Generate runner
    gen_RFAA_runner(workpath)
    