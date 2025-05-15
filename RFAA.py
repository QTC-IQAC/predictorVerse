import sys
import os
import utils as uu
from utils import System, Predictor, RunnerParams

from rdkit import Chem # main tools
from rdkit.Chem import AllChem # additional tools, including 3D



rfaa_prot_template = """protein_inputs:
  A:
    fasta_file: ./{predictor.inputs_unix}/fastas/{system.name}_prot.fasta"""

rfaa_lig_template = """sm_inputs:
  B:
    input: ./{predictor.inputs_unix}/fastas/{system.name}_lig.sdf
    input_type: 'sdf'"""

rfaa_yaml_template ="""defaults:
- base
- _self_
job_name: {system.name}
output_path: ./{predictor.outputs_unix}
{input}
"""
#.format(system.name, fastas_dir, output_dir)

extra_cmds = "cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed"

exec_command = """# Change the relative path for absolute ones (will need to run this from the folder before Galaxy42_inputs)
sed -i "s|\./|$cwd/|g" $file


# Split input file into file name and folder
filename=$(basename "$file")
foldername=$(dirname "$file")

python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername

"""


# def gen_RFAA_input(system: System, predictor:Predictor, paths:list) -> None:
#     """
#     Generate RFAA yaml inputs


#     # save input as predictor+"_inputs/RFAA/system.name"
#     # save output in predictor+"_outputs/RFAA/

#     Also needs a FUCKING fasta for each protein and ligand
#     """
    
#     # Handle input/output directories
#     fastas_dir, fastas_dir_unix, output_dir_unix = paths

#     # Generate a fasta for protein and and sdf for the ligand (in fasta folder in RFAA)
#     uu.gen_fasta(system, fastas_dir, mode="protein")
#     uu.lig_smiles_to_sdf(system, fastas_dir)

#     # Generate yaml
#     yaml_file = os.path.join(predictor.inputs_predictor, system.name+".yaml")
#     yaml_str = rfaa_yaml_template.format(system.name, fastas_dir_unix, output_dir_unix)
    
#     with open(yaml_file,"w") as yy:
#         yy.write(yaml_str)


# # TODO make a gestor for slurm options
# # TODO make it a job array
# runner_temp ="""#!/bin/bash
# #SBATCH -J RFAA
# #SBATCH -e %J.err
# #SBATCH -o %J.out
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --partition=long
# #SBATCH --gres=gpu:1
# ##SBATCH --constraint=gpu8

# cwd=$(pwd)
# echo $cwd
# inputs_dir=$cwd/{0}

# echo "Starting RFAA..."

# cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed

# for input_file in $inputs_dir/*.yaml; do 
#     # Change the relative path for absolute ones (will need to run this from the folder before Galaxy42_inputs)
#     sed -i "s|\./|$cwd/|g" $file


#     # Split input file into file name and folder
#     filename=$(basename "$input_file")
#     foldername=$(dirname "$input_file")

#     python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername
# done

# cd $cwd

# """
# #.format(predictor)

# def gen_RFAA_runner(predictor:Predictor) -> None:
#     runner_file = os.path.join(predictor.runners,"RFAA_runner.sh")

#     runner_str = runner_temp.format(predictor.inputs_predictor_unix)
    
#     with open(runner_file,"w") as rr:
#         rr.write(runner_str)


# def main(system_list:list[System], predictor:Predictor):
#     # Change current predictor in in Predictor
#     predictor.predictor = "RFAA"

#     # Create directories
#     os.makedirs(predictor.inputs_predictor,exist_ok=True)
#     os.makedirs(predictor.outputs_predictor,exist_ok=True)

#     fastas_dir = os.path.join(".",predictor.inputs_predictor,"fastas")
#     fastas_dir_unix = "./"+predictor.inputs_predictor_unix+"/fastas"
#     output_dir_unix = "./"+predictor.outputs_predictor_unix

#     path_list = [fastas_dir,fastas_dir_unix,output_dir_unix]
    
#     os.makedirs(fastas_dir,exist_ok=True)

#     # Generate input
#     for system in system_list:
#       gen_RFAA_input(system, predictor, path_list)

#     # Generate runner
#     gen_RFAA_runner(predictor)


def gen_prot_fasta(system: System, predictor: Predictor):
    """
    Generate fasta file for the protein
    """
    # Create fastas folder
    fastas_dir = os.path.join(predictor.inputs,"fastas")
    os.makedirs(fastas_dir,exist_ok=True)

    # Generate protein fasta file
    fasta_file = os.path.join(predictor.inputs,"fastas",f"{system.name}_prot.fasta")
    
    with open(fasta_file,"w") as fastaff:
        fastaff.write(f">protein|{system.name}_prot\n{system.seq}\n")


def lig_smiles_to_sdf(system:System, predictor: Predictor):
    """
    From ligand smiles generate a sdf file 
    """
    sdf_file = os.path.join(predictor.inputs,"fastas",f"{system.name}_lig.sdf")
    
    mol = Chem.MolFromSmiles(system.smiles) # initialize molecule

    mol = Chem.AddHs(mol) # adding explicit Hs for 3D generation
    cid = AllChem.EmbedMolecule(mol) # returns the id of the generated conformer,
                                    # and -1 if no conformers were generated

    AllChem.MMFFOptimizeMolecule(mol) # optimize molecule with MMFF94
    writer = Chem.SDWriter(sdf_file) # Write in .sdf file
    writer.write(mol)

runner_params = RunnerParams(header="clusteriqac",
                             extra_cmds=True,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)

    
rfaa_data = Predictor(name= "RFAA",
            prot_temp= rfaa_prot_template,
            lig_temp= rfaa_lig_template,
            joiner="\n",
            prot_lig_temp= rfaa_yaml_template,
            input_extension= ".yaml",
            other_funcs=[gen_prot_fasta, lig_smiles_to_sdf],
            extra_cmds=extra_cmds,
            main_cmds=exec_command,
            runner_params=runner_params
)
#             "runner_params":{header = False # depends on local or cluster. Cluster should have an option to choose cluster header
# extra_cmds = False
# extra_inputs = False # Only True for job arrays
# looper = True # True if looper indicat, false if not (in certain cases and in jobarrays)
# jobarray = False
# }
