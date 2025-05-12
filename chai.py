"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
import os
from utils import System, Workspace, gen_fasta

chai_prot_fasta = """>protein|{system.name}_prot
{system.seq}
"""

chai_lig_fasta = """>ligand|{system.name}_lig
{system.smiles}
"""

chai_fasta_template = chai_prot_fasta + chai_lig_fasta

chai_runner_template = """ #!/bin/bash

inputs_folder={0}
outputs_folder={1}

for file in $inputs_folder/*.fasta; 
    do
    # Get name of system
    system=$(basename $file | sed "s/.fasta//g")
    mkdir -p $outputs_folder/$system
    chai-lab fold --seed 42 --num-diffn-samples 10 $file $outputs_folder/$system
done

"""
# .format(workspace.inputs_predictor_unix, workspace.outputs_predictor_unix )

def gen_chai_runner(workspace:Workspace) -> None:
    runner_file = os.path.join(workspace.runners,"Chai_runner.sh")

    runner_str = chai_runner_template.format(workspace.inputs_predictor_unix, workspace.outputs_predictor_unix)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)


def main(system_list:System, workspace:Workspace)-> None:
    # Change current predictor in in Workspace
    workspace.predictor = "Chai"

    # Create directories
    os.makedirs(workspace.inputs_predictor,exist_ok=True)
    os.makedirs(workspace.outputs_predictor,exist_ok=True)

    # Generate input (the default fasta file)
    for system in system_list:
        gen_fasta(system, workspace.inputs_predictor,mode=None)

    # Generate runner
    gen_chai_runner(workspace)

chai_data = {"name": "Chai-1",
            "prot_temp": chai_prot_fasta
            "lig_temp": chai_lig_fasta
            "prot_lig_temp": chai_fasta_template,
            "input_extension": ".fasta",
            # "runner_temp":
}