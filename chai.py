"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
import os
from utils import System, Workpath, gen_fasta


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
# .format(workpath.inputs_predictor_unix, workpath.outputs_predictor_unix )

def gen_chai_runner(workpath:Workpath) -> None:
    runner_file = os.path.join(workpath.runners,"Chai_runner.sh")

    runner_str = chai_runner_template.format(workpath.inputs_predictor_unix, workpath.outputs_predictor_unix)
    
    with open(runner_file,"w") as rr:
        rr.write(runner_str)


def main(system_list:System, workpath:Workpath)-> None:
    # Change current predictor in in Workpath
    workpath.predictor = "Chai"

    # Create directories
    os.makedirs(workpath.inputs_predictor,exist_ok=True)
    os.makedirs(workpath.outputs_predictor,exist_ok=True)

    # Generate input (the default fasta file)
    for system in system_list:
        gen_fasta(system, workpath.inputs_predictor,mode=None)

    # Generate runner
    gen_chai_runner(workpath)
