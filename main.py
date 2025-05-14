"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


You will execute this via command line and give
- csv path with system name and sequences and smiles
- the predictors to use
- A place to dump the outputs
"""


from utils import System, Workspace
from utils import read_input_csv, check_predictor_exists, gen_input
import sys
import os
from info import predictors_library

# Input arguments
input_csv = sys.argv[1]
workspace_name = "GalaxyTEST" # TODO This will be a defalt with argparse
input_predictors = ["AF3"]
# ["AF3",
#               "RFAA",
#               "Chai",
#               "Boltz",
#               "OF"
# ]

# Read inputs
system_list = read_input_csv(input_csv)
workspace = Workspace(workspace_name)
predictors_list = check_predictor_exists(input_predictors, predictors_library)

# Start doing things (TODO: put this in funcs)
for predictor in predictors_list:
    print(f"------{predictor}------")
    
    # Get predictor data
    data = predictors_library[predictor]

    # Change current predictor in in Workspace
    workspace.predictor = predictor

    # Create directories
    os.makedirs(workspace.inputs_predictor,exist_ok=True)
    os.makedirs(workspace.outputs_predictor,exist_ok=True)

    # Generate input (the default fasta file)
    for system in system_list:
        gen_input(system, workspace, data,mode="prot")

    # # Generate runner
    # gen_runner(workspace)
