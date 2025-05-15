"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


You will execute this via command line and give
- csv path with system name and sequences and smiles
- the predictors to use
- A place to dump the outputs
"""


from utils import System, Predictor
from utils import read_input_csv, check_predictor_exists, gen_input
import sys
import os
from info import predictors_library
from jobscripts import gen_runner

# Input arguments
input_csv = sys.argv[1]
input_predictors_name = ["AF3",
              "RFAA",
              "Chai",
              "Boltz",
              "OF"
]

# Read inputs
system_list = read_input_csv(input_csv)

predictors_name_list = check_predictor_exists(input_predictors_name, predictors_library)

# Start doing things (TODO: put this in funcs)
for predictor_name in predictors_name_list:
    print(f"------{predictor_name}------")
    
    # Get predictor data
    predictor = predictors_library[predictor_name]


    # Create directories
    predictor.create_folders()

    # Generate input (the default fasta file)
    for system in system_list:
        gen_input(system, predictor, mode="prot")

    # # Generate runner
    gen_runner(system_list,predictor)
