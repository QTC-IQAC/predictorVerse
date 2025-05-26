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
import argparse
import os
from info import predictors_library
from jobscripts import gen_runner

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Generate input files for the given sequences in a .csv. The structure predictors to be used can be specified with -p')

    # Add the required argument for the CSV file
    parser.add_argument('input_csv', type=str, help='Path to the input .csv file')

    # Add the optional argument for a list of strings
    parser.add_argument('--predictors','-p', type=str, nargs='*' ,help='Optional list of predictors to use. Default all predictors ')
    parser.add_argument("--only_prot", action="store_true", help="Optional key to only generate inputs of the protein part.")
    
    # Parse the arguments
    args = parser.parse_args()

    # Read inputs
    system_list = read_input_csv(args.input_csv)

    if args.predictors:
        predictors_name_list = check_predictor_exists(args.predictors, predictors_library)
    else:
        predictors_name_list = predictors_library.keys()


    # Start doing things (TODO: put this in funcs)
    print("Generating inputs for the following predictors:")
    for predictor_name in predictors_name_list:
        print(f"------{predictor_name}------")
        
        # Get predictor data
        predictor = predictors_library[predictor_name]

        # Create directories
        predictor.create_folders()

        # Generate input (the default fasta file)
        for system in system_list:
            gen_input(system, predictor, only_prot=args.only_prot)

        # Generate runner
        gen_runner(system_list,predictor)


if __name__ == '__main__':
    main()