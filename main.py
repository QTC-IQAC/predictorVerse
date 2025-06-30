"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


You will execute this via command line and give
- json path with system name and sequences and smiles
- the predictors to use
"""


from utils import System, Predictor
from utils import read_input_json, check_predictor_exists, gen_input
import argparse
import os
from info import predictors_library
from jobscripts import gen_runner


def main(args):

    # Read inputs
    system_list = read_input_json(args.input_json)

    if args.predictors:
        predictors_name_list = check_predictor_exists(args.predictors, predictors_library)
    else:
        predictors_name_list = predictors_library.keys()


    # Start doing things
    print("Generating inputs for the following predictors:")
    for predictor_name in predictors_name_list:
        print(f"------{predictor_name}------")
        
        # Get predictor data
        predictor = predictors_library[predictor_name]

        # Create directories
        predictor.create_folders()

        # Generate input (the default fasta file)
        for system in system_list:
            gen_input(system, predictor)

        # Generate runner
        gen_runner(system_list,predictor)


if __name__ == '__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description='Generate input files for the given sequences in a .json. The structure predictors to be used can be specified with -p')
    # Add the required argument for the JSON file
    parser.add_argument('input_json', type=str, help='Path to the input .json file')

    # Add the optional argument for a list of strings
    parser.add_argument('--predictors','-p', type=str, nargs='*' ,
                        help=f'Optional list of predictors to use. Default all predictors. Predictors available: {list(predictors_library.keys())}')

    # Parse the arguments
    args = parser.parse_args()
    
    main(args)