import os 
from utils import System, Predictor


general_temp ="""#!/bin/bash
{header}

cwd=$(pwd)

{inputs_section}

{extras}

#### Main program ####

start_date=$(date)

{execution_section}

final_date=$(date)

cd $cwd

echo "START DATE: $start_date"
echo "FINAL DATE: $final_date"

"""

# flag for runner in jobarray
extra_inputs_txt = """
input_name=$1
input=$inputs_dir/$input_name

# Check if each argument is provided
if [ -z "$1" ] ; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <input.json> "
    exit 1
fi
"""

# Basic inputs
basic_input_temp = """
inputs_dir=$cwd/{predictor.inputs_unix}
outputs_dir=$cwd/{predictor.outputs_unix}
"""


# Looper
looper_temp = """
for file in $inputs_dir/*{file_extension}; do
    {execution_command}
done
"""

header_clusteriqac="""#SBATCH --job-name={predictor.name}
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
"""

header_csuc="""#SBATCH --job-name={predictor.name}
#SBATCH -e %j.err
#SBATCH -o %j.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --reservation=cuda_12.6
#SBATCH --gres=gpu:1
#SBATCH -n 32
"""

headers_dict = {"clusteriqac": header_clusteriqac,
                "csuc": header_csuc,
                "default": header_clusteriqac}

# Line for jobarrays
jobarr_temp = """
#SBATCH --array=1-{0}%{1}
"""

def get_header(header_key: bool | str) -> str:
    if isinstance(header_key, str):
        try:
            return headers_dict[header_key]
        except:
            print(f"{header_key} header is not found. Returning default one")
            return headers_dict["default"]
    elif isinstance(header_key, bool) and header_key is True:
        print("No header selected. Returning default one")
        return headers_dict["default"]



def gen_jobarray(system_list:list[System], predictor:Predictor, max_cap_jobs=None | int) -> None:
    """
    
    max_cap_jobs : max number of jobs that can run at the same time.
                  If None, set to number of jobs
    """
    # Set variables
    num_jobs = len(system_list)
    num_cap_jobs = num_jobs if max_cap_jobs is None else max_cap_jobs
    header_txt = get_header(predictor.runner_params.header)
    jobarr_header = header_txt.format(predictor=predictor) + jobarr_temp.format(num_jobs, num_cap_jobs)

    # Runner file relative path
    runner_file_rel_path = os.path.join(".",predictor.runners,f"{predictor.name}_runner.sh")

    # Open jobarr file
    jobarr_file = os.path.join(predictor.runners,f"{predictor.name}_jobarray.sh")
    with open(jobarr_file, "w") as jobarr:
        # Write header with correct parameters
        jobarr.write(jobarr_header)

        # Write case start
        jobarr.write("case $SLURM_ARRAY_TASK_ID in\n")

        # Start looping through systems
        for ii, system in enumerate(system_list):
            jobarr.write(f"{ii+1}) {runner_file_rel_path} {system.name}{predictor.input_extension} ;;\n")
        
        jobarr.write("esac")


def fix_commands_for_loop(txt:str) -> str:
    """
    Adds \t to the heading of each line, because I like indentations
    Probably in operating systems outside of Linux does not work
    """
    cmds_list = txt.split("\n")
    new_txt = "\n\t".join(cmds_list)
    return new_txt


###### Params ######
"""
header = False | True # depends on local or cluster. Cluster should have an option to choose cluster header
extra_cmds = False | True # If there are extra commands
extra_inputs = False | True # Only True for job arrays
looper = False | True # True if looper indicat, false if not (in certain cases and in jobarrays)
jobarray = False | True # True only if jobarray
"""


def gen_runner(system_list:list[System], predictor: Predictor):
    # Getting header. TODO: when we have params for this, we will have a function here
    if not predictor.runner_params.header:
        header_txt = ""
    else:
        header_txt = get_header(predictor.runner_params.header).format(predictor=predictor) # check which header to use

    
    # Processing inputs
    basic_input = basic_input_temp.format(predictor=predictor)
    
    if predictor.runner_params.extra_inputs:
        inputs_txt = "\n".join([basic_input_temp, extra_inputs_txt])
    else:
        inputs_txt = basic_input
    
    # Processing extra cmds
    if predictor.extra_cmds is not None or predictor.extra_cmds != "":
        extra_cmds_txt = predictor.extra_cmds
    else:
        extra_cmds_txt = ""

    # Processing main commands 
    if predictor.runner_params.looper:
        new_cmds = fix_commands_for_loop(predictor.main_cmds)
        main_cmds_txt = looper_temp.format(file_extension=predictor.input_extension,
                                           execution_command=new_cmds)
    else:
        main_cmds_txt = predictor.main_cmds
    
    # Assemble everything
    runner_txt = general_temp.format(header=header_txt,
                                    inputs_section=inputs_txt,
                                    extras=extra_cmds_txt,
                                    execution_section=main_cmds_txt)
    runner_file = os.path.join(predictor.runners,f"{predictor.name}_runner.sh")

    
    with open(runner_file,"w") as rr:
        rr.write(runner_txt)
    
    # Generate jobarray if needed
    if predictor.runner_params.jobarray:
        gen_jobarray(system_list,predictor, max_cap_jobs=2)