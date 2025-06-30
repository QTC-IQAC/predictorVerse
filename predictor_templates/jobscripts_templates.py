"""
Berta Bori Bru - IQAC-CSIC
Spring 2025

Text templates for creation of runner scripts
"""


# General template of a runner file
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


# Line for jobarrays
jobarr_temp = """#SBATCH --array=1-{0}%{1}

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


### HEADERS OF OPTIONS FOR CLUSTERS ###
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
#SBATCH --gres=gpu:1
#SBATCH -n 32
"""

headers_dict = {"clusteriqac": header_clusteriqac,
                "csuc": header_csuc,
                "default": header_clusteriqac}




###### Params ######
"""
header = False | True # depends on local or cluster. Cluster should have an option to choose cluster header
extra_cmds = False | True # If there are extra commands
extra_inputs = False | True # Only True for job arrays
looper = False | True # True if looper indicat, false if not (in certain cases and in jobarrays)
jobarray = False | True # True only if jobarray
"""