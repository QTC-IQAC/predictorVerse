general_template ="""#!/bin/bash
{header}

cwd=$(pwd)

{inputs_section}

{extras}

start_date=$(date)
{execution_section}
final_date=$(date)

echo "START DATE: $start_date"
echo "FINAL DATE: $final_date"

"""

# flag for runner in jobarray
"""
# Check if each argument is provided
if [ -z "$1" ] ; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <input.json> "
    exit 1
fi
"""




# Basic inputs
"""
inputs_dir=$cwd/{workspace.inputs_predictor_unix}
outputs_dir=$cwd/{workspace.outputs_predictor_unix}
"""


# Addition to inputs in jobarray
"""
input_name=$1
input=$inputs_dir/$input_name
"""


# Looper
"""
for input_file in $inputs_dir/*{file_extension}; do
    {execution_command}
done
"""

header_clusteriqac="""#SBATCH -J boltz
#SBATCH -e %J.err
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
"""

header_csuc="""#SBATCH --job-name=AF3
#SBATCH -e AF3.err
#SBATCH -o AF3.log
#SBATCH -t 00-01:00
#SBATCH -p gpu
#SBATCH --reservation=cuda_12.6
#SBATCH --gres=gpu:1
#SBATCH -n 32
"""

# Line for jobarrays
"""
#SBATCH --array=1-10%2
"""

def gen_af3_jobarray(system_list:list[System], workspace:Workspace, max_cap_jobs=None | int) -> None:
    """
    
    max_cap_jobs : max number of jobs that can run at the same time.
                  If None, set to number of jobs
    """
    # Set variables
    num_jobs = len(system_list)
    num_cap_jobs = num_jobs if max_cap_jobs is None else max_cap_jobs
    af3_jobarr_header = af3_jobarray_template.format(num_jobs, num_cap_jobs)

    # Runner file relative path
    runner_file_rel_path = os.path.join(".",workspace.runners,"AF3_runner.sh")

    # Open jobarr file
    jobarr_file = os.path.join(workspace.runners,"AF3_jobarray.sh")
    with open(jobarr_file, "w") as jobarr:
        # Write header with correct parameters
        jobarr.write(af3_jobarr_header)

        # Write case start
        jobarr.write("case $SLURM_ARRAY_TASK_ID in\n")

        # Start looping through systems
        for ii, system in enumerate(system_list):
            jobarr.write(f"{ii+1}) {runner_file_rel_path} {system.name}.json ;;\n")
        
        jobarr.write("esac")