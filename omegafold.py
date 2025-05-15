import os
from utils import System, Predictor, RunnerParams

"""
TODO: OF only produces 1 output. Make it so:
    - We produce fasta for each system
    - runner has loop like
    for input file
        <generate subfolder with system name>
        for ii in {1..10}
            <run omegafold>
"""
of_prot_fasta = """>{system.name}_prot\n{system.seq}"""
of_fasta_template = """{input}"""

exec_command = "omegafold $file $output_dir"


runner_params = RunnerParams(header=True,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)


of_data = Predictor(name= "OF",
            prot_temp= of_fasta_template,
            lig_temp= "",
            joiner="\n",
            prot_lig_temp= of_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)


# def gen_of_input(system_list:list[System], predictor:Predictor)-> None:
#     """
#     Generate a fasta for all the proteins in system.
#     OF is only a predictor for proteins. 
#     """

#     fasta_file = os.path.join(predictor.inputs_predictor,"OF_input.fasta")

#     with open(fasta_file, "w") as oo:
#         for system in system_list:
#             oo.write(f">{system.name}_prot\n")
#             oo.write(f"{system.seq}\n")

# of_runner_temp = """#!/bin/bash
# #SBATCH -J OF
# #SBATCH -e %J.err
# #SBATCH -o %J.out
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --partition=long
# #SBATCH --gres=gpu:1
# ##SBATCH --constraint=gpu8

# input={0}/OF_input.fasta
# output_dir={1}

# omegafold $input $output_dir

# """
# #.format(inputs_dir,outputs_dir)


# def gen_of_runner(predictor:Predictor) -> None:
    
#     runner_file = os.path.join(predictor.runners,"OF_runner.sh")


#     runner_str = of_runner_temp.format(predictor.inputs_predictor_unix,
#                                        predictor.outputs_predictor_unix)
    
#     with open(runner_file,"w") as rr:
#         rr.write(runner_str)

# def main(system_list:System, predictor:Predictor):
#     # Change current predictor in in Predictor
#     predictor.predictor = "OmegaFold"

#     # Create directories
#     os.makedirs(predictor.inputs_predictor,exist_ok=True)
#     os.makedirs(predictor.outputs_predictor,exist_ok=True)

#     # Generate input
#     gen_of_input(system_list, predictor)

#     # Generate runner
#     gen_of_runner(predictor)

of_data = {"name": "OF",
            "prot_temp": of_fasta_template,
            "lig_temp": "",
            "joiner":"\n",
            "prot_lig_temp": of_fasta_template,
            "input_extension": ".fasta",
            # "runner_temp":
}