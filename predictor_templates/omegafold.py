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
of_prot_fasta = """>{system.name}_prot\n{seq}"""
of_fasta_template = """{input}"""

exec_command = "omegafold $file $output_dir"


runner_params = RunnerParams(header="clusteriqac",
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)


of_data = Predictor(name= "OF",
            prot_temp= of_prot_fasta,
            lig_temp= "",
            joiner="\n",
            prot_lig_temp= of_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)


