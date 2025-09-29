"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from predictorVerse.utils import Predictor, RunnerParams

chai_prot_fasta = """>protein|{system.name}_prot_{letter}\n{seq}"""
chai_lig_fasta = """>ligand|{system.name}_lig_{letter}\n{seq}"""
chai_fasta_template = """{input}"""

exec_command = """# Get name of system

system=$(basename $file | sed "s/.fasta//g")
mkdir -p $outputs_dir/$system
chai-lab fold --seed 13 --num-diffn-samples 10 --num-trunk-recycles 3 $file $outputs_dir/$system
"""

runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)

chai_data = Predictor(name = "Chai",
            prot_temp= chai_prot_fasta,
            lig_temp= chai_lig_fasta,
            joiner="\n",
            prot_lig_temp= chai_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)


