"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from utils import System, Predictor, RunnerParams
import os
from boltz1 import boltz_prot_fasta, boltz_lig_fasta, boltz_fasta_template

# boltz_prot_fasta = """>A|protein||{system.name}_prot\n{system.seq}"""
# boltz_lig_fasta = """>B|smiles||{system.name}_lig\n{system.smiles}"""
# boltz_fasta_template = """{input}"""

exec_command = "boltz predict $inputs_dir --use_msa_server --diffusion_samples 10 --out_dir $outputs_dir --output_format pdb"

runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=False,
                             jobarray=False)


boltz1x_data = Predictor(name= "Boltz1x",
            prot_temp= boltz_prot_fasta,
            lig_temp= boltz_lig_fasta,
            joiner="\n",
            prot_lig_temp= boltz_fasta_template,
            input_extension= ".fasta",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)