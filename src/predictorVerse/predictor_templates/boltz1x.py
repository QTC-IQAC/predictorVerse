"""
Berta Bori Bru - IQAC-CSIC
Spring 2025
"""
from predictorVerse.utils import Predictor, RunnerParams
from predictorVerse.predictor_templates.boltz1 import boltz_prot_fasta, boltz_lig_fasta, boltz_fasta_template

exec_command = "boltz predict $inputs_dir --use_msa_server --diffusion_samples {samples} --recycling_steps {recycles} --out_dir $outputs_dir --output_format pdb"

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