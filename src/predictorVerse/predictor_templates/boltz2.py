"""
Berta Bori Bru - IQAC-CSIC
Winter 2026
"""
from predictorVerse.utils import Predictor, RunnerParams


boltz_prot_chunk = """\
    - protein:
        id: {letter}
        sequence: {seq}
"""
boltz_lig_chunk = """\
    - ligand:
        id: {letter}
        smiles: '{seq}'
"""
boltz_rna_chunk = """\
    - rna:
        id: {letter}
        sequence: {seq}
"""
boltz_dna_chunk = """\
    - dna:
        id: {letter}
        sequence: {seq}
"""


boltz_all_template = """\
sequences:
{input}

"""


exec_command = "boltz predict $inputs_dir --use_msa_server --use_potentials --diffusion_samples {samples} --recycling_steps {recycles} --out_dir $outputs_dir --output_format pdb"



runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=False,
                             jobarray=False)

boltz2_data = Predictor(name = "Boltz2",
            prot_temp= boltz_prot_chunk,
            lig_temp= boltz_lig_chunk,
            rna_temp= boltz_rna_chunk,
            dna_temp= boltz_dna_chunk,
            joiner="\n",
            prot_lig_temp= boltz_all_template,
            input_extension= ".yaml",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)