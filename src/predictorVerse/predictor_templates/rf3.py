"""
Berta Bori Bru - IQAC-CSIC
Winter 2026
"""
from predictorVerse.utils import Predictor, RunnerParams

rf3_prot_chunk = """            {{
                "seq": "{seq}",
                "chain_id": "{letter}"
            }}
"""
rf3_lig_chunk = """            {{
                "smiles": "{seq}",
                "chain_id": "{letter}"
            }}
"""
rf3_rna_chunk = """            {{
                "seq": "{seq}",
                "chain_id": "{letter}"
            }}
"""
rf3_dna_chunk = """            {{
                "seq": "{seq}",
                "chain_id": "{letter}"
            }}
"""


rf3_all_template = """{{
    "name": "{system.name}",
    "components":[
{input}
    ]
}}"""


exec_command = "rf3 fold inputs=$inputs_dir ckpt_path=~/.foundry/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt diffusion_batch_size={samples} n_recycles={recycles} seed=42 out_dir=$outputs_dir "



runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=False,
                             jobarray=False)

rf3_data = Predictor(name = "RF3",
            prot_temp= rf3_prot_chunk,
            lig_temp= rf3_lig_chunk,
            rna_temp= rf3_rna_chunk,
            dna_temp= rf3_dna_chunk,
            joiner=",\n",
            prot_lig_temp= rf3_all_template,
            input_extension= ".json",
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)