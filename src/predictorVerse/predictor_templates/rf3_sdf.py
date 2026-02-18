"""
Berta Bori Bru - IQAC-CSIC
Winter 2026
"""
from predictorVerse.utils import System, Predictor, RunnerParams
from rdkit import Chem # main tools
from rdkit.Chem import AllChem # additional tools, including 3D
import os

rf3_prot_chunk = """            {{
                "seq": "{seq}",
                "chain_id": "{letter}"
            }}
"""
rf3_lig_chunk = """            {{
                "path": "./{predictor.inputs_unix}/sdfs/{system.name}_lig{ii}.sdf",
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


def lig_smiles_to_sdf(system:System, predictor: Predictor, *args):
    """
    From ligand smiles generate a sdf file 
    """
    try: 
       system.ligand
    except:
       return
    
    sdf_dir = os.path.join(predictor.inputs,"sdfs")
    os.makedirs(sdf_dir,exist_ok=True)

    
    for ii,lig_smile in enumerate(system.ligand):
      sdf_file = os.path.join(sdf_dir,f"{system.name}_lig{ii}.sdf")
      
      mol = Chem.MolFromSmiles(lig_smile) # initialize molecule
      mol = Chem.AddHs(mol) # adding explicit Hs for 3D generation
      # cid = AllChem.EmbedMolecule(mol) # returns the id of the generated conformer,
      #                                 # and -1 if no conformers were generated

      # AllChem.MMFFOptimizeMolecule(mol) # optimize molecule with MMFF94
      AllChem.EmbedMolecule(mol, AllChem.ETKDG())
      AllChem.UFFOptimizeMolecule(mol)
      
      writer = Chem.SDWriter(sdf_file) # Write in .sdf file
      writer.write(mol)


exec_command = "rf3 fold inputs=$inputs_dir ckpt_path=~/.foundry/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt diffusion_batch_size={samples} n_recycles={recycles} seed=42 out_dir=$outputs_dir "



runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=False,
                             jobarray=False)

rf3_sdf_data = Predictor(name = "RF3_sdf",
            prot_temp= rf3_prot_chunk,
            lig_temp= rf3_lig_chunk,
            rna_temp= rf3_rna_chunk,
            dna_temp= rf3_dna_chunk,
            joiner=",\n",
            prot_lig_temp= rf3_all_template,
            input_extension= ".json",
            other_funcs=[lig_smiles_to_sdf],
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)