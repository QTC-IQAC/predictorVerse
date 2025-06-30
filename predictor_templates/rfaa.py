import os
from utils import System, Predictor, RunnerParams

from rdkit import Chem # main tools
from rdkit.Chem import AllChem # additional tools, including 3D



rfaa_prot_template = """protein_inputs:
  {letter}:
    fasta_file: ./{predictor.inputs_unix}/fastas/{system.name}_prot{ii}.fasta"""

rfaa_lig_template = """sm_inputs:
  {letter}:
    input: ./{predictor.inputs_unix}/fastas/{system.name}_lig{ii}.sdf
    input_type: 'sdf'"""

rfaa_yaml_template ="""defaults:
- base
- _self_
job_name: {system.name}
output_path: ./{predictor.outputs_unix}
{input}
"""

extra_cmds = "cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed"

exec_command = """# Change the relative path for absolute ones
sed -i "s|\\./|$cwd/|g" $file


# Split input file into file name and folder
filename=$(basename "$file")
foldername=$(dirname "$file")

python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername

"""

def gen_prot_fasta(system: System, predictor: Predictor):
    """
    Generate fasta file for the protein
    """
    try: 
       system.protein
    except:
       return
    
    # Create fastas folder
    fastas_dir = os.path.join(predictor.inputs,"fastas")
    os.makedirs(fastas_dir,exist_ok=True)

    # Generate protein fasta files
    for ii,seq in enumerate(system.protein):
      fasta_file = os.path.join(predictor.inputs,"fastas",f"{system.name}_prot{ii}.fasta")
      
      with open(fasta_file,"w") as fastaff:
          fastaff.write(f">protein|{system.name}_prot\n{seq}\n")


def lig_smiles_to_sdf(system:System, predictor: Predictor):
    """
    From ligand smiles generate a sdf file 
    """
    try: 
       system.ligand
    except:
       return
    
    for ii,lig_smile in enumerate(system.ligand):
      sdf_file = os.path.join(predictor.inputs,"fastas",f"{system.name}_lig{ii}.sdf")
      
      mol = Chem.MolFromSmiles(lig_smile) # initialize molecule
      mol = Chem.AddHs(mol) # adding explicit Hs for 3D generation
      # cid = AllChem.EmbedMolecule(mol) # returns the id of the generated conformer,
      #                                 # and -1 if no conformers were generated

      AllChem.MMFFOptimizeMolecule(mol) # optimize molecule with MMFF94
      writer = Chem.SDWriter(sdf_file) # Write in .sdf file
      writer.write(mol)


runner_params = RunnerParams(header="clusteriqac",
                             extra_cmds=True,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)

    
rfaa_data = Predictor(name= "RFAA",
            prot_temp= rfaa_prot_template,
            lig_temp= rfaa_lig_template,
            joiner="\n",
            prot_lig_temp= rfaa_yaml_template,
            input_extension= ".yaml",
            other_funcs=[gen_prot_fasta, lig_smiles_to_sdf],
            extra_cmds=extra_cmds,
            main_cmds=exec_command,
            runner_params=runner_params
)
