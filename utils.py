"""
Berta Bori Bru

Utility functions
"""
import os
from rdkit import Chem # main tools
from rdkit.Chem import AllChem # additional tools, including 3D



class System:
    def __init__(self, name, seq, smiles ):
        self.name = name #system name
        self.seq = seq #protein sequence
        self.smiles = smiles #smiles sequence

class Workspace:
    def __init__(self, name, pred=""):
        self.name = name
        self.inputs = self.name+"_inputs"
        self.outputs = self.name+"_outputs"
        self.runners = os.path.join(self.name+"_inputs","runners")

        os.makedirs(self.inputs,exist_ok=True)
        os.makedirs(self.runners,exist_ok=True)
        os.makedirs(self.outputs,exist_ok=True)
        
        self.predictor = pred

        # @property
        # def predictor(self):
        #     return self._predictor

    @property
    def inputs_predictor(self):
        return os.path.join(self.inputs,self.predictor)
    
    @property
    def inputs_predictor_unix(self):
        return self.inputs+"/"+self.predictor
    
    @property
    def outputs_predictor(self):
        return os.path.join(self.outputs,self.predictor)
    
    @property
    def outputs_predictor_unix(self):
        return self.outputs+"/"+self.predictor
        


def read_input_csv(csv_file: str, sep=",") -> list:
    """
    Read a .csv to extract info of the systems to write
    
    """
    
    system_list = []
    with open(csv_file,"r") as ff:
        header  = ff.readline() # separate header
        
        for sys in ff.readlines():
            name,seq,smiles = sys.split(sep)
            smiles = smiles.strip()
            system_list.append(System(name,seq,smiles))
    
    return system_list



def gen_fasta(system:System ,out_path:str, mode=None|str)-> None:
    """
    Generate fasta for system. Default structure is Chai-1 type
    out_path: path were to write the fastas
    mode = None (makes fasta with all the things in system). 
           protein (makes fasta of the protein only)
           ligand (makes fasta only of the ligand only)

    """
    if mode is None:                suffix = ""
    elif mode.lower() == "protein": suffix = "_prot"
    elif mode.lower() == "ligand":  suffix = "_lig"
    else: raise ValueError("Mode is not correct. Should be None, protein or ligand.")

    fasta_file = os.path.join(out_path,f"{system.name}{suffix}.fasta")
    
    fff = open(fasta_file,"w")
    
    if mode is None or mode.lower() == "protein" :
        fff.write(f">protein|{system.name}_prot\n{system.seq}\n")

    if mode is None or mode.lower() == "ligand":
        fff.write(f">ligand|{system.name}_lig\n{system.smiles}\n")

    fff.close()

def lig_smiles_to_sdf(system:System, out_path:str):
    """
    
    """
    sdf_file = os.path.join(out_path,f"{system.name}_lig.sdf")
    
    mol = Chem.MolFromSmiles(system.smiles) # initialize molecule

    mol = Chem.AddHs(mol) # adding explicit Hs for 3D generation
    cid = AllChem.EmbedMolecule(mol) # returns the id of the generated conformer,
                                    # and -1 if no conformers were generated

    AllChem.MMFFOptimizeMolecule(mol) # optimize molecule with MMFF94
    writer = Chem.SDWriter(sdf_file) # Write in .sdf file
    writer.write(mol)

def gen_input(system:System, workspace: Workspace, input_txt: str, input_extension: str) -> None:
    """
    input_text: text of the input. formated with system and workspace attributes
    input_extension: extension of the file. Ex: ".json", ".fasta"
    """

    input_file = os.path.join(workspace.inputs_predictor,system.name+input_extension)
    input_str = af3_json_template.format(system=system, workspace=workspace)
    
    with open(input_file, "w") as file:
        file.write(input_str)

