"""
Berta Bori Bru

Utility functions
"""
import os


class System:
    def __init__(self, name, seq, smiles ):
        self.name = name #system name
        self.seq = seq #protein sequence
        self.smiles = smiles #smiles sequence

def read_fasta(fasta_file: str) -> tuple:
    pass

def read_input_csv(csv_file: str, sep=",") -> list:
    """
    Read a .csv to extract info of the systems to write
    
    """
    
    system_list = []
    with open(csv_file,"r") as ff:
        header  = ff.readline() # separate header
        
        for sys in ff.readlines():
            name,seq,smiles = sys.split(sep)
            system_list.append(System(name,seq,smiles))
    
    return system_list

def make_folder_structure(workspace_name="Galaxy42", workspace_path=".") -> None : #TODO: migrate defaults to argparse
    """
    Make folders to store inputs, outputs and runners for the predictors
    """
    inputs_dir = os.path.join(workspace_path, workspace_name+"_inputs")
    runner_dir = os.path.join(workspace_path, workspace_name+"_inputs","runners")
    outputs_dir = os.path.join(workspace_path, workspace_name+"_outputs")

    os.makedirs(inputs_dir,exist_ok=True)
    os.makedirs(runner_dir,exist_ok=True)
    os.makedirs(outputs_dir,exist_ok=True)

    #TODO: would be cool if it created subfolders only for the predictors you use
            
    return workspace_name


# def checkSeqType(seq:str) -> str:
#     """
#     Checks is the sequence corresponds to RNA, DNA, protein or ligand (SMILES) type
#     This program is only able to identify basic types of sequences and large SMILES.
#     Check afterwards if the program has identified the type correctly.

#     Parameters:
#         seq: string with the sequence to identify

#     Returns
#         mode: identified type of the sequence
#     """
#     # Convert sequence to list
#     seq_list = list(seq)

#     # Definition of characteristic characters of each sequence
#     smiles = list("=#:+-[]()/\@.%1234567890")
#     rna=list("AUGC")
#     dna=list("ATGC")

#     # Check if the sequence is RNA, DNA or ligand (SMILES)
#     if all(char1 in rna for char1 in seq_list):
#         mode="rna"
#     elif all(char1 in dna for char1 in seq_list):
#         mode="dna"
#     elif any(char1 in smiles for char1 in seq_list) or len(seq_list)<50:
#         mode="ligand"
#     else:
#         mode="protein"

#     return mode


def gen_fasta(system:System ,write_path:str, mode=None|str)-> None:
    """
    Generate fasta for system. Default structure is Chai-1 type
    write_path: path were to write the fastas
    mode = None (makes fasta with all the things in system). 
           protein (makes fasta of the protein only)
           ligand (makes fasta only of the ligand only)

    """
    if mode is None:                suffix = ""
    elif mode.lower() == "protein": suffix = "_prot"
    elif mode.lower() == "ligand":  suffix = "_lig"
    else: raise ValueError("Mode is not correct. Should be None, protein or ligand.")



    fasta_file = os.path.join(write_path,f"{system.name}{suffix}.fasta")
    
    fff = open(fasta_file,"w")
    
    if mode is None or mode.lower() == "protein" :
        fff.write(f">protein|{system.name}_prot\n{system.seq}\n")

    if mode is None or mode.lower() == "ligand":
        fff.write(f">ligand|{system.name}_lig\n{system.smiles}\n")

    fff.close()


    
