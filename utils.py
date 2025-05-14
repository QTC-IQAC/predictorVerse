"""
Berta Bori Bru

Utility functions
"""
import os
# from info import predictors_library



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



def gen_input(system:System, workspace: Workspace, predictor_data: dict, mode="all") -> None:
    """
    input_template: text of the input. formated with system and workspace attributes
    input_extension: extension of the file. Ex: ".json", ".fasta"
    """
    input_extension = predictor_data["input_extension"]
    general_template = predictor_data["prot_lig_temp"]
    prot_template = predictor_data["prot_temp"]
    lig_template = predictor_data["lig_temp"]
    joiner = predictor_data["joiner"]

    input_file = os.path.join(workspace.inputs_predictor,system.name+input_extension)

    if mode == "all":
        prot_text = prot_template.format(system=system,workspace=workspace)
        lig_text = lig_template.format(system=system,workspace=workspace)
        inputs_to_join = [prot_text,lig_text]
    elif mode == "prot":
        prot_text = prot_template.format(system=system,workspace=workspace)
        inputs_to_join = [prot_text]
    
    input_text = joiner.join(inputs_to_join)
    final_input = general_template.format(input=input_text, system=system, workspace=workspace)
    
    with open(input_file, "w") as file:
        file.write(final_input)

    # If there are other functions (in predictor_data["other_funcs"]), execute them
    try:
        if predictor_data["other_funcs"] is not None:
            for func in predictor_data["other_funcs"]:
                func(system, workspace)
    except:
        print(f"WARNING: No other funcs were found for {predictor_data["name"]} or is not iterable. Skipping")



def check_predictor_exists(predictor_id_list: list[str], predictors_library: dict) -> list[str]: 
    
    pred_to_remove = []
    
    for predictor in predictor_id_list:
        try:
            aa = predictors_library[predictor]
        except:
            print(f"WARNING: {predictor} id has no data associated. Removing it from predictors list.")
            pred_to_remove.appen(predictor)
    
    checked_predictors = [x for x in predictor_id_list if x not in pred_to_remove]

    # Check if there is some predictor left
    if len(checked_predictors) == 0: 
        raise ValueError("There are no predictors. Cannot continue.")
    
    return checked_predictors
