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

class RunnerParams:
    def __init__(self,
                 header=False,
                 extra_cmds=False,
                 extra_inputs=False,
                 looper=True,
                 jobarray=False):
        self._header = header
        self._extra_cmds = extra_cmds
        self._extra_inputs = extra_inputs
        self._looper = looper
        self._jobarray = jobarray

    @property
    def header(self):
        return self._header
    @header.setter
    def header(self, value):
        if isinstance(value, (bool,str)):
            self._header = value
        else:
            raise ValueError("header must be a boolean or a string value.")
   
    @property
    def extra_cmds(self):
        return self._extra_cmds
    @extra_cmds.setter
    def extra_cmds(self, value):
        if isinstance(value, bool):
            self._extra_cmds = value
        else:
            raise ValueError("extra_cmds must be a boolean value.")
    
    @property
    def extra_inputs(self):
        return self._extra_inputs
    @extra_inputs.setter
    def extra_inputs(self, value):
        if isinstance(value, bool):
            self._extra_inputs = value
        else:
            raise ValueError("extra_inputs must be a boolean value.")
    
    @property
    def looper(self):
        return self._looper
    @looper.setter
    def looper(self, value):
        if isinstance(value, bool):
            self._looper = value
        else:
            raise ValueError("looper must be a boolean value.")
    
    @property
    def jobarray(self):
        return self._jobarray
    @jobarray.setter
    def jobarray(self, value):
        if isinstance(value, bool):
            self._jobarray = value
        else:
            raise ValueError("jobarray must be a boolean value.")



class Predictor:
    def __init__(self, name,
                 prot_temp,
                 lig_temp,
                 prot_lig_temp,
                 input_extension,
                 main_cmds,
                 joiner="\n",
                 other_funcs=None,
                 extra_cmds=None,
                 runner_params=RunnerParams()):
        
        self.name = name
        self.inputs = os.path.join(self.name, "inputs")
        self.outputs = os.path.join(self.name, "outputs")
        self.runners = os.path.join(self.name, "runners")
        self.inputs_unix = self.name+"/"+"inputs"
        self.outputs_unix = self.name+"/"+"outputs"

        self.prot_temp = prot_temp
        self.lig_temp = lig_temp
        self.prot_lig_temp = prot_lig_temp
        self.input_extension = input_extension
        self.main_cmds = main_cmds
        self.joiner = joiner
        self.other_funcs = other_funcs
        self.extra_cmds = extra_cmds
        self.runner_params = runner_params

    @property
    def prot_temp(self):
        return self._prot_temp

    @prot_temp.setter
    def prot_temp(self, value):
        if isinstance(value, str):
            self._prot_temp = value
        else:
            raise ValueError("prot_temp must be a string.")

    @property
    def lig_temp(self):
        return self._lig_temp

    @lig_temp.setter
    def lig_temp(self, value):
        if isinstance(value, str):
            self._lig_temp = value
        else:
            raise ValueError("lig_temp must be a string.")

    @property
    def prot_lig_temp(self):
        return self._prot_lig_temp

    @prot_lig_temp.setter
    def prot_lig_temp(self, value):
        if isinstance(value, str):
            self._prot_lig_temp = value
        else:
            raise ValueError("prot_lig_temp must be a string.")

    @property
    def input_extension(self):
        return self._input_extension

    @input_extension.setter
    def input_extension(self, value):
        if isinstance(value, str) or value is None:
            self._input_extension = value
        else:
            raise ValueError("input_extension must be a string or None.")

    @property
    def main_cmds(self):
        return self._main_cmds

    @main_cmds.setter
    def main_cmds(self, value):
        if isinstance(value, str):
            self._main_cmds = value
        else:
            raise ValueError("main_cmds must be a string.")

    @property
    def joiner(self):
        return self._joiner

    @joiner.setter
    def joiner(self, value):
        if isinstance(value, str):
            self._joiner = value
        else:
            raise ValueError("joiner must be a string.")

    @property
    def other_funcs(self):
        return self._other_funcs

    @other_funcs.setter
    def other_funcs(self, value):
        if isinstance(value, list) or value is None:
            self._other_funcs = value
        else:
            raise ValueError("other_funcs must be a list of functions or None.")

    @property
    def extra_cmds(self):
        return self._extra_cmds

    @extra_cmds.setter
    def extra_cmds(self, value):
        if isinstance(value, str) or value is None:
            self._extra_cmds = value
        else:
            raise ValueError("extra_cmds must be a string or None.")

    @property
    def runner_params(self):
        return self._runner_params

    @runner_params.setter
    def runner_params(self, value):
        if isinstance(value, RunnerParams):
            self._runner_params = value
        else:
            raise ValueError("runner_params must be an instance of RunnerParams.")

    def create_folders(self):
        os.makedirs(self.inputs, exist_ok=True)
        os.makedirs(self.runners, exist_ok=True)
        os.makedirs(self.outputs, exist_ok=True)
        

##### Functions #####

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



def gen_input(system:System, predictor: Predictor, only_prot=False) -> None:
    """
    input_template: text of the input. formated with system and predictor attributes
    input_extension: extension of the file. Ex: ".json", ".fasta"
    """
    input_extension = predictor.input_extension
    general_template = predictor.prot_lig_temp
    prot_template = predictor.prot_temp
    lig_template = predictor.lig_temp
    joiner = predictor.joiner

    input_file = os.path.join(predictor.inputs,system.name+input_extension)

    if not only_prot:
        prot_text = prot_template.format(system=system,predictor=predictor)
        lig_text = lig_template.format(system=system,predictor=predictor)
        inputs_to_join = [prot_text,lig_text]
    
    elif only_prot:
        prot_text = prot_template.format(system=system,predictor=predictor)
        inputs_to_join = [prot_text]
    
    input_text = joiner.join(inputs_to_join)
    final_input = general_template.format(input=input_text, system=system, predictor=predictor)
    
    with open(input_file, "w") as file:
        file.write(final_input)

    # If there are other functions (in predictor_data["other_funcs"]), execute them
    try:
        if predictor.other_funcs is not None:
            for func in predictor.other_funcs:
                func(system, predictor)
    except:
        print(f"WARNING: No other funcs were found for {predictor.name} or is not iterable. Skipping")



def check_predictor_exists(predictor_id_list: list[str], predictors_library: dict) -> list[str]: 
    
    pred_to_remove = []
    
    for predictor in predictor_id_list:
        try:
            aa = predictors_library[predictor]
        except:
            print(f"WARNING: {predictor} id has no data associated. Removing it from predictors list.")
            pred_to_remove.append(predictor)
    
    checked_predictors = [x for x in predictor_id_list if x not in pred_to_remove]

    # Check if there is some predictor left
    if len(checked_predictors) == 0: 
        raise ValueError("There are no predictors. Cannot continue.")
    
    return checked_predictors
