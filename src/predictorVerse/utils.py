"""
Berta Bori Bru - IQAC-CSIC
Spring 2025

Utility functions
"""
import os
import json


class System:
    def __init__(self, name: str, data_dict: dict):
        self.name = name #system name
        self.inner_dct = dict()
        
        for kk, vv in data_dict.items():
            if kk not in ["protein","ligand"]: # Handle other types of names
                raise ValueError(f"Incorrect data field for '{kk}': can only be protein or ligand")
            
            if type(vv) is not list and type(vv) is not str: # Handle other types of data 
                raise TypeError(f"Incorrect data type for '{vv}' in '{self.name}, {kk}': can only be str or list")
            
            # Handle the passing of list and str differently
            if type(vv) is list:
                self.inner_dct[kk] = vv
                setattr(self, kk, vv)
            
            elif type(vv) is str:
                self.inner_dct[kk] = [vv]
                setattr(self, kk, [vv])
    


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

def read_input_json(json_file: str) -> list:
    """
    Read a .json to extract info of the systems to write
    
    """
    input_info = json.load(open(json_file))
    system_list = [System(kk,vv) for kk,vv in input_info.items()]
    print(f"Found {len(system_list)} systems")
    
    return system_list


def alphabet_generator():
    """
    An alphabet generator to put correct chain letters in AF3, RFAA and the like
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for letter in alphabet:
        yield letter


def combinations(field:str, system:System, predictor:Predictor): #predictor tmb aniria aqui
    """ 
    Returns the combinations of system sequences in a determinate field and its corresponding predictor input template
    For now only handles protein and ligand fields
    """

    if field.lower() == "protein":
        return system.protein, predictor.prot_temp

    
    elif field.lower() == "ligand":
        return system.ligand, predictor.lig_temp


def gen_subinputs(system:System, predictor:Predictor): # Aquí també va predictor
    """ 
    Generates the input text that contain information about the sequences
    """
    alphabet = alphabet_generator()

    list_subinputs = []
    for field in ["protein","ligand"]:
        try:
            data,txt = combinations(field, system, predictor)
        except:
            continue
        list_subinputs.append( predictor.joiner.join(
            [txt.format(system=system, predictor=predictor, seq=dd, ii=ii,letter=next(alphabet)) 
            for ii,dd in enumerate(data)])   
            )
    
    return list_subinputs



def gen_input(system:System, predictor: Predictor) -> None:
    """
    input_template: text of the input. formated with system and predictor attributes
    input_extension: extension of the file. Ex: ".json", ".fasta"
    """
    input_extension = predictor.input_extension
    general_template = predictor.prot_lig_temp
    joiner = predictor.joiner

    input_file = os.path.join(predictor.inputs,system.name+input_extension)

    subinputs_to_join = gen_subinputs(system,predictor)
    
    input_text = joiner.join(subinputs_to_join)
    final_input = general_template.format(input=input_text, system=system, predictor=predictor)
    
    with open(input_file, "w") as file:
        file.write(final_input)

    # If there are other functions (in predictor_data["other_funcs"]), execute them
    if predictor.other_funcs is not None:
        try:
            for func in predictor.other_funcs:
                func(system, predictor)
        except:
            print(f"WARNING: other_funcs value in {predictor.name} is not iterable. Skipping")
    


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
