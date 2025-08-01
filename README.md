
# Welcome to predictorVerse!

The universe of protein structure predictors.
This is a highly customizable package that generates several input and runner scripts for different protein-ligand systems. 

You give:
- A .json with protein sequence and/or ligand SMILES for different systems
- Which predictor(s) do you want to use

You get:
- A series of organized inputs and running scripts to make the predictions
- Peace of mind

# Requirements, intallation and stuff
## Installation
```bash
git clone blablabalba
cd predictorVerse
pip install .
```

## Execution
predictorVerse can be run through CLI:
```bash
cd test
predictorVerse test.json        # Generates inputs for all the available predictors
predictorVerse test.json -p AF3 # Only generates inputs for AF3
```

or through module import 
```python
from predictorVerse.main import gen_predictor_inputs
gen_predictor_inputs("test/test.json")  # Generates inputs for all the available predictors
gen_predictor_inputs("test/test.json",["AF3"])
```

Note that the predictor keys passed to predictorVerse have to be available. To check that, see the help message of the program available using
```bash
python main.py -h
```


## Input file
The .json input file accepts any number of systems, each containing any number (max 27 in total) of proteins and ligands (dna and rna might come someday). The input should look like:
```json
{
    "name_system01":
        {"protein":["AMINO ACID SEQ 1",
                    "AMINO ACID SEQ 2",
                    "..."]
        ,
        "ligand":"SMILES SEQ 1"

        },

    "name_system02":
        {"protein":["..."],
        "ligand":["..."]
        },
    
    "name_system03":
        {"protein":["..."]
        }

}
```
There can be a different number of ligands and protein in the defined systems. If the "ligand" or "protein" properties are not defined, there is no problem. If both of them are not defined, the program generates an empty input file.

## Output
predictorVerse generates a series of folders in the **current folder** you are executing it in:
```
.
└── current_folder/
    ├── AF3/
    │   ├── inputs
    │   ├── outputs
    │   └── runners
    ├── Chai
    ├── Boltz
    └── ...
```
- **inputs**: Here the inputs for all your systems will be generated, for each predictor
- **outputs**: Here the outputs are expected to be saved after you make the predictions
- **runners**: Here are generated the runner files to make the predictions

## Execution of the runner files
First, ensure you have all the requirements and activated environments needed to run the predictor.
Then, in `current_folder` execute the runner file. For example:
```bash
source AF3/runners/AF3_runner.sh
```

# predictorVerse customisation
As previously mentioned, this is a highly customizable package (and thus might not be the most efficient), so you are expected to tweek some parts of it. 
predictorVerse works by filling in the gaps in a bunch of template chunks and then assembling them together. In this way it is easy to be adapted to the preferences of the user. The templates are stored in the `predictor_templates` folder in a bunch of `.py` files. 

The predictor files contain the necessary templates to generate the inputs and extra functions to generate extra files, if needed (this is the case of RFAA, which needs to generate .sdf for the ligands to run). It also contains the exact commands needed to run the predictor and any necessary "extra commands" (for example, loading a module in a cluster).

The `jobscripts_templates.py` file contains the necessary templates to build the runner files. From the general structure of a runner to cluster header chunks. The cluster specifications use SLURM.

Finally, the `info.py` file contains a dictionary that tells the main program the name of a predictor and its parameters.


## I want to change the runner generation, I get some cluster flags that I don't want
That's probably the most important customization, as you will run the predictors in different ways that I might do. In any predictor file in `predictor_templates` you will find a `runner_params` variable. This variable of class `RunnerParams` has some arguments that control what is written on the runner file:
- **header** (str or bool): Which header to use for the runner. Used to select the cluster flags. 
    - If you pass False, no header is selected. 
    - If True, gets the default one (defined in `predictor_templates/jobscripts_templates.py`). 
    - If you pass a string it will get the one defined under that name in the dictionary `headers_dict` in `predictor_templates/jobscripts_templates`.
- **extra_cmds** (bool): If there are any extra commands you want to add to the runner that need to run before calling the main ones. The actual commands will need to be passed through the `Predictor` class in the same file.
- **extra_inputs** (bool): If True, adds option to pass the input file to the runner through command line. For example: `source AF3/AF3_runner.sh AF3/inputs/file.json`. I only use it when I need to generate jobarrays.
- **looper** (bool): If you want to run the input files in a sequential manner. Usage is recommended if the predictor has no form of batch processing. 
- **jobarray** (bool): Adds special line to the file for jobarrays. It indicates the number of jobs to expect in the array.
                             
Let's run through some examples. I run AF3 in CSUC cluster (`header="csuc"`) with jobarrays (`jobarray=True` and `extra_inputs=True`). To do that I need to load the AF3 module in the cluster (`extra_cmds=True`), so I will define that command and pass it to the `Predictor` class. Finally, AF3 has batch processing (`looper=False`). The options would look like this:

```python 
runner_params = RunnerParams(header="csuc",
                            extra_cmds=True,
                            extra_inputs=True,
                            looper=False,
                            jobarray=True)
```


For Chai-1, I have it installed locally (`header=False` and `jobarray=False`). It does not need any extra command (`extra_cmds=False`) and I want it to run all the structures with 1 script (`extra_inputs=False`). Chai-1 has no batch processing that I know of (`looper=True`). All together they are:

```python 
runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=True,
                             jobarray=False)
```

## I want to add another predictor to predictorVerse! How do I do it?
Well, there are a lot of things to do!
First, create a `.py` file in `predictor_templates`. Then add:
```python
from utils import Predictor, RunnerParams
```
Now, you should define
1. A string with the format to pass a protein sequence
2. A string with the format to pass a ligand sequence
3. The general structure of the input template
4. The command to execute the predictor from a bash script
5. Any extra commands needed for execution


For steps 1 and 2 you can use the following formatting parameters:
- `{system.name}`: the system name (for example the pdb name) for identification
- `{predictor.inputs_unix}` and `{predictor.inputs_unix}`: the relative path to the inputs and outputs folder, respectively
- `{seq}`: the protein amino acid sequence or the ligand SMILEs
- `{ii}`: the number of protein and ligand sequence from the ones defined in the system
- `{letter}`: chain letter for the pdb or cif

For step 3, the formatting parameters are:
- `{system.name}`: the system name (for example the pdb name) for identification
- `{predictor.inputs_unix}` and `{predictor.inputs_unix}`: the relative path to the inputs and outputs folder, respectively
- `{input}`: where to put the information of the proteins and ligands

If any other operations need to be made to generate the inputs, you should define here the corresponding functions that do that. The functions must have as arguments `system: System, predictor: Predictor` (even if Predictor is not used). Take into account that the protein and ligand sequences are stored in a list under `system.protein` and `system.ligand`.

After all that, you need to make a `RunnerParams` object (see previous section) and a `Predictor` object.
The prediction object takes as arguments:
- **name** (str): The name of the predictor.
- **prot_temp** (str): The formatted string for protein sequences (step 1)
- **lig_temp** (str): The formatted string for ligand sequences (step 2)
- **joiner** (str): The string by which to join the protein and ligand templates. The usual is `"\n"`.
- **prot_lig_temp** (str): The general template (step 3)
- **input_extension** (str): The extension for the input. Example: `".fasta"`
- **other_funcs** (list(function)): A list containing the auxiliary functions to run for input generation.
- **extra_cmds** (str): The string containing extra execution commands to run before the main execution commands.
- **main_cmds** (str): The string containing the main execution commands
- **runner_params** (RunnerParams): The object containing the options for runner generation

As a final step, you should add the predictor in `info.py`:
1. Import the `Predictor` object that you just created
2. Add it to the dictionary with the predictor id (the name by which to call it through the command line) as key and the `Predictor` object as value. 

Now, the only thing that's left is try if it works!
    
## I want to modify the strings in jobscripts_templates. Can I get a guide into what's in there?
Sure! Here is what's in there:
- **general_temp**: General template of a runner.     Contains the formatting strings:
    - `{header}`: place for cluster headers and such
    - `{inputs_section}`: here the input and output paths are defined
    - `{extras}`: here go the extra commands (to run before the main ones)
    - `{execution_section}`: here go the main commands and looper (if specified)
- **jobarr_temp**: extra line to add in header for jobarrays in SLURM
- **extra_inputs_txt**: bash exerpt that make the runner take the input file as command line argument
- **basic_input_temp**: Basic definition of input and output directories
- **looper_temp**: Template to run the main execution commands in a loop
- **headers_dict**: dictionary containing the available headers for SLURM clusters. 

