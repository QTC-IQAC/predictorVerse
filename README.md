
# Welcome to predictorVerse!

The universe of protein structure predictors.
This package will serve to generate several input and runner scripts for different protein-ligand pairs. 

You give:
- A csv with sequence (and SMILES) for different compounds
- Which predictor(s) do you want to use
- Others

You get:
- A series of organized inputs and running scripts to make the predictions
- Peace of mind

# Requirements, intallation and stuff
TO DO xd

# I want to add another predictor to predictorVerse! How do I do it?
I'm glad you want to contribute to this project. You need to:
- Create a .py file with all the stuff necessary for you predictor. It should contain:
    - A formated string with the necessary input content
    - A formated string with the execution command for the predictor
    - Any additional functions that create files or do whatever that is necessary for the predictor to work. Note that they should only have as input the classes System and Worspace (defined in utils).
    - A dictionary with all the necessary information. Example:
    ```python
    rfaa_data = {"name": "RFAA",
            # "prot_temp": 
            # "lig_temp": 
            "prot_lig_temp": rfaa_yaml_template,
            "input_extension": ".yaml",
            # "runner_temp":
            "other_funcs":[gen_prot_fasta, lig_smiles_to_sdf]
    ```
- Add in info file the following information:
    - Import the data dictionary in the .py file you've created
    - Add to the dictionary the name by which you want to call the predictor (key) and the data dictionary (value)
    

