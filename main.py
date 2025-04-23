"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


You will execute this via command line and give
- csv path with system name and sequences and smiles
- the predictors to use
- A place to dump the outputs
"""

import utils as uu
import sys
import os
from RFAA import gen_RFAA_input, gen_RFAA_runner
from boltz1 import gen_boltz_input, gen_boltz_runner

input_csv = sys.argv[1]
# provide output path input, if not, default to current path
#workspace_name = str: default Galaxy_42

workspace_name = "Galaxy42"
workspace_path = "."
workspace_prefx = os.path.join(workspace_path, workspace_name) # TODO: fix this ugly thing with the paths (make them global ??)


uu.make_folder_structure(workspace_name=workspace_name)
system_list = uu.read_input_csv(input_csv)

for system in system_list[:]:
    print(system.name)
    # uu.gen_fasta(system,".",mode="protein")
    try:
        gen_RFAA_input(system,workspace_prefx)
    except: 
        continue
gen_RFAA_runner(workspace_name)

    # uu.lig_smiles_to_sdf(system,".")

#     gen_boltz_input(system, workspace_prefx)
# gen_boltz_runner(workspace_name)



print(len(system_list))
