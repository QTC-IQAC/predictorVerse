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
# from RFAA import main
from boltz1 import main
# from omegafold import main
# from chai import main
# from af3 import main

input_csv = sys.argv[1]
# provide output path input, if not, default to current path
#workspace_name = str: default Galaxy_42

workspace_name = "Galaxy42"

system_list = uu.read_input_csv(input_csv)

ww = uu.Workpath(workspace_name)

main(system_list,ww)

# print(len(system_list))

# gen_of_input(system_list,workspace_prefx)
# gen_of_runner(workspace_name)

# for system in system_list[:]:
#     print(system.name)
#     # uu.gen_fasta(system,".",mode="protein")
#     try:
#         gen_RFAA_input(system,workspace_prefx)
#     except: 
#         continue
# gen_RFAA_runner(workspace_name)

    # uu.lig_smiles_to_sdf(system,".")

#     gen_boltz_input(system, workspace_prefx)
# gen_boltz_runner(workspace_name)
