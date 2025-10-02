"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


Functions to generate runner files for the predictors
"""



import os 
from predictorVerse.utils import System, Predictor
import predictorVerse.predictor_templates.jobscripts_templates as jt



def get_header(header_key: bool | str) -> str:
    if isinstance(header_key, str):
        try:
            return jt.headers_dict[header_key]
        except:
            print(f"{header_key} header is not found. Returning default one")
            return jt.headers_dict["default"]
    elif isinstance(header_key, bool) and header_key is True:
        print("No header selected. Returning default one")
        return jt.headers_dict["default"]



def gen_jobarray(system_list:list[System], predictor:Predictor, max_cap_jobs=None | int) -> None:
    """
    
    max_cap_jobs : max number of jobs that can run at the same time.
                  If None, set to number of jobs
    """
    # Set variables
    num_jobs = len(system_list)
    num_cap_jobs = num_jobs if max_cap_jobs is None else max_cap_jobs
    header_txt = get_header(predictor.runner_params.header)
    jobarr_header = "#!/bin/bash\n" + header_txt.format(predictor=predictor) + jt.jobarr_temp.format(num_jobs, num_cap_jobs)

    # Runner file relative path
    runner_file_rel_path = os.path.join(".",predictor.runners,f"{predictor.name}_runner.sh")

    # Open jobarr file
    jobarr_file = os.path.join(predictor.runners,f"{predictor.name}_jobarray.sh")
    with open(jobarr_file, "w") as jobarr:
        # Write header with correct parameters
        jobarr.write(jobarr_header)

        # Write case start
        jobarr.write("case $SLURM_ARRAY_TASK_ID in\n")

        # Start looping through systems
        for ii, system in enumerate(system_list):
            jobarr.write(f"{ii+1}) {runner_file_rel_path} {system.name}{predictor.input_extension} ;;\n")
        
        jobarr.write("esac")


def fix_commands_for_loop(txt:str) -> str:
    """
    Adds \t to the heading of each line, because I like indentations
    Probably in operating systems outside of Linux does not work
    """
    cmds_list = txt.split("\n")
    new_txt = "\n\t".join(cmds_list)
    return new_txt


def gen_runner(system_list:list[System], predictor: Predictor, samples:int, recycles:int):
    # Getting header. TODO: when we have params for this, we will have a function here
    if not predictor.runner_params.header:
        header_txt = ""
    else:
        header_txt = get_header(predictor.runner_params.header).format(predictor=predictor) # check which header to use

    
    # Processing inputs
    basic_input = jt.basic_input_temp.format(predictor=predictor)
    
    if predictor.runner_params.extra_inputs:
        inputs_txt = "\n".join([basic_input, jt.extra_inputs_txt])
    else:
        inputs_txt = basic_input
    
    # Processing extra cmds
    if predictor.extra_cmds is not None or predictor.extra_cmds != "":
        extra_cmds_txt = predictor.extra_cmds
    else:
        extra_cmds_txt = ""

    # Processing main commands 
    ## Add sample and recycles info
    main_cmd = predictor.main_cmds.format(samples=samples, recycles=recycles)

    if predictor.runner_params.looper:
        new_cmds = fix_commands_for_loop(main_cmd)
        main_cmds_txt = jt.looper_temp.format(file_extension=predictor.input_extension,
                                           execution_command=new_cmds)
    else:
        main_cmds_txt = main_cmd
    
    # Assemble everything
    runner_txt = jt.general_temp.format(header=header_txt,
                                    inputs_section=inputs_txt,
                                    extras=extra_cmds_txt,
                                    execution_section=main_cmds_txt)
    runner_file = os.path.join(predictor.runners,f"{predictor.name}_runner.sh")

    
    with open(runner_file,"w") as rr:
        rr.write(runner_txt)
    
    # Generate jobarray if needed
    if predictor.runner_params.jobarray:
        gen_jobarray(system_list,predictor, max_cap_jobs=2)