from predictorVerse.utils import System, Predictor, RunnerParams


# Double braces to not confuse with {0}

px_prot_json = """      {{
        "proteinChain": {{
            "count": 1,
            "id": ["{letter}"],
            "sequence": "{seq}"
        }}
      }}"""



px_lig_json = """      {{
        "ligand": {{
            "count": 1,
            "id": ["{letter}"],
            "ligand": "{seq}"
        }}
      }}"""

px_rna_json = """ {{
        "rnaSequence": {{
            "count": 1,
            "id": ["{letter}"],
            "sequence": "{seq}"
        }}
      }}"""

px_dna_json = """ {{
        "dnaSequence": {{
            "count": 1,
            "id": ["{letter}"],
            "sequence": "{seq}"
          }}
      }}"""

px_json_template = """[{{
    "name": "{system.name}",
    "sequences": [
      {input}
    ]
}}]
"""

exec_command = "protenix pred -i $inputs_dir -s SEEDS -c {recycles} -o $outputs_dir --use_msa true --need_atom_confidence true"



def calc_number_of_seeds_from_samples(system:System, predictor: Predictor, samples:int, *args):
  """
  From the number of samples calculate how many seeds are needed to achieve that many samples
  and write them in the template
  """
  # Assert multiple of 5
  # If not, get the nearest highest multiple of 5 and divide by 5
  seeds_num = samples//5 if samples%5 == 0 else samples//5 + 1

  # Get a list of the seeds
  seed_list = list(map(str,range(42,42+seeds_num)))
  seed_str =",".join(seed_list)
  
  # write it in prot_lig_temp
  predictor.main_cmds = predictor.main_cmds.replace("SEEDS",seed_str)


runner_params = RunnerParams(header=False,
                             extra_cmds=False,
                             extra_inputs=False,
                             looper=False,
                             jobarray=False)


protenix_data = Predictor(name= "Protenix",
            prot_temp= px_prot_json,
            lig_temp= px_lig_json,
            rna_temp= px_rna_json,
            dna_temp= px_dna_json,
            joiner=",\n",
            prot_lig_temp= px_json_template,
            input_extension= ".json",
            other_funcs=[calc_number_of_seeds_from_samples],
            extra_cmds="",
            main_cmds=exec_command,
            runner_params=runner_params
)