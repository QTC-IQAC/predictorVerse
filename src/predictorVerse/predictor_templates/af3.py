from predictorVerse.utils import System, Predictor, RunnerParams


# Double braces to not confuse with {0}

af3_prot_json = """      {{
        "protein": {{
          "id": "{letter}",
          "sequence": "{seq}"
        }}
      }}"""



af3_lig_json = """      {{
        "ligand": {{
          "id": "{letter}",
          "smiles": "{seq}"
        }}
      }}"""

af3_rna_json = """ {{
        "rna": {{
          "id": "{letter}",
          "sequence": "{seq}"
        }}
      }}"""

af3_dna_json = """ {{
        "dna": {{
          "id": "{letter}",
          "sequence": "{seq}"
        }}
      }}"""

af3_json_template = """{{
    "name": "{system.name}",
    "modelSeeds": SEEDS,
    "sequences": [
      {input}
    ],
    "dialect": "alphafold3",
    "version": 2

}}
"""

extra_cmds = "module load alphafold/3"

exec_command = """singularity exec  \\
      \\
     --bind $inputs_dir:/root/$inputs_dir \\
     --bind $outputs_dir:/root/$outputs_dir \\
     --bind $cwd/models:/root/models \\
     --bind /data/ddbb/alphafold3:/root/public_databases \\
           /prod/container/alphafold3/alphafold3.sif \\
     python run_alphafold.py --model_dir=/root/models --db_dir=/root/public_databases \\
              --json_path=/root/$input \\
              --output_dir=/root/$outputs_dir --num_recycles {recycles}
"""

def calc_number_of_seeds_from_samples(system:System, predictor: Predictor, samples:int, *args):
  """
  From the number of samples calculate how many seeds are needed to achieve that many samples
  and write them in the template
  """
  # Assert multiple of 5
  # If not, get the nearest highest multiple of 5 and divide by 5
  seeds_num = samples//5 if samples%5 == 0 else samples//5 + 1

  # Get a list of the seeds
  seed_list =list(range(42,42+seeds_num))
  
  # write it in prot_lig_temp
  predictor.prot_lig_temp = predictor.prot_lig_temp.replace("SEEDS",str(seed_list))


runner_params = RunnerParams(header="csuc",
                             extra_cmds=True,
                             extra_inputs=True,
                             looper=False,
                             jobarray=True)


af3_data = Predictor(name= "AF3",
            prot_temp= af3_prot_json,
            lig_temp= af3_lig_json,
            rna_temp= af3_rna_json,
            dna_temp= af3_dna_json,
            joiner=",\n",
            prot_lig_temp= af3_json_template,
            input_extension= ".json",
            other_funcs=[calc_number_of_seeds_from_samples],
            extra_cmds=extra_cmds,
            main_cmds=exec_command,
            runner_params=runner_params
)



