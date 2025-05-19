from af3 import af3_data
from chai import chai_data
from boltz1 import boltz_data
from omegafold import of_data
from RFAA import rfaa_data
from boltz1x import boltz1x_data

# TODO: Make a help command that indicates all correspondences between full name and short name

predictors_library = { "AF3": af3_data, # AlphaFold3
                    "Chai": chai_data, # Chai-1
                    "Boltz": boltz_data, # Boltz-1
                    "OF": of_data, # OmegaFold
                    "RFAA": rfaa_data, # RosettaFold All-Atom
                    "Boltz1x": boltz1x_data # Boltz-1x

}

