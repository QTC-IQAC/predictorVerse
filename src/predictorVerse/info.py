"""
Berta Bori Bru - IQAC-CSIC
Spring 2025


Dictionary with the information of the predictors available for the program
"""

from predictor_templates.af3 import af3_data
from predictor_templates.chai import chai_data
from predictor_templates.boltz1 import boltz_data
from predictor_templates.omegafold import of_data
from predictor_templates.rfaa import rfaa_data
from predictor_templates.boltz1x import boltz1x_data


predictors_library = { "AF3": af3_data, # AlphaFold3
                    "Chai": chai_data, # Chai-1
                    "Boltz": boltz_data, # Boltz-1
                    "OF": of_data, # OmegaFold
                    "RFAA": rfaa_data, # RosettaFold All-Atom
                    "Boltz1x": boltz1x_data # Boltz-1x

}

