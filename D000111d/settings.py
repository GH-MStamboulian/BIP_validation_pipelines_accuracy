import os
import yaml
import pandas as pd
from D000111d import __path__ as report_namespace_path
import shutil


#: The path to the folder for this study report
REPORT_DIR = report_namespace_path[0]
path = REPORT_DIR +"/joblib.cache/"

if (os.path.isdir(path)):
    shutil.rmtree(os.path.join(REPORT_DIR, 'joblib.cache'))

#: The directory of CSV and Excel files within the repository
DATA_DIR = os.path.join(REPORT_DIR, 'data')

#: Where to dump LaTeX tables
LATEX_DIR = os.path.join(REPORT_DIR, 'latex')

#: The joblib cache
CACHE_DIR = os.path.join(REPORT_DIR, 'joblib.cache')

#: Display names for the ldt, cdx, and mdl assays
LONG_ASSAY_NAMES = {'ldt': 'G360 LDT', 'mdl': 'LBP70', 'cdx': 'G360 CDx'}

#: The bed files that define the reportable range
PARAM_DIR = '/ghds/ivd/analytical_validation/ghpipeline-3.5.3-8857b98/parameter_sets/'
bed_file_names = {
    'GH2.10': 'G360/v2.10/GH2.10_annotated_reportable_regions.bed',
    'GH2.10.1': 'G360/v2.10.1/GH2.10.1_reportable_regions.bed',
    'GH2.11': 'G360/v2.11/GH2.11_reportable_regions.bed',
    'MDACC_v1.0': 'MDACC/v1.0/MDACC_v1.0_annotated_reportable_regions.bed'
}
blacklist_beds = {
    'GH2.10': 'G360/v2.10/GH2.10_blacklist.bed',
    'GH2.10.1': 'G360/v2.10.1/GH2.10.1_blacklist.bed',
    'GH2.11': 'G360/v2.11/GH2.11_blacklist.bed',
    'MDACC_v1.0': 'MDACC/v1.0/MDACC_v1.0_blacklist.bed'
}
BED_FILE_PATHS = {panel: os.path.join(PARAM_DIR, file_name)
                  for panel, file_name in bed_file_names.items()}
BLACKLIST_BED_FILE_PATHS = {panel: os.path.join(PARAM_DIR, file_name)
                            for panel, file_name in blacklist_beds.items()}

# Load the priors from the yaml
with open(os.path.join(DATA_DIR, 'priors.yaml'), 'r') as file_pointer:
    PRIORS = yaml.load(file_pointer)

# Loads the "priors" for the restricted reportable range
with open(os.path.join(DATA_DIR, 'priors_rrr.yaml'), 'r') as file_pointer:
    PRIORS_RRR = yaml.load(file_pointer)


LDT_PRIORS = pd.Series(PRIORS['category']['prior'])
LDT_PRIORS.name = 'prior'
N_VAR_SPACE = pd.Series(PRIORS['category']['n_var'])
N_VAR_SPACE.name = 'n_var'

#: LLCI AC for PPA
PPA_AC = pd.Series({('cs', 'snv'): 0.72,
                    ('panel', 'snv'): 0.74,
                    ('cs', 'indel'): 0.72,
                    ('panel', 'indel'): 0.76,
                    ('cs', 'cnv'): 0.60,
                    ('cs', 'fusion'): 0.70})
PPA_AC.index.names = ['panel_type', 'variant_type']

#: LLCI AC for NPA
NPA_AC = pd.Series({('cs', 'snv'): 0.996,
                    ('panel', 'snv'): 0.999,
                    ('cs', 'indel'): 0.998,
                    ('panel', 'indel'): 0.999,
                    ('cs', 'cnv'): 0.88,
                    ('cs', 'fusion'): 0.98})
NPA_AC.index.names = ['panel_type', 'variant_type']
