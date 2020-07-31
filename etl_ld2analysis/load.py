"""
Load the line data from the excel spreadsheet into data frames that are ready to go. The locations
of the files are in the settings file for this analysis.
"""
from joblib import Memory
from D000111d.etl_ld2analysis.extract import extract
from D000111d.etl_ld2analysis.transform import transform
from D000111d.settings import CACHE_DIR


memory = Memory(CACHE_DIR, verbose=False)
c_extract = memory.cache(extract)


def load(line_data_path=None, use_rrr=False):
    """
    Load the line data for analysis
    """
    print("HIIIIIIIIIIIIIIIIIIIIIIIII")
    print(line_data_path)
    extracted_data = c_extract(line_data_path=line_data_path)
    transformed_data = transform(extracted_data, use_rrr=use_rrr)
    return transformed_data
