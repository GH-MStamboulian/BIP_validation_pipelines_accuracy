# first line: 9
def extract(line_data_path=None):
    """
    Load the data from the line data spreadsheet into memory.
    """
    if  line_data_path is None:
        file_path = os.path.join(DATA_DIR, 'line_data.xlsx')
    else:
        file_path = line_data_path
    key_sheet = {
        'manifest': 'Sample Manifest',
        'sample_qc': 'Sample QC',
        'calls': 'Variant Call Data',
        'cnv_calls': 'CNV Call Data',
        'reviewed_calls': 'Reviewed Calls'
    }
    extracted_data = dict()
    for key, sheet_name in key_sheet.items():
        extracted_data[key] = pd.read_excel(file_path, sheet_name=sheet_name)
    return extracted_data
