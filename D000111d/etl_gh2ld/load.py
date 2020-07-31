"""
Load all of the transformed data into the excel spreadsheet as the line data.
"""
import os
import xlsxwriter
import numpy as np
import pandas as pd
from av360.backend.find_backend import start_backend
from D000111d.etl_gh2ld.column_definitions import COLUMN_KEYS
from D000111d.etl_gh2ld.extract import extract
from D000111d.etl_gh2ld.transform import transform
from D000111d.settings import DATA_DIR


def load(transformed_data, upload_calls=False, file_name='line_data.xlsx'):
    """
    Load the transformed data into an excel spreadsheet.
    """
    line_file = os.path.join(DATA_DIR, file_name)
    workbook = xlsxwriter.Workbook(line_file)

    key_sheet = {
        'manifest': 'Sample Manifest',
        'flowcell_qc': 'Flowcell QC',
        'control_qc': 'Control QC',
        'sample_qc': 'Sample QC',
        'calls': 'Variant Call Data',
        'cnv_calls': 'CNV Call Data',
        'reviewed_calls': 'Reviewed Calls'
    }

    colnames = set()
    for key in key_sheet:
        colnames = colnames.union(transformed_data[key].columns)
    def_names = sorted(list(colnames))

    # First save the sheet describing all of the columns
    df_terms = pd.DataFrame(pd.Series(COLUMN_KEYS)).loc[def_names].reset_index()
    df_terms.columns = ['Term', 'Definition']
    df_to_worksheet(df_terms, 'Column Definitions', workbook, column_keys=None)

    # and then the sheets that actually have the data
    for key, sheet_name in key_sheet.items():
        df_to_worksheet(transformed_data[key], sheet_name, workbook)
    workbook.close()

    # Load set of calls to server in a CSV format so that other studies can use it
    if upload_calls:
        df_dump = transformed_data['call_dump']
        sftp_client = start_backend('sftp')
        # server_location = '/ghds/ivd/analytical_validation/D_000111/call_dump2.csv'
        # Updated to "call_dump3" for BIP 3.5.3 output (Jan 14, 2020)
        # Updated to "call_dump_bip3p5p3_74genes" and "call_dump_bip3p5p3_55genes" (Jan 30, 2020)
        server_location = '/ghds/ivd/analytical_validation/D_000111/call_dump_bip3p5p3_74genes.csv'
        with sftp_client.open(server_location, 'w') as file_pointer:
            df_dump.to_csv(file_pointer, index=False)
        server_location = '/ghds/ivd/analytical_validation/D_000111/call_dump_bip3p5p3_55genes.csv'
        df_dump = df_dump[df_dump['in_rrr']]
        with sftp_client.open(server_location, 'w') as file_pointer:
            df_dump.to_csv(file_pointer, index=False)


def df_to_worksheet(df_str, sheet_name, workbook, column_keys=COLUMN_KEYS):
    """
    Write a data frame to a new sheet of an existing excel workbook

    Args:
        df_str (pandas.DataFrame): the data
        sheet_name (str): name of the sheet
        workbook (xlsxwriter.workbook.Workbook): the existing excel workbook object

    Returns:
        xlsxwriter.workbook.Workbook: updated with the new sheet
    """
    df_str = df_str.astype(str)
    # Make sure that every column has a definition
    if column_keys is not None:
        assert set(df_str.columns).intersection(list(COLUMN_KEYS.keys())) == set(df_str.columns)
    worksheet = workbook.add_worksheet(sheet_name)

    # Write the header row
    bold = workbook.add_format({'bold': True})
    worksheet.write_row(0, 0, df_str.columns, bold)

    # Write the data rows
    row_counter = 1
    for _, row in df_str.iterrows():
        worksheet.write_row(row_counter, 0, row)
        row_counter += 1

    # Adjust all of the column widths
    header_chars = [len(str(col_name)) for col_name in df_str.columns]
    max_data_chars = df_str.applymap(len).apply(max).values
    max_widths = np.maximum(header_chars, max_data_chars)
    for i_col, width in enumerate(max_widths):
        worksheet.set_column(i_col, i_col, min(width * 1.1, 50))
    return workbook
