"""
Run all of the steps to generate the line data file.
"""
from D000111d.etl_gh2ld.extract import extract
from D000111d.etl_gh2ld.transform import transform
from D000111d.etl_gh2ld.load import load


extracted_data = extract(use_old_bip=False)
transformed_data = transform(extracted_data)
load(transformed_data, upload_calls=False, file_name='line_data_new.xlsx')

extracted_data = extract(use_old_bip=True)
transformed_data = transform(extracted_data)
load(transformed_data, upload_calls=False, file_name='line_data_old.xlsx')
