import pandas as pd


def country_code_dict(conv_table_path):
    return pd.read_csv(conv_table_path, index_col='country_name').to_dict()['country_code'])
