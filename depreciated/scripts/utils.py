import pandas as pd


def country_name2code_dict(conv_table_path):
    return pd.read_csv(conv_table_path, index_col="country_name").to_dict()[
        "country_code"
    ]


def country_code2name_dict(conv_table_path):
    return pd.read_csv(conv_table_path, index_col="country_code").to_dict()[
        "country_name"
    ]


def convert_namecol_to_code(a_series, conv_table_path):
    d = country_name2code_dict(conv_table_path)
    return a_series.apply(lambda x: d[x])


def convert_codecol_to_name(a_series, conv_table_path):
    d = country_code2name_dict(conv_table_path)
    return a_series.apply(lambda x: d[x])
