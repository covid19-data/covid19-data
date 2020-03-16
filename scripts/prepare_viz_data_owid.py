"""Prepare the CFR data visualization data from OWID (WHO)."""
import json

import pandas as pd

import utils


def prepare_data_structure(df, meta_dict, code2name_dict):
    data = []
    for g in df.groupby(["country_code"]):
        code = g[0]
        cntry_df = g[1]
        try:
            country_data = {
                "country_code": code,
                "country_name": code2name_dict[code],
                "population": meta_dict["population"][code],
                "region": meta_dict["region"][code],
                "income": meta_dict["income"][code],
                "confirmed": list(zip(cntry_df.date, cntry_df.total_cases)),
                "deaths": list(zip(cntry_df.date, cntry_df.total_deaths)),
            }
            data.append(country_data)
        except KeyError:
            print("metadata doesn't exist for: ", code, code2name_dict[code])
            continue
    return data


df = pd.read_csv(snakemake.input[0])
meta_dict = pd.read_csv(snakemake.input[1]).set_index("country_code").to_dict()
data = prepare_data_structure(
    df, meta_dict, utils.country_code2name_dict(snakemake.input[2])
)
open(snakemake.output[0], "w").write(json.dumps(data, indent=2))
