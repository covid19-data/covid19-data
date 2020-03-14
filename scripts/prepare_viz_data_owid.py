"""Prepare the CFR data visualization data from OWID (WHO)."""
import json

import pandas as pd

import utils


def load_OWID(fname):
    return (
        pd.read_csv(fname, index_col=0, parse_dates=True)
        .rename(columns={"location": "country_name"})
        .fillna(0, downcast="infer")
    )


def load_country_metadata(fname):
    return pd.read_csv(fname).set_index("country_code").to_dict()


def prepare_data_structure(df, date_idx, meta_dict, code2name_dict):
    all_dates = date_idx.strftime("%Y-%m-%d")

    data = []
    for group in df.groupby(["country_code"]):
        cntry_df = group[1][["total_cases", "total_deaths"]].reindex(
            date_idx, fill_value=0
        )
        if group[0] in {"---", "WWW"}:
            continue
        try:
            country_data = {
                "country_code": group[0],
                "country_name": code2name_dict[group[0]],
                "population": meta_dict["population"][group[0]],
                "region": meta_dict["region"][group[0]],
                "income": meta_dict["income"][group[0]],
                "confirmed": list(zip(all_dates, cntry_df["total_cases"])),
                "deaths": list(zip(all_dates, cntry_df["total_deaths"])),
            }
            data.append(country_data)
        except KeyError:
            print("key error!", group[0])
            continue
    return data


df = load_OWID(snakemake.input[0])
idx = pd.date_range(df.index.min(), df.index.max())
meta_dict = load_country_metadata(snakemake.input[1])

df.loc[:, "country_code"] = utils.convert_namecol_to_code(
    df.country_name, snakemake.input[2]
)

data = prepare_data_structure(
    df, idx, meta_dict, utils.country_code2name_dict(snakemake.input[3])
)
open(snakemake.output[0], "w").write(json.dumps(data))
