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


def extract_cntry_dfs(df, date_idx, code2name_dict):
    all_dates = date_idx.strftime("%Y-%m-%d")

    dfs = []
    for group in df.groupby(["country_code"]):
        if group[0] in {"---"}:
            continue
        cntry_df = group[1][["total_cases", "total_deaths"]].reindex(
            date_idx, fill_value=0
        )
        cntry_df.loc[:, "country_code"] = group[0]
        cntry_df.loc[:, "country_name"] = group[1]["country_name"][-1]
        dfs.append(
            cntry_df[["country_code", "country_name", "total_cases", "total_deaths"]]
        )
    return dfs


df = load_OWID(snakemake.input[0])
idx = pd.date_range(df.index.min(), df.index.max())

df.loc[:, "country_code"] = utils.convert_namecol_to_code(
    df.country_name, snakemake.input[1]
)

pd.concat(
    extract_cntry_dfs(df, idx, utils.country_code2name_dict(snakemake.input[2]))
).reset_index().rename(columns={"index": "date"}).to_csv(
    snakemake.output[0], index=False
)
