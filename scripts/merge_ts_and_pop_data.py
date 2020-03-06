"""Merge time series data with worldbank's population and income dataset.

input[0]: combined time series data
input[1]: population, region, and income data.
"""
import json

import pandas as pd


def build_json(df, meta_dict):
    """Build the final json file.

    :df: the molten, merged df that contains all time series data.
    :pop_df: population, region, and income
    :returns: a list of objects that will be saved as a json file.

    """
    data = []
    for group in df.groupby(["country"]):
        if group[0] == "Others":
            continue

        # ignore regions within a country. TODO: incorporate regions.
        temp_df = (
            group[1]
            .groupby(["date"])[["confirmed", "deaths", "recovered"]]
            .sum()
            .reset_index()
        )
        country_data = {
            "country": group[0],
            "population": meta_dict["population"][group[0]],
            "region": meta_dict["region"][group[0]],
            "income": meta_dict["income"][group[0]],
            "confirmed": temp_df[["date", "confirmed"]].values.tolist(),
            "deaths": temp_df[["date", "deaths"]].values.tolist(),
            "recovered": temp_df[["date", "recovered"]].values.tolist(),
        }
        data.append(country_data)
    return data


data = build_json(
    pd.read_csv(snakemake.input[0]),
    pd.read_csv(snakemake.input[1]).set_index("country").to_dict(),
)

open(snakemake.output[0], "w").write(json.dumps(data))
