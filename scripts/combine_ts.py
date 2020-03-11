"""Convert CSSE COVID-19 dataset into a json file for visualization."""
import datetime
import json

import pandas as pd

DATA_TYPES = ("confirmed", "deaths", "recovered")

# {"Republic of Ireland" } this issue was fixed here:
# https://github.com/CSSEGISandData/COVID-19/issues/316
TO_BE_DROPPED = set([])
TO_BE_REPLACED = {
    "Vatican City": "Holy See",
    "Hong Kong SAR": "Hong Kong",
    "Iran (Islamic Republic of)": "Iran",
    "Macao SAR": "Macau",
    "Republic of Korea": "South Korea",
    "Republic of Moldova": "Moldova",
    "Russian Federation": "Russia",
    "Saint Martin": "St. Martin",
    "Taipei and environs": "Taiwan",
    "Viet Nam": "Vietnam",
    "occupied Palestinian territory": "Palestine",
}


def convert_date(x):
    """Convert the date ("month/day/year") into ISO format."""
    month, day, year = map(int, x.split("/"))
    year = 2000 + year
    return datetime.date(year, month, day).isoformat()


def load_data(input_f):
    """Load time series data and return dataframe.

    TODO: incorporate locations later.
    """
    return (
        pd.read_csv(input_f)
        .drop(columns=["Lat", "Long"])
        .rename(columns={"Province/State": "state", "Country/Region": "country"})
    )


def melt_df(df, value_name):
    """Convert wide-formatted df to long, making it tidy."""
    return df.melt(id_vars=["country", "state"], var_name="date", value_name=value_name)


def melt_and_merge(dfs, values=DATA_TYPES, merge_on=("country", "state", "date")):
    """Melt the three dfs and then merge them together."""
    molten_dfs = [melt_df(df, vname) for df, vname in zip(dfs, values)]
    df = (
        molten_dfs[0]
        .merge(molten_dfs[1], on=merge_on)
        .merge(molten_dfs[2], on=merge_on)
    )
    df.loc[:, "date"] = df.date.apply(convert_date)
    df.loc[:, "country"] = df.country.replace(TO_BE_REPLACED)
    return df[~df.country.isin(TO_BE_DROPPED)]


melt_and_merge(list(map(load_data, snakemake.input[0:3]))).to_csv(
    snakemake.output[0], index=False
)
