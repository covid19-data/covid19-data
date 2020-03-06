"""Extract national-level population and other metadata (from Worldbank).

country name, country code, region, income level, and population (2018).
"""

import pandas as pd

df = pd.read_csv(snakemake.input[0], skiprows=4)[
    ["Country Name", "Country Code", "2018"]
].rename(
    columns={
        "Country Name": "country",
        "Country Code": "country_code",
        "2018": "population",
    }
)

mdf = pd.read_csv(snakemake.input[1])[["Country Code", "Region", "IncomeGroup"]].rename(
    columns={
        "Country Code": "country_code",
        "Region": "region",
        "IncomeGroup": "income",
    }
)
df.merge(mdf, on="country_code", how="left").to_csv(snakemake.output[0], index=False)
