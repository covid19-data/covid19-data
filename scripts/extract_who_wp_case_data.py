import pandas as pd

pd.concat(
    [pd.read_csv(snakemake.input[1]), pd.read_csv(snakemake.input[0])]
).drop_duplicates(subset=["date", "country_code"], keep="first").sort_values(
    ["country_code", "date"]
).to_csv(
    snakemake.output[0], index=False
)
