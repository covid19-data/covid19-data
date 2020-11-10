"""Extract country-level metadata (population, income, etc.) from Worldbank data."""
import pandas as pd

columns = ["country_code", "population", "region", "income"]

pd.concat(
    [
        pd.read_csv(snakemake.input[0], usecols=columns),
        pd.read_csv(snakemake.input[1], usecols=columns),
    ]
).to_csv(snakemake.output[0], index=False)
