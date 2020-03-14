"""Extract country-level metadata (population, income, etc.) from Worldbank data."""
import pandas as pd

pd.read_csv(
    snakemake.input[0], usecols=["country_code", "population", "region", "income"]
).to_csv(snakemake.output[0], index=False)
