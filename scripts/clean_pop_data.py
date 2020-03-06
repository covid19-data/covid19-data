"""Clean the population data by adding Taiwan and converting country names.

input[0]: population CSV extracted from worldbank data.
input[1]: addition table. Now it just contains Taiwan.
input[2]: country name conversion table from worldbank's name to JHU's name.

"""
import pandas as pd

df = pd.concat([pd.read_csv(snakemake.input[0]), pd.read_csv(snakemake.input[1])])
conv = pd.read_csv(snakemake.input[2]).set_index("worldbank").to_dict()["covid19"]
df.replace({"country": conv}).to_csv(snakemake.output[0], index=False)
