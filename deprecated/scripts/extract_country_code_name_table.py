import pandas as pd

df = pd.read_csv(snakemake.input[0])
df.loc[:, "country_name_len"] = df.country_name.str.len()
df.sort_values(["country_code", "country_name_len"]).drop_duplicates(
    "country_code"
).drop("country_name_len", axis="columns").to_csv(snakemake.output[0], index=False)
