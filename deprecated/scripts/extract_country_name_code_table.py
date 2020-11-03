"""Create a table with country names and country code (ISO3166 Alpha 3)."""
import pandas as pd

pd.concat(
    [
        pd.read_csv(snakemake.input[0], usecols=["name", "ISO3166-1_alpha-3"]).rename(
            columns={"name": "country_name", "ISO3166-1_alpha-3": "country_code"}
        ),
        pd.read_csv(
            snakemake.input[0], usecols=["official_state_name", "ISO3166-1_alpha-3"]
        ).rename(
            columns={
                "official_state_name": "country_name",
                "ISO3166-1_alpha-3": "country_code",
            }
        ),
    ]
    + [
        pd.read_csv(x, usecols=["country_name", "country_code"])
        for x in snakemake.input[1:]
    ]
)[["country_code", "country_name"]].drop_duplicates().sort_values(
    by=["country_code", "country_name"]
).to_csv(
    snakemake.output[0], index=False
)
