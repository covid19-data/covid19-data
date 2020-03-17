"""Extract lat,lng coordinates from a time series file."""
import pandas as pd

df = (
    pd.read_csv(snakemake.input[0])
    .rename(
        columns={
            "Province_State": "state",
            "Country_Region": "country_name",
            "Lat": "lat",
            "Long": "lng",
        }
    )[["country_name", "state", "lat", "lng"]]
    .drop_duplicates()
)

df.to_csv(snakemake.output[0], index=False)
