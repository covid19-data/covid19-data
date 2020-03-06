"""Extract lat,lng coordinates from a time series file."""
import pandas as pd

df = pd.read_csv(snakemake.input[0]).rename(
    columns={
        "Province/State": "state",
        "Country/Region": "country",
        "Lat": "lat",
        "Long": "lng",
    }
)[["country", "state", "lat", "lng"]]

df.to_csv(snakemake.output[0], index=False)
