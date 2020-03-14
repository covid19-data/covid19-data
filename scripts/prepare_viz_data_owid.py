"""Prepare the CFR data visualization data from OWID (WHO)."""
import json

import pandas as pd

# Load the OWID data
df = (
    pd.read_csv(snakemake.input[0], index_col=0, parse_dates=True)
    .rename(columns={"location": "country_name"})
    .fillna(0, downcast="infer")
)

# full date range
idx = pd.date_range(df.index.min(), df.index.max())
all_dates = idx.strftime("%Y-%m-%d")

# Load the worldbank data.
country_dict = (
    pd.read_csv(
        snakemake.input[1], usecols=["country_name", "population", "region", "income"]
    )
    .replace(
        {
            "Brunei Darussalam": "Brunei",
            "Congo, Dem. Rep.": "Democratic Republic of Congo",
            "Egypt, Arab Rep.": "Egypt",
            "Iran, Islamic Rep.": "Iran",
            "Korea, Rep.": "South Korea",
        }
    )
    .set_index("country_name")
    .to_dict()
)

data = []
for group in df.groupby(["country_name"]):
    cntry_df = group[1][["total_cases", "total_deaths"]].reindex(idx, fill_value=0)
    try:
        country_data = {
            "country": group[0],
            "population": country_dict["population"][group[0]],
            "region": country_dict["region"][group[0]],
            "income": country_dict["income"][group[0]],
            "confirmed": list(zip(all_dates, cntry_df["total_cases"])),
            "deaths": list(zip(all_dates, cntry_df["total_deaths"])),
        }
        data.append(country_data)
    except KeyError:
        print("key error!", group[0])
        continue

open(snakemake.output[0], "w").write(json.dumps(data))
