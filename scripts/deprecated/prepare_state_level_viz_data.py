"""Prepare the CFR data visualization data from OWID (WHO)."""
import json
import sys

import pandas as pd

CASE_DATA = "/Users/Tal/Desktop/covid19-data/output/cases/cases_us_states_nyt.csv"
OUT_FNAME = "/Users/Tal/Desktop/covid19-data/output/us_state_nyt.json"

def prepare_data_structure(df, gby="state_code"):
    data = []
    for g in df.groupby([gby]):
        code = g[0]
        temp_df = g[1]
        try:
            country_data = {
                "code": code,
                "name": temp_df["state_name"].iloc[0],
                "confirmed": list(zip(temp_df.date, temp_df.total_cases)),
                "deaths": list(zip(temp_df.date, temp_df.total_deaths)),
            }
            data.append(country_data)
        except KeyError:
            print("metadata doesn't exist for: ", code)
            continue
    return data


df = pd.read_csv(CASE_DATA)
data = prepare_data_structure(df)
open(OUT_FNAME, "w").write(json.dumps(data, separators=(",", ":")))
