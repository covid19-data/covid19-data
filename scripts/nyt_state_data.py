"""attach state_code to the NYT data and sort/normalize."""
import sys

import pandas as pd

NYT_STATE_DATA = sys.argv[1]
USA_STATE_CODE_DATA = sys.argv[2]
OUT_FNAME = sys.argv[3]

state_code_dict = pd.read_csv(USA_STATE_CODE_DATA, index_col="state_name").to_dict()[
    "state_code"
]
df = pd.read_csv(NYT_STATE_DATA).rename(
    columns={"cases": "total_cases", "deaths": "total_deaths", "state": "state_name"}
)
df["state_code"] = df.state_name.apply(lambda x: state_code_dict[x])
df.sort_values(["state_code", "date"])[
    ["date", "state_code", "fips", "state_name", "total_cases", "total_deaths"]
].to_csv(OUT_FNAME, index=False)
