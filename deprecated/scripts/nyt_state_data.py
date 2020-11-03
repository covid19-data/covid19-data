"""attach state_code to the NYT data and sort/normalize."""
import sys

import pandas as pd

NYT_STATE_DATA = sys.argv[1]
USA_STATE_CODE_DATA = sys.argv[2]
OUT_FNAME = sys.argv[3]


def extract_normalized_state_dfs(df, date_idx):
    """Normalize each state data by using the same date range."""
    dfs = []
    for group in df.groupby(["state_code"]):
        state_df = (
            group[1][["total_cases", "total_deaths"]]
            .reindex(date_idx, method="pad")
            .fillna(0)
        )
        state_df["total_cases"] = state_df["total_cases"].astype("int64")
        state_df["total_deaths"] = state_df["total_deaths"].astype("int64")
        state_df.loc[:, "state_code"] = group[0]
        state_df.loc[:, "state_name"] = group[1]["state_name"][-1]
        state_df.loc[:, "fips"] = group[1]["fips"][-1]
        dfs.append(
            state_df[
                ["state_code", "fips", "state_name", "total_cases", "total_deaths"]
            ]
        )
    return dfs


state_code_dict = pd.read_csv(USA_STATE_CODE_DATA, index_col="state_name").to_dict()[
    "state_code"
]
df = pd.read_csv(NYT_STATE_DATA, index_col=0, parse_dates=True).rename(
    columns={"cases": "total_cases", "deaths": "total_deaths", "state": "state_name"}
)
df["state_code"] = df.state_name.apply(lambda x: state_code_dict[x])
idx = pd.date_range(df.index.min(), df.index.max())

pd.concat(extract_normalized_state_dfs(df, idx)).reset_index().rename(
    columns={"index": "date"}
).to_csv(OUT_FNAME, index=False)

# df.sort_values(["state_code", "date"])[
# ["date", "state_code", "fips", "state_name", "total_cases", "total_deaths"]
# ].to_csv(OUT_FNAME, index=False)
