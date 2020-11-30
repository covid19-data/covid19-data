"""attach state_code to the NYT data and sort/normalize."""
import sys
import http.client as http
http.HTTPConnection._http_vsn = 10
http.HTTPConnection._http_vsn_str = 'HTTP/1.0'

import pandas as pd

NYT_STATE_DATA = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"
OUT_FNAME = "/Users/Tal/Desktop/covid19-data/output/cases/cases_us_states_nyt.csv"

# Source: https://gist.github.com/rogerallen/1583593
us_state_abbrev = {
    'Alabama': 'AL',
    'Alaska': 'AK',
    'American Samoa': 'AS',
    'Arizona': 'AZ',
    'Arkansas': 'AR',
    'California': 'CA',
    'Colorado': 'CO',
    'Connecticut': 'CT',
    'Delaware': 'DE',
    'District of Columbia': 'DC',
    'Florida': 'FL',
    'Georgia': 'GA',
    'Guam': 'GU',
    'Hawaii': 'HI',
    'Idaho': 'ID',
    'Illinois': 'IL',
    'Indiana': 'IN',
    'Iowa': 'IA',
    'Kansas': 'KS',
    'Kentucky': 'KY',
    'Louisiana': 'LA',
    'Maine': 'ME',
    'Maryland': 'MD',
    'Massachusetts': 'MA',
    'Michigan': 'MI',
    'Minnesota': 'MN',
    'Mississippi': 'MS',
    'Missouri': 'MO',
    'Montana': 'MT',
    'Nebraska': 'NE',
    'Nevada': 'NV',
    'New Hampshire': 'NH',
    'New Jersey': 'NJ',
    'New Mexico': 'NM',
    'New York': 'NY',
    'North Carolina': 'NC',
    'North Dakota': 'ND',
    'Northern Mariana Islands':'MP',
    'Ohio': 'OH',
    'Oklahoma': 'OK',
    'Oregon': 'OR',
    'Pennsylvania': 'PA',
    'Puerto Rico': 'PR',
    'Rhode Island': 'RI',
    'South Carolina': 'SC',
    'South Dakota': 'SD',
    'Tennessee': 'TN',
    'Texas': 'TX',
    'Utah': 'UT',
    'Vermont': 'VT',
    'Virgin Islands': 'VI',
    'Virginia': 'VA',
    'Washington': 'WA',
    'West Virginia': 'WV',
    'Wisconsin': 'WI',
    'Wyoming': 'WY'
}

# thank you to @kinghelix and @trevormarburger for this idea
abbrev_us_state = dict(map(reversed, us_state_abbrev.items()))

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

df = pd.read_csv(NYT_STATE_DATA, index_col=0, parse_dates=True).rename(
    columns={"cases": "total_cases", "deaths": "total_deaths", "state": "state_name"}
)
df['state_code'] = df.state_name.apply(lambda x: us_state_abbrev[x])
idx = pd.date_range('2019-12-31', df.index.max())

pd.concat(extract_normalized_state_dfs(df, idx)).reset_index().rename(
    columns={"index": "date"}
).to_csv(OUT_FNAME, index=False)
