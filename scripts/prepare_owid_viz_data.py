import pandas as pd
import numpy as np
import json
import math

# prepare owid raw data:
df = pd.read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv", 
                     index_col = ["date"],
                     parse_dates = True, 
                     usecols=["date", "iso_code", "location", 
                              "population", "continent",
                              "total_cases", "total_deaths"]).rename(
    columns = {'iso_code': 'country_code', 'location': 'country_name'})

# prepare metadata from World Bank
metadata = pd.read_csv("https://raw.githubusercontent.com/hongtaoh/covid19-data/master/data_sources/metadata/worldbank/country_metadata.csv",
                      usecols=["Country Code", "Region"]).rename(
    columns = {'Country Code': 'country_code', 'Region': 'world_region'}
)

# get the full date range:
all_dates = pd.date_range(df.index.min(), df.index.max())

def extract_cntry_dfs(df): # input is df
    dfs = []
    for group in df.groupby("country_code"):
        cntry_df = (
            group[1].reindex(all_dates, method="pad")
        )
        cntry_df.loc[:,"total_cases"] = cntry_df.loc[:,"total_cases"]
        cntry_df.loc[:,"total_deaths"] = cntry_df.loc[:,"total_deaths"]
        cntry_df.loc[:,"country_code"] = group[0]
        cntry_df.loc[:,"country_name"] = group[1].loc[:,"country_name"][-1]
        cntry_df.loc[:,"continent"] = group[1].loc[:,"continent"][-1]
        cntry_df.loc[:,"population"] = group[1].loc[:,"population"][-1]
        dfs.append(
            cntry_df
        )
    return dfs

def fill_first_case_death_with_zero(df): # input is dfs
    for i in np.arange(0, len(dfs)):
        if math.isnan(dfs[i].loc[:, "total_cases"].iloc[0]):
            dfs[i].loc[:,"total_cases"].iloc[0] = 0
        if math.isnan(dfs[i].loc[:, "total_deaths"].iloc[0]):
            dfs[i].loc[:, "total_deaths"].iloc[0] = 0
    return dfs

dfs = extract_cntry_dfs(df)

dfs_first_zero_filled = fill_first_case_death_with_zero(dfs)

def merge_with_meta(df): #input should be dfs_first_zero_filled
    concat_df = pd.concat(dfs).fillna(method="ffill").reset_index().rename(
        columns={"index": "date"})
    #To change the original country codes of "KOS" and World to match metadata from WB:
    concat_df.loc[(concat_df.country_code == "OWID_KOS"), ('country_code')] = "XKX"
    concat_df.loc[(concat_df.country_code == "OWID_WRL"), ('country_code')] = "WLD"
    # To get the column of "world_region" in concat_df by merging with WB metadata
    left_join_df = pd.merge(concat_df, metadata, on = "country_code", how = "left")
    left_join_df.loc[:,'date'] = left_join_df.loc[:,'date'].dt.strftime('%Y-%m-%d')
    return left_join_df

left_join_df = merge_with_meta(dfs_first_zero_filled)

def prepare_data_structure(df, gby="country_code"): # input should be left_join_df
    data = []
    for g in df.groupby([gby]):
        code = g[0]
        cntry_df = g[1]
        try:
            country_data = {
                "country_code": code,
                "country_name": cntry_df.loc[:,"country_name"].iloc[0],
                "population": cntry_df.loc[:,"population"].iloc[0],
                "region": cntry_df.loc[:,"world_region"].iloc[0],
                "confirmed": list(zip(cntry_df.date, cntry_df.total_cases)),
                "deaths": list(zip(cntry_df.date, cntry_df.total_deaths)),
            }
            data.append(country_data)
        except KeyError:
            print("metadata doesn't exist for: ", code)
            continue
    return data

data = prepare_data_structure(left_join_df)

open("/Users/Tal/Desktop/covid19-data/output/cntry_stat_owid.json", "w").write(
    json.dumps(data, separators=(",", ":")))
