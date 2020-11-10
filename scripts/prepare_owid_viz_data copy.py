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
# input is df
def extract_cntry_dfs(df): 
    dfs = []
    for group in df.groupby("country_code"):
        cntry_df = (
            group[1].reindex(all_dates, method="pad")
        )
        cntry_df["total_cases"] = cntry_df["total_cases"]
        cntry_df["total_deaths"] = cntry_df["total_deaths"]
        cntry_df["country_code"] = group[0]
        cntry_df["country_name"] = group[1]["country_name"][-1]
        cntry_df["continent"] = group[1]["continent"][-1]
        cntry_df["population"] = group[1]["population"][-1]
        dfs.append(
            cntry_df
        )
    return dfs

def fill_first_case_death_with_zero(df): # input is dfs
    for i in np.arange(0, len(df)):
        if (df[i].head(1).total_cases.isnull()[0] & df[i].head(1).total_deaths.isnull()[0]):
            df[i][0:1] = [df[i].country_code[-1],
            df[i].continent[-1],
            df[i].country_name[-1],
            0,
            0, 
            df[i].population[-1]
            ]
        if (df[i].head(1).total_cases.isnull()[0] & df[i].head(1).total_deaths.notnull()[0]):
            df[i][0:1] = [df[i].country_code[-1],
            df[i].continent[-1],
            df[i].country_name[-1],
            0,
            df[i].total_deaths[0],
            df[i].population[-1]
            ]
        if (df[i].head(1).total_cases.notnull()[0] & df[i].head(1).total_deaths.isnull()[0]):
            df[i][0:1] = [df[i].country_code[-1],
            df[i].continent[-1],
            df[i].country_name[-1],
            df[i].total_cases[0],
            0,
            df[i].population[-1]
            ]
    return df # output is the dfs with first case & death conditionally filled with zero. 
              # Later, I name this output to be "dfs_first_zero_filled" 

dfs = extract_cntry_dfs(df)

dfs_first_zero_filled = fill_first_case_death_with_zero(dfs)

def merge_with_meta(df): #input should be dfs_first_zero_filled
    concat_df = pd.concat(df).fillna(method="ffill").reset_index().rename(
        columns={"index": "date"})
    #To change the original country codes of "KOS" and World to match metadata from WB:
    concat_df.loc[(concat_df.country_code == "OWID_KOS"), ('country_code')] = "XKX"
    concat_df.loc[(concat_df.country_code == "OWID_WRL"), ('country_code')] = "WLD"
    # To get the column of "world_region" in concat_df by merging with WB metadata
    left_join_df = pd.merge(concat_df, metadata, on = "country_code", how = "left")
    left_join_df.loc[:,'date'] = left_join_df.loc[:,'date'].dt.strftime('%Y-%m-%d')
    left_join_df.loc[(left_join_df.country_code == "AIA"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "BES"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "ESH"), ('world_region')] = "Middle East & North Africa"
    left_join_df.loc[(left_join_df.country_code == "FLK"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "GGY"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "JEY"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "MSR"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "TWN"), ('world_region')] = "East Asia & Pacific"
    left_join_df.loc[(left_join_df.country_code == "VAT"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "WLF"), ('world_region')] = "East Asia & Pacific"
    left_join_df.loc[(left_join_df.country_code == "WLD"), ('world_region')] = "World"
    return left_join_df

left_join_df = merge_with_meta(dfs_first_zero_filled)

def fallBehind_zero_to_nan (df): # input should be left_join_df
    left_join_copy_group1 = []
    for group in df.groupby('country_code'):
        for i in np.arange(0, len(fallBehind_list)):
            if group[1].tail(1).country_code.iloc[0] == fallBehind_list.iloc[i, 0]:
                group[1].tail(len(all_dates) - 1 - fallBehind_list.iloc[i, 2]).total_cases = np.nan
                group[1].tail(len(all_dates) - 1 - fallBehind_list.iloc[i, 2]).total_deaths = np.nan
        left_join_copy_group1.append(group[1])
    return left_join_copy_group1 # I will name the output later to be fallBehind_nan_changed

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
