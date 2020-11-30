import pandas as pd
import numpy as np
import json
import os
import http.client as http
http.HTTPConnection._http_vsn = 10
http.HTTPConnection._http_vsn_str = 'HTTP/1.0'


# prepare owid raw data:
df = pd.read_csv('https://covid.ourworldindata.org/data/owid-covid-data.csv',
    index_col = ['date'], parse_dates = True,
    usecols=["date", "iso_code", "location", "population", 
                      "continent", "total_cases", "total_deaths"]).rename(
    columns = {'iso_code': 'country_code', 'location': 'country_name'}
)

# prepare metadata from World Bank
metadata = pd.read_csv("https://raw.githubusercontent.com/hongtaoh/covid19-data/master/data_sources/metadata/worldbank/country_metadata.csv",
                      usecols=["Country Code", "Region"]).rename(
    columns = {'Country Code': 'country_code', 'Region': 'world_region'}
)

# get the full date range:
all_dates = pd.date_range(df.index.min(), df.index.max())

# The following codes are to find out which places do not have the most up-to-date data. 
# All 'nan' will be filled by forward filling in steps below. 
# So I need change the cells that shouldn't have been filled to 'nan' before exporting data. 
# To do that, I need to find out which places fall behind and which dates sholud be changed to nan. 

# I changed theose cells to 'null' instead of np.nan, so that Observablehq can read the final output.
# That's something I'll explain later. 

def get_fallBehind_place_date_dateframe(df): # input should be df
    fallBehind_place_list = []
    fallBehind_last_date_available_list = []
    fallBehind_last_date_index_list = []
    for group in df.groupby("country_code"):
        if group[1].tail(1).index != df.index.max():
            fallBehind_place_list.append(group[0])
    for p in fallBehind_place_list:
        fallBehind_last_date_index_list.append(all_dates.get_loc(
            df[df.loc[:, 'country_code']== p].index[-1].strftime("%Y-%m-%d")))
        fallBehind_last_date_available_list.append(
            df[df.loc[:, 'country_code']== p].index[-1].strftime("%Y-%m-%d"))
    d = {'fallBehind_place': fallBehind_place_list, 
         'fallBehind_last_date_available': fallBehind_last_date_available_list,
         'fallBehind_last_date_index': fallBehind_last_date_index_list}
    fallBehind_list = pd.DataFrame(data = d)
    return fallBehind_list
# output is called fallBehind_list

fallBehind_list = get_fallBehind_place_date_dateframe(df)

def extract_cntry_dfs(df): # input is df
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
# output is literally dfs

dfs = extract_cntry_dfs(df)

# What I intend to do in the following part of the script is to change the data for cases and deaths
# on Dec. 31st, 2019 for each country / area to be zero. Otherwise, I cannot use forward filling 
# to replace nan later. 

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


# In the following step, I converted to nan the data that shouldn't have been filled with forward filling for places 
# that do not have the updated data. Then I converted these NaNs to null because NaN is not recognized by JSON. 

def fallBehind_filled_to_null (df): # input should be left_join_df
    left_join_copy_group1_with_nan = []
    for group in df.groupby('country_code'):
        for i in np.arange(0, len(fallBehind_list)):
            if group[1].tail(1).country_code.iloc[0] == fallBehind_list.iloc[i, 0]:
                group[1].tail(len(all_dates) - 1 - fallBehind_list.iloc[i, 2]).total_cases = np.nan
                group[1].tail(len(all_dates) - 1 - fallBehind_list.iloc[i, 2]).total_deaths = np.nan
        left_join_copy_group1_with_nan.append(group[1])
    left_join_copy_group1_with_nan_concated = pd.concat(left_join_copy_group1_with_nan)
    left_join_copy_group1_concated_with_null = left_join_copy_group1_with_nan_concated
    left_join_copy_group1_concated_with_null.replace(np.nan, 'null', inplace=True)
    return left_join_copy_group1_concated_with_null 
# Later, I'll name the output to be fallBehind_with_null

fallBehind_with_null = fallBehind_filled_to_null(left_join_df)

def prepare_data_structure(df, gby="country_code"): # input should be fallBehind_with_null
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

data = prepare_data_structure(fallBehind_with_null)

os.chdir("../output")

fallBehind_with_null.to_csv("race_chart_data.csv", index=False)

open("cntry_stat_owid.json", "w").write(json.dumps(data, separators=(",", ":")))
