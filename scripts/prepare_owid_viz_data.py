import pandas as pd
import numpy as np
import json
import math

# prepare owid raw data:
df = pd.read_csv('https://covid.ourworldindata.org/data/owid-covid-data.csv',
                 index_col = ['date'], 
                 parse_dates = True,
                 usecols=["date", "iso_code", "location", 
                          "population", "continent", 
                          "total_cases", "total_deaths"]).rename(
    columns = {'iso_code': 'country_code', 'location': 'country_name'}
)

# prepare metadata from World Bank
metadata = pd.read_csv("https://raw.githubusercontent.com/hongtaoh/covid19-data/master/data_sources/metadata/worldbank/country_metadata.csv",
                      usecols=["Country Code", "Region"]).rename(
    columns = {'Country Code': 'country_code', 'Region': 'world_region'}
)

# get the full date range:
all_dates = pd.date_range(df.index.min(), df.index.max())
# all_dates = pd.date_range('22/01/2020', df.index.max())

def extract_cntry_dfs(df): # input is df
    dfs = [] # Initiating dfs, which is a list
    for group in df.groupby("country_code"):
        # Normalize countries with all_dates. Using forward filling for rows of newly 
        # produced dates:
        cntry_df = (
            group[1].reindex(all_dates, method="ffill")
        )
        # Forward filling makes sense for total_cases & total_deaths, but not for
        # all other variables that are constant. 
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

# I used forward filling in `def extract_cntry_dfs`. The problem is that most countries/regions do not have data for the first date, i.e., January 1st, 2020. This means the total_cases and total_deaths on 2021-01-01 for most countries/regions are `nan`. Since it's forward filling, data for following dates that do not have available data will be fillled with `nan`. 
#
# To better visualize the trend, we'd better change unavailable data in the beginning to be ZERO. (**I din't fully understand why. Couldn't just leave them as `nan`? Why do we have to change them to zero?**)

# +
# def fill_first_case_death_with_zero(df): # input should be dfs, NOT df!
#     for i in np.arange(0, len(df)):
# # dfs[i] will be a dataframe for each country/region
# # dfs[i].head(1).total_cases.isnull()[0] will return either True or False

# # There are three if statements. The 1st one: if both the first total_cases and the first total_deaths are null
# # The 2nd one: if the first total_cases is null and the first total_deaths is NOT null
# # The 3rd one: if the first total_cases is NOT null and the first total_deaths is null
#         if (df[i].head(1).total_cases.isnull()[0] & df[i].head(1).total_deaths.isnull()[0]):
#         # dfs[i].iloc[0:0] is the variable names: country_code, continent, etc.
#         # dfs[i].iloc[0:1] is the variable names + the first row. 
#         # Try `dfs[1].iloc[0:1] = [1,'a','c',2,4,'t']` and you'll see that this changes the data of first row 
#         #   but not the variable names. 
#             df[i].iloc[0:1] = [df[i].country_code[-1],
#             df[i].continent[-1],
#             df[i].country_name[-1],
#             0,
#             0, 
#             df[i].population[-1]
#             ]
#         if (df[i].head(1).total_cases.isnull()[0] & df[i].head(1).total_deaths.notnull()[0]):
#             df[i].iloc[0:1] = [df[i].country_code[-1],
#             df[i].continent[-1],
#             df[i].country_name[-1],
#             0,
#             df[i].total_deaths[0],
#             df[i].population[-1]
#             ]
#         if (df[i].head(1).total_cases.notnull()[0] & df[i].head(1).total_deaths.isnull()[0]):
#             df[i].iloc[0:1] = [df[i].country_code[-1],
#             df[i].continent[-1],
#             df[i].country_name[-1],
#             df[i].total_cases[0],
#             0,
#             df[i].population[-1]
#             ]
#     return df # output is the dfs with first case & death conditionally filled with zero. 
#               # Later, I name this output to be "dfs_first_zero_filled" 
# -

# dfs_first_zero_filled = fill_first_case_death_with_zero(dfs)
dfs_first_zero_filled = dfs

def merge_with_meta(df): #input should be dfs_first_zero_filled
    # pd.concat() will stack up each group's data and produce a DataFrame containing all data
    # reset_index() will reset the index and use the default one. Basically,
    #  reset_index() will treat the "date" as a column and index each row using integers starting from 0
    concat_df = pd.concat(df).fillna(method="ffill").reset_index().rename(
        columns={"index": "date"})
    #To change the original country codes of "KOS" and World to match metadata from WB:
    concat_df.loc[(concat_df.country_code == "OWID_KOS"), ('country_code')] = "XKX"
    concat_df.loc[(concat_df.country_code == "OWID_WRL"), ('country_code')] = "WLD"
    # To get the column of "world_region" in concat_df by merging with WB metadata
    left_join_df = pd.merge(concat_df, metadata, on = "country_code", how = "left")
    # to convert datetime to strings in the form of 2001-01-01
    left_join_df.loc[:,'date'] = left_join_df.loc[:,'date'].dt.strftime('%Y-%m-%d')
    # Some countries/regions are in the covid19 data but not in the metadata from WB, so I need to manually change 
    #  these countries/regions' world_region value:
    left_join_df.loc[(left_join_df.country_code == "AIA"), ('world_region')] = "Latin America & Caribbean"
#     left_join_df.loc[(left_join_df.country_code == "BES"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "ESH"), ('world_region')] = "Middle East & North Africa"
#     left_join_df.loc[(left_join_df.country_code == "FLK"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "GGY"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "JEY"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "SHN"), ('world_region')] = "Sub-Saharan Africa"
#     left_join_df.loc[(left_join_df.country_code == "MSR"), ('world_region')] = "Latin America & Caribbean"
    left_join_df.loc[(left_join_df.country_code == "TWN"), ('world_region')] = "East Asia & Pacific"
    left_join_df.loc[(left_join_df.country_code == "VAT"), ('world_region')] = "Europe & Central Asia"
    left_join_df.loc[(left_join_df.country_code == "FLK"), ('world_region')] = "Latin America & Caribbean"
#     left_join_df.loc[(left_join_df.country_code == "WLF"), ('world_region')] = "East Asia & Pacific"
    left_join_df.loc[(left_join_df.country_code == "WLD"), ('world_region')] = "World"
    return left_join_df

left_join_df = merge_with_meta(dfs_first_zero_filled)


# Forward filling works for most of the missing data. For example, if there is no data for Ghana from March 3rd, 2020 to March 17th, 2020, it's safe and ideal to fill these date's data with that of March 2nd, 2020. However, if, say, Ghana's most updated data is on June 2nd, 2021 whereas the latest date in `all_dates` is June 22nd, 2021, then we shouldn't fill data for 06/03/2021 - 06/22/2021 with data for 06/02/2021. 
#
# The following codes are to find out which places do not have the most up-to-date data. 
# All 'nan' will be filled by forward filling in steps below. 
# So I need to change the cells that shouldn't have been filled to 'nan' before exporting data. 
# To do that, I need to find out which places fall behind and which dates sholud be changed to nan. 

def get_fallBehind_place_date_dateframe(df): # input should be df
    fallBehind_place_list = [] #initiating
    fallBehind_last_date_available_list = [] #initiating
    fallBehind_last_date_index_list = [] #initiating
    for group in df.groupby("country_code"): 
# What the following lines of codes does is to find out which group does not have the 
# the most updated data, then add this group's country_code to FallBehind_place_list
        if group[1].tail(1).index != df.index.max():
            fallBehind_place_list.append(group[0])
    for p in fallBehind_place_list:
# For each country/region in fallBehind_place_list, find out what the most updated date is and 
# add this date and it's index in all_dates
        fallBehind_last_date_available_list.append(
            df[df.loc[:, 'country_code']== p].index[-1].strftime("%Y-%m-%d"))
        fallBehind_last_date_index_list.append(all_dates.get_loc(
            df[df.loc[:, 'country_code']== p].index[-1].strftime("%Y-%m-%d")))
# Create a dictionary and then create a DataFrame from this dictionary
    d = {'fallBehind_place': fallBehind_place_list, 
         'fallBehind_last_date_available': fallBehind_last_date_available_list,
         'fallBehind_last_date_index': fallBehind_last_date_index_list}
    fallBehind_list = pd.DataFrame(data = d)
    return fallBehind_list
# output is called fallBehind_list

fallBehind_list = get_fallBehind_place_date_dateframe(df)

# In the following step, for places that do not have the most updated data, I converted data that shouldn't have been filled with forward filling to `nan` (by `np.nan`). Then I converted these NaNs to null because NaN is not recognized by JSON. 

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

# os.chdir("../output")

fallBehind_with_null.to_csv("/Users/Tal/Desktop/covid19-data/output/race_chart_data.csv", index=False)

open("/Users/Tal/Desktop/covid19-data/output/cntry_stat_owid.json", "w").write(json.dumps(data, separators=(",", ":")))
