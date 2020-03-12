from os.path import join as j

# Time series raw data
# TIME_SERIES_PATH = "COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
# TIME_SERIES_FNAME = j(TIME_SERIES_PATH, "time_series_19-covid-{data_type}.csv")
# DATA_TYPES = ["Confirmed", "Deaths", "Recovered"]
# CONFIRMED_TS, DEATHS_TS, RECOVERED_TS = expand(TIME_SERIES_FNAME, data_type=DATA_TYPES)
# COMBINED_TS = 'combined_ts.csv'
COORDINATES = 'coordinates.csv'

TABLEAU_TS = 'tableau_ts.csv'

# Data file from Our World in Data
OWID_TS = 'data_sources/our_world_in_data/owid_ts.csv'

# Population data from Worldbank
WB_DATA_DIR = "data_sources/worldbank_population_data"
WB_POP_RAW = j(WB_DATA_DIR, 'API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
WB_META = j(WB_DATA_DIR, 'Metadata_Country_API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
WB_RAW = j(WB_DATA_DIR, 'wb_raw.csv')

# Population data cleaning and matching with JHU data
POP_ADDITION = j(WB_DATA_DIR, 'pop_addition.csv')
POP_CONVERSION = 'pop_conversion.csv'
POP_CLEANED_CSV = 'pop.csv'


# Final json file for the visualization
CNTRY_STAT_JSON = 'cntry_stat.json'

include: "data_sources/tableau/Snakefile"
include: "data_sources/our_world_in_data/Snakefile"
include: "data_sources/worldbank_population_data/Snakefile"

rule all:
    input: CNTRY_STAT_JSON #, COORDINATES

rule extract_coordinates:
    input: TABLEAU_TS
    output: COORDINATES
    script: "scripts/extract_coordinates.py"

# rule combine_ts:
    # input: CONFIRMED_TS, DEATHS_TS, RECOVERED_TS
    # output: COMBINED_TS
    # script: "scripts/combine_ts.py"

# rule merge_ts_and_pop_data:
    # input: COMBINED, POP_CLEANED_CSV
    # output: CNTRY_STAT_JSON
    # script: "scripts/merge_ts_and_pop_data.py"
