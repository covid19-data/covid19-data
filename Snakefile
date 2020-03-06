from os.path import join as j

# Time series raw data
TIME_SERIES_PATH = "COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
TIME_SERIES_FNAME = j(TIME_SERIES_PATH, "time_series_19-covid-{data_type}.csv")
DATA_TYPES = ["Confirmed", "Deaths", "Recovered"]
CONFIRMED_TS, DEATHS_TS, RECOVERED_TS = expand(TIME_SERIES_FNAME, data_type=DATA_TYPES)

COMBINED_TS = 'combined_ts.csv'
COORDINATES = 'coordinates.csv'

# Population data from Worldbank
POP_RAW_DATA_DIR = "worldbank_population_data"
POP_RAW_DATA = j(POP_RAW_DATA_DIR, 'API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
CNTRY_META = j(POP_RAW_DATA_DIR, 'Metadata_Country_API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')

POP_CSV = 'pop_raw.csv'

# Population data cleaning and matching with JHU data
POP_ADDITION = 'pop_addition.csv'
POP_CONVERSION = 'pop_conversion.csv'
POP_CLEANED_CSV = 'pop.csv'


# Final json file for the visualization
CNTRY_STAT_JSON = 'cntry_stat.json'


rule all:
    input: CNTRY_STAT_JSON, COORDINATES

rule extract_coordinates:
    input: CONFIRMED_TS
    output: COORDINATES
    script: "scripts/extract_coordinates.py"

rule combine_ts:
    input: CONFIRMED_TS, DEATHS_TS, RECOVERED_TS
    output: COMBINED_TS
    script: "scripts/combine_ts.py"

rule clean_pop_data:
    input: POP_CSV, POP_ADDITION, POP_CONVERSION
    output: POP_CLEANED_CSV
    script: "scripts/clean_pop_data.py"

rule world_pop_data:
    input: POP_RAW_DATA, CNTRY_META
    output: POP_CSV
    script: "scripts/extract_pop_data.py"

rule merge_ts_and_pop_data:
    input: COMBINED_TS, POP_CLEANED_CSV
    output: CNTRY_STAT_JSON
    script: "scripts/merge_ts_and_pop_data.py"
