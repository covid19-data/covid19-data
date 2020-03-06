from os.path import join as j

# Time series raw data
TIME_SERIES_PATH = "COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
TIME_SERIES_FNAME = j(TIME_SERIES_PATH, "time_series_19-covid-{data_type}.csv")
DATA_TYPES = ["Confirmed", "Deaths", "Recovered"]
CONFIRMED_TS, DEATHS_TS, RECOVERED_TS = expand(TIME_SERIES_FNAME, data_type=DATA_TYPES)

COMBINED_TS = 'combined_ts.csv'
COORDINATES = 'coordinates.csv'

# Population data from Worldbank
POP_DATA_DIR = "worldbank_population_data"
POP_DATA = j(POP_DATA_DIR, 'API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
CNTRY_META = j(POP_DATA_DIR, 'Metadata_Country_API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')

POP_CSV = 'pop.csv'


CNTRY_STAT_JSON = 'cntry_stat.json'

rule all:
    input: POP_CSV, COMBINED_TS, COORDINATES

rule extract_coordinates:
    input: CONFIRMED_TS
    output: COORDINATES
    script: "extract_coordinates.py"

rule combine_ts:
    input: CONFIRMED_TS, DEATHS_TS, RECOVERED_TS
    output: COMBINED_TS
    script: "combine_ts.py"

rule world_pop_data:
    input: POP_DATA, CNTRY_META
    output: POP_CSV
    script: "extract_pop_data.py"
