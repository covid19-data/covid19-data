from os.path import join as j

###############################################################################
# Raw datasets
###############################################################################

# Data file from Our World in Data. Directly from WHO
OWID_TS = 'data_sources/our_world_in_data/owid_ts.csv'

# Data feed from Tableau. Raw data is from JHU CSSE GitHub. currently not used.
TABLEAU_TS = 'data_sources/tableau/tableau_ts.csv'

# Population data from Worldbank
WB_DATA_DIR = "data_sources/worldbank_population_data"
WB_POP_RAW = j(WB_DATA_DIR, 'API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
WB_META = j(WB_DATA_DIR, 'Metadata_Country_API_SP.POP.TOTL_DS2_en_csv_v2_821007.csv')
WB_RAW = j(WB_DATA_DIR, 'wb_raw.csv')

# Population data cleaning and matching with JHU data
# POP_ADDITION = j(WB_DATA_DIR, 'pop_addition.csv')
# POP_CONVERSION = j(WB_DATA_DIR, 'pop_conversion.csv')
# POP_CLEANED_CSV = j(WB_DATA_DIR, 'pop.csv')

###############################################################################
# Final outputs
###############################################################################

# Location data pulled from tableau/jhu dataset.
COORDINATES = 'output/location/coordinates.csv'

# json file for the visualization: http://yyahn.com/covid19
CNTRY_STAT_JSON = 'output/cntry_stat.json'
CNTRY_STAT_JSON_FROM_OWID = 'output/cntry_stat_owid.json'


###############################################################################
# Workflows
###############################################################################

rule all:
    input: CNTRY_STAT_JSON_FROM_OWID, COORDINATES

rule extract_coordinates:
    input: TABLEAU_TS
    output: COORDINATES
    script: "scripts/extract_coordinates.py"

rule prepare_viz_data:
    input: OWID_TS, WB_RAW
    output: CNTRY_STAT_JSON_FROM_OWID
    script: "scripts/prepare_viz_data_owid.py"
