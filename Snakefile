from os.path import join as j

###############################################################################
# Raw datasets
###############################################################################

# Data file from Our World in Data. Directly from WHO
OWID_TS = 'data_sources/our_world_in_data/owid_ts.csv'

# Data feed from Tableau. Raw data is from JHU CSSE GitHub. currently not used.
TABLEAU_TS = 'data_sources/tableau/tableau_ts.csv'

# Population data from Worldbank
WB_RAW = "data_sources/worldbank_population_data/wb_raw.csv"

# country code data from wikipedia
WP_CNTRY_RAW =  'data_sources/wikipedia/ISO3166_country_code/iso3166_country_code.csv'

# Population data cleaning and matching with JHU data
# POP_CONVERSION = j(WB_DATA_DIR, 'pop_conversion.csv')
# POP_CLEANED_CSV = j(WB_DATA_DIR, 'pop.csv')


###############################################################################
# Manually curated datasets
###############################################################################

WB_ADDED = 'curation_data/country/extra_country_metadata.csv'
EXTRA_CNTRY_NAME_CODE = 'curation_data/country/extra_country_name_code.csv'

###############################################################################
# Final outputs
###############################################################################

OUT_META_DIR = 'output/metadata'

# Location data pulled from tableau/jhu dataset.
COORDINATES = 'output/location/coordinates.csv'

# Country name and code conversion table
CNTRY_NAME_CODE_TABLE = j(OUT_META_DIR, 'country/country_name_code.csv')
CNTRY_CODE_NAME_TABLE = j(OUT_META_DIR, 'country/country_code_name.csv')

# Country-level metadata
CNTRY_META = 'output/metadata/country/country_metadata.csv'

# json file for the visualization: http://yyahn.com/covid19
CNTRY_STAT_JSON = 'output/cntry_stat.json'
CNTRY_STAT_JSON_FROM_OWID = 'output/cntry_stat_owid.json'


###############################################################################
# Workflows
###############################################################################

rule all:
    input: CNTRY_STAT_JSON_FROM_OWID, COORDINATES, CNTRY_META, CNTRY_CODE_NAME_TABLE

rule extract_country_code_name_table:
    input: CNTRY_NAME_CODE_TABLE
    output: CNTRY_CODE_NAME_TABLE
    script: "scripts/extract_country_code_name_table.py"

rule extract_country_name_code_table:
    input: WP_CNTRY_RAW, WB_RAW, WB_ADDED, EXTRA_CNTRY_NAME_CODE
    output: CNTRY_NAME_CODE_TABLE
    script: "scripts/extract_country_name_code_table.py"

rule extract_coordinates:
    input: TABLEAU_TS
    output: COORDINATES
    script: "scripts/extract_coordinates.py"

rule prepare_viz_data:
    input: OWID_TS, WB_RAW
    output: CNTRY_STAT_JSON_FROM_OWID
    script: "scripts/prepare_viz_data_owid.py"

rule extract_country_metadata:
    input: WB_RAW, WB_ADDED
    output: CNTRY_META
    script: "scripts/extract_country_metadata.py"
