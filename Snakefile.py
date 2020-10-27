from os.path import join as j

import pandas as pd

###############################################################################
# Raw datasets
###############################################################################

# NYT data
NYT_STATE_TS_RAW = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'

###############################################################################
# Final outputs
###############################################################################

# USA state data from NYT
NYT_STATE_TS = 'output/cases/cases_us_states_nyt.csv'

## JSON file for visualizations

USA_STAT_JSON_FROM_NYT = 'output/us_state_nyt.json'

###############################################################################
# Workflows
###############################################################################

rule nyt_state_ts:
    input: NYT_STATE_TS_RAW
    output: NYT_STATE_TS
    shell: "python scripts/nyt_state_data.py {input} {output}"

rule prepare_us_states_viz_data:
    input: NYT_STATE_TS
    output: USA_STAT_JSON_FROM_NYT
    shell: "python scripts/prepare_state_level_viz_data.py {input} {output}"