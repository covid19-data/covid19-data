# COVID-19 Data Processing Pipelines and datasets

[![Join the chat at https://gitter.im/covid19-data/community](https://badges.gitter.im/covid19-data/community.svg)](https://gitter.im/covid19-data/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This repository hosts workflows to process several data sources and cleaned datasets for COVID-19 cases across the world.

## Datasets 

### Historical (daily) case data 

- [`owid-covid-data.json`](https://covid.ourworldindata.org/data/owid-covid-data.json): European Centre for Disease Prevention and Control (ECDC) historical world-wide case data (currently through [Our World in Data](https://ourworldindata.org/coronavirus-source-data)).

- [`output/cases/cases_us_states_nyt.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_us_states_nyt.csv): US state-level historical case data from [New York Times](https://github.com/nytimes/covid-19-data).

### Country metadata 

- World Bank Data: Classification of **world region** (Latin America & Caribben, South Asia, Sub-Saharan African, Europea & Central Asia, Middle East & North Africa, East Asia & Pacific, North America) for each country or area is based on [`data_source/metadata/worldbank/country_metadata.csv`](https://github.com/covid19-data/covid19-data/blob/master/data_source/metadata/worldbank/country_metadata.csv) from [World Bank](https://data.worldbank.org/indicator/SP.POP.TOTL).
- ISO-3166-Countries-with-Regional-Codes: For countries or areas not found in World Bank data, their **world region** is found in [`ISO-3166-Countries-with-Regional-Codes`](https://github.com/hongtaoh/covid19-data/blob/master/data_source/metadata/ISO-3166-Countries-with-Regional-Codes/all/all.csv).


### Historical case data for visualizations

- [`cntry_stat_owid.json`](https://github.com/honghaoh/covid19-data/blob/master/output/cntry_stat_owid.json): ECDC historical data merged with Worldbank's classification of world regions. Used in:
  - [an interactive visualization of case fatality rate of COVID-19](http://yyahn.com/covid19)
    - Website source code: https://github.com/covid19-data/covid19-dashboard
    - visualization source code on ObservableHQ: https://observablehq.com/@yy/covid-19-fatality-rate and https://observablehq.com/@yy/covid-19-spreading-trends
  - [An example to create case time series charts in ObservableHQ](https://observablehq.com/@benjyz/covid-chart-alpha) by [benjyz](https://github.com/benjyz)
  - You can find the data manipulation process of [`cntry_stat_owid.json`](https://github.com/hongtaoh/covid19-data/blob/master/output/cntry_stat_owid.json) [here](https://observablehq.com/@hongtaoh/day-46-2020-10-08).
- [`us_state_nyt.json`](https://github.com/hongtaoh/covid19-data/blob/master/output/us_state_nyt.json): New York Time historical data. Used in:
    - [New cases vs. all confirmed interactive visualization](https://observablehq.com/@yy/covid-19-confirmed-vs-new-cases).

### Deprecated

WHO dataset is deprecated. See Our World in Data's announcement: [Why we stopped relying on data from the World Health Organization](https://ourworldindata.org/coronavirus#why-we-stopped-relying-on-data-from-the-world-health-organization)

- [`output/cases/cases_WHO.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_WHO.csv)
- https://www.worldometers.info/coronavirus/
- [`coordinates.csv`](https://github.com/yy/covid19-data/blob/master/output/location/coordinates.csv): Lat Lng location data from JHU dataset (Unreliable).

- ISO 3166-1 Alpha-3 country code conversion table. 
    - [`output/metadata/country/country_name_code.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_name_code.csv): a conversion table from country name to code (ISO 3166 Alpha 3). Note that multiple names point to the same code.
    - [`output/metadata/country/country_code_name.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_code_name.csv): a conversion table from country code (ISO 3166 Alpha 3) to country name. The shortest country names are picked from the above dataset.

## Usage

Install [pandas](https://pandas.pydata.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/) using `conda`.

```sh
conda install -c bioconda -c conda-forge snakemake pandas numpy
```

or `pip`:

```sh
pip install pandas snakemake numpy
```

