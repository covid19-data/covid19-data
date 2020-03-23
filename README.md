# 2019-nCoV Data Processing Pipelines and datasets

We are looking for volunteers who want to contribute to the cleaning of the raw datasets!

[![Join the chat at https://gitter.im/covid19-data/community](https://badges.gitter.im/covid19-data/community.svg)](https://gitter.im/covid19-data/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This repository hosts workflows to process several data sources and cleaned datasets for COVID-19 cases across the world.

## Data sources

### Currently used

- European Centre for Disease Prevention and Control, ECDC
  - Processed by Our World in Data: https://ourworldindata.org/coronavirus-source-data

- previously WHO
  - Processed by Our World in Data: https://ourworldindata.org/coronavirus-source-data
  - Tableau: Tableau cleans the JHU CSSE dataset and provides a tidy-formatted dataset. However, as of now, it does not
 address the data consistency issues in the raw dataset.
- [Worldbank country population data and country metadata](https://data.worldbank.org/indicator/SP.POP.TOTL)
- Wikipedia ISO3166 Country code data
- https://www.worldometers.info/coronavirus/

### Not used at the moment

- [COVID-19 daily report by JHU](https://github.com/CSSEGISandData/COVID-19): This has many consistency issues regarding country names and aggregation of US data. Aggregation mechanism is not so transparent.

## Data Outputs and Usages

### Daily case data based on WHO report, cleaned by Our World in Data

- [`output/cases/cases_WHO.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_WHO.csv): This converts the CSV dataset cleaned by Our World in Data team by using ISO 3166 Alpha-3 country code. It also fills up non-existing dates so that for every country, the dataset starts from the same date (Jan. 21st). One may want to combine this with country-level metadata or alternative country names [here](https://github.com/covid19-data/covid19-data/tree/master/output/metadata/country).

### Daily case data (WHO + US case data from Wikipedia)

- [`output/cases/cases_WHO_WP.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_WHO_WP.csv): similar to `cases_WHO.csv`, but US data is overrided by the data from Wikipedia.

### Daily case data in JSON

- [`cntry_stat_owid.json`](https://github.com/yy/covid19-data/blob/master/output/cntry_stat_owid.json)
  - Used in [an interactive visualization of case fatality rate of COVID-19](http://yyahn.com/covid19)
    - Website source code: https://github.com/covid19-data/covid19-dashboard
    - visualization source code on ObservableHQ: https://observablehq.com/@yy/covid-19-fatality-rate and https://observablehq.com/@yy/covid-19-trends
  - [An example to create case time series charts in ObservableHQ](https://observablehq.com/@benjyz/covid-chart-alpha) by [benjyz](https://github.com/benjyz)


### Country name conversion table

- [`output/metadata/country/country_name_code.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_name_code.csv): a conversion table from country name to code (ISO 3166 Alpha 3). Note that multiple names point to the same code.
- [`output/metadata/country/country_code_name.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_code_name.csv): a conversion table from country code (ISO 3166 Alpha 3) to country name. The shortest country names are picked from the above dataset.

### Country metadata

- [`output/metadata/country/country_metadata.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_metadata.csv): Country metadata, such as population, region, and income group, indexed by the ISO 3166 Alpha 3 codes.

### Location data

- [`coordinates.csv`](https://github.com/yy/covid19-data/blob/master/output/location/coordinates.csv): Lat Lng location data from JHU dataset (Unreliable).


## Usage

Install [pandas](https://pandas.pydata.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/) using `conda`.

```sh
conda install -c bioconda -c conda-forge snakemake pandas numpy
```

or `pip`:

```sh
pip install pandas snakemake numpy
```

