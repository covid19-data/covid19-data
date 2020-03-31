# COVID-19 Data Processing Pipelines and datasets

[![Join the chat at https://gitter.im/covid19-data/community](https://badges.gitter.im/covid19-data/community.svg)](https://gitter.im/covid19-data/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This repository hosts workflows to process several data sources and cleaned datasets for COVID-19 cases across the world.

## Datasets 

### Historical (daily) case data 

- [`output/cases/cases_ECDC.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_ECDC.csv): European Centre for Disease Prevention and Control (ECDC) historical world-wide case data (currently through [Our World in Data](https://ourworldindata.org/coronavirus-source-data)).

### Country metadata 

- Country metadata: [`output/metadata/country/country_metadata.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/metadata/country/country_metadata.csv) from [Worldbank](https://data.worldbank.org/indicator/SP.POP.TOTL).

- ISO 3166-1 Alpha-3 country code conversion table. 
    - [`output/metadata/country/country_name_code.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_name_code.csv): a conversion table from country name to code (ISO 3166 Alpha 3). Note that multiple names point to the same code.
    - [`output/metadata/country/country_code_name.csv`](https://github.com/yy/covid19-data/blob/master/output/metadata/country/country_code_name.csv): a conversion table from country code (ISO 3166 Alpha 3) to country name. The shortest country names are picked from the above dataset.


### Historical case data for visualizations

- [`cntry_stat_owid.json`](https://github.com/yy/covid19-data/blob/master/output/cntry_stat_owid.json): ECDC historical data merged with Worldbank's country metadata and ISO 3166-1 Alpha-3 data. 
  - Used in [an interactive visualization of case fatality rate of COVID-19](http://yyahn.com/covid19)
    - Website source code: https://github.com/covid19-data/covid19-dashboard
    - visualization source code on ObservableHQ: https://observablehq.com/@yy/covid-19-fatality-rate and https://observablehq.com/@yy/covid-19-trends
  - [An example to create case time series charts in ObservableHQ](https://observablehq.com/@benjyz/covid-chart-alpha) by [benjyz](https://github.com/benjyz)


### Deprecated

WHO dataset is deprecated. See Our World in Data's announcement: [Why we stopped relying on data from the World Health Organization](https://ourworldindata.org/coronavirus#why-we-stopped-relying-on-data-from-the-world-health-organization)

- [`output/cases/cases_WHO.csv`](https://github.com/covid19-data/covid19-data/blob/master/output/cases/cases_WHO.csv)

- https://www.worldometers.info/coronavirus/

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

