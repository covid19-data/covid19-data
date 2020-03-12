# 2019-nCoV Data Processing Pipelines

[![Join the chat at https://gitter.im/covid19-data/community](https://badges.gitter.im/covid19-data/community.svg)](https://gitter.im/covid19-data/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This repository hosts workflows to process several data sources. 

## Data sources

- WHO
  - Processed by Our World in Data: https://ourworldindata.org/coronavirus-source-data
- [COVID-19 daily report by JHU](https://github.com/CSSEGISandData/COVID-19): Note that this data source has many consistency issues regarding country names and aggregation of US data. 
  - Tableau: Tableau cleans the JHU CSSE dataset and provides a tidy-formatted dataset. However, as of now, it does not 
 address the data consistency issues in the raw dataset. 
- [Worldbank country population data and country metadata](https://data.worldbank.org/indicator/SP.POP.TOTL)

# Outputs

- [`cntry_stat.json`](https://github.com/yy/covid19-data/blob/master/output/cntry_stat.json): this is used in [an interactive visualization of case fatality rate of COVID-19](http://yyahn.com/covid19) ([source code on ObservableHQ](https://observablehq.com/@yy/covid-19-fatality-rate)).


# Usage

After cloning the repository, pull the submodule (JHU CSSE dataset). 

```sh
git submodule init
git submodule
```

Install [pandas](https://pandas.pydata.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/) using `conda` 

```sh
conda install -c bioconda -c conda-forge snakemake pandas
```

or `pip`:

```sh
pip install pandas snakemake
```

Run `snakemake`.

```sh
snakemake
```
