# covid19-data

A data workflow to combine [COVID-19 daily report by JHU](https://github.com/CSSEGISandData/COVID-19) with [Worldbank country metadata](https://data.worldbank.org/indicator/SP.POP.TOTL) to produce [`cntry_stat.json`](https://github.com/yy/covid19-data/blob/master/cntry_stat.json) this is used in an interactive visualization of case fatality rate of COVID-19.

- Visualization: http://yyahn.com/covid19
- Visualization code: https://observablehq.com/@yy/covid-19-fatality-rate

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

Run `make`.

```sh
make
```
