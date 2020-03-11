# covid19-data

A data workflow to combine [COVID-19 daily report by JHU](https://github.com/CSSEGISandData/COVID-19) with [Worldbank country metadata](https://data.worldbank.org/indicator/SP.POP.TOTL), 
for an interactive visualization of case fatality rate of COVID-19.

- Visualization: http://yyahn.com/covid19
- Visualization code: https://observablehq.com/@yy/covid-19-fatality-rate

# Usage

After cloning the repository, pull the submodule (JHU CSSE dataset). 

```sh
git submodule init
git submodule
```

Install [pandas](https://pandas.pydata.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/) using `conda` or `pip`.

```sh
conda install -c bioconda -c conda-forge snakemake pandas
```

```sh
pip install pandas snakemake
```

Run 

```sh
make
```
