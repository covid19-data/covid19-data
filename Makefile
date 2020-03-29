.PHONY: all watch

all:
	snakemake -R `snakemake --lc --li --lp`

owid:
	cd data_sources/our_world_in_data; ./download.sh

nyt:
	cd data_sources/nyt/covid-19-data; git pull origin master

wiki:
	cd data_sources/wikipedia/cases; snakemake -R download_and_parse_countries download_and_parse_states

tableau:
	cd data_sources/tableau; ./download.sh
