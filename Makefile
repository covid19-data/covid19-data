.PHONY: all watch

all:
	cd data_sources/COVID-19; git pull origin master
	cd data_sources/our_world_in_data; ./download.sh
	cd data_sources/tableau; ./download.sh
	snakemake -R `snakemake --lc --li --lp`

