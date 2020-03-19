.PHONY: all watch

all:
	snakemake -R `snakemake --lc --li --lp`

owid:
	cd data_sources/our_world_in_data; ./download.sh

tableau:
	cd data_sources/tableau; ./download.sh
