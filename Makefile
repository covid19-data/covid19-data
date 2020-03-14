.PHONY: all watch

all:
	git pull
	cd data_sources/COVID-19; git pull origin master
	cd data_sources/our_world_in_data; ./download.sh
	snakemake -R `snakemake --lc --li --lp`

watch:
	watch -n 100 "python check_new_commit.py"
