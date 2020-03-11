.PHONY: all watch

all:
	git pull
	cd COVID-19; git pull origin master
	snakemake -R `snakemake --lc --li --lp`

watch:
	watch -n 100 "python check_new_commit.py"
