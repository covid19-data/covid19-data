.PHONY: all watch

all:
	git pull
	cd COVID-19; git pull origin master
	snakemake -R `snakemake --lc --li --lp`
	git commit -a -m "data update"

watch:
	watch -n 300 "python check_new_commit.py"
