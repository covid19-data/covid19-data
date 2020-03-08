all:
	cd COVID-19; git pull origin master
	snakemake -R `snakemake --lc --li --lp`
