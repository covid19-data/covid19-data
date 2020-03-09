all:
	cd COVID-19; git pull origin master
	snakemake -R `snakemake --lc --li --lp`
	git commit -a -m "data update"
	git pull
	git push
