all:
	snakemake -R `snakemake --list-code-changes`
