
run:
	nextflow run workflows/run_bam_readcount.nf \
		-profile test,singularity \
		--input tests/files/manifest.csv \
		--sites tests/files/sites.bed \
		--fasta https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta \
		--outdir results
