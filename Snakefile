

# go in each folder, zip up data into a thing


account="dhm267"

rule all:
  input: expand("barseqSortedOnGAP1GFP_M{mismatches}.counts",mismatches=[0])

# make tmp folder


rule gunzip_fastq:
  input:
    gzipped="/scratch/cgsb/gencore/out/Gresham/2017-04-14_AVA3A/1/000000000-AVA3A_l01n01.3310000008a5cf.fastq.gz"
  output: "{basename}.fastq"
  shell: 
    """
    mkdir -p data
    sbatch --mail-user={account}@nyu.edu --mail-type ALL \
      --job-name "gunzip" --output all_%j.out \
      --nodes 1 --tasks-per-node 1 \
      --mem=16GB -t 00:30:00 \
      --wrap 'zcat {input.gzipped} > {output}'
    """

rule barnone:
  input:
    gzipped="{basename}.fastq",
    sample_index="sampleBarcodesRobinson2014.txt",
    strain_index="nislow_revised.txt",
  output:
    counts="{basename}_M{mismatches}.counts",
    mismatch="{basename}_M{mismatches}.mismatch",
    revised="{basename}_M{mismatches}.revised",
  shell: 
    """
    sbatch --mail-user={account}@nyu.edu --mail-type ALL \
      --job-name 'barnone' --output all_%j.out \
      --nodes 1 --tasks-per-node 1 \
      --mem=30GB -t 00:36:00 \
      --wrap \
      '\
      module purge ;\
      module barnone/intel/20170501 ;\
      BarNone -f fastq \
        --multiplexfile {input.sample_index} \
        --multiplexstart 1 --multiplexlength 5 \
        --tagstart 18 --taglength 3 --start 9 \
        --mismatches {wildcards.mismatches} \
        --mismatchfile {output.mismatch} \
        --revisedcatalog {output.revised} \
        -p 500000 \
        {input.gzipped} \
        {output.counts} \
        {input.strain_index} ; \
      '
    """

grep 'Note=Dubious' ../ref/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff > tmp/dubious_orfs.gff
grep -v '#' ../ref/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff | grep -v 'Note=Dubious' | grep 'gene' > tmp/notdubious_orfs.gff
bedtools slop -b 200 -i tmp/notdubious_orfs.gff -g ../ref/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.fasta.fai > tmp/notdubious_expanded_orfs.gff 
bedtools subtract -a tmp/dubious_orfs.gff -b tmp/notdubious_expanded_orfs.gff > tmp/dorfs_no_overlap.gff
grep -o 'ID=[^;]\+;' tmp/dorfs_no_overlap.gff | sed 's/ID=//g' | sed 's/;//' > tmp/dorfs_no_overlap_list.txt


archive_examples_data.zip : examples/data/
	zip -r $@ $<

.PHONY: unpack
unpack:
	unzip archive_examples_data.zip
Rscript -e rmarkdown::render('term_enrichment_clusterProfiler.Rmd')
