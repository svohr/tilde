tilde (Testing Identity with Linkage DisEquilibrium)
====================================================


Introduction 
------------

This repository contains scripts and additional files related to the 
manuscript:

> Samuel H. Vohr, Carlos Fernando Buen Abad Najar, Beth Shapiro, Richard E.
> Green. **A method for positive forensic identification of samples 
> from extremely low-coverage sequence data.** *BMC Genomics* 2015, **16**:1034
> \[[link](http://www.biomedcentral.com/1471-2164/16/1034)\]

This manuscript presents a method for determining whether two samples are from
the same individual from extremely low-coverage, high-throughput sequencing
data. In the cases we consider, coverage of the nuclear genome is so low that 
there is essentially no overlap between samples. Bases can be observed 
for many potentially informative positions (single nucleotide polymorphims, 
SNPs), but the positions observed in one sample are completely different
from the positions observed in another. This is common in highly-degraded 
forensic and ancient DNA samples and makes finding a "match" by directly 
comparing samples impossible. Instead, we propose an approach that exploits 
patterns of linkage disequilibrium in human populations and our growing 
knowledge of human genetic diversity to ask whether observations made in two 
samples are consistent with originating in a single individual.

The approach described in the manuscript is implemented as a series of 
Python scripts. This pipeline was designed to make it easy to check 
intermediate results and to adjust the procedures and parameters of our 
method. The pipeline takes mapped reads in BAM format and a reference panel 
of haplotypes in [IMPUTE hap/legend/sample format](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample) and ultimately produces
aggregated log-likelihood ratios that can be used to assess (through plotting,
or statistical analysis) whether the observations made in two samples are 
consistent with originating in a single individual from the population 
described by the reference panel.

Our approach requires a reference panel of phased haplotypes to model the
allele and haplotype frequencies of the population from which our samples 
originated. We used statistically-phased SNP data from the 
[1000 Genomes Project](http://www.1000genomes.org/), which includes a large
number of individuals from a diverse range of populations. To build the
reference panels, we used [VCFtools](https://vcftools.github.io/) to apply
basic filters and to generate the haplotype tables in [IMPUTE format](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample). 
We observe bases for these SNP positions from mapped reads using 
[Samtools](http://www.htslib.org/). 


Scripts
-------

### make\_obs\_tab
Takes an [IMPUTE format](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample) legend file and one or more BAM files and 
produces a table of SNP positions with single-base observations from the BAM 
files. Base observations from BAM files are made using output from calls to 
`samtools mpileup`. Mapping quality, base quality, and coverage filters are 
applied to reduce the chance of errors due to mismapped reads. 

### pairs\_in\_range
This script takes as input a tab-delimited file of SNP positions and base 
observations produced by `make_obs_tab`. From this table, the script finds 
all pairs of SNPs with an observation from each sample that are within a 
specified range (in base pairs) and pass filtering criteria.

### indv\_test
Reads in a tab-delimited table of pairs of SNP positions and observations
from 2 samples as produced by `pairs_in_range` and calculates the 
probabilities of observing those two bases under a single individual model or 
two unrelated individuals model. The log-likelihood ratio comparing the two 
is written out with the original SNP positions.

### sample\_pairs
Downsamples output from `indv_test` run(s) and aggregates log-likelihood 
ratios by sliding window. Resample-aggregate steps are repeated to assess the 
mean and variance of the genome-wide aggregated log-likelihood ratio. Unlike
the other scripts, input can be made up of results from multiple chromosomes.


Example
-------

### 1. Build the reference panel.
   ```
   vcftools --gzvcf 1000g/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 
            --out chr21.CEU_panel 
            --remove-indels 
            --min-alleles 2 
            --max-alleles 2 
            --mac 10 
            --keep CEU_indv.txt 
            --bed mappable_regions.bed 
            --IMPUTE
   ```
   VCFTools is used to filter for bi-allelic SNPs with a minor allele counts of
   at least 10. We also filter SNPs based on a BED file of mappable regions. 
   This can be generated using [mappability annotation](http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=454283543_xadktKfUePJnQNSDZdmRKBlu4Tvg&c=chr21&g=wgEncodeMapability) from the [UCSC Genome Browser](http://www.genome.ucsc.edu/) or 
   Heng Li's [SNPable](http://lh3lh3.users.sourceforge.net/snpable.shtml) 
   program, for example. The file `CEU_indv.txt` contains the identifiers
   for the individuals from the CEU (Utah residents with Northern and Western 
   European ancestry). 

### 2. Observe single bases from sequence data
   ```
   make_obs_tab -q 30 -Q 30 chr21 chr21.CEU_panel.impute.legend sample1.bam sample2.bam > chr21.CEU_panel.obs.tab
   ```
   This builds a table containing the SNP positions and base observations 
   for all positions that have an overlapping read in any of the samples
   that pass quality filters. Here we set the minimum mapping quality and 
   minimum base quality to 30. Positions that are not observed in any sample
   are not included in final table.

### 3. Find pairs of observed SNPs in close proximity. 
   ```
   pairs_in_range -r 100-50000 -d 4 5 chr21.CEU_panel.obs.tab > chr21.CEU_panel.s1_s2.pairs.tab
   ```
   This invocation finds all pairs of SNPs with 100 bp and 50,000 bp where
   one SNP has observation from Sample 1 and the other SNP has an observation
   from Sample 2. The `-d` flag excludes SNPs with observations that may be
   the result of damage (e.g. where a T is observed and the reference or
   alternate base is a C). The first two positional arguments give the
   column in the table from which to draw observations (where the first 
   column of the input table is 0). The same column can be specified for 
   both arguments.

### 4. Calculate the log-likelihood ratios for all pairs
   ```
   indv_test chr21.CEU_panel.impute.hap chr21.CEU_panel.impute.legend chr21.CEU_panel.s1_s2.pairs.tab > chr21.CEU_panel.s1_s2.probs.tab
   ```
   This script calculates the log-likelihood ratios and reports other 
   statistics for all pairs of SNPs.

   Alternatively, steps 2 through 4 can be run together by piping the output of
   one script to the next.
   ```
   make_obs_tab -q 30 -Q 30 chr21 chr21.CEU_panel.impute.legend sample1.bam sample2.bam 
   | pairs_in_range -r 100-50000 -d 4 5 
   | indv_test chr21.CEU_panel.impute.hap chr21.CEU_panel.impute.legend
   > chr21.CEU_panel.s1_s2.pairs.tab
   ```

### 5. Aggregate log-likelihood ratios by sampling with sliding windows.
   ```
   samples_pairs chr21.CEU_panel.s1_s2.probs.tab > chr21.CEU_panel.s1_s2.samp.tab
   ```
   Pairs of SNPs are sampled from sliding windows to mitigate the affects of 
   linkage of the results. The log-likelihood ratios are summed together to
   summarize the results. This is repeated (default=1,000 times) to build
   an empirical distribution of the aggregated log-likelihood ratio. This 
   can be run on the combined results from multiple chromosomes.


Simulations
-----------
Scripts that were used to convert genotypes from coalescent simulations
to simulated low-coverage libraries and reference panels are included in
the directory `simulations`.

