---
title: Variant Calling Workflow
teaching: 35
exercises: 25
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand the steps involved in variant calling.
- Describe the types of data formats encountered during variant calling.
- Use command line tools to perform variant calling.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I find sequence variants between my sample and a reference genome?

::::::::::::::::::::::::::::::::::::::::::::::::::

We mentioned before that we are working with files from a long-term evolution study of an *E. coli* population (designated Ara-3). Now that we have looked at our data to make sure that it is high quality, and removed low-quality base calls, we can perform variant calling to see how the population changed over time. We care how this population changed relative to the original population, *E. coli* strain REL606. Therefore, we will align each of our samples to the *E. coli* REL606 reference genome, and see what differences exist in our reads versus the genome.

## Alignment to a reference genome

![](fig/variant_calling_workflow_align.png){alt='workflow\_align'}

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be
using the [Burrows Wheeler Aligner (BWA)](https://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent
sequences against a large reference genome.

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

## Setting up

First we download the reference genome for *E. coli* REL606. Although we could copy or move the file with `cp` or `mv`, most genomics workflows begin with a download step, so we will practice that here.

```bash
$ cd ~/dc_workshop_YEAST
$ mkdir -p data/ref_genome
$ cd data/ref_genome
$ ln -s /home/dc_workshop_YEAST/data/ref_genome/SACE_S288C_v1_allChr.fasta .
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise

We saved this file as `data/ref_genome/SACE_S288C_v1_allChr.fasta`.
How many chromosomes does have this genome?

:::::::::::::::  solution

### Solution

```bash
$ head data/ref_genome/SACE_S288C_v1_allChr.fasta
```
```bash
>SACE_S288C_v1_chr_01
>SACE_S288C_v1_chr_02
>SACE_S288C_v1_chr_03
>SACE_S288C_v1_chr_04
>SACE_S288C_v1_chr_05
>SACE_S288C_v1_chr_06
>SACE_S288C_v1_chr_07
>SACE_S288C_v1_chr_08
>SACE_S288C_v1_chr_09
>SACE_S288C_v1_chr_10
>SACE_S288C_v1_chr_11
>SACE_S288C_v1_chr_12
>SACE_S288C_v1_chr_13
>SACE_S288C_v1_chr_14
>SACE_S288C_v1_chr_15
>SACE_S288C_v1_chr_16
>SACE_S288C_v1_chr_mt
>SACE_S288C_v1_chr_2m
```


:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

We will also subset trimmed FASTQ files with [seqtk](https://www.biostars.org/p/6544/#42695) to work with. These are small subsets of our real trimmed data,
and will enable us to run our variant calling workflow quite quickly.

```bash
$ mkdir  ~/dc_workshop_YEAST/data/trimmed_fastq_small
$ cd  ~/dc_workshop_YEAST/data/trimmed_fastq_small
$ seqtk sample -s100 ~/dc_workshop_YEAST/data/trimmed_fastq/YMX005645_R1.trimmed.fastq 0.1 > YMX005645_R1_small.fastq
$ seqtk sample -s100 ~/dc_workshop_YEAST/data/trimmed_fastq/YMX005645_R2.trimmed.fastq 0.1 > YMX005645_R2_small.fastq
```

You will also need to create directories for the results that will be generated as part of this workflow. We can do this in a single
line of code, because `mkdir` can accept multiple new directory
names as input.

```bash
$ cd ~/dc_workshop_YEAST
$ mkdir -p results/sam results/bam results/bcf results/vcf
```

#### Index the reference genome

Our first step is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

```bash
$ bwa index data/ref_genome/SACE_S288C_v1_allChr.fasta
```

While the index is created, you will see output that looks something like this:

```output
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.08 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.75 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index data/ref_genome/SACE_S288C_v1_allChr.fasta
[main] Real time: 3.307 sec; CPU: 2.942 sec
```

#### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an
aligner. We will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it
is faster and more accurate.

An example of what a `bwa` command looks like is below. This command will not run, as we do not have the files `ref_genome.fa`, `input_file_R1.fastq`, or `input_file_R2.fastq`.

```bash
$ bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sam
```

Have a look at the [bwa options page](https://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default
parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using
and make sure the options you use are appropriate for your data.*

We are going to start by aligning the reads from just one of the
samples in our dataset (`YMX005645`). Later, we will be
iterating this whole process on all of our sample files.

```bash
$ bwa mem data/ref_genome/SACE_S288C_v1_allChr.fasta data/trimmed_fastq_small/YMX005645_R1_small.fastq data/trimmed_fastq_small/YMX005645_R2_small.fastq > results/sam/YMX005645.aligned.sam
```

You will see output that starts like this:

```output
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 67454 sequences (10000030 bp)...
[M::process] read 67342 sequences (10000276 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (2, 29284, 46, 0)
[M::mem_pestat] skip orientation FF as there are not enough pairs
[M::mem_pestat] analyzing insert size distribution for orientation FR...
[M::mem_pestat] (25, 50, 75) percentile: (252, 291, 333)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (90, 495)
[M::mem_pestat] mean and std.dev: (293.03, 60.86)
[M::mem_pestat] low and high boundaries for proper pairs: (9, 576)
[M::mem_pestat] analyzing insert size distribution for orientation RF...
[M::mem_pestat] (25, 50, 75) percentile: (56, 141, 515)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 1433)
[M::mem_pestat] mean and std.dev: (270.44, 234.23)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 1892)
[M::mem_pestat] skip orientation RR as there are not enough pairs
[M::mem_pestat] skip orientation RF
```

##### SAM/BAM format

The [SAM file](https://genome.sph.umich.edu/wiki/SAM),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not
have time to go into detail about the features of the SAM format, the paper by
[Heng Li et al.](https://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. The header is used to describe the source of data, reference sequence, method of
alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line
that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential
mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is
displayed below with the different fields highlighted.

![](fig/sam_bam.png){alt='sam\_bam1'}

![](fig/sam_bam3.png){alt='sam\_bam2'}

Reduce size of sam file
```bash
 head -n20 results/sam/YMX005645.aligned.sam > results/sam/output_header.sam 
 shuf -n 100000 results/sam/YMX005645.aligned.sam | head -n 100000 > results/sam/output_random.sam 
 cat results/sam/output_random.sam >> results/sam/output_header.sam
```

We will convert the SAM file to BAM format using the `samtools` program with the `view` command and tell this command that the input is in SAM format (`-S`) and to output BAM format (`-b`):

```bash
$ samtools view -S -b results/sam/output_header.sam > results/bam/output.bam #convert to bam
```

```output
[samopen] SAM header is present: 1 sequences.
```

#### Sort BAM file by coordinates

Next we sort the BAM file using the `sort` command from `samtools`. `-o` tells the command where to write the output.

```bash
$ samtools sort -o results/bam/output.sorted.bam results/bam/output.bam  # sort el bam
```

Our files are pretty small, so we will not see this output. If you run the workflow with larger files, you will see something like this:

```output
[bam_sort_core] merging from 2 files...
```

SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.

You can use samtools to learn more about this bam file as well.

```bash
samtools flagstat results/bam/output.sorted.bam 
```

This will give you the following statistics about your sorted bam file:

```output
100001 + 0 in total (QC-passed reads + QC-failed reads)
99545 + 0 primary
0 + 0 secondary
456 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
98872 + 0 mapped (98.87% : N/A)
98416 + 0 primary mapped (98.87% : N/A)
99545 + 0 paired in sequencing
49783 + 0 read1
49762 + 0 read2
97163 + 0 properly paired (97.61% : N/A)
98300 + 0 with itself and mate mapped
116 + 0 singletons (0.12% : N/A)
793 + 0 with mate mapped to a different chr
544 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Variant calling

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given position in an individual genome
or transcriptome, often referred to as a Single Nucleotide Variant (SNV). The call is usually accompanied by an estimate of
variant frequency and some measure of confidence. Similar to other steps in this workflow, there are a number of tools available for
variant calling. In this workshop we will be using `bcftools`, but there are a few things we need to do before actually calling the
variants.

![](fig/variant_calling_workflow.png){alt='workflow'}

#### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with
[bcftools](https://samtools.github.io/bcftools/bcftools.html). We will
use the command `mpileup`. The flag `-O b` tells bcftools to generate a
bcf format output file, `-o` specifies where to write the output file, and `-f` flags the path to the reference genome:

```bash
$ bcftools mpileup -O b -o results/bcf/output_raw.bcf -f data/ref_genome/SACE_S288C_v1_allChr.fasta results/bam/output.sorted.bam 
```

```output
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 250
```

We have now generated a file with coverage information for every base.

#### Step 2: Detect the single nucleotide variants (SNVs)

Identify SNVs using bcftools `call`. We have to specify ploidy with the flag `--ploidy`, which is one for the haploid *E. coli*. `-m` allows for multiallelic and rare-variant calling, `-v` tells the program to output variant sites only (not every site in the genome), and `-o` specifies where to write the output file:

```bash
$ bcftools call --ploidy 2 -m -v -o results/vcf/output_variants.vcf results/bcf/output_raw.bcf
```

#### Step 3: Filter and report the SNV variants in variant calling format (VCF)

Filter the SNVs for the final output in VCF format, using `vcfutils.pl`:

```bash
$ vcfutils.pl varFilter results/vcf/output_variants.vcf  > results/vcf/output_final_variants.vcf
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Filtering

The `vcfutils.pl varFilter` call filters out variants that do not meet minimum quality default criteria, which can be changed through
its options. Using `bcftools` we can verify that the quality of the variant call set has improved after this filtering step by
calculating the ratio of [transitions(TS)](https://en.wikipedia.org/wiki/Transition_%28genetics%29) to
[transversions (TV)](https://en.wikipedia.org/wiki/Transversion) ratio (TS/TV),
where transitions should be more likely to occur than transversions:

```bash
$ bcftools stats results/bcf/output_variants.vcf | grep TSTV
# TSTV, transitions/transversions:
# TSTV	[2]id	[3]ts	[4]tv	[5]ts/tv	[6]ts (1st ALT)	[7]tv (1st ALT)	[8]ts/tv (1st ALT)
TSTV	0	628	58	10.83	628	58	10.83
$ bcftools stats results/vcf/output_final_variants.vcf | grep TSTV
# TSTV, transitions/transversions:
# TSTV	[2]id	[3]ts	[4]tv	[5]ts/tv	[6]ts (1st ALT)	[7]tv (1st ALT)	[8]ts/tv (1st ALT)
TSTV	0	621	54	11.50	621	54	11.50
```

::::::::::::::::::::::::::::::::::::::::::::::::::

### Explore the VCF format:

```bash
$ less -S results/vcf/output_final_variants.vcf
```

You will see the header (which describes the format), the time and date the file was
created, the version of bcftools that was used, the command line parameters used, and
some additional information:

```output##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13+ds
##bcftoolsCommand=mpileup -O b -o results/bcf/output_raw.bcf -f data/ref_genome/SACE_S288C_v1_allChr.>
##reference=file://data/ref_genome/SACE_S288C_v1_allChr.fasta
##contig=<ID=SACE_S288C_v1_chr_01,length=230218>
##contig=<ID=SACE_S288C_v1_chr_02,length=813184>
##contig=<ID=SACE_S288C_v1_chr_03,length=316620>
##contig=<ID=SACE_S288C_v1_chr_04,length=1531933>
##contig=<ID=SACE_S288C_v1_chr_05,length=576874>
##contig=<ID=SACE_S288C_v1_chr_06,length=270161>
##contig=<ID=SACE_S288C_v1_chr_07,length=1090940>
##contig=<ID=SACE_S288C_v1_chr_08,length=562643>
##contig=<ID=SACE_S288C_v1_chr_09,length=439888>
##contig=<ID=SACE_S288C_v1_chr_10,length=745751>
##contig=<ID=SACE_S288C_v1_chr_11,length=666816>
##contig=<ID=SACE_S288C_v1_chr_12,length=1078177>
##contig=<ID=SACE_S288C_v1_chr_13,length=924431>
##contig=<ID=SACE_S288C_v1_chr_14,length=784333>
##contig=<ID=SACE_S288C_v1_chr_15,length=1091291>
##contig=<ID=SACE_S288C_v1_chr_16,length=948066>
##contig=<ID=SACE_S288C_v1_chr_mt,length=85779>
##contig=<ID=SACE_S288C_v1_chr_2m,length=6318>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
```

Followed by information on each of the variations observed:

```output
##bcftools_callVersion=1.13+htslib-1.13+ds
##bcftools_callCommand=call --ploidy 2 -m -v -o results/vcf/output_variants.vcf results/bcf/output_raw.bcf; Date=Thu Oct  5 02:16:48 2023
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  results/bam/output.sorted.bam
SACE_S288C_v1_chr_01    537     .       T       C       14.5687 .       DP=2;VDB=0.96;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=23     GT:PL   1/1:44,6,0
SACE_S288C_v1_chr_01    1352    .       A       G       24.9996 .       DP=4;VDB=0.36;SGB=-0.453602;RPBZ=-0.774597;MQBZ=0.408248;MQSBZ=-0.942809;BQBZ=1.54919;SCBZ=1;FS=0;MQ0F=0.25;AC=1;>
SACE_S288C_v1_chr_01    1658    .       C       T       22.4391 .       DP=2;VDB=0.84;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=30     GT:PL   1/1:52,6,0
SACE_S288C_v1_chr_01    1774    .       A       G       30.4038 .       DP=5;VDB=0.52;SGB=-0.453602;RPBZ=0;MQBZ=1.77705;MQSBZ=0.296174;BQBZ=1.82574;FS=0;MQ0F=0;AC=1;AN=2;DP4=1,2,1,1;MQ=>
SACE_S288C_v1_chr_01    2706    .       A       G       38.415  .       DP=2;VDB=0.02;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60     GT:PL   1/1:68,6,0
SACE_S288C_v1_chr_01    2748    .       A       T       35.4156 .       DP=2;VDB=0.06;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60     GT:PL   1/1:65,6,0
SACE_S288C_v1_chr_01    2790    .       C       A       32.4168 .       DP=2;VDB=0.06;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60     GT:PL   1/1:62,6,0
SACE_S288C_v1_chr_01    2796    .       G       A       34.4159 .       DP=2;VDB=0.06;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60     GT:PL   1/1:64,6,0
SACE_S288C_v1_chr_01    2797    .       A       G       38.415  .       DP=2;VDB=0.06;SGB=-0.453602;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60     GT:PL   1/1:68,6,0
```

This is a lot of information, so let's take some time to make sure we understand our output.

The first few columns represent the information we have about a predicted variation.

| column  | info                                                                                                                                                                                                                                                                                                                                    | 
| ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CHROM   | contig location where the variation occurs                                                                                                                                                                                                                                                                                              | 
| POS     | position within the contig where the variation occurs                                                                                                                                                                                                                                                                                   | 
| ID      | a `.` until we add annotation information                                                                                                                                                                                                                                                                                                                                      | 
| REF     | reference genotype (forward strand)                                                                                                                                                                                                                                                                                                     | 
| ALT     | sample genotype (forward strand)                                                                                                                                                                                                                                                                                                        | 
| QUAL    | Phred-scaled probability that the observed variant exists at this site (higher is better)                                                                                                                                                                                                                                               | 
| FILTER  | a `.` if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed                                                                                                                                                                                                                                                                                                                                      | 

In an ideal world, the information in the `QUAL` column would be all we needed to filter out bad variant calls.
However, in reality we need to filter on multiple other metrics.

The last two columns contain the genotypes and can be tricky to decode.

| column  | info                                                                                                                                                                                                                                                                                                                                    | 
| ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FORMAT  | lists in order the metrics presented in the final column                                                                                                                                                                                                                                                                                | 
| results | lists the values associated with those metrics in order                                                                                                                                                                                                                                                                                 | 

For our file, the metrics presented are GT:PL:GQ.

| metric  | definition                                                                                                                                                                                                                                                                                                                              | 
| ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| AD, DP  | the depth per allele by sample and coverage                                                                                                                                                                                                                                                                                             | 
| GT      | the genotype for the sample at this loci. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. A 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. | 
| PL      | the likelihoods of the given genotypes                                                                                                                                                                                                                                                                                                  | 
| GQ      | the Phred-scaled confidence for the genotype                                                                                                                                                                                                                                                                                            | 

The Broad Institute's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place
to learn more about the VCF file format.

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise

Use the `grep` and `wc` commands you have learned to assess how many variants are in the vcf file.

:::::::::::::::  solution

### Solution

```bash
$ grep -v "#" results/vcf/output_final_variants.vcf | wc -l
```

```output
20336
```

There are 20336 variants in this file.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Assess the alignment (visualization) - optional step

It is often instructive to look at your data in a genome browser. Visualization will allow you to get a "feel" for
the data, as well as detecting abnormalities and problems. Also, exploring the data in such a way may give you
ideas for further analyses.  As such, visualization tools are useful for exploratory analysis. In this lesson we
will describe two different tools for visualization: a light-weight command-line based one and the Broad
Institute's Integrative Genomics Viewer (IGV) which requires
software installation and transfer of files.

In order for us to visualize the alignment files, we will need to index the BAM file using `samtools`:

```bash
$ samtools index results/bam/output.sorted.bam
```

#### Viewing with `tview`

[Samtools](https://www.htslib.org/) implements a very simple text alignment viewer based on the GNU
`ncurses` library, called `tview`. This alignment viewer works with short indels and shows [MAQ](https://maq.sourceforge.net/) consensus.
It uses different colors to display mapping quality or base quality, subjected to users' choice. Samtools viewer is known to work with a 130 GB alignment swiftly. Due to its text interface, displaying alignments over network is also very fast.

In order to visualize our mapped reads, we use `tview`, giving it the sorted bam file and the reference file:

```bash
$ samtools tview -p 1:10000 results/bam/output.sorted.bam data/ref_genome/SACE_S288C_v1_allChr.fasta
```

```output
1         11        21        31        41        51        61        71        81        91        101       111       121
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATAC
..................................................................................................................................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,.............................
...................................,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................   ................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................................   ....................................      ,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................  ,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,     .......
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .............................  ,,,,,,,,,,,,,,,,,g,,,,,    ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ...........................T.......   ,,,,,,,,,,,,,,,,,,,,,,,c,          ......
......................... ................................   ,g,,,,,,,,,,,,,,,,,,,      ...........................
,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,       ..........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ................................T..  ..............................   ,,,,,,
...........................       ,,,,,,g,,,,,,,,,,,,,,,,,   ....................................         ,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................  ...................................        ....
....................................  ........................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
........................            .................................. .............................     ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ....................................        ..........................
...............................       ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................
...................................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ..................................
.................................... ,,,,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,        ,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ............................ ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
```

The first line of output shows the genome coordinates in our reference genome. The second line shows the reference
genome sequence. The third line shows the consensus sequence determined from the sequence reads. A `.` indicates
a match to the reference sequence, so we can see that the consensus from our sample matches the reference in most
locations. That is good! If that was not the case, we should probably reconsider our choice of reference.

Below the horizontal line, we can see all of the reads in our sample aligned with the reference genome. Only
positions where the called base differs from the reference are shown. You can use the arrow keys on your keyboard
to scroll or type `?` for a help menu. To navigate to a specific position, type `g`. A dialogue box will appear. In
this box, type the name of the "chromosome" followed by a colon and the position of the variant you would like to view
(e.g. for this sample, type `CP000819.1:50` to view the 50th base. Type `Ctrl^C` or `q` to exit `tview`.

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise

Visualize the alignment of the reads for our `SRR2584866` sample. What variant is present at
position 4377265? What is the canonical nucleotide in that position?

:::::::::::::::  solution

### Solution

```bash
$ samtools tview ~/dc_workshop_YEAST/results/bam/SRR2584866.aligned.sorted.bam ~/dc_workshop_YEAST/data/ref_genome/SACE_S288C_v1_allChr.fasta
```

Then type `g`. In the dialogue box, type `CP000819.1:4377265`.
`G` is the variant. `A` is canonical. This variant possibly changes the phenotype of this sample to hypermutable. It occurs
in the gene *mutL*, which controls DNA mismatch repair.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

#### Viewing with IGV

[IGV](https://www.broadinstitute.org/igv/) is a stand-alone browser, which has the advantage of being installed locally and providing fast access. Web-based genome browsers, like [Ensembl](https://www.ensembl.org/index.html) or the [UCSC browser](https://genome.ucsc.edu/), are slower, but provide more functionality. They not only allow for more polished and flexible visualization, but also provide easy access to a wealth of annotations and external data sources. This makes it straightforward to relate your data with information about repeat regions, known genes, epigenetic features or areas of cross-species conservation, to name just a few.

In order to use IGV, we will need to transfer some files to our local machine. We know how to do this with `scp`.
Open a new tab in your terminal window and create a new folder. We will put this folder on our Desktop for
demonstration purposes, but in general you should avoide proliferating folders and files on your Desktop and
instead organize files within a directory structure like we have been using in our `dc_workshop` directory.

```bash
$ mkdir ~/Desktop/files_for_igv
$ cd ~/Desktop/files_for_igv
```

Now we will transfer our files to that new directory. Remember to replace the text between the `@` and the `:`
with your AWS instance number. The commands to `scp` always go in the terminal window that is connected to your
local computer (not your AWS instance).

```bash
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop_YEAST/results/bam/SRR2584866.aligned.sorted.bam ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop_YEAST/results/bam/SRR2584866.aligned.sorted.bam.bai ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop_YEAST/data/ref_genome/SACE_S288C_v1_allChr.fasta ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop_YEAST/results/vcf/SRR2584866_final_variants.vcf ~/Desktop/files_for_igv
```

You will need to type the password for your AWS instance each time you call `scp`.

Next, we need to open the IGV software. If you have not done so already, you can download IGV from the [Broad Institute's software page](https://www.broadinstitute.org/software/igv/download), double-click the `.zip` file
to unzip it, and then drag the program into your Applications folder.

1. Open IGV.
2. Load our reference genome file (`SACE_S288C_v1_allChr.fasta`) into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
3. Load our BAM file (`SRR2584866.aligned.sorted.bam`) using the **"Load from File..."** option under the **"File"** pull-down menu.
4. Do the same with our VCF file (`SRR2584866_final_variants.vcf`).

Your IGV browser should look like the screenshot below:

![](fig/igv-screenshot.png){alt='IGV'}

There should be two tracks: one coresponding to our BAM file and the other for our VCF file.

In the **VCF track**, each bar across the top of the plot shows the allele fraction for a single locus. The second bar shows
the genotypes for each locus in each *sample*. We only have one sample called here, so we only see a single line. Dark blue =
heterozygous, Cyan = homozygous variant, Grey = reference.  Filtered entries are transparent.

Zoom in to inspect variants you see in your filtered VCF file to become more familiar with IGV. See how quality information
corresponds to alignment information at those loci.
Use [this website](https://software.broadinstitute.org/software/igv/AlignmentData) and the links therein to understand how IGV colors the alignments.

Now that we have run through our workflow for a single sample, we want to repeat this workflow for our other five
samples. However, we do not want to type each of these individual steps again five more times. That would be very
time consuming and error-prone, and would become impossible as we gathered more and more samples. Luckily, we
already know the tools we need to use to automate this workflow and run it on as many files as we want using a
single line of code. Those tools are: wildcards, for loops, and bash scripts. We will use all three in the next
lesson.

:::::::::::::::::::::::::::::::::::::::::  callout

### Installing software

It is worth noting that all of the software we are using for
this workshop has been pre-installed on our remote computer.
This saves us a lot of time - installing software can be a
time-consuming and frustrating task - however, this does mean that
you will not be able to walk out the door and start doing these
analyses on your own computer. You will need to install
the software first. Look at the [setup instructions](https://www.datacarpentry.org/wrangling-genomics/setup.html) for more information
on installing these software packages.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

### BWA alignment options

BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence
reads up to 100bp, while the other two are for sequences ranging from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such
as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it
is faster and more accurate.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- Bioinformatic command line tools are collections of commands that can be used to carry out bioinformatic analyses.
- To use most powerful bioinformatic tools, you will need to use the command line.
- There are many different file formats for storing genomics data. It is important to understand what type of information is contained in each file, and how it was derived.

::::::::::::::::::::::::::::::::::::::::::::::::::


