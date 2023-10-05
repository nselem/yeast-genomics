---
title: Trimming and Filtering
teaching: 30
exercises: 25
---

::::::::::::::::::::::::::::::::::::::: objectives

- Clean FASTQ reads using Trimmomatic.
- Select and set multiple options for command-line bioinformatic tools.
- Write `for` loops with two variables.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How can I get rid of sequence data that does not meet my quality standards?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Cleaning reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This does not mean,
though, that our samples should be thrown out! It is very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called
[Trimmomatic](https://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic options

Trimmomatic has a variety of options to trim your reads. If we run the following command, we can see some of our options.

```bash
$ trimmomatic
```

Which will give you the following output:

```output
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
```

This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads.
Next, we specify what flag we would like to run. For example, you can specify `threads` to indicate the number of
processors on your computer that you want Trimmomatic to use. In most cases using multiple threads (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important.
In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

| option         | meaning                                                                                                      | 
| -------------- | ------------------------------------------------------------------------------------------------------------ |
| \<inputFile1>   | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.                                          | 
| \<inputFile2>   | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.                                          | 
| \<outputFile1P> | Output file that contains surviving pairs from the `_1` file.                                                          | 
| \<outputFile1U> | Output file that contains orphaned reads from the `_1` file.                                                           | 
| \<outputFile2P> | Output file that contains surviving pairs from the `_2` file.                                                          | 
| \<outputFile2U> | Output file that contains orphaned reads from the `_2` file.                                                           | 

The last thing trimmomatic expects to see is the trimming parameters:

| step           | meaning                                                                                                      | 
| -------------- | ------------------------------------------------------------------------------------------------------------ |
| `ILLUMINACLIP`               | Perform adapter removal.                                                                                     | 
| `SLIDINGWINDOW`               | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. | 
| `LEADING`               | Cut bases off the start of a read, if below a threshold quality.                                             | 
| `TRAILING`               | Cut bases off the end of a read, if below a threshold quality.                                               | 
| `CROP`               | Cut the read to a specified length.                                                                          | 
| `HEADCROP`               | Cut the specified number of bases from the start of the read.                                                | 
| `MINLEN`               | Drop an entire read if it is below a specified length.                                                       | 
| `TOPHRED33`               | Convert quality scores to Phred-33.                                                                          | 
| `TOPHRED64`               | Convert quality scores to Phred-64.                                                                          | 

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](https://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

However, a complete command for Trimmomatic will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

```bash
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
```

In this example, we have told Trimmomatic:

| code           | meaning                                                                                                      | 
| -------------- | ------------------------------------------------------------------------------------------------------------ |
| `PE`               | that it will be taking a paired end file as input                                                            | 
| `-threads 4`               | to use four computing threads to run (this will speed up our run)                                            | 
| `SRR_1056_1.fastq`               | the first input file name                                                                                    | 
| `SRR_1056_2.fastq`               | the second input file name                                                                                   | 
| `SRR_1056_1.trimmed.fastq`               | the output file for surviving pairs from the `_1` file                                                                | 
| `SRR_1056_1un.trimmed.fastq`               | the output file for orphaned reads from the `_1` file                                                                 | 
| `SRR_1056_2.trimmed.fastq`               | the output file for surviving pairs from the `_2` file                                                                | 
| `SRR_1056_2un.trimmed.fastq`               | the output file for orphaned reads from the `_2` file                                                                 | 
| `ILLUMINACLIP:SRR_adapters.fa`               | to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa`                     | 
| `SLIDINGWINDOW:4:20`               | to use a sliding window of size 4 that will remove bases if their phred score is below 20                    | 

:::::::::::::::::::::::::::::::::::::::::  callout

## Multi-line commands

Some of the commands we ran in this lesson are long! When typing a long
command into your terminal, you can use the `\` character
to separate code chunks onto separate lines. This can make your code more readable.


::::::::::::::::::::::::::::::::::::::::::::::::::

## Running Trimmomatic

Now we will run Trimmomatic on our data. To begin, navigate to your `untrimmed_fastq` data directory:

```bash
$ cd ~/dc_workshop_YEAST/data/untrimmed_fastq
```

We are going to run Trimmomatic on one of our paired-end samples.
While using FastQC we did not saw that Nextera adapters were present in our samples.
If they were, the adapter sequences came with the installation of trimmomatic, so we would have needed to copy these sequences into our current directory.

```bash
# Do not run
$ cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
```

We will  use a sliding window of size 4 that will remove bases if their
phred score is below 20 (like in our example above). We will also
discard any reads that do not have at least 100 bases remaining after
this trimming step. This command will take a few minutes to run.

```bash
$ trimmomatic PE -threads 4 YMX005645_R1.fastq YMX005645_R2.fastq \
 YMX005645_R1.trimmed.fastq YMX005645_R1un.trimmed.fastq \
 YMX005645_R2.trimmed.fastq YMX005645_R2un.trimmed.fastq \
 SLIDINGWINDOW:4:20 MINLEN:100 LEADING:5 TRAILING:5
```

```output
TrimmomaticPE: Started with arguments:
 -threads 4 YMX005645_R1.fastq YMX005645_R2.fastq YMX005645_R1.trimmed.fastq YMX005645_R1un.trimmed.fastq YMX005645_R2.trimmed.fastq YMX005645_R2un.trimmed.fastq SLIDINGWINDOW:4:20 MINLEN:100 LEADING:5 TRAILING:5
Quality encoding detected as phred64
Input Read Pairs: 8020412 Both Surviving: 7029362 (87.64%) Forward Only Surviving: 433562 (5.41%) Reverse Only Surviving: 369853 (4.61%) Dropped: 187635 (2.34%)
TrimmomaticPE: Completed successfully
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Exercise

Use the output from your Trimmomatic command to answer the
following questions.

1) What percent of reads did we discard from our sample?
2) What percent of reads did we keep both pairs?

:::::::::::::::  solution

## Solution

1) 2\.34%
2) 87\.64%
  
  

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

You may have noticed that Trimmomatic automatically detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output files:

```bash
$ ls  YMX005645_R1*
```

```output
 YMX005645_R1.fastq.gz        YMX005645_R1un.trim.fastq.gz   YMX005645_R2.trim.fastq.gz
 YMX005645_R1.trim.fastq.gz   YMX005645_R2.fastq.gz          YMX005645_R2un.trim.fastq.gz
```

The output files are also FASTQ files. It should be smaller than our
input file, because we have removed reads. We can confirm this:

```bash
$ ls  YMX005645* -lL -h
```

```output
-rw-r--r-- 1 root    root          2.6G Oct  4 22:00 YMX005645_R1.fastq
-rw-r--r-- 1 user rstudio-users 2.2G Oct  5 00:08 YMX005645_R1.trimmed.fastq
-rw-r--r-- 1 user rstudio-users 135M Oct  5 00:08 YMX005645_R1un.trimmed.fastq
-rw-r--r-- 1 root    root          2.6G Oct  4 23:17 YMX005645_R2.fastq
-rw-r--r-- 1 user rstudio-users 2.2G Oct  5 00:08 YMX005645_R2.trimmed.fastq
-rw-r--r-- 1 user rstudio-users 113M Oct  5 00:08 YMX005645_R2un.trimmed.fastq
```

We have just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that If we have more samples we can use a `for` loop to iterate through our sample files
quickly!

We unzipped one of our files before to work with it, let's compress it again before we run our for loop.

```bash
gzip YMX005645_R1.fastq 
```

```bash
$ for infile in *_1.fastq.gz
> do
>   base=$(basename ${infile} _1.fastq.gz)
>   trimmomatic PE ${infile} ${base}_2.fastq.gz \
>                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
>                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
>                SLIDINGWINDOW:4:20 MINLEN:100 LEADING:5 TRAILING:5 
> done
```
For example if you have another sample SRR2589044 _1.fastq.gz.
Go ahead and run the for loop. It should take a few minutes for
Trimmomatic to run for each of our six input files. Once it is done
running, take a look at your directory contents. You will notice that even though we ran Trimmomatic on file `SRR2589044` before running the for loop, there is only one set of files for it. Because we matched the ending `_1.fastq.gz`, we re-ran Trimmomatic on this file, overwriting our first results. That is ok, but it is good to be aware that it happened.

```bash
$ ls
```

```output
YMX005645_R1.fastq          YMX005645_R1un.trimmed.fastq  YMX005645_R2.trimmed.fastq
YMX005645_R1.trimmed.fastq  YMX005645_R2.fastq            YMX005645_R2un.trimmed.fastq
```

:::::::::::::::::::::::::::::::::::::::  challenge


We have now completed the trimming and filtering steps of our quality
control process! Before we move on, let's move our trimmed FASTQ files
to a new subdirectory within our `data/` directory.

```bash
$ cd ~/dc_workshop_YEAST/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ mv *.trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
```

```output
YMX005645_R1.trimmed.fastq  YMX005645_R1un.trimmed.fastq  YMX005645_R2.trimmed.fastq  YMX005645_R2un.trimmed.fastq
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Bonus exercise (advanced)

Now that our samples have gone through quality control, they should perform
better on the quality tests run by FastQC. Go ahead and re-run
FastQC on your trimmed FASTQ files and visualize the HTML files
to see whether your per base sequence quality is higher after
trimming.

:::::::::::::::  solution

## Solution

In your AWS terminal window do:

```bash
$ fastqc ~/dc_workshop_YEAST/data/trimmed_fastq/*.fastq*
```

In a new tab in your terminal do:

```bash
$ mkdir ~/Desktop/fastqc_html/trimmed
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop_YEAST/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
```

Then take a look at the html files in your browser.

Remember to replace everything between the `@` and `:` in your scp
command with your AWS instance number.

After trimming and filtering, our overall quality is much higher,
we have a distribution of sequence lengths, and more samples pass
adapter content. However, quality trimming is not perfect, and some
programs are better at removing some sequences than others. Because our
sequences still contain 3' adapters, it could be important to explore
other trimming tools like [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove these, depending on your
downstream application. Trimmomatic did pretty well though, and its performance
is good enough for our workflow.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- The options you set for the command-line tools you use are important!
- Data cleaning is an essential step in a genomics workflow.

::::::::::::::::::::::::::::::::::::::::::::::::::


