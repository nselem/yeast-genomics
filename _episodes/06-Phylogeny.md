---
title: Running a Phylogenetic tree
teaching: 30
exercises: 15
---

::::::::::::::::::::::::::::::::::::::: objectives

- Run a phylogenetic tree.
- Visualize the result with associated metadata.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How are the sequenced isolates related?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Run the phylogenetic tree

You got vcf files with variants called for each one of the samples using a bash script in the previous lesson.

Due to time and space constrains, we will jump from individual .vcf files to multiple .gvcf files. The latter, are a specific form of vcf that stores information not only about variant sites, but also about invariant blocks of sequence. Also, they contain information about many samples and therefore are a good starting point for comparison of variant sites among the dataset.

Now we are going to 

```bash
cd dc_workshop_YEAST/multvcf/

less  MATRIX_24SACE.vcf

```

Navigate the file until you see how the gvcf file is organized.

To run the phylogeny in RaXML, we need to convert the multivcf

```bash
python vcf2phylip.py -i dc_workshop_YEAST/multvcf/MATRIX_24SACE.vcf 
```

Finally, we will run the maximum likelihood tree

```bash
raxmlHPC-PTHREADS-AVX2 -f a -x 12345 -p 12345 -N 100  -T 20 -m GTRGAMMA -s dc_workshop_YEAST/multvcf/MATRIX_24SACE.vcf.min4.phy -n Phylogeny_Yeast.tree 
```

To visualize the tree, we can use [Microreact](https://microreact.org/)  
  
  
load the .tree output along with the metadata to visualize the results.


:::::::::::::::::::::::::::::::::::::::: keypoints

- We can build Phylip alignments from multi gvcf.
- RaXML runs maximum likelihood phylogenetic trees.

::::::::::::::::::::::::::::::::::::::::::::::::::


