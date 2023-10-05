---
title: Background and Metadata
teaching: 10
exercises: 5
---

::::::::::::::::::::::::::::::::::::::: objectives

- Why study *S. cerevisiae*?
- Understand the data set.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What data are we using?
- Why is this experiment important?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Background

We are going to use a country-wide sequencing dataset from *Saccharomyces cerevisiae* collected from agave fermentation facilities.

- **What is *S. cerevisiae*?**
  - *S. cerevisiae* are unicellular fungus also known as budding yeasts. It is commonly found in fermentations of dairy, spirits, and other food. 


<!-- https://es.wikipedia.org/wiki/Saccharomyces_cerevisiae#/media/Archivo:S_cerevisiae_under_DIC_microscopy.jpg -->

- **Why is *S. cerevisiae* important?**
  - *S. cerevisiae* are one of the most well-studied model organisms in science. As a single-celled organism, *S. cerevisiae* typically doubling its population every two hours in the lab, which means it can be manipulated easily in experiments. In addition, most naturally occurring strains of *S. cerevisiae* are harmless. Finally, it is a great model for eukaryotic genetics and it is easily manipulated to study molecular biology, biochemistry, adaptation, and evolution. Still, its ecological roles are currently being studied.

## The data

- The data we will use is part of the Yeast Genomes MX Initiative led by [Lucía Morales](https://liigh.unam.mx/profile/dra-lucia-morales/), [Alexander DeLuna](https://langebio.cinvestav.mx/Dr-Alexander-de-Luna-Fors) and [Eugenio Mancera](https://ira.cinvestav.mx/ingenieriagenetica/dr-eugenio-mancera/) at UNAM and CINVESTAV.

- The experiment was designed to explore the genomic diversity of Saccharomyces yeasts across a megadiverse context that surrounds the culturally and biotechnologically relevant traditional agave fermentations, which are carried on in regions that go from highland forests to lowland deserts. During the sampling, we collected samples from fermentations of other species being carried out in the same distilleries.

- To read the first publications of the project, check out this [article describing a nationwide collection of isolates by Gallegos-Casillas et al. 2023](https://www.biorxiv.org/content/10.1101/2023.07.02.547337v1), and this [paper with the first genomic draft of *Kz. humilis* by Garcia-Ortega et al. 2022](https://journals.asm.org/doi/10.1128/mra.01154-21).

### View the metadata

We will be working on sequences from isolates sampled from fermentation tanks used to produce different spirits at the seven producing regions. The metadata file associated with this lesson can be [downloaded directly here]([files/Ecoli_metadata_composite.csv](https://docs.google.com/spreadsheets/d/1mQFl-YEGzwSK77qzT8YL1FSzJ5kfe2ujvHBwST1JY5s/edit?usp=sharing) 

| Column           | Description                                     | 
| ---------------- | ----------------------------------------------- |
| ID               | strain name                                     | 
| Origin           | Where did the strain come from                  | 
| Spirit           | What was being produced with this strain        | 
| Read_type        | library type of reads                           | 


:::::::::::::::::::::::::::::::::::::::  challenge

### Challenge

Based on the metadata, can you answer the following questions?

1. How many different types of spirits exist in the data?
2. How many states are represented in this data?
3. Why do you think are there isolates from China, USA and the French Guiana?

:::::::::::::::  solution

### Solution

1. 8 different types of spirits
2. 11 different states
3. They will serve as reference

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::::: keypoints

- It is important to record and understand your experiment's metadata.

::::::::::::::::::::::::::::::::::::::::::::::::::




