---
title: Analysis of metabarcoding data with QIIME2
author: Sean Harrington, Josh Harrison, & Félix Brédoire
date: April 30, 2022
---


## Table of Contents


- [1. Introduction](#introduction)
- [2. Installation and setup](#installation-and-setup)
- [3. Read and process the data](#read-and-process-the-data)
- [4. Phylogenetic diversity analyses](#phylogenetic-diversity-analyses)


<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>
<br><br><br><br><br><br>


## 1. Introduction

Qiime2 is a program for metabarcoding analysis. SOME MORE STUFF HERE!!!

This tutorial is adapted from the following tutorials in the QIIME2 documentation: [“Atacama soil microbiome” tutorial](https://docs.qiime2.org/2022.2/tutorials/atacama-soils/), [“Moving Pictures” tutorial](https://docs.qiime2.org/2022.2/tutorials/moving-pictures/).

describe the dataset, etc. 


<br><br><br>


## 2. Installation and setup

The QIIME2 software is built from a lot of different pieces and dependencies, and so installation is done using conda or virtual machines. I find conda installation to be the easiest (both here and for many other programs). conda allows users to create separate environments that programs are installed into, which allows users to have multiple versions of the same program all installed on the same machine, each within a different conda environment. This can be very useful with QIIME because some versions of QIIME have different features. Full details of QIIME2 installation options can be found [here](https://docs.qiime2.org/2022.2/install/), but we'll stick with conda:

If you do not already have a conda install on your system, install miniconda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).


Then run the following to download the conda environment file and install QIIME from file into that environment at the same time. Note that most conda installs of other programs can be done from a remote source without first downloading a file.

```
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-osx-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-osx-conda.yml
rm qiime2-2022.2-py38-osx-conda.yml
```

To use QIIME now or at any point in the future, you will need to activate the `qiime2-2022.2` environment. If the environment is not active, your system will not know where to look for the QIIME2 software. 

Activate the environment and test the install:

```
conda activate qiime2-2022.2
qiime --help
```

This should spit out the help menu for QIIME if everything worked as expected.


You can deactivate the environment by running `conda deactivate` -- if you run it now, you will need to activate it again, though.

You can always check what conda environments you have and which you are currently in by running `conda env list`. The active environemnt will have `*` next to it.



Now that we have QIIME installed, let's download the data that we'll be working with. Create a new directory, move into that directory, and then download the data there. We'll work from within this directory for the duration of our QIIME analyses.

```
mkdir qiime2-atacama-tutorial
cd qiime2-atacama-tutorial

conda install wget
wget \
  -O "sample-metadata.tsv" \
  "https://data.qiime2.org/2022.2/tutorials/atacama-soils/sample_metadata.tsv"
  
mkdir emp-paired-end-sequences
wget \
  -O "emp-paired-end-sequences/forward.fastq.gz" \
  "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/forward.fastq.gz"
wget \
  -O "emp-paired-end-sequences/reverse.fastq.gz" \
  "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/reverse.fastq.gz"
wget \
  -O "emp-paired-end-sequences/barcodes.fastq.gz" \
  "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/barcodes.fastq.gz"
  
```

Note that the slashes are just escaping line endings so that `wget` commands are interpreted as a single-line commands while allowing them to easily fit on this page.

We should now have the sample metadata, raw reads, and barcodes for the samples.

Take a look at the sample metadata in the `sample_metadata.tsv` file to get familiar with what we're working with here. We have a mixture of categorical variables, such as site-name, and continuous variables, such as elevation.




## 3. Read and process the data

###########
########### I need to add in more stuff about each of these steps below
###########     qza vs qzv files - qiime tools view for qzv files

Now we can import the data and prepare it for analysis in Qiime. 


```
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path emp-paired-end-sequences \
   --output-path emp-paired-end-sequences.qza
```

This creates a QIIME artifact with the sequences. 

We now need to demultiplex the sequences, which sorts the reads into samples based on barcodes added in the library preparation.

```
qiime demux emp-paired \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs emp-paired-end-sequences.qza \
  --o-per-sample-sequences demux-full.qza \
  --o-error-correction-details demux-details.qza
```


To make things run faster, we'll subsample only 30% of the reads. *This is not a general step for most pipelines*. You will typically want to use all of your data unless you have some specific reason for subsampling.

```
qiime demux subsample-paired \
  --i-sequences demux-full.qza \
  --p-fraction 0.3 \
  --o-subsampled-sequences demux-subsample.qza

qiime demux summarize \
  --i-data demux-subsample.qza \
  --o-visualization demux-subsample.qzv
```

#### Bring the qvz into qiime view

We have a lot of samples that have fewer than 100 reads in them. This isn't really enough reads for meaningful analysis, so let's filter these out. You do not necessarily need to run through this step in other analyses. Think carefully about what thresholds you may want to use when filtering out any data.

```
qiime tools export \
  --input-path demux-subsample.qzv \
  --output-path ./demux-subsample/

qiime demux filter-samples \
  --i-demux demux-subsample.qza \
  --m-metadata-file ./demux-subsample/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux demux.qza
```

Denoising the data: - add more stuff here, error correction - note there are other options

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```


Generate summaries of the output stats from that:


```
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```


## 4. Phylogenetic diversity analyses

We need to start by making a tree:

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```


Next, we can use that to calculate some diversity metrics. We have to specify a minimum sampling depth (`--p-sampling-depth`). This threshold is important: any samples with more than the specified depth will be randomly sampled to contain only that many features, whereas samples with fewer features will be discarded from analysis. This is necessary because using different numbers of features across samples can bias estimates of diversity. There is no easy rule for selecting this threshold, other than that we want to try to select a threshold that maximizes the number of features without dropping out too many samples. Let's look at the *Interactive Sample Detail* section of `table.qzv` to help us figure out what threshold to use.

```
qiime tools view table.qzv
```

In this document, click on the *Interactive Sample Detail*  tab up top. As you move the *Sampling Depth* slider around you can see a visual representation of how many samples will be retained. You can also scroll down the table to see if there is a point at which there is a sharp decline in the Feature Counts, the value just prior to such a drop off can be good to use.

What do you think is a reasonable value to use?


```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 733 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
```

This generates a lot of output that is worth looking through. 

```
qiime tools view core-metrics-results/bray_curtis_emperor.qzv
```
For continuous variables, try selecting a color scheme from the sequential or diverging sets of colors, these should make it easier to identify trends.

Which variables seem to be most strongly associated with beta diversity?


### 4.1 Categorical tests

Now that we have computed the diversity metrics and done some qualitative exploration of the emperor plots, we can test for associations between these diversity metrics and the sample metadata. We can test for differences among categorical groups (e.g., using the `vegetation` column, which is yes/no) or correlations with continuous variables, such as `percentcover`. 

Let's start with a categorical test with the evenness and Faith Phylogenetic Diversity metrics. 


```
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
  

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv  
```

Take a look at these results:

```
qiime tools view core-metrics-results/evenness-group-significance.qzv

qiime tools view core-metrics-results/faith-pd-group-significance.qzv
```

What associations are statistically significant?


### I need to put in more interpretation here

Now let's look at how beta diversity composition varies across categorical variables using PERMANOVA. This tests if distances between samples within a group are more similar to each other than to samples from other groups. Because it uses permutation to assess significance, this command can be slow, and so we will only run it on the vegetation column of our metadata right now. The addition of the `--p-pairwise` option will perform pairwise tests to determine which groups are significantly different from each other--note that this is redundant here because we only have two groups for vegetation.

```
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization core-metrics-results/unweighted-unifrac-vegetation-significance.qzv \
  --p-pairwise
```

Take a look:

```
qiime tools view core-metrics-results/unweighted-unifrac-vegetation-significance.qzv
```


Are there any other categorical comparisons that are worth making?


### 4.2 Continuous tests

We can also make correlations between diversity and continuous variables from the metadata. As above, we'll start with alpha diversity:

```
qiime diversity alpha-correlation \
	--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-method spearman \
	--o-visualization core-metrics-results/faith-pd-correlation.qzv  
```


View it:

```
qiime tools view core-metrics-results/faith-pd-correlation.qzv
```

What relationships are significant? Note that we can also run this correlation using evenness or with a Pearson correlation test instead of Spearman.


Let's now look for significant correlations using beta diversity.

To do this, we need to first calculate a distance matrix from the column of interest from the metadata. Here we'll use just elevation, but you are free to explore other variables as well. 


```
qiime metadata distance-matrix \
	--m-metadata-file sample-metadata.tsv \
	--m-metadata-column elevation \
	--o-distance-matrix core-metrics-results/elevation-dist-mat.qza
```


Then we can use a Mantel test to test for an association between this distance matrix and one of our metrics of beta diversity. We'll use just unweighted unifrac distance, but we can do this with any.


```
qiime diversity mantel \
	--i-dm1 core-metrics-results/elevation-dist-mat.qza \
	--i-dm2 core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--p-method spearman \
	--p-label1 elevation \
	--p-label2 unweighted_unifrac \
	--p-intersect-ids \
	--o-visualization core-metrics-results/unweight_unifrac_elevation_mantel.qzv
```

View the resulting visualization:


```
qiime tools view core-metrics-results/unweight_unifrac_elevation_mantel.qzv
```

look into the `qiime diversity bioenv` - don't know what that does


### 4.3 Alpha rarefaction





**Stuff to look up later**


Figure out why intersect IDs is necessary:

	Plugin error from diversity:

	  The following ID(s) are not contained in both distance matrices. This sometimes 		occurs when mismatched files are passed. If this is not the case, you can use `intersect_ids` to discard these mismatches and apply the Mantel test to only those IDs that are found in both distance matrices.

	  BAQ1370.1.2, BAQ1370.1.3, BAQ1370.3, BAQ1552.1.1, BAQ1552.2, BAQ2838.1, BAQ2838.2, BAQ2838.3, BAQ895.2, BAQ895.3, YUN1005.1.3, YUN1005.2, YUN1242.2, YUN1609.3, YUN2029.1, YUN2029.3, YUN3008.1.1, YUN3008.1.2, YUN3008.1.3, YUN3008.2, YUN3008.3, YUN3153.1, YUN3184.1, YUN3184.2, YUN3259.1.1, YUN3259.1.3, YUN3346.2

