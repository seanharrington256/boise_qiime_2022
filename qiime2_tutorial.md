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
- [5. Taxonomic analysis](#taxonomic-analysis)
- [6. Differential abundance testing with ANCOM](#differential-abundance-testing-with-ancom)



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

QIIME2 is a pipeline for metabarcoding analysis described in [this paper](https://www.nature.com/articles/s41587-019-0209-9). It is a piece of software that wraps around several other programs via "plugins" in QIIME2 lingo. For example, in the error correction or "denoising" step, one option is to use the `dada2` plugin to run [DADA2](https://benjjneb.github.io/dada2/). DADA2 is incorporated into the QIIME2 pipeline, but is not developed or maintained by the same group that develops and maintains QIIME2. Users can develop their own plugins and contribute them to the QIIME2 pipeline to extend functionality as new methods are developed.

The greatest advantage of QIIME2 is that it aggregates various tools into a single pipeline that uses a common grammar. This includes unifying data formatting so that users do not need to worry about complex file conversions to prep output from one program for input into another, which I personally find to be one of the most obnoxious hassles in bioinformatic analysis.

Another advantage of QIIME2 is that the files track the provenance of the data--i.e., QIIME2 tracks exactly how the data has been processed at every step up to any given analysis. We will explore this further as we start looking at QIIME2 files.


Note that hroughout this tutorial, I may simply refer to QIIME2 as QIIME, but in all cases, we are referring to QIIME2


Today, we will use QIIME2 to explore soil microbiome data from the Atacama desert in Chile. The data are presented in [this paper](https://journals.asm.org/doi/full/10.1128/mSystems.00195-16).

This tutorial is adapted from the following tutorials in the QIIME2 documentation: [“Atacama soil microbiome” tutorial](https://docs.qiime2.org/2022.2/tutorials/atacama-soils/), [“Moving Pictures” tutorial](https://docs.qiime2.org/2022.2/tutorials/moving-pictures/), [Training feature classifiers](https://docs.qiime2.org/2022.2/tutorials/feature-classifier/). The QIIME2 documentation is extensive and features several tutorials. More detail on all of the concepts that we will cover today can be found in this documentation. 

Blocks of code to be entered into your command line terminal will look like this:

```
# This is a code block
```

The `#` denotes a comment in the bash language, and anything following that will not be interpreted by the system.



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



GO OVER what plugins are and the fact that a lot of this wraps other pieces of software (https://docs.qiime2.org/2022.2/citation/)

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

Alpha rarefaction plotting will compute alpha diversity at multiple sampling depths and then visualize how increased sampling affects estimates of alpha diversity. 



```
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1175 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```


The parameter `--p-max-depth` should be chosen from the information in `table.qzv` that we generated previously. Using roughly the median frequency is generally recommended, but you may want to increase the value if your rarefaction plot does not level out or decrease it if you are losing many of your samples at this depth. We've chosen the median above, let's see how it looks:


```
qiime tools view alpha-rarefaction.qzv
```


The top plot shows us how much diversity we detect at varying sequencing depths. If the plot has leveled off, we can be reasonably confident that we are accurately characterizing the diversity in our samples. If the plots do not level off, that suggests that further sequencing may be needed to detect additional features in the samples and accurately characterize diversity.

How does this plot look to you? Do you think that the sequencing depth is adequate?

The bottom plot shows the number of samples that remain in each category when grouping by metadata columns. This is important to look at because if the diversity metric in the top plot is calculated from very few samples at a given sampling depth, that estimate of diversity may be unreliable.



## 5. Taxonomic analysis


Next up, we'll start to explore the taxonomic composition of our samples. We will start by training a taxonomic classifier. There are several existing, pre-trained classifiers that exist for QIIME, but the developers recommend training your own classifier, as they perform best when trained on your specific data.

We'll follow the documentation [here](https://docs.qiime2.org/2022.2/tutorials/feature-classifier/) to train our classifier. As stated in that documentation, we'll use the Greengenes 13_8 85% OTU data set -- **note that is not recommended for real data** -- we are using it here for the sake of time efficiency.


Start by creating a new directory, moving into it, and then downloading the data that we will use to train the classifier.

```
mkdir training-feature-classifiers
cd training-feature-classifiers
wget \
  -O "85_otus.fasta" \
  "https://data.qiime2.org/2022.2/tutorials/training-feature-classifiers/85_otus.fasta"

wget \
  -O "85_otu_taxonomy.txt" \  "https://data.qiime2.org/2022.2/tutorials/training-feature-classifiers/85_otu_taxonomy.txt"
```

If that all ran successfully, you should now have three files in your current directory, a fasta, txt, and qza file.


Then we need to import these data as QIIME artifact files.

```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 85_otus.fasta \
  --output-path 85_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 85_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
```


PUT IN MORE DETAIL HERE ABOUT THE CLASSIFIER:

Next up, we'll extract the reads from the reference sequences to match our data:


```
qiime feature-classifier extract-reads \
  --i-sequences 85_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza
```


Then we can run the Naive Bayes classifier.

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza
```

Now we can run the classifier on our data and then generate a visualization.


```
cd ..  # We first need to move up one level to get out of the classifier directory

qiime feature-classifier classify-sklearn \
  --i-classifier training-feature-classifiers/classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```


Let's take a look.


```
qiime tools view taxonomy.qzv
```

This isn't too informative or at all pretty to look at, so let's make some bar plots:

```
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv


qiime tools view taxa-bar-plots.qzv
```


Set the `Taxonomic level` to *Level 2* and then `Sort Samples By` *vegetation*. Are there any noticeable differences in the phyla that are represented by samples from each of the vegetation categories?



## 6. Differential abundance testing with ANCOM

We can also perform explicit tests for features that are differentially abundant across sample groups using [ANCOM](https://pubmed.ncbi.nlm.nih.gov/26028277/). 


To do this, we need to use `add-pseudocount` to generate a `FeatureTable[Composition]` QIIME artifact. This is necessary because ANCOM cannot handle frequencies of zero, and so these values must be imputed prior to analysis with ANCOM.

```
qiime composition add-pseudocount \
  --i-table table.qza \
  --o-composition-table comp-table.qza
```

Then we'll run ANCOM on the `vegetation` column to see if we have any significant differential abundance across vegetation categories.

```
qiime composition ancom \
  --i-table comp-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization ancom-vegetation.qzv
```

Visualize this:

```
qiime tools view ancom-vegetation.qzv
```

NEED a lot more interpretation here 

We can also run this at the genus level by collapsing down the taxonomy THIS explanation sucks


``` 
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-l6.qza

qiime composition add-pseudocount \
  --i-table table-l6.qza \
  --o-composition-table comp-table-l6.qza

qiime composition ancom \
  --i-table comp-table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization l6-ancom-vegetation.qzv
```

Take a look.


```
qiime tools view l6-ancom-vegetation.qzv
```


We can also do the same for phyla:

``` 
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table table-l2.qza

qiime composition add-pseudocount \
  --i-table table-l2.qza \
  --o-composition-table comp-table-l2.qza

qiime composition ancom \
  --i-table comp-table-l2.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization l2-ancom-vegetation.qzv
  
  
  
qiime tools view l2-ancom-vegetation.qzv
```


What phyla differ among the vegetated and non-vegetated sites.


**Stuff to look up later**


May want to add in Gneiss - at least reference it


Figure out why intersect IDs is necessary:

	Plugin error from diversity:

	  The following ID(s) are not contained in both distance matrices. This sometimes 		occurs when mismatched files are passed. If this is not the case, you can use `intersect_ids` to discard these mismatches and apply the Mantel test to only those IDs that are found in both distance matrices.

	  BAQ1370.1.2, BAQ1370.1.3, BAQ1370.3, BAQ1552.1.1, BAQ1552.2, BAQ2838.1, BAQ2838.2, BAQ2838.3, BAQ895.2, BAQ895.3, YUN1005.1.3, YUN1005.2, YUN1242.2, YUN1609.3, YUN2029.1, YUN2029.3, YUN3008.1.1, YUN3008.1.2, YUN3008.1.3, YUN3008.2, YUN3008.3, YUN3153.1, YUN3184.1, YUN3184.2, YUN3259.1.1, YUN3259.1.3, YUN3346.2

