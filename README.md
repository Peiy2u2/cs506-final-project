## Project Description
Patch-seq is a powerful technique that integrates single-cell transcriptomics with neuronal morphology and electrophysiology, enabling a comprehensive understanding of neuronal diversity. While single-cell RNA sequencing has revolutionized neuroscience by uncovering gene expression patterns across individual neurons, the complexity of neuronal function often requires a multi-modal approach. 

By combining Patch-seq data with single-cell RNA sequencing and embedding them in a two-dimensional space, we bridge the gap between electrophysiological measurements—traditionally used by neuroscientists to study neuronal activity—and transcriptomic profiles, which provide insights into gene expression. This integration allows for a more holistic characterization of neuronal subtypes, linking their electrical properties to molecular identities and ultimately advancing our understanding of brain function.

## Motivations
The primary motivation is  to determine what types of cells those 36 are. In a bigger scheme, this would help bridge the gap between the electrophysiological measurements that neuroscientists use with the transcriptomic measurements that we rely on in genomics. 
For cells that are biophysically and transcriptionally homogeneous, there may be a developmental change in their transcriptomic state. And transcriptomically similar cells may also differ in terms of morphology and electrophysiology.  ([Reference](https://www.jneurosci.org/content/41/5/937))

We probably can find something intriguing from our 36 cells data when embedded into the space of ~10k single cells. It would be interesting if two cells that belong to two totally different cell types A and B are very close to each other in the embedding space. 
And likewise, if the distance of Cell A - another Cell A is further than the distance of Cell A - Cell B, this may indicate those cells A are in different transcriptomic states though they still have similar biophysical properties. 


## Data Collection
Data were collected from different regions of brains (including prefrontal cortex) of Macaca Mulatta (let's say monkeys). 
* Data availability
	* RNA-seq counts data of 36 cells from different regions of 7 (4F + 3M) individuals, with counts data of 25,432 genes: `data/featureCounts_matrix.tsv`
	* RNA-seq counts data of 27,964 single cells from 2 individuals, with  of 19,605 genes (download the file from this [link](https://drive.google.com/file/d/17wdUNMDrk_AGIX44MOa9itZobufrxgy8/view?usp=drive_link)): `data/seurat_merge_processed.csv`

## Data Cleaning
- **Dimension allignment**: 
	- 19,605 (geneids + genenames) in snRNA-seq counts data
	- 35,432 geneids in patch-seq counts data
    - Align the dimensions of the two datasets by removing unique genes/geneids in each dataset and keeping the common ones. 

- **Normalization**: 
    - After dimension alignment, normalize the data by summing up gene expression counts within each cell, and divide each count data count by the sum.
    - This step is to make the data comparable across different cells.

- **Batch effect checking**: 
    - Since we have two different batches of snRNA-seq, we need to check batch effect and if necessary, do batch corrections
	- After normalization, perform PCA on the normalized gene expression data, and draw the PC1-PC2 figure. Color the PCA plot by batch. If the batches cluster separately, it suggests a batch effect (different batches with different colors). 

## Feature Extraction
- K-means clustering
    - Perform PCA on the normalized gene expression data.
    - Perform K-means clustering on the PCA results.
- Graph embedding method (e.g. deepwalk, node2vec)
    - Perform graph embedding method on the normalized gene expression data.
    - Perform PCA on the vector data from the graph embedding method.

## Data Visulization
- K-means clustering plot
    - Draw the K-means clustering plot with PC1 and PC2 as the x and y axes.
- Graph embedding plot
    - Draw the graph embedding plot with with PC1 and PC2 as the x and y axes.

## Model Training
- Train the K-means model to achive the optimal number of clusters by referring to the regions and types of neuron cells. 
- Likewise, train the graph embedding model to achive the optimal layers or imporant hyperparameters. 


