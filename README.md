We have taken the following steps to clarify our methodology:
1.	Data acquisition: Acquired the datasets from the GEO repository, using their respective accession numbers.
2.	Preprocessing:  The raw microarray data were processed using the robust multi-array average (RMA) method for normalization and background correction.
3.	Selection of Genes: We have now explicitly stated the conditions between which differential expression was analyzed (e.g., latent tuberculosis infection versus active infection). The method used for identifying differentially expressed genes (DEGs) has been clearly outlined, specifying that we used [e.g., limma, DESeq2] for statistical analysis. We applied a threshold of p < 0.05 and a fold-change > 2 to select significant DEGs, resulting in the identification of 12,256 DEGs.
4.	Machine Learning methodology: 
3.1.	Clustering Methodology: We clarified the clustering process, explaining the choice of the initial 20 clusters and the subsequent reduction to 8. This was done using [e.g., k-means or hierarchical clustering] based on similarities in gene expression patterns. The reduction to 8 clusters was determined by optimizing cluster cohesion and separation metrics. . this process of refining the DEGs from 12,256 to 7,610 is now thoroughly explained. This reduction was achieved by applying additional filters such as removing genes with low expression variability across the datasets, as well as further refining based on specific biological relevance to latent tuberculosis infection.
3.2.	Gene classification (Network analysis): We have outlined the steps taken to derive the final set of 305 genes, including any additional filters or criteria used in this process, as well as how this number was further reduced to 200 genes.
5.	Generating direct miRNA: We described the process of generating miRNA clusters and their relevance to LTBI through network and functional analyses.
6.	Hypothesis validation: We validated the identified gene clusters by integrating them into functional hypotheses relevant to latent tuberculosis infection.


By addressing these points, we aim to provide a comprehensive and clear description of our analysis pipeline, enhancing the reproducibility and understanding of our methods. 
