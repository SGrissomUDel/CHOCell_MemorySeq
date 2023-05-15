# CHOCell_MemorySeq
RNA-Seq data processing pipeline and MemorySeq analysis for identifying stress-resistance associated heritable gene states including heritable gene state identification, GO enrichment analysis, and differential gene expression analysis

Folder structure is broken down into 3 main folders:
1. Scripts used in HCP Linux cluster for processing raw RNA-Seq data to produce a combined count table for batches of RNA-Seq data. For example, used to produce the combined tables for MemorySeq samples, Noise Control samples, and Fed-Batch samples
2. R scripts used for MemorySeq analysis. Scripts are written to complete multiple tasks in the same file. Main scripts include:
    a. MemorySeq Analysis
        1. Melt and filter MemorySeq and Noise control reads for each samples
        2. Determine summary parameters for each sample and each gene
        3. Perform Poisson regression fitting to genes and determine top 2% in residuals for heritable gene identification
        4. Plot data such as Coefficient of Variation (CV) of transcripts per million (TPM) vs log2 mean TPM, histogram residuals, ratio of CVs, histogram of             gene expresison.
        5. Calculating Pearson pairwise correlation coefficiencys and developing a heatmap
        6. Performing Cooks Distance analysis for detecting significant outliers in the replicates
        7. Performing k-clique community network identification
    b. GO Enrichment Analysis
        1. Perform GO Enrichment for 1 of 3 different inputs (List of gene IDs, list of gene names, and output from differential gene expression analysis)
        2. Produce csv and pdf files describing the Fisher/Kolmogorov-smirnov (classic, weighted, and elim) for the enrichment
        3. Allows searching genes for those with specific GO Terms
     c. Differential gene expression analysis
        1. Load data in and define the experimental design space
        2. Transform data and visualize the transformation
        3. Run principal component analysis for the conditions tested
        4. Run differential gene expression analysis and isolate these genes
        5. Demonstrate overlap of differentially expressed genes between conditions and between heritable genes
        6. Produce transcriptomic heat maps
3. Sundry R scripts for plotting
    
