# inflammatome_package_sandbox

R package accompanying the manuscript "The inflammatome: A meta-analysis of human genes regulated during inflammation"

inflammatomeR is a collection of functions that can be used to assess the presence of inflammation in transcriptomics or proteomics data. 

There are two main use case scenarios, requiring different input data:
- Gene set enrichment analysis: inflammatomeR can be used to assess the enrichment of a core set of inflammation-related genes within a list of results from differential expression analysis
- Volcano plot visualisation of the results list, marking the inflammatome genes; if many of the inflammatome genes are among the top hits, one can choose to interpret the results as being driven 
by inflammation; if the focus of the study is to understand more disease-specific processes, the inflammatome genes can be filtered out e.g. when selecting genes for further validation, biomarkers, drug targets, etc  

- Inflammatome score: inflammatomeR can be used to calculate an inflammation score for each sample in a dataset, which represents the level of inflammation in that sample, relative to the other samples;  
