# Pan-Can Metabolism Analysis

Input data and informatic pipeline, analysis, visualization R code for the analysis of metabolomics data from over 900 tissue samples spanning 7 cancer types. 

# Code Structure
* analysis: R code for analyses performed in the project. Scripts prefixed by Run* are primary scripts used in the project analysis. Some of the key scripts are described below:
    * RunVisualizeClinical.R: Visualize clinical data
    * RunSummary: Generate a summary of the merged metabolomics data
    * RunKeggCoverage: Calculate the coverage of KEGG represented by the project dataset 
    * RunPathwayCommonsCoverage: Calculate the coverage of Pathway Commons represented by the project dataset 
    * RunPharmaCoverage: Determine enzymes making use of metabolites (as substrates or products) in the study are targetable by drugs
* data: Primary data that was processed as part of the project 
    * merged_metabolomics: Merged datasets after informatic standardization (e.g. name mapping across studies) pipeline, including metabolic profiling values, clinical features, and paired tumor-normal fold changes
    * studies: Primary data collected from studies (often supplementary tables) both metabolomic profiling and clinical variable data
    * pharmacology_drugbank: CHEBI IDs from Pathway Commons (pathwaycommons.org) dataset categorized using DrugBank drug categories
* import: Scripts for importing data (e.g. KEGG pathway data, metabolite ID mapping)
	* RunMergeMetabolomics: Merges the individual metabolomics files, 
	* RunManualImport: Imports metabolomics data from a few studies that do not have normalized data
    * RunMap2Kegg: Maps our metabolite data to KEGG IDs
* results: Results from analysis 
* shinyapp: Data used for the R Shiny web application 
