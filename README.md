# BioID-Network-Analysis
The following contains the code used to develop the networks models for the following publication: "A protein-proximity screen reveals that Ebola virus co-opts the mRNA decapping complex through the scaffold protein EDC4"
![Representation of the processing steps used to assemble the virus-host protein interaction networks used in this study](/BioID-Network-Analysis/Pipeline Inputs/
# Methods Overview
BioID2 was used to identify putative interacting host proteins of six Ebola virus structural proteins. These host proteins were assembled into multiple network models, including a traditional hub and spoke protein-protein interaction network and Prize-Collecting Steiner Forest (PCSF) network, assembled using the Human Integrated Protein-Protein Interaction rEfernece (HIPPIE) interactome. The PCSF network was used to identify protein complexes present in the BioID hit lists and guide selection of hits for follow-up.
## Software dependencies
Network analysis and gene-set enrichment analyses were performed locally and through APIs to Cytoscape and EnrichR
- R Version 4.3.1
- BiocManager Version 1.30.22
- devtools Version 2.4.5
- igraph Version 1.5.1
- RCy3 Version 2.20.1
- PCSF Version 0.99.1
- ndexr Version 1.22.0
- readxl Version 1.4.3
- tidyverse Version 2.0.0
- tidyr Version 1.3.0
- enrichR Version 3.2
- xlsx Version 0.6.5
- Cytoscape Version 3.10
## Hub and Spoke Network assembly
Using the script "HnS_for_cytoscape_script.R" an edge list from Extended data table 2 was uploaded to R. Using iGraph, a network was assembled from these edges and exported to Cytoscape using the RCy3 package. Further aesthetic changes were then made through the Cytoscape GUI
## PCSF Network assembly
PCSF was installed in R from the github repository "IOR-BioinformaticsPCSF" found here: https://github.com/IOR-Bioinformatics/PCSF using devtools
PCSF network mapping, network annotation, and gene-set enrichment analysis were performed according to the script "PCSF Network Analysis.R". The HIPPIE interactome was pulled from the NdexBio Repository using ndexr, and assembled into an interactome using PCSF. A prize file of scored BioID hits was then uploaded to R from Extended data table 2 and used to create a Prize File for PCSF analysis. PCSF was run as described in the script. 
### Network annotation
Viral protein interaction rescores were uploaded from Extended data table 3, and added to the PCSF final network by adding the rescores as numerical data attributes to the genes within the PCSF network using the igraph package. This attributed network was then exported to Cytoscape using the RCy3 package. Network annotation was then completed using the "Style" section and Image/Chart tab.
### Gene-set enrichment analysis
Network gene-set enrichment was performed using the code listed in "PCSF Network Analysis.R" script under the "Perform gene set enrichment for network analysis" section. This code was adapted from the source code for enrichment_analysis() function in the PCSF package (function is now deprecated). The network was clustered using edge_betweenness community detection in iGraph. Genelists were then submitted for gene-set enrichment against Gene Ontology, KEGG, and Reactome databases using the "enrichR" package. Gene-set enrichments were exported to an excel table (see Extended data table 4). Network enrichment was manually annotated onto the PCSF network using the "Style" section in Cytoscape.  


