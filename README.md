# Supporting code for manuscript:
### "dsRNA-based viromics: A novel tool unveiled the hidden soil viral communities"

#### This github repository contains the following files:

1- [Taxonomy Extraction](https://github.com/poursalavati/SOVAP_Soil_dsRNA/blob/main/Taxonomy.md): Extract taxonomy information for IMG/VR Blast result   
  
2- [Prepare data for analysis in R](https://github.com/poursalavati/SOVAP_Soil_dsRNA/blob/main/Preparation.R): Manipulate and prepare necessary data and tables to plot taxonomy results  
  
3- [Plotting Stacked Bar Charts and HeatMaps](https://github.com/poursalavati/SOVAP_Soil_dsRNA/blob/main/StackedBarcharts_HeatMaps.R): Creates two functions that generate stacked bar chart and heatmaps for each taxonomic level  
  
4- [Diversity Indices Calculation](https://github.com/poursalavati/SOVAP_Soil_dsRNA/blob/main/Diversity_Indices.R): Calculates richness and diversity indices (shannon and simpson) for each taxonomic level and Violin Plot   
  
5- [Beta diversity and PCoA](https://github.com/poursalavati/SOVAP_Soil_dsRNA/blob/main/Beta_and_PCoA.R): Calculate Bray-Curtis dissimilarity matrix and PCoA plot  
  
#### Folders:  
1- [Abundance_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/Abundance_Results): This folder contains the abundance results of the SOVAP pipeline for all datasets  
  
2- [RdRp_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/RdRp_Results): This folder contains the RdRp Sequences from all datasets  
- 2.1- [RdRp_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/RdRp_Results/RdRp-Scan): This folder contains the RdRp Analysis Results with RdRp-Scan approach: alignments and Newick files
- 2.2- [RdRp_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/RdRp_Results/WSH): This folder contains the RdRp Analysis Results with WSH approach: alignments and Newick files  

3- [Plasmid_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/geNomad_Plasmid): This folder contains the plasmid results of the geNomad tool for all datasets  (final.contigs_plasmid_summary.tsv)
  
4- [Virus_dir](https://github.com/poursalavati/SOVAP_Soil_dsRNA/tree/main/geNomad_Virus): This folder contains the plasmid results of the geNomad tool for all datasets  (final.contigs_virus_summary.tsv)

#### If you have any other questions or concerns related to the code or analysis, please create a new issue on Github.
#### If you find this code helpful, please cite:  
- Manuscript DOI  
- Poursalavati A. (2023). SOVAP v.1.3 : Soil Virome Analysis Pipeline (1.3). Zenodo. https://doi.org/10.5281/zenodo.7700081
