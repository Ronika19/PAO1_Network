################################## Readme File ##############################

Requirements:
1. Ubuntu OS
2. Python version 3.x 


To run the python scripts and generate output you need to run the shell script (Run_Codes.sh) by the following command:
sh Run_Codes.sh

Input Files for Acute and Chronic Wound Infections:
1. 'KEGG_PAO1_AllPathways.txt' file contains the pathway annotations for PAO1
2. 'Module_Assignment_PAO1_WGCNA.tsv' file contains the module assignment of PAO1 genes involved in the network.
3. 'Murine_Dataset_UpRegulatedDEGs.txt' and 'Murine_Dataset_DownRegulatedDEGs.txt' files contain the upregulated and downregulated DEGs respectively.
4. 'PAO1_Annotations.tsv' file contain the Locustags and their corresponding Transcript ids.
5. 'PAO1_Nodes.txt' file contain the nodes (transcript ids) of PAO1 involved in the network.
6. 'pgen.1004518.s005.txt' file contains the fold changes and p-values for all PAO1 genes.
7. 'Pseudomonas_aeruginosa_PAO1.gene_info' file containd the Entrez ids and Locustags for PAO1 genes.

Input Files for Antibiotic Resistance:
1. 'KEGG_PAO1_AllPathways.txt' file contains the pathway annotations for PAO1
2. 'Module_Assignment_PAO1_WGCNA.tsv' file contains the module assignment of PAO1 genes involved in the network.
3. 'UpRegulated_AZM_Vs_SERUM.txt' and 'DownRegulated_AZM_Vs_SERUM.txt' files contain the upregulated and downregulated DEGs respectively.
4. 'PAO1_Annotations.tsv' file contain the Locustags and their corresponding Transcript ids.
5. 'PAO1_Nodes.txt' file contain the nodes (transcript ids) of PAO1 involved in the network.
6. 'AZM_SERUM_Paper.txt' file contains the fold changes and p-values for all PAO1 genes.
7. 'Pseudomonas_aeruginosa_PAO1.gene_info' file containd the Entrez ids and Locustags for PAO1 genes.

Input Files for Cystic Fibrosis:
1. 'KEGG_PAO1_AllPathways.txt' file contains the pathway annotations for PAO1
2. 'Module_Assignment_PAO1_WGCNA.tsv' file contains the module assignment of PAO1 genes involved in the network.
3. 'CF_Vs_LB_Significant.txt' file contain the upregulated and downregulated DEGs respectively.
4. 'PAO1_Annotations.tsv' file contain the Locustags and their corresponding Transcript ids.
5. 'PAO1_Nodes.txt' file contain the nodes (transcript ids) of PAO1 involved in the network.
6. 'CF_Vs_LB.txt' file contains the fold changes and p-values for all PAO1 genes.
7. 'Pseudomonas_aeruginosa_PAO1.gene_info' file containd the Entrez ids and Locustags for PAO1 genes.

Output Files:
The 'DEG_NonDEG_Module_Pathway.txt' file is the final output file containing the DEG count, LDEG count, and the percentages of DEGs and DEGs+LDEGs in each of the Modules.

For Acute and Chronic Wound Infections:
'UpRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the elevated LDEGs identified in Modules enriched in DEGs.
'DownRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the suppressed LDEGs identified in Modules enriched in DEGs.
'PAO1_DEGMod_Enrichment.txt': This file contains the Fold Change and P-value (Fisher test) for identification of Modules enriched in DEGs.

For Antibiotic Resistance and Cystic Fibrosis:
'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the elevated and suppressed LDEGs identified in Modules enriched in DEGs.
'PAO1_DEGMod_Enrichment.txt': This file contains the Fold Change and P-value (Fisher test) for identification of Modules enriched in DEGs.
