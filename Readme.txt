###################################################################### Readme File ###################################################################

Requirements:
1. Ubuntu OS
2. Python version 3.x 


To run the python scripts and generate output you need to run the shell script (Run_Codes.sh) by the following command:
sh Run_Codes.sh

The 'DEG_NonDEG_Module_Pathway.txt' file is the final output file containing the DEG count, LDEG count, and the percentages of DEGs and DEGs+LDEGs in each of the Modules.

For Acute and Chronic Wound Infections:
'UpRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the elevated LDEGs identified in Modules enriched in DEGs.
'DownRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the suppressed LDEGs identified in Modules enriched in DEGs.
'PAO1_DEGMod_Enrichment.txt': This file contains the Fold Change and P-value (Fisher test) for identification of Modules enriched in DEGs.

For Antibiotic Resistance:
'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the elevated and suppressed LDEGs identified in Modules enriched in DEGs.
'PAO1_DEGMod_Enrichment.txt': This file contains the Fold Change and P-value (Fisher test) for identification of Modules enriched in DEGs.
