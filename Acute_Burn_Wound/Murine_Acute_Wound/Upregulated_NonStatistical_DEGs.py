import re
import sys
from .DataProcess import Data_Processer
from .DataProcess import *
from .DEG_Module_Enrichment import DEG_Enrichment

class NonStatistical_DEG:
	def UpDown_NonStatistical_DEGs(self, infile1, infile2, infile3, infile4, infile5, infile6, out_file1, out_file2, out_file3):
		DEG_Enrichment().Enrichment('Input_Files/PAO1_Nodes.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/PAO1_DEG_Modules.txt', 'Output_Files/PAO1_DEGMod_Enrichment.txt')
		dicts1 = Data_Processer().data_extract(infile1, 1)
		module, observed, expected, foldchange, fishertest = [],[],[],[],[]
		stats_module, stats_foldchange, stats_fishertest = [],[],[]
		modules, observed_degs, expected_degs, foldschange, fisher_test = dicts1['arr_0'], dicts1['arr_1'], dicts1['arr_2'], dicts1['arr_3'], dicts1['arr_4']
		for i in range(len(foldschange)):		
			if (float(foldschange[i]) > 1): # DEG Enriched Modules with Fold Change is greater than 1
				module.append(modules[i])
				observed.append(observed_degs[i])
				expected.append(expected_degs[i])
				foldchange.append(foldschange[i])
				fishertest.append(fisher_test[i])
				if (float(fisher_test[i]) <= 0.05): # DEG Enriched Modules with p-value is less than equal to 0.05
					stats_module.append(modules[i])
					stats_foldchange.append(foldschange[i])
					stats_fishertest.append(fisher_test[i])

		f = open(infile2,'r')
		m=0; mod, transcript = [],[];
		for line in f:
			split_line = (line.rstrip()).split('\t')
			mod.append(str(m))
			transcript.append(split_line)
			m += 1

		dicts2 = Data_Processer().data_extract(infile3, 0)
		locusid_down, transcriptid_down, genename_down, foldchanges_down, pvalue_down, koid_down, cogid_down, cogcategory_down =  dicts2['arr_0'], dicts2['arr_1'], dicts2['arr_2'], dicts2['arr_3'], dicts2['arr_4'],  dicts2['arr_5'], dicts2['arr_6'], dicts2['arr_7']

		dicts3 = Data_Processer().data_extract(infile4, 1)
		locusid, transcriptid, genename, foldchanges, pvalue, koid, cogid, cogcategory = dicts3['arr_0'], dicts3['arr_1'], dicts3['arr_2'], dicts3['arr_3'], dicts3['arr_4'],  dicts3['arr_5'], dicts3['arr_6'], dicts3['arr_7']

		dicts4 = Data_Processer().data_extract(infile5, 1)
		locus_id, gene_name, fold_change, p_value, ko_id, cog_id, cog_category = dicts4['arr_0'], dicts4['arr_1'], dicts4['arr_2'], dicts4['arr_3'], dicts4['arr_6'], dicts4['arr_7'], dicts4['arr_8']

		dicts5 = Data_Processer().data_extract(infile6, 1)
		locus_tag, transcript_tag = dicts5['arr_6'], dicts5['arr_7']

		outfile1 = open(out_file1,'w')
		outfile1.write('Module'+'\t'+'Transcript_Id'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
		outfile2 = open(out_file2,'w')
		outfile2.write('Module'+'\t'+'Transcript_Id'+'\t'+'Locus_Tag'+'\t'+'GeneName'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'KO_Id'+'\t'+'COG'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\t'+'COG_Category'+'\n')
		outfile3 = open(out_file3,'w')
		outfile3.write('Module'+'\t'+'Transcript_Id'+'\t'+'Locus_Tag'+'\t'+'GeneName'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'KO_Id'+'\t'+'COG'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\t'+'COG_Category'+'\n')
		modid, transids = [],[];
		for x in range(len(module)): # DEG Enriched Modules
			if str(module[x]) in mod:
				index1 = int(mod.index(str(module[x])))
				mod_transcripts = transcript[index1]; 
				for t in range(len(mod_transcripts)): # All transcripts in DEG Enriched Modules
					outfile1.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
					modid.append(str(module[x]))
					transids.append(mod_transcripts[t])
					if (mod_transcripts[t] in transcript_tag):
						index2 = int(transcript_tag.index(mod_transcripts[t]))
						locus_gene = locus_tag[index2]
						if (locus_gene in locus_id):
							index3 = int(locus_id.index(locus_gene))
							if (float(fold_change[index3])>=4): # Up-regulated DEGs. You may need to change to >=4 instead of >2 
								if (locus_gene not in locusid): # Up-regulated DEGs not in paper.
									outfile2.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+locus_gene+'\t'+gene_name[index3]+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+ko_id[index3]+'\t'+cog_id[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\t'+cog_category[index3]+'\n')
							if ((float(fold_change[index3])<=0.25) and (float(fold_change[index3]) > 0)): # Down-regulated DEGs 
								if (locus_gene not in locusid_down): # Up-regulated DEGs not in paper.
									outfile3.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+locus_gene+'\t'+gene_name[index3]+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+ko_id[index3]+'\t'+cog_id[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\t'+cog_category[index3]+'\n')

		f.close(); outfile1.close(); outfile2.close(); outfile3.close();
		
if __name__ == "__main__":
	NonStatistical_DEG().UpDown_NonStatistical_DEGs('Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Input_Files/Murine_Dataset_DownRegulatedDEGs.txt', 'Input_Files/Murine_Dataset_UpRegulatedDEGs.txt', 'Input_Files/pgen.1004518.s005.txt', 'Input_Files/PAO1_Annotations.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DownRegulatedNonDEGs_In_DEGenrichedModules.txt')



