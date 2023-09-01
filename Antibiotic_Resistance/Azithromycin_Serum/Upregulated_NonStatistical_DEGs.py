import re
from .DataProcess import Data_Processer
from .DEG_Module_Enrichment import DEG_Enrichment

class NonStatistical_DEG:
	def UpDown_NonStatistical_DEGs(self, infile1, infile2, infile3, infile4, infile5, out_file1, out_file2):
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

		dicts2 = Data_Processer().data_extract(infile3, 1)
		locusid, foldchanges, pvalue =  dicts2['arr_0'], dicts2['arr_1'], dicts2['arr_2']

		dicts3 = Data_Processer().data_extract(infile4, 1)
		locus_id, fold_change, p_value = dicts3['arr_0'], dicts3['arr_6'], dicts3['arr_4']
		
		dicts4 = Data_Processer().data_extract(infile5, 1)
		locus_tag, transcript_tag = dicts4['arr_0'], dicts4['arr_1']
		
		outfile1 = open(out_file1,'w')
		outfile1.write('Module'+'\t'+'Transcript_Id'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
		outfile2 = open(out_file2,'w')
		outfile2.write('Module'+'\t'+'Transcript_Id'+'\t'+'Locus_Tag'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
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
						if (locus_gene in locusid):
							if (locus_gene in locus_id): # Up & Down Regulated DEGs in paper
								index3 = int(locus_id.index(locus_gene))
								#outfile2.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+locus_gene+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
							elif (locus_gene not in locus_id): # Up & Down Regulated DEGs not in paper
								index4 = int(locusid.index(locus_gene))
								if ((float(foldchanges[index4]) > 0.58) and (float(pvalue[index4]) > 0.05)): # UpRegulated Non-DEGs
									outfile2.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+locus_gene+'\t'+foldchanges[index4]+'\t'+pvalue[index4]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
								elif ((float(foldchanges[index4]) < -0.58) and (float(pvalue[index4]) > 0.05)): # DownRegulated Non-DEGs
									outfile2.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+locus_gene+'\t'+foldchanges[index4]+'\t'+pvalue[index4]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')


if __name__ == "__main__":
	NonStatistical_DEG().UpDown_NonStatistical_DEGs('Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/AZM_SERUM_BiologicalReplicate.txt', 'Input_Files/AZM_SERUM_Paper.txt', 'Input_Files/PAO1_AMR_Annotations.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt')




