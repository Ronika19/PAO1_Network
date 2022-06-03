import re
from .DEG_Module_Enrichment import DEG_Enrichment

class NonStatistical_DEG:
	def UpDown_NonStatistical_DEGs(self, infile1, infile2, infile3, infile4, infile5, out_file1, out_file2):
		DEG_Enrichment().Enrichment('PAO1_Nodes.txt', 'WGCNA_PAO1_Modules_Clusters.txt', 'WGCNA_PAO1_Modules_Clusters.txt', 'PAO1_DEG_Modules.txt', 'PAO1_DEGMod_Enrichment.txt')

		file1 = open( infile1,'r')
		for i in range(2):
			line1 = file1.readline()
		module, observed, expected, foldchange, fishertest = [],[],[],[],[]
		stats_module, stats_foldchange, stats_fishertest = [],[],[]
		while line1:
			split_line1 = (line1.rstrip()).split('\t')
			if (float(split_line1[3]) > 1): # DEG Enriched Modules with Fold Change is greater than 1
				module.append(split_line1[0])
				observed.append(split_line1[1])
				expected.append(split_line1[2])
				foldchange.append(split_line1[3])
				fishertest.append(split_line1[4])
				if (float(split_line1[4]) <= 0.05): # DEG Enriched Modules with p-value is less than equal to 0.05
					stats_module.append(split_line1[0])
					stats_foldchange.append(split_line1[3])
					stats_fishertest.append(split_line1[4])
			line1 = file1.readline()

		file5 = open( infile2,'r')
		m=0; mod, transcript = [],[];
		for line5 in file5:
			split_line5 = (line5.rstrip()).split('\t')
			mod.append(str(m))
			transcript.append(split_line5)
			m += 1

		file3 = open( infile3,'r')
		for i in range(2):
			line3 = file3.readline()
		locusid, foldchanges, pvalue = [],[],[];
		while line3:
			split_line3 = (line3.rstrip()).split('\t')
			locusid.append(split_line3[0])
			foldchanges.append(split_line3[1])
			pvalue.append(split_line3[2])
			line3 = file3.readline()

		file4 = open( infile4,'r')
		for i in range(2):
			line4 = file4.readline()
		locus_id, fold_change, p_value = [],[],[];
		while line4:
			split_line4 = (line4.rstrip()).split('\t')
			locus_id.append(split_line4[0])
			#gene_name.append(split_line4[1])
			fold_change.append(split_line4[-1])
			p_value.append(split_line4[4])
			line4 = file4.readline()

		file6 = open( infile5,'r')
		for i in range(2):
			line6 = file6.readline()
		locus_tag, transcript_tag = [],[];
		while line6:
			split_line6 = (line6.rstrip()).split('\t')
			locus_tag.append(split_line6[0])
			transcript_tag.append(split_line6[1]); #print(split_line6[6], split_line6[7]);
			line6 = file6.readline()

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
	NonStatistical_DEG().UpDown_NonStatistical_DEGs('PAO1_DEGMod_Enrichment.txt', 'WGCNA_PAO1_Modules_Clusters.txt', 'AZM_SERUM_BiologicalReplicate.txt', 'AZM_SERUM_Paper.txt', 'PAO1_AMR_Annotations.tsv', 'All_Transcripts_In_DEGenrichedModules.txt', 'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt')




