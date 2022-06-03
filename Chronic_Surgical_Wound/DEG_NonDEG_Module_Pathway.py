import re
import sys
from Murine_Chronic_Wound.Upregulated_NonStatistical_DEGs import NonStatistical_DEG

class DEG_NonDEG_Module:
	def DEG_NonDEG_Count_Module(self, infile1, infile2, infile3, infile4, out_file):
		file_list = ['PAO1_DEGMod_Enrichment.txt', 'WGCNA_PAO1_Modules_Clusters.txt', 'Murine_Dataset_DownRegulatedDEGs.txt', 'Murine_Dataset_UpRegulatedDEGs.txt', 'pgen.1004518.s005.txt', 'PAO1_Annotations.tsv', 'All_Transcripts_In_DEGenrichedModules.txt', 'UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'DownRegulatedNonDEGs_In_DEGenrichedModules.txt']
		NonStatistical_DEG().UpDown_NonStatistical_DEGs(file_list[0], file_list[1], file_list[2], file_list[3], file_list[4], file_list[5], file_list[6], file_list[7], file_list[8])

		file1 = open(infile1,'r')
		mod_count = 0; module, transcripts, module_length = [],[],[];
		for line1 in file1:
			split_line1 = (line1.rstrip()).split('\t')
			module.append(mod_count)
			transcripts.append(split_line1)
			module_length.append(len(split_line1))
			mod_count += 1

		file2 = open(infile2,'r')
		for i in range(2):
			line2 = file2.readline()
		transcript_id, locus_tag, fold_change, gene_name, ko_id, pvalue, cog = [],[],[],[],[],[],[];
		while line2:
			split_line2 = (line2.rstrip()).split('\t')
			transcript_id.append(split_line2[0])
			locus_tag.append(split_line2[1])
			fold_change.append(split_line2[2])
			gene_name.append(split_line2[3])
			ko_id.append(split_line2[4])
			pvalue.append(split_line2[5])
			cog.append(split_line2[6])
			line2 = file2.readline()

		file3 = open(infile3,'r')
		for i in range(2):
			line3 = file3.readline()
		modules, transcriptid, locustag, genename, expression_foldchange, expression_pvalue = [],[],[],[],[],[];
		koid, cogid, module_foldchange, module_fishertest, cog_category = [],[],[],[],[];
		while line3:
			split_line3 = (line3.rstrip()).split('\t')
			modules.append(split_line3[0])
			transcriptid.append(split_line3[1]); #print(transcriptid);
			locustag.append(split_line3[2])
			genename.append(split_line3[3])
			expression_foldchange.append(split_line3[4])
			expression_pvalue.append(split_line3[5])
			koid.append(split_line3[6])
			cogid.append(split_line3[7])
			module_foldchange.append(split_line3[8])
			module_fishertest.append(split_line3[9])
			cog_category.append(split_line3[10])
			line3 = file3.readline()

		file4 = open(infile4,'r')
		for i in range(2):
			line4 = file4.readline()
		while line4:
			split_line4 = (line4.rstrip()).split('\t')
			modules.append(split_line4[0])
			transcriptid.append(split_line4[1]); #print(transcriptid);
			locustag.append(split_line4[2])
			genename.append(split_line4[3])
			expression_foldchange.append(split_line4[4])
			expression_pvalue.append(split_line4[5])
			koid.append(split_line4[6])
			cogid.append(split_line4[7])
			module_foldchange.append(split_line4[8])
			module_fishertest.append(split_line4[9])
			cog_category.append(split_line4[10])
			line4 = file4.readline()

		outfile = open(out_file,'w')
		outfile.write('Module'+'\t'+'DEG_Count'+'\t'+'LDEG_Count'+'\t'+'Total_DEG_Count'+'\t'+'Module_Gene_Count'+'\t'+'DEG_Percent_In_Module'+'\t'+'LDEG_Percent_In_Module'+'\n')
		for i in range(len(transcripts)):
			transcript_deg_count = 0; transcript_nondeg_count = 0;
			for j in range(len(transcripts[i])):
				if ((transcripts[i])[j] in transcript_id):
					transcript_deg_count += 1
				if ((transcripts[i])[j] in transcriptid):
					transcript_nondeg_count += 1; print((transcripts[i])[j]);
			transcript_count = transcript_deg_count+transcript_nondeg_count
			#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i])
			#if (float(transcript_count/int(module_length[i])) >= 0.9):
			#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i], str(float(transcript_count/int(module_length[i]))))
			outfile.write(str(module[i])+'\t'+str(transcript_deg_count)+'\t'+str(transcript_nondeg_count)+'\t'+str(transcript_count)+'\t'+str(module_length[i])+'\t'+str(float(transcript_count/int(module_length[i])))+'\t'+str(float(transcript_deg_count/int(module_length[i])))+'\n') # LDEG_Percent_In_Module: These are the DEGs that show Fold change >= 4 and, Pvalue < 0.01.


DEG_NonDEG_Module().DEG_NonDEG_Count_Module('WGCNA_PAO1_Modules_Clusters.txt', 'Pathogenesis_Vs_BiologicalReplicate.txt', 'UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'DownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'DEG_NonDEG_Module_Pathway.txt')






