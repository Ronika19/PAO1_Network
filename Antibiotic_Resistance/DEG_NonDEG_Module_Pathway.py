import re
import sys
from Azithromycin_Serum.Upregulated_NonStatistical_DEGs import NonStatistical_DEG

class DEG_NonDEG_Module:
	def DEG_NonDEG_Count_Module(self, infile1, infile2, infile3, infile4, infile5, out_file):
		file_list = ['PAO1_DEGMod_Enrichment.txt', 'WGCNA_PAO1_Modules_Clusters.txt', 'AZM_SERUM_BiologicalReplicate.txt', 'AZM_SERUM_Paper.txt', 'PAO1_AMR_Annotations.tsv', 'All_Transcripts_In_DEGenrichedModules.txt', 'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt']
		NonStatistical_DEG().UpDown_NonStatistical_DEGs(file_list[0], file_list[1], file_list[2], file_list[3], file_list[4], file_list[5], file_list[6])

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
		transcript_id, locus_tag, fold_change, pvalue = [],[],[],[];
		while line2:
			split_line2 = (line2.rstrip()).split('\t')
			transcript_id.append(split_line2[3])
			locus_tag.append(split_line2[0])
			fold_change.append(split_line2[1])
			pvalue.append(split_line2[2])
			line2 = file2.readline()

		file4 = open(infile3,'r')
		ltag, tid = [],[];
		for line4 in file4:
			split_line4 = (line4.rstrip()).split('\t')
			ltag.append(split_line4[0])
			tid.append(split_line4[1])

		file5 = open(infile4,'r')
		for i in range(2):
			line5 = file5.readline()
		locusid, transcript_tag = [],[]
		while line5:
			split_line5 = (line5.rstrip()).split('\t')
			locusid.append(split_line5[0])
			if (split_line5[0] in ltag):
				indexes = int(ltag.index(split_line5[0]))
				transcript_tag.append(tid[indexes])
			line5 = file5.readline()

		file3 = open(infile5,'r')
		for i in range(2):
			line3 = file3.readline()
		modules, transcriptid, locustag, expression_foldchange, expression_pvalue = [],[],[],[],[];
		module_foldchange, module_fishertest = [],[];
		while line3:
			split_line3 = (line3.rstrip()).split('\t')
			modules.append(split_line3[0])
			transcriptid.append(split_line3[1]); #print(transcriptid);
			locustag.append(split_line3[2])
			expression_foldchange.append(split_line3[3])
			expression_pvalue.append(split_line3[4])
			module_foldchange.append(split_line3[5])
			module_fishertest.append(split_line3[6])
			line3 = file3.readline()

		outfile = open(out_file,'w')
		counter=0;
		outfile.write('Module'+'\t'+'DEG_Count'+'\t'+'LDEG_Count'+'\t'+'Total_DEG_Count'+'\t'+'Module_Gene_Count'+'\t'+'DEG_Percent_In_Module'+'\t'+'LDEG_Percent_In_Module'+'\n')
		for i in range(len(transcripts)):
			transcript_deg_count = 0; transcript_nondeg_count = 0; non_author_deg = 0;
			for j in range(len(transcripts[i])):
				if ((transcripts[i])[j] in transcript_id):
					ind = int(transcript_id.index((transcripts[i])[j]))
					if (pvalue[ind] != 'NA'):
						if ((transcripts[i])[j] in transcript_tag):
							transcript_deg_count += 1
						elif (((transcripts[i])[j] not in transcript_tag) and (float(pvalue[ind]) > 0.05)): #Insignificant Up&Down Genes
							transcript_nondeg_count += 1; print(transcript_nondeg_count)
						elif (((transcripts[i])[j] not in transcript_tag) and (float(pvalue[ind]) <= 0.05)): # Statistically significant DEGs not found by Authors
							non_author_deg += 1; print(non_author_deg)
							counter += 1
					elif (pvalue[ind] == 'NA'):
						transcript_nondeg_count += 1
				elif (((transcripts[i])[j] in transcriptid) and (pvalue[ind] != 'NA')):
					transcript_nondeg_count += 1
			transcript_count = transcript_deg_count+transcript_nondeg_count
			print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i])
			#if (float(transcript_count/int(module_length[i])) >= 0.9):
			print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i], str(float(transcript_count/int(module_length[i]))))
			outfile.write(str(module[i])+'\t'+str(transcript_deg_count)+'\t'+str(transcript_nondeg_count)+'\t'+str(transcript_count)+'\t'+str(module_length[i])+'\t'+str(float(transcript_count/int(module_length[i])))+'\t'+str(float(transcript_deg_count/int(module_length[i])))+'\n') # LDEG_Percent_In_Module: These are the DEGs that show Fold change > 1.5 or Fold change < -1.5, and Pvalue <= 0.05.


DEG_NonDEG_Module().DEG_NonDEG_Count_Module('WGCNA_PAO1_Modules_Clusters.txt', 'AZM_SERUM_BiologicalReplicate.txt', 'PAO1_AMR_Annotations.tsv', 'AZM_SERUM_Paper.txt', 'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'DEG_NonDEG_Module_Pathway.txt')





