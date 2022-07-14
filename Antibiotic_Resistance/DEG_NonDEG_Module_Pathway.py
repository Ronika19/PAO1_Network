import re
import sys
import glob
import pathlib
from pathlib import Path
from Azithromycin_Serum.DataProcess import Data_Processer
from Azithromycin_Serum.Upregulated_NonStatistical_DEGs import NonStatistical_DEG

class DEG_NonDEG_Module:
	def DEG_NonDEG_Count_Module(self, infile1, infile2, infile3, infile4, infile5, out_file):
		file_list = ['Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/AZM_SERUM_BiologicalReplicate.txt', 'Input_Files/AZM_SERUM_Paper.txt', 'Input_Files/PAO1_AMR_Annotations.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt']
		NonStatistical_DEG().UpDown_NonStatistical_DEGs(file_list[0], file_list[1], file_list[2], file_list[3], file_list[4], file_list[5], file_list[6])

		file1 = open(infile1,'r')
		mod_count = 0; module, transcripts, module_length = [],[],[];
		for line1 in file1:
			split_line1 = (line1.rstrip()).split('\t')
			module.append(mod_count)
			transcripts.append(split_line1)
			module_length.append(len(split_line1))
			mod_count += 1

		dict_1 = Data_Processer().data_extract(infile2, 1)
		transcript_id, locus_tag, fold_change, pvalue = dict_1['arr_3'], dict_1['arr_0'], dict_1['arr_1'], dict_1['arr_2']
		
		dict_2 = Data_Processer().data_extract(infile3, 0)
		ltag, tid = dict_2['arr_0'], dict_2['arr_1']

		dict_3 = Data_Processer().data_extract(infile4, 1)
		locusid, transcript_tag = dict_3['arr_0'], []
		for i in range(len(locusid)):
			if (locusid[i] in ltag):
				indexes = int(ltag.index(locusid[i]))
				transcript_tag.append(tid[indexes])
		
		dict_4 = Data_Processer().data_extract(infile5, 1)
		modules, transcriptid, locustag, expression_foldchange, expression_pvalue = dict_4['arr_0'], dict_4['arr_1'], dict_4['arr_2'], dict_4['arr_3'], dict_4['arr_4']
		module_foldchange, module_fishertest = dict_4['arr_5'], dict_4['arr_6']
		
		outfile = open(out_file,'w')
		counter=0;
		outfile.write('Module'+'\t'+'DEG_Count'+'\t'+'LDEG_Count'+'\t'+'Total_DEG_Count'+'\t'+'Module_Gene_Count'+'\t'+'DEG+LDEG_Percent_In_Module'+'\t'+'DEG_Percent_In_Module'+'\n')
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
				#elif (((transcripts[i])[j] in transcriptid) and (pvalue[ind] != 'NA')):
				#	transcript_nondeg_count += 1
			transcript_count = transcript_deg_count+transcript_nondeg_count
			print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i])
			#if (float(transcript_count/int(module_length[i])) >= 0.9):
			print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i], str(float(transcript_count/int(module_length[i]))))
			outfile.write(str(module[i])+'\t'+str(transcript_deg_count)+'\t'+str(transcript_nondeg_count)+'\t'+str(transcript_count)+'\t'+str(module_length[i])+'\t'+str(float(transcript_count/int(module_length[i])))+'\t'+str(float(transcript_deg_count/int(module_length[i])))+'\n') # LDEG_Percent_In_Module: These are the DEGs that show Fold change > 1.5 or Fold change < -1.5, and Pvalue <= 0.05.

		output_files = glob.glob('Output_Files/*'); print(output_files);
		keep_files = ['Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DEG_NonDEG_Module_Pathway.txt']
		for f in output_files:
			if (f not in keep_files):
				path = pathlib.Path(f)
				path.unlink()


DEG_NonDEG_Module().DEG_NonDEG_Count_Module('Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/AZM_SERUM_BiologicalReplicate.txt', 'Input_Files/PAO1_AMR_Annotations.tsv', 'Input_Files/AZM_SERUM_Paper.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DEG_NonDEG_Module_Pathway.txt')





