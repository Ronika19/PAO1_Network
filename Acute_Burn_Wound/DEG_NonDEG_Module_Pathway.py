import re
import sys
import glob
import pathlib
from pathlib import Path
from Murine_Acute_Wound.DataProcess import Data_Processer
from Murine_Acute_Wound.Upregulated_NonStatistical_DEGs import NonStatistical_DEG

class DEG_NonDEG_Module:
	def DEG_NonDEG_Count_Module(self, infile1, infile2, infile3, infile4, out_file):
		file_list = ['Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Input_Files/Murine_Dataset_DownRegulatedDEGs.txt', 'Input_Files/Murine_Dataset_UpRegulatedDEGs.txt', 'Input_Files/pgen.1004518.s005.txt', 'Input_Files/PAO1_Annotations.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DownRegulatedNonDEGs_In_DEGenrichedModules.txt']
		NonStatistical_DEG().UpDown_NonStatistical_DEGs(file_list[0], file_list[1], file_list[2], file_list[3], file_list[4], file_list[5], file_list[6], file_list[7], file_list[8])

		file1 = open(infile1,'r')
		mod_count = 0; module, transcripts, module_length = [],[],[];
		for line1 in file1:
			split_line1 = (line1.rstrip()).split('\t')
			module.append(mod_count)
			transcripts.append(split_line1)
			module_length.append(len(split_line1))
			mod_count += 1

		dict_1 = Data_Processer().data_extract(infile2, 1)
		transcript_id, locus_tag, fold_change, gene_name, ko_id, pvalue, cog = dict_1['arr_0'], dict_1['arr_1'], dict_1['arr_2'], dict_1['arr_3'], dict_1['arr_4'],  dict_1['arr_5'], dict_1['arr_6']

		dict_2 = Data_Processer().data_extract(infile3, 1)
		modules, transcriptid, locustag, genename, expression_foldchange, expression_pvalue = dict_2['arr_0'], dict_2['arr_1'], dict_2['arr_2'], dict_2['arr_3'], dict_2['arr_4'],  dict_2['arr_5']
		koid, cogid, module_foldchange, module_fishertest, cog_category = dict_2['arr_6'], dict_2['arr_7'], dict_2['arr_8'], dict_2['arr_9'], dict_2['arr_10']

		dict_3 = Data_Processer().data_extract(infile4, 1)
		modules, transcriptid, locustag, genename, expression_foldchange, expression_pvalue = modules+dict_3['arr_0'], transcriptid+dict_3['arr_1'], locustag+dict_3['arr_2'], genename+dict_3['arr_3'], expression_foldchange+dict_3['arr_4'],  expression_pvalue+dict_3['arr_5']
		koid, cogid, module_foldchange, module_fishertest, cog_category = koid+dict_3['arr_6'], cogid+dict_3['arr_7'], module_foldchange+dict_3['arr_8'], module_fishertest+dict_3['arr_9'], cog_category+dict_3['arr_10']

		outfile = open(out_file,'w')
		outfile.write('Module'+'\t'+'DEG_Count'+'\t'+'LDEG_Count'+'\t'+'Total_DEG_Count'+'\t'+'Module_Gene_Count'+'\t'+'DEG_Percent_In_Module'+'\t'+'LDEG_Percent_In_Module'+'\n')
		for i in range(len(transcripts)):
			transcript_deg_count = 0; transcript_nondeg_count = 0;
			for j in range(len(transcripts[i])):
				if ((transcripts[i])[j] in transcript_id):
					transcript_deg_count += 1
				if ((transcripts[i])[j] in transcriptid):
					transcript_nondeg_count += 1; #print((transcripts[i])[j]);
			transcript_count = transcript_deg_count+transcript_nondeg_count
			#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i])
			#if (float(transcript_count/int(module_length[i])) >= 0.9):
			#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i], str(float(transcript_count/int(module_length[i]))))
			outfile.write(str(module[i])+'\t'+str(transcript_deg_count)+'\t'+str(transcript_nondeg_count)+'\t'+str(transcript_count)+'\t'+str(module_length[i])+'\t'+str(float(transcript_count/int(module_length[i])))+'\t'+str(float(transcript_deg_count/int(module_length[i])))+'\n') # LDEG_Percent_In_Module: These are the DEGs that show Fold change >= 4 and, Pvalue < 0.01.

		output_files = glob.glob('Output_Files/*'); #print(output_files);
		keep_files = ['Output_Files/PAO1_DEGMod_Enrichment.txt', 'Output_Files/DownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DEG_NonDEG_Module_Pathway.txt']
		for f in output_files:
			if (f not in keep_files):
				path = pathlib.Path(f)
				path.unlink()


DEG_NonDEG_Module().DEG_NonDEG_Count_Module('Output_Files/WGCNA_PAO1_Modules_Clusters.txt', 'Output_Files/Pathogenesis_Vs_BiologicalReplicate.txt', 'Output_Files/UpRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DownRegulatedNonDEGs_In_DEGenrichedModules.txt', 'Output_Files/DEG_NonDEG_Module_Pathway.txt')



