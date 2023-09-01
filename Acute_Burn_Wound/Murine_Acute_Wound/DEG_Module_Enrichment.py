import scipy.stats as stats
from .DataProcess import Data_Processer

class DEG_Enrichment:
	def DataPreprocess(self, outfile):
		outfile = open(outfile,'w')
		outfile.write('Transcript_Id'+'\t'+'Locus_Tag'+'\t'+'Fold_Change'+'\t'+'Gene_Name'+'\t'+'KO_Id'+'\t'+'Pvalue'+'\t'+'COG'+'\n')
		infile_list = ['Input_Files/Murine_Dataset_UpRegulatedDEGs.txt', 'Input_Files/Murine_Dataset_DownRegulatedDEGs.txt']
		for infl in infile_list:
			Format = Data_Processer().Pathogenesis_Format(infl)
			for f in Format:
				outfile.write(f)
		outfile.close()

		fl_list = ['Output_Files/Pathogenesis_Vs_BiologicalReplicate.txt', 'Output_Files/PAO1_Upregulated_DEGs.txt', 'Output_Files/PAO1_Downregulated_DEGs.txt', 'Output_Files/PAO1_Cluster_DEG_BiologicalReplicates.txt']
		Cluster = Data_Processer().DEG_Cluster(fl_list[0], fl_list[1], fl_list[2], fl_list[3])
	
		mod_clustlist = ['Input_Files/Module_Assignment_PAO1_WGCNA.tsv', 'Output_Files/PAO1_Modules_Clusters.txt', 'Output_Files/WGCNA_PAO1_Modules_Clusters.txt']
		mod_cluster = Data_Processer().Modules_Cluster(mod_clustlist[0], mod_clustlist[1], mod_clustlist[2])

		clus_mod_list = ['Output_Files/PAO1_Cluster_DEG_BiologicalReplicates.txt', 'Input_Files/Module_Assignment_PAO1_WGCNA.tsv', 'Output_Files/PAO1_DEG_Modules.txt', 'Output_Files/PAO1_DEG_Modules.txt', 'Output_Files/PAO1_DEG_Genes_Modules.txt']
		Cluster_Modules = Data_Processer().DEGClusters_2_WGCNAModules(clus_mod_list[0], clus_mod_list[1], clus_mod_list[2], clus_mod_list[3], clus_mod_list[4])

	def Enrichment(self, infile1, infile2, infile3, infile4, outfile):
		self.DataPreprocess('Output_Files/Pathogenesis_Vs_BiologicalReplicate.txt')

		nodes_file = open(infile1,'r')
		line = nodes_file.readlines()
		len_nodes = len(line[1:]); #print(len_nodes);	# Total no of Genes in WGCNA Network
	
		file1_1 = open(infile2,'r')
		Total_Modules = len(file1_1.readlines()); # print(Total_Modules)
		file1_1.close()
		file1 = open(infile3,'r')
		line1 = file1.readlines()
		Mod_Size, Mod_Size_Percent = [],[]; 

		zero_modgenes = [value for x, value in enumerate (((line1[0]).rstrip()).split('\t'))]; print(zero_modgenes);

		for l1 in line1:
			l1 = l1.rstrip()
			split_line1 = l1.split('\t'); #print(len(split_line1), (len(split_line1)/len_nodes)*100);
			Mod_Size.append(len(split_line1))	# No of Genes or Gene Count in each WGCNA Modules = Mod_Size[1:]
			Mod_Size_Percent.append((len(split_line1)/(len_nodes-len(zero_modgenes))*100))	# Percentage of Genes in each WGCNA Modules = Mod_Size_Percent[1:]

		file2 = open(infile4,'r')
		line2 = file2.readline()
		PAO1_dict = {}; len_PAO1_degs = 0
		while line2:
			line2 = line2.rstrip()
			split_line2 = line2.split('\t')
			PAO1_deg_mods = int(split_line2[1]); #print(PAO1_deg_mods);
			for i in range(Total_Modules):
				if (PAO1_deg_mods == i):
					PAO1_dict[split_line2[0]] = i
			len_PAO1_degs += 1	# No of PAO1 DEGs
			line2 = file2.readline()
		#print(PAO1_dict); print(len_PAO1_degs); print(Mod_Size_Percent);

		dict_genes, zero_moddeg = [],[]
		for i in PAO1_dict.keys():
			dict_genes.append(i)

		for i in dict_genes:
			if i in zero_modgenes:
				zero_moddeg.append(i)
			
		zero_modnondeg = len(zero_modgenes)-len(zero_moddeg); print(len(zero_moddeg));

		outfile2 = open(outfile,'w')
		outfile2.write('Module'+'\t'+'Observed_DEGs'+'\t'+'Expected_DEGs'+'\t'+'Fold_Change'+'\t'+'Fisher_Test'+'\n')
		PAO1_degs_count = {}; PAO1_degs_counter = []
		for j in range(Total_Modules):
			PAO1_count = sum(x == j for x in PAO1_dict.values()); #print(PAO1_count);
			PAO1_degs_count[j] = PAO1_count	# dictionary where key = modules & values = DEG count for Pathogenesis
			PAO1_degs_counter.append(PAO1_count)	# List of DEG counts in each module for Pathogenesis
			PAO1_deg_mod = float(PAO1_count)	# DEG in each Module
			PAO1_nondeg_mod = float(Mod_Size[j])-float(PAO1_count)	# Non-DEG in each Module
			PAO1_deg_nonmod = float(len_PAO1_degs)-float(PAO1_count)-float(len(zero_moddeg))	# DEG not in each Module
			PAO1_nondeg_nonmod = float(len_nodes)-float(PAO1_nondeg_mod)-float(zero_modnondeg)	# Non-DEG not in each Module
			oddsratio, pvalue = stats.fisher_exact([[PAO1_deg_mod, PAO1_nondeg_mod], [PAO1_deg_nonmod, PAO1_nondeg_nonmod]]); # p-values of DEGs in each Module
			PAO1_observed_degs_mod = PAO1_count;	# DEGs observed in each Module
			PAO1_expected_degs_mod = (len_PAO1_degs*Mod_Size_Percent[j])/100	# DEGs in each Module, expected by random chance
			PAO1_fold_change = float(PAO1_observed_degs_mod)/float(PAO1_expected_degs_mod)	# Fold Change of DEGs in each Module
			#print("P-value = ", pvalue);
			outfile2.write(str(j)+'\t'+str(PAO1_observed_degs_mod)+'\t'+str(PAO1_expected_degs_mod)+'\t'+str(PAO1_fold_change)+'\t'+str(pvalue)+'\n')
		nodes_file.close(); file1.close(); file2.close(); outfile2.close()

#if __name__ == "__main__":
#	DEG_Enrichment().Enrichment()


