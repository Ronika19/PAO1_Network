import re

class Data_Processer:
	def data_extract(self, infile, deletes):
		f = open(infile, 'r')
		lines = f.readlines()
		dict_array = {}
		split_l = ((lines[5]).rstrip()).split('\t')
		for i in range(len(split_l)):
			dict_array['arr_'+str(i)] = []
		if (deletes == 0):
			for l in lines:
				split_l = (l.rstrip()).split('\t')
				for x in range(len(split_l)):
					dict_array['arr_'+str(x)].append(split_l[x])
		elif (deletes > 0):
			for l in lines[deletes:]:
				split_l = (l.rstrip()).split('\t')
				for x in range(len(split_l)):
					dict_array['arr_'+str(x)].append(split_l[x])
		f.close(); #print(dict_array);
		return dict_array

	def DEG_Cluster(self, infile, outfile):
		f = open(outfile,'w')
		dict1 = self.data_extract(infile,1)
		gene_ids, fold_change, pval = dict1['arr_0'], dict1['arr_2'], dict1['arr_6']
		for i in range(len(fold_change)):
			if ((float(fold_change[i]) >= 1) and (float(pval[i]) <= 0.05)):
				f.write(gene_ids[i]+"\t")
			if ((float(fold_change[i]) <= -1) and (float(pval[i]) <= 0.05)):
				f.write(gene_ids[i]+"\t")
		f.write("\n"); f.close();

	def Module_Clusters(self, infile, outfile1, outfile2):
		dict2 = self.data_extract(infile,1)
		genes_ids, modules = dict2['arr_0'], dict2['arr_1']
		genes = [ids.replace('"','') for i, ids in enumerate(genes_ids)]

		mods = []
		for module in modules:
			if int(module) not in mods:
				mods.append(int(module))
		mods.sort()
		x = 0; gene=[];
		file2 = open(outfile1,'w'); file3 = open(outfile2,'w');
		for mod in mods:
			file2.write(str(mod)+"\t")
			for x in range(len(modules)):
				if (int(modules[x]) == int(mod)):
					gene.append(genes[x])
					file2.write(genes[x]+"\t"); file3.write(genes[x]+"\t");
			file2.write("\n"); file3.write("\n");
		file2.close(); file3.close();

	def DEGClusters_2_WGCNAModules(self, infile1, infile2, outfile1, infile3, outfile2):
		# All DEGs
		file1 = open(infile1,'r'); file3 = open(outfile1,'w'); file5 = open(outfile2,'w');
		line1 = file1.readline()
		DEG_transcripts = []; indices = [];
		while line1:
			line1 = line1.rstrip()
			split_line1 = line1.split('\t')
			for tids in split_line1:
				DEG_transcripts.append(tids)
			DEG_transcripts.append('\n')
			line1 = file1.readline()

		for items in range(len(DEG_transcripts)):
			if DEG_transcripts[items] == '\n':
				indices.append(items); print(items)

		# Modules using WGCNA
		dict3 = self.data_extract(infile2,1) 
		gene_ids, modules = dict3['arr_0'], dict3['arr_1']
		gene_id = [val.replace('"','') for i, val in enumerate(gene_ids)]
		
		k = 0; DEG1_Module = [];
		for genes in DEG_transcripts:
			if (k < int(indices[0])):	# Cluster of DEGs
				if (genes in gene_id):
					indexes = gene_id.index(genes); print(indexes, genes, gene_id[indexes]);
					DEG1_Module.append(modules[indexes])
					file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
			k += 1

		file3.close(); file1.close();

		dict4 = self.data_extract(infile3,0) 
		geneid, mod = dict4['arr_0'], dict4['arr_1']

		modset = sorted(set(mod), reverse=False); print(modset);
		for i in modset:
			file5.write(str(i)+'\t')
			indices = [index for index, element in enumerate(mod) if element == i]; #print(i, indices);
			m = 0
			for j in indices:
				m += 1; print(geneid[j], mod[j])
				if m < len(indices):
					file5.write(geneid[j]+',')
				elif m == len(indices):
					file5.write(geneid[j]+'\n')

		file5.close();

if __name__ == "__main__":
	Data_Processer().DEG_Cluster('CF_Vs_LB_Significant.txt', 'PAO1_Cluster_DEG_BiologicalReplicates.txt')
	Data_Processer().Module_Clusters('Module_Assignment_PAO1_WGCNA.tsv', 'PAO1_Modules_Clusters.txt', 'WGCNA_PAO1_Modules_Clusters.txt')
	Data_Processer().DEGClusters_2_WGCNAModules('PAO1_Cluster_DEG_BiologicalReplicates.txt', 'Module_Assignment_PAO1_WGCNA.tsv', 'PAO1_DEG_Modules.txt', 'PAO1_DEG_Modules.txt', 'PAO1_DEG_Genes_Modules.txt')



