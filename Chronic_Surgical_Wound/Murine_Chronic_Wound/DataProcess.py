import re

class Data_Processer:
	def data_extract(self, infile, deletes):
		f = open(infile, 'r')
		lines = f.readlines()
		dict_array = {}
		split_l = ((lines[1]).rstrip()).split('\t')
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

	def Pathogenesis_Format(self, infile): 
		dict1 = self.data_extract(infile, 1)
		transcript_id, locus_tag, fold_change, gene_name, ko_id, pvalue, cog =  dict1['arr_1'], dict1['arr_0'], dict1['arr_3'], dict1['arr_2'], dict1['arr_5'], dict1['arr_4'], dict1['arr_6']
		outfile_list = []
		for i in range(len(transcript_id)):
			outfile_list.append(transcript_id[i]+'\t'+locus_tag[i]+'\t'+fold_change[i]+'\t'+gene_name[i]+'\t'+ko_id[i]+'\t'+pvalue[i]+'\t'+cog[i]+'\n')
		return outfile_list

	def DEG_Cluster(self, infile, outfile1, outfile2, outfile3):
		dict2 = self.data_extract(infile, 1) 
		genes, fold_change, pvalue = dict2['arr_0'], dict2['arr_2'], dict2['arr_5']
		file2 = open(outfile1,'w'); file3 = open(outfile2,'w'); file4 = open(outfile3,'w')
		for i in range(len(genes)):
			if ((float(fold_change[i]) >= 4) and (float(pvalue[i]) < 0.01)):
				file4.write(genes[i]+"\t"); file2.write(genes[i]+"\n");
			if ((float(fold_change[i]) <= 0.25) and (float(pvalue[i]) < 0.01)):
				file4.write(genes[i]+"\t"); file3.write(genes[i]+"\n");
		file4.write("\n")
		file2.close(); file3.close(); file4.close();

	def Modules_Cluster(self, infile, outfile1, outfile2):
		dict3 = self.data_extract(infile, 1) 
		genes, modules = dict3['arr_0'], dict3['arr_1']
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
					#print(mod, modules[x], genes[x])
					gene.append(genes[x])
					file2.write(genes[x]+"\t")
					file3.write(genes[x]+"\t")
			file2.write("\n"); file3.write("\n")
		file2.close(); file3.close();

	def DEGClusters_2_WGCNAModules(self, infile1, infile2, infile3, outfile1, outfile2):
		# All DEGs
		file1 = open(infile1,'r'); file3 = open(outfile1,'w'); file5 = open(outfile2,'w')
		line = file1.readline()
		DEG_transcripts = []; indices = [];
		while line:
			split_line = (line.rstrip()).split('\t')
			for tids in split_line:
				DEG_transcripts.append(tids)
			DEG_transcripts.append('\n')
			line = file1.readline()

		for items in range(len(DEG_transcripts)):
			if DEG_transcripts[items] == '\n':
				indices.append(items); #print(items)

		# Modules using WGCNA
		dict4 = self.data_extract(infile2,1) 
		gene_ids, modules = dict4['arr_0'], dict4['arr_1']
		gene_id = [val.replace('"','') for i, val in enumerate(gene_ids)]
		
		k = 0; DEG1_Module = [];
		for genes in DEG_transcripts:
			if (k < int(indices[0])):	# Cluster of DEGs
				if (genes in gene_id):
					indexes = gene_id.index(genes); #print(indexes, genes, gene_id[indexes]);
					DEG1_Module.append(modules[indexes])
					file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
			k += 1
		file3.close(); file1.close(); 

		dict5 = self.data_extract(infile3,0) 
		geneid, mod = dict5['arr_0'], dict5['arr_1']
		modset = sorted(set(mod), reverse=False); #print(modset);
		for i in modset:
			file5.write(str(i)+'\t')
			indices = [index for index, element in enumerate(mod) if element == i]; #print(i, indices);
			m = 0
			for j in indices:
				m += 1; #print(geneid[j], mod[j])
				if m < len(indices):
					file5.write(geneid[j]+',')
				elif m == len(indices):
					file5.write(geneid[j]+'\n')
		file5.close()


