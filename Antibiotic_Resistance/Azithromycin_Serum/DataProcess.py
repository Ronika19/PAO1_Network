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

	def AMR_Format(self, infile1, infile2, infile3, infile4, out_file): 
		dict1 = self.data_extract(infile1, 1)
		gene, log2FC, padj = dict1['arr_0'], dict1['arr_6'], dict1['arr_4']
		
		dict2 = self.data_extract(infile2, 1)
		genes, log2fc, padjv = dict2['arr_0'], dict2['arr_2'], dict2['arr_6']
		
		dict3 = self.data_extract(infile3, 1)
		genes, log2fc, padjv = genes+dict3['arr_0'], log2fc+dict3['arr_2'], padjv+dict3['arr_6']
		
		dict4 = self.data_extract(infile4, 0)
		gene_id, transcript_id = dict4['arr_0'], dict4['arr_1']
		
		outfile = open(out_file,'w')
		outfile.write('Gene'+'\t'+'Log2FC'+'\t'+'Padj'+'\t'+'Transcript'+'\n')
		for g in genes:
			if g in gene:
				index1 = int(gene.index(g))
				if g in gene_id:
					index3 = int(gene_id.index(g))
					outfile.write(g+'\t'+log2FC[index1]+'\t'+padj[index1]+'\t'+transcript_id[index3]+'\n')
			elif g not in gene:
				index2 = int(genes.index(g))
				if g in gene_id:
					index4 = int(gene_id.index(g))
					outfile.write(g+'\t'+log2fc[index2]+'\t'+padjv[index2]+'\t'+transcript_id[index4]+'\n')

		for gn in gene:
			if gn not in genes:
				index6 = int(gene.index(gn)); #print(gn); print('\n\n')
				if gn in gene_id:
					index5 = int(gene_id.index(gn)); #print(gn)
					outfile.write(gn+'\t'+log2FC[index6]+'\t'+padj[index6]+'\t'+transcript_id[index5]+'\n')
		outfile.close()
	
			
	def DEG_Cluster(self, infile1, infile2, infile3, outfile):
		f = open(outfile,'w')

		dict5 = self.data_extract(infile2, 0)
		gene_id, transcript_id = dict5['arr_0'], dict5['arr_1']
		
		dict6 = self.data_extract(infile3, 1)
		gene, log2FC, padj = dict6['arr_0'], dict6['arr_6'], dict6['arr_4']

		dict7 = self.data_extract(infile1, 1)
		genes = dict7['arr_0']
		for i in range(len(genes)):
			if ((genes[i] in gene_id) and (genes[i] in gene)):
				indexes = int(gene_id.index(genes[i]))
				f.write(transcript_id[indexes]+"\t")
		f.write("\n")
		f.close()

	def Modules_Cluster(self, infile, outfile1, outfile2):
		dict8 = self.data_extract(infile, 1)
		genes, modules = dict8['arr_0'], dict8['arr_1']
		
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
					file2.write(genes[x]+"\t")
					file3.write(genes[x]+"\t")
			file2.write("\n"); file3.write("\n")
		file2.close(); file3.close()


	def DEGClusters_2_WGCNAModules(self, infile1, infile2, infile3, outfile1, outfile2):
		# All DEGs
		file1 = open(infile1,'r'); file3 = open(outfile1,'w'); file5 = open(outfile2,'w')
		line = file1.readline()
		DEG_transcripts = []; indices = [];
		while line:
			line = line.rstrip()
			split_line = line.split('\t')
			for tids in split_line:
				DEG_transcripts.append(tids)
			DEG_transcripts.append('\n')
			line = file1.readline()

		for items in range(len(DEG_transcripts)):
			if DEG_transcripts[items] == '\n':
				indices.append(items); print(items)

		# Modules using WGCNA
		dict9 = self.data_extract(infile2,1) 
		gene_ids, modules = dict9['arr_0'], dict9['arr_1']
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

		dict5 = self.data_extract(infile3,0) 
		geneid, mod = dict5['arr_0'], dict5['arr_1']

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

		file5.close()



