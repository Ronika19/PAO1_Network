import re

class Data_Processer:
	def AMR_Format(infile1, infile2, infile3, infile4, out_file): 
		file1 = open(infile1,'r')
		for i in range(2):
			line1 = file1.readline()
		gene, log2FC, padj = [],[],[];
		while line1:
			split_line1 = (line1.rstrip()).split('\t'); #print(split_line1)
			gene.append(split_line1[0])
			log2FC.append(split_line1[-1])
			padj.append(split_line1[4])
			line1 = file1.readline()

		file2 = open(infile2,'r')
		for i in range(2):
			line2 = file2.readline()
		genes, log2fc, padjv = [],[],[];
		while line2:
			split_line2 = (line2.rstrip()).split('\t')
			genes.append(split_line2[0])
			log2fc.append(split_line2[2])
			padjv.append(split_line2[-1])
			line2 = file2.readline()

		file3 = open(infile3,'r')
		for i in range(2):
			line3 = file3.readline()
		while line3:
			split_line3 = (line3.rstrip()).split('\t')
			genes.append(split_line3[0])
			log2fc.append(split_line3[2])
			padjv.append(split_line3[-1])
			line3 = file3.readline()

		file4 = open(infile4,'r')
		#for i in range(2):
		line4 = file4.readline()
		gene_id, transcript_id = [],[];
		while line4:
			split_line4 = (line4.rstrip()).split('\t')
			gene_id.append(split_line4[0])
			transcript_id.append(split_line4[1])
			line4 = file4.readline()

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
	
			
	def DEG_Cluster(infile1, infile2, infile3, outfile):
		file1 = open(infile1,'r')
		file4 = open(outfile,'w')

		file2 = open(infile2,'r')
		gene_id, transcript_id = [],[];
		for line2 in file2:
			split_line2 = (line2.rstrip()).split('\t')
			gene_id.append(split_line2[0])
			transcript_id.append(split_line2[1])

		file3 = open(infile3,'r')
		for i in range(2):
			line3 = file3.readline()
		gene, log2FC, padj = [],[],[];
		while line3:
			split_line3 = (line3.rstrip()).split('\t'); #print(split_line3)
			gene.append(split_line3[0])
			log2FC.append(split_line3[-1])
			padj.append(split_line3[4])
			line3 = file3.readline()

		for i in range(2):
			line1 = file1.readline()
		while line1:
			line1 = line1.rstrip()
			split_line1 = line1.split('\t'); print(split_line1)
			genes1 = split_line1[0]
			if ((genes1 in gene_id) and (genes1 in gene)):
				indexes = int(gene_id.index(genes1))
				file4.write(transcript_id[indexes]+"\t")
			line1 = file1.readline()
		file4.write("\n")
		file4.close()


	def Modules_Cluster(infile, outfile1, outfile2):
		file1 = open(infile,'r')
		for i in range(2):
			line1 = file1.readline()
		genes, modules = [],[];
		while line1:
			line1 = (line1.rstrip()).replace('"','')
			split_line1 = line1.split()
			genes.append(split_line1[0])
			modules.append(int(split_line1[1]))
			line1 = file1.readline()
		#print(len(genes),len(modules))

		mods = []
		for module in modules:
			if module not in mods:
				mods.append(int(module))
		mods.sort()
		x = 0; gene=[];
		file2 = open(outfile1,'w')
		file3 = open(outfile2,'w')
		for mod in mods:
			file2.write(str(mod)+"\t")
			for x in range(len(modules)):
				if (modules[x] == mod):
					#print(mod, modules[x], genes[x])
					gene.append(genes[x])
					file2.write(genes[x]+"\t")
					file3.write(genes[x]+"\t")
					#file2.write(str(mod)+"\t"+str(modules[x])+"\t"+genes[x]+"\n")
			file2.write("\n"); file3.write("\n")
		file1.close(); file2.close(); file3.close()


	def DEGClusters_2_WGCNAModules(infile1, infile2, infile3, outfile1, outfile2):
		
		# All DEGs
		file1 = open(infile1,'r')
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
		file2 = open(infile2,'r')
		for i in range(2):
			line2 = file2.readline()
		modules, gene_id = [],[];
		while line2:
			line2 = line2.rstrip()
			split_line2 = line2.split('\t')
			gene_id.append(split_line2[0].replace('"',''))
			modules.append(split_line2[1])
			line2 = file2.readline()
		print(DEG_transcripts,"\n",gene_id)

		file3 = open(outfile1,'w')
		k = 0; DEG1_Module = [];
		for genes in DEG_transcripts:
			if (k < int(indices[0])):	# Cluster of Pathogenesis DEGs
				if (genes in gene_id):
					indexes = gene_id.index(genes); print(indexes, genes, gene_id[indexes]);
					DEG1_Module.append(modules[indexes])
					file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
			k += 1

		file3.close(); file1.close(); file2.close();

		file4 = open(infile3,'r')
		line4 = file4.readline()
		geneid, mod = [],[];
		while line4:
			line4 = line4.rstrip()
			split_line4 = line4.split('\t')
			geneid.append(split_line4[0])
			mod.append(int(split_line4[1]))
			line4 = file4.readline()

		modset = sorted(set(mod), reverse=False); print(modset);

		file5 = open(outfile2,'w')
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





