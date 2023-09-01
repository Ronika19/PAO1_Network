# Load Library and Library-Specific Params
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
  
# Data Cleaning
datExpr0 = t(as.data.frame(read.table("PAO1_Network.tsv",sep="\t",header=TRUE,row.names=1)));
# head(datExpr0); 
length(datExpr0); cat('\n\n');

# Pull "bad" genes/experiments
gsg = goodSamplesGenes(datExpr0,verbose=3)
datExpr0 = datExpr0[gsg$goodSamples,gsg$goodGenes]; 
 
# # Check the results
# head(datExpr0); 
dim (datExpr0); length(datExpr0); cat('\n\n');
 
# Cluster experiments to detect outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex=0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 3e+05, col = "red");

# Remove outliers
remove = c("SRR1654901","SRR9720457", "ERR2275088")
datExpr0 = datExpr0[!row.names(datExpr0) %in% remove,]
dim (datExpr0); length(datExpr0); cat('\n\n');

# Select power range for WGCNA tests
powers = c(c(1:10),seq(from=12,to=20,by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
 
# Plot results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
 
# Use established soft-thresh power to build network with adjusted settings
net = blockwiseModules(datExpr0, power = 5,TOMType = "signed", minModuleSize = 10,reassignThreshold = 0, mergeCutHeight = 0.1, deepSplit = 2,  numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "PAO1_WGCNA",verbose = 3,maxBlockSize=50000)

# Write modules to file
write.table(file="Module_Assignment_PAO1_WGCNA.tsv",net$colors,sep="\t")
 
# Reload the TOM
load("PAO1_WGCNA-block.1.RData")
 
# Export *most edges from network for further study
cyt = exportNetworkToCytoscape(TOM,edgeFile="PAO1_WGCNA_Edges.txt",nodeFile="PAO1_WGCNA_Nodes.txt",weighted=TRUE,threshold=0.04,nodeNames=colnames(datExpr0))
 


