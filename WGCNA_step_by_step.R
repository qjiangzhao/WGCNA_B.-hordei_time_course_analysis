#=====================================================================================
#
#  Code chunk 1 import gene nomolized table and feature table 
#
#=====================================================================================

setwd("G:/WGCNA")
# Load the WGCNA package
library(WGCNA)
library(tidyverse)


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the RNA seq data set
TPM_Data = read_csv("bgh_dh14_v4-2_mRNA_and_lncRNA1_2_count_table_with_lengh_TPM_filtered_by_TEseq_v3", col_names = F);

### read in the trail table 
feature_table <- read_csv("design.tc-table.txt")

# Take a quick look at what is in the data set:
dim(TPM_Data);

### convert column to rownames 
TPM_Data <- column_to_rownames(TPM_Data, var = "X1")

### assign feature table to column names 
names(TPM_Data) <- feature_table$sampleid

#=====================================================================================
#
#  Code chunk 2 calculate distances among samples and eliminate samples that far away with other samples
#
#=====================================================================================

### t(TPM_data) will convert column and rows 
### dist will calculate the distance for each rows, this is the reason why the conversion of column 
### and rows is required

TPM_Data <- t(TPM_Data)

sampleTree = hclust(dist(TPM_Data), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
### this function adds one or more straight lines through the current plot. 
### h: the y-values for horizontal lines.
abline(h = 15000, col = "red"); 

# Determine cluster under the line
### cutHeight: height at which branches are to be cut.
### minSize: minimum number of object on a branch to be considered a cluster.
clust = cutreeStatic(sampleTree, cutHeight = 15000, minSize = 10)
table(clust)

# clust 0 contains the samples we want to keep.
keepSamples = (clust==0) ### notice the way how to filter cluster number 
datExpr = TPM_Data[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

### convey feature_table to datTraits
datTraits <- feature_table
datTraits <- column_to_rownames(datTraits, var = "sampleid")

#=====================================================================================
#
#  Code chunk 3 plot tree of samples again 
#
#=====================================================================================

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
png("Sample_cluster.png", width = 600, height = 600)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 4 choose a soft-threshold for building network  
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
### verbose simply defines how "talky" a function is so how much status messages it prints to screen.
### This has no effect on performance or output.
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
png("Choose_soft_threshold.png", width = 800, height = 400)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 5 manually network construction 
#
#=====================================================================================

### datExpr include all TPM values 
### calculate pcc number of gene based on TPM value
### apply function might be better for the step, change it when you have time 
adjacency_neg <- cor(datExpr)

nrow(adjacency_neg)

for (i in 1:4836 ){
  for(j in 1:4836){
    if (adjacency_neg[i,j] > 0  && adjacency_neg[i,j] != 1){
      adjacency_neg[i,j] = 0.0001^8 } else {
        adjacency_neg[i,j] = adjacency_neg[i,j]^8
      }
    }
  }

### alternatively the built in function adjacency can be used
### not only pearson but also spearman can be used for adjaceny construction.
### the problem of spearman is the calculation time is long. 
softPower = 8;
adjacency = adjacency(datExpr,
                      power = softPower,
                      type = "unsigned",
                      corFnc = "cor",
                      #corOptions = list(use = 'p', method = 'spearman')
                      );

save(adjacency, file = "WGCNA_step_by_step_pearson_adjacency")

#=====================================================================================
#
#  Code chunk 6 calculate TOM value 
#
#=====================================================================================

### TOM value will be used to build the cluster tree
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 7 plot cluster tree
#
#=====================================================================================

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=====================================================================================
#
#  Code chunk 8 set module size and use dynamic tree cut 
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
### from the table you can see how many genes for each clusters
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 9 convert cluster number to color and plot clusters
#
#=====================================================================================

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 10 cluster tree construction, which can guide us cut the cluster tree later 
#
#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 11 cut cluster tree
#
#=====================================================================================


MEDissThres = 0.25 ### need to be changed 
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors; ### re-assign color to genes
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 12 show the new clusters 
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 13 save data 
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "WGCNA_stepByStep.RData")

#############################################################################################################################
#=====================================================================================
#
#  Code chunk 14 correlate modules with trail features
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


########################################################################################

sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 13, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
ggsave("Cluster_trail_correlation.png", width = 1000, height = 800)

##############################################################################################################################

#=====================================================================================
#
#  Code chunk 15 connect module with gene annotation file 
#
#=====================================================================================
names(datExpr)


### improt annotation file 

annot <- read_csv("out.emapper.csv", skip = 2)
dim(annot)
probes = colnames(datExpr)
probes2annot = match(probes, annot$query)

### indicate genes without annotation
sum(is.na(probes2annot))

### generate a dataframe contains all annotation information 
geneInfo0 <- data.frame(substanceBXH = probes, 
                       geneSymbol = annot$PFAMs[probes2annot],
                       geneDiscirption = annot$Description[probes2annot],
                      # LocusLinkID = annot$LocusLinkID[probes2annot], 
                      # moduleColor = moduleColors 
                      # geneTraitSignificance, GSPvalue
                      )
### import effector annotation file

eff_annot <- read_delim("in silico Bgh proteins +SP-TM_PGI_vers_1.csv", skip = 1, delim = ",", col_names = F )
eff_annot <- eff_annot%>%
  mutate(ID = paste(eff_annot$X1,"-mRNA-1", sep = ""))
probes = colnames(datExpr)
probes2eff_annot = match(probes, eff_annot$ID)
geneInfo0 <- geneInfo0%>%
  mutate(effector = eff_annot$X2[probes2eff_annot] )%>%
  mutate(bind_c = paste(geneInfo0$geneSymbol, geneInfo0$effector, sep = ""))


##############################################################################################################################

#=====================================================================================
#
#  Code chunk 16 export to cytoscape
#
#=====================================================================================

### export all module to cytoscape
modules = unique(moduleColors)


### exprot one module to cytoscape
modules = c("blue");

# Select module probes
#probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = colnames(datExpr)[inModule]
modProbes_tibble <- as_tibble(modProbes)
#write_csv(modProbes_tibble, "WGCNA_pearson_unsigned_powder_8_threshold_0.3_darkgreen_module_gene_ID.csv")


#modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.3,
                               nodeNames = modProbes,
                               altNodeNames = geneInfo0$bind_c[inModule],
                               nodeAttr = moduleColors[inModule]);

cyt$edgeData

#write_csv(cyt$edgeData, "WGCNA_darkgreen_06h_time_points_pearson_unsigned_powder_8_threshold_0.3_cytoscape.csv")

### bash coding for file modification
cat CytoscapeInput-edges-blue.txt | sed 's/NANA/-/g' | sed 's/NA//g' | sed 's/-effector/effector/g' > CytoscapeInput-edges-blue_m1.txt
cat CytoscapeInput-nodes-blue.txt | sed 's/NANA/-/g' | sed 's/NA//g' | sed 's/-effector/effector/g' > CytoscapeInput-nodes-blue_m1.txt


##############################################################################################################################
#=====================================================================================
#
#  Code chunk 17 visulization TOM_plot
#
#=====================================================================================

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)

TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")



