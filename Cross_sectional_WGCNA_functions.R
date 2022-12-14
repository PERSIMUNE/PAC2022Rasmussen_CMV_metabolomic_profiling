#------------------------------------------------------------------------------
# Script contains custom functions for data cleaning, WGCNA functions, and
# figure making.
#
# Author: Kirstine K. Rasmussen
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Data loading and formatting
#------------------------------------------------------------------------------

data_init <- function(data, IDs){
  # Checks data for non-variating variables and removes them.

  datExpr0 <- as.data.frame(data)
  rownames(datExpr0) <- IDs
  
  # Check for observations with no variation across patients
  # Should not do much as near-zero variation was removed in the pre-processing
  gsg = goodSamplesGenes(datExpr0, verbose = 0);
  if (gsg$allOK){
    print("No samples removed")
  }
  
  # Remove non-variated observations
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  return(datExpr0)
}

#------------------------------------------------------------------------------
# Outlier detection of samples
#------------------------------------------------------------------------------

check_outliers <- function(datExpr0, cutoff = 50){
  # Make hierarchical clustering of data
  sampleTree = hclust(dist(datExpr0), method = "average");
  
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6) # text size
  par(mar = c(0,4,2,0)) # margins of plot
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  # Plot a line to show the cut
  abline(h = cutoff, col = "red");
  
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = cutoff)
  print('Clusters detected:')
  print(table(clust))
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)
  datExpr = datExpr0[keepSamples, ]
  nMetabolites = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  return(datExpr)
}

#------------------------------------------------------------------------------
# Include trait data
#------------------------------------------------------------------------------

importClinical <- function(clinicalData, datExpr){
  # Load trait file and inspect
  #traitData = read.csv("Data/metabolites_clinical_closestto30_realcases.csv")
  allTraits = as.data.frame(read_tsv(clinicalData))
  dim(allTraits)
  names(allTraits)
  
  # Select columns for analysis
  allTraits = allTraits[, c("PATIENT","Group","CMV_risk","ttype","Sex","age_at_T")]
  names(allTraits) <- c("Patient","CMV infection","CMV risk score","Conditioning regimen","Sex","Age at aHSCT")
  printFlush(paste("\nClinical data included:", paste(names(allTraits), collapse = ", ")))
  
  # Chech if all are integers/floats/doubles
  # If characters exist, they need to be transformed
  for(i in 1:length(names(allTraits))){
    if(typeof(allTraits[1,i]) != "integer"){
      if(typeof(allTraits[1,i]) != "double"){
        print(names(allTraits[i]))
        print('Is not an integer or double. Please change.')
        print(typeof(allTraits[1,i]))}
    }
  }
  
  # Form a data frame analogous to expression data that will hold the selected clinical traits
  patientSamples = rownames(datExpr);
  traitRows = match(patientSamples, allTraits$Patient);
  datTraits = allTraits[traitRows, -1]; # Remove "Row.names" in first column
  rownames(datTraits) = allTraits[traitRows, 1]
  
  return(datTraits)
}

sampleTraitDend <- function(datExpr, datTraits){
  # Create sample-trait correlation dendrogram plot.
  
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  name <- paste('Figures/SampleTraitDend_',dataType,'.png', sep='')
  print(paste('Figure made:',name))
  png(name, width=700, height=400)
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits), hang = 0.02,
                      #main = "Sample dendrogram and trait heatmap",
                      main = "",
                      cex.colorLabels = 0.8, cex.dendroLabels = 0.8,
                      marAll = c(0,5,2,0), cex.main = 1.5)
  invisible(dev.off())
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits), hang = 0.02,
                      #main = "Sample dendrogram and trait heatmap",
                      main = "",
                      cex.colorLabels = 0.8, cex.dendroLabels = 0.8,
                      marAll = c(0,5,2,0), cex.main = 1.5)
}

#------------------------------------------------------------------------------
# Determine power parameters
#------------------------------------------------------------------------------

determine_power <- function(datExpr){
  # Makes plots of scale-free topology model fits for different
  # values of power. 
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 10, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = 'signed')
  #print(sft$fitIndices[,1:2])
  
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  # You want the lowest power that results in a scale-free topology network
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  # You do not want too low of a mean connectivity
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  powerList <- as.data.frame(cbind("power"=sft$fitIndices[,1],
                               "Fit"=-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  
  power <- ifelse(is.na(powerList[powerList$Fit > 0.9, ][1,1]),
                  30, # If model fit > 0.9 is not available
                  powerList[powerList$Fit > 0.9, ][1,1])
  
  return(power)
  }

scale_free_topology <- function(datExpr, power){
  ADJ = adjacency(datExpr, type = 'signed', power = power)

  k = as.vector(apply(ADJ, 2, sum, na.rm = T))
  
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  hist(k, main = "Adjacency counts in data", xlab = "Adjacencies")
  scaleFreePlot(k, main="Check scale free topology\n")
}

#------------------------------------------------------------------------------
# Network construction and module detection
#------------------------------------------------------------------------------

plotHist <- function(dissTOM, dataType, dend){
  vect <- as.vector(dissTOM)
  
  data_frame(val = vect) %>%
    ggplot(., aes(val)) + 
    geom_histogram(binwidth = 0.002, fill = 'darkmagenta') +
    ggtitle(dend) +
    xlim(c(0.9,1)) +
    xlab('TOM distances') +
    theme(plot.title = element_text(22), axis.title = element_text(20)) +
    theme_light()
  
  name <- paste('Figures/Histogram_',dataType,'_',dend,'.png', sep=''); print(name)
  ggsave(name, height = 5, width = 7)
}

plotTree <- function(dissTOM, dataType, dend){
  
  geneTree <- fastcluster::hclust(as.dist(dissTOM), method = 'average')
  
  name <- paste('Figures/Tree_',dataType,'_',dend,'.png', sep=''); print(name)
  png(name, height = 500, width = 800)
  
  plot(geneTree, xlab = '', sub = '', labels = F, hang = 0.04, main = dend)
  
  dev.off()
  
  return(geneTree)
}

plotdissTOM <- function(geneTree, dissTOM, dataType, dend){
  
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                              deepSplit = 2, cutHeight = 0.995, minClusterSize = 20,
                              pamRespectsDendro = FALSE)
  
  # Modules do not fit with blockwiseModules
  moduleColors = labels2colors(dynamicMods)
  
  # Raise matrix to value for better visualization
  plotTOM = dissTOM^7;
  
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
  
  # Save the plot
  name <- paste('Figures/Heatmap_',dataType,'_',dend,'.png', sep=''); print(name)
  png(name, height = 800, width = 800)
  
  TOMplot(plotTOM, geneTree, moduleColors, main = dend)
  
  dev.off()
}

network_construction <- function(datExpr, corType, power, minMod, dataType){
  print('---------------------------------')
  x <- sprintf("corType: %s, power: %d, minModuleSize: %d",corType,power,minMod)
  print('A network is constructed with the following parameter settings:')
  print(x)
  
  net <- suppressWarnings(blockwiseModules(datExpr, corType = corType, #correlation matrix
                         networkType = "signed", power = power, # adjacency matrix
                         minModuleSize = minMod, TOMType = 'signed', # dendrogram
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = TRUE))
  
  # Print how many modules are produced and how many analytes are in each
  modules <- length(table(net$colors))
  
  # Save information from network
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  
  print(table(moduleColors))
  
  moduleContent = cbind(moduleLabels,moduleColors)
  #MEs = net$MEs;
  # Calculates module eigengenes (MEs) (1st principle component) of modules in data set
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  geneTree = net$dendrograms[[1]]
  
  y <- sprintf("corType: %s, power: %d, minModuleSize: %d\n Modules: %d", 
               corType,power,minMod,modules)
  
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                      main = y)
  name <- sprintf("Figures/04_%s_corType_%s_power_%d_minModuleSize_%d.png",
                  dataType, corType, power, minMod)
  
  png(name, width = 900, height = 550)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module\ncolors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                      main = y,
                      cex.colorLabels = 1.5,
                      cex.axis = 1.5,
                      cex.lab = 1.5)

  dev.off()
  
  return(list(moduleColors, moduleContent, MEs, geneTree))
}

#------------------------------------------------------------------------------
# Module-trait correlation
#------------------------------------------------------------------------------

traitHeatmap <- function(datExpr, datTraits, MEs, moduleColors, minMod, dataType){
  # Define numbers of samples
  nSamples = nrow(datExpr);

  # Correlations between all modules and all traits
  # use = "p" indicates pairwise.complete.obs
  # Pick correlation method (pearson, spearman, kendall)
  moduleTraitCor = cor(MEs, datTraits, use = "p", method = "spearman");
  
  # Calculate p-values 
  moduleTraitPvalue = corPvalueFisher(moduleTraitCor, nSamples);

  # Calculate false discovery rate (FDR) adjusted p-values
  moduleTraitPvalue.FDR <- matrix(p.adjust(moduleTraitPvalue, "fdr"), 
                                  ncol = length(names(datTraits)))
  
  # Matrix of p-value significance asterisks
  asterisk = symnum(moduleTraitPvalue.FDR, corr = FALSE, na = FALSE, 
                    cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                    symbols = c("***", "**", "*", ""))
  textMatrix = as.matrix(paste(asterisk, sep=""))
  dim(textMatrix) = dim(moduleTraitCor)
  
  mods <- length(unique(moduleColors))
  
  name = paste('Figures/04_ModuleTraitHeatmap_',dataType,'_',minMod,'.png',sep='')
  print(name)
  
  png(name, width = 6, height = mods/3+1.5, units="in", res = 1000)
  par(mar = c(6, 7.5, 1, 1));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 textAdj = c(0.5, 0.8),
                 zlim = c(-1,1),
                 cex.lab = 0.9,
                 # lines between cells
                 horizontalSeparator.y = c(0:ncol(MEs)),
                 horizontalSeparator.col = "white",
                 horizontalSeparator.lwd = 2,
                 horizontalSeparator.ext = 0,
                 verticalSeparator.x = c(0:ncol(datTraits)),
                 verticalSeparator.col = "white",
                 verticalSeparator.lwd = 2,
                 verticalSeparator.ext = 0
                 )
                 #main = paste(dataType, ", minModuleSize = ", minMod, sep=''))
  dev.off()
}

interactiveHeatmap <- function(datExpr, datTraits, MEs, minMod, dataType){
  moduleTraitCor = cor(MEs, datTraits, use = "p", method = "spearman");
  nSamples = nrow(datExpr)
  
  # Calculate p-values 
  moduleTraitPvalue = corPvalueFisher(moduleTraitCor, nSamples);
  
  # Calculate false discovery rate (FDR) adjusted p-values
  moduleTraitPvalue.FDR <- matrix(p.adjust(moduleTraitPvalue, "fdr"), 
                                  ncol = length(names(datTraits)))
  
  # Matrix of p-value significance asterisks
  asterisk = symnum(moduleTraitPvalue.FDR, corr = FALSE, na = FALSE, 
                    cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                    symbols = c("***", "**", "*", ""))
  textMatrix = as.matrix(paste(asterisk, sep=""))
  dim(textMatrix) = dim(moduleTraitCor)
  
  # Matrix for interactive hovering
  hoverMatrix = paste("q-value: ",signif(moduleTraitPvalue.FDR,2), sep="")
  dim(hoverMatrix) = dim(moduleTraitCor)
  
  # List of module colors
  moduleNames0 <- substring(names(MEs), 3)
  moduleNames <- data.frame(" "=moduleNames0,check.names=F)
  
  # Make interactive heatmap
  interactiveHeatmap <- heatmaply(moduleTraitCor,
                                  label_names = c("Module","Clinical trait","Correlation"),
                                  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                                    low = "dodgerblue", 
                                    high = "red",
                                    limits = c(-1,1)),
                                  RowSideColors = moduleNames, 
                                  Colv = F,
                                  cellnote = textMatrix,
                                  cellnote_textposition = "bottom center",
                                  custom_hovertext = hoverMatrix,
                                  grid_color = 'white',
                                  file = paste('Figures/ModuleTraitHeatmap_interactive_',dataType,'_',minMod,'.png',sep=''),
                                  margins = c(t = 5, r = 1, b = 2, l = 2)
  ) %>% layout(width=700, height=135+((ncol(MEs)/3)*90))
  #interactiveHeatmap
  
  name = paste('Figures/04_ModuleTraitHeatmap_interactive_',dataType,'_',minMod,'.html',sep='')
  saveWidget(interactiveHeatmap, file = name, knitrOptions = list(width=700, height=135+((ncol(MEs)/3)*90)))
  
  return(interactiveHeatmap)
}

#------------------------------------------------------------------------------
# Module membership plots
#------------------------------------------------------------------------------

moduleMembershipPlot <- function(datExpr, MEs, moduleColors, minMod, dataType, plot=F){
  # Make matrix of correlations between all analytes and the module eigengene of each module
  # A measure of intramodular connectivity
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p", method = "spearman")) %>%
    cbind(., moduleColors)
  
  # Make list of module names
  modNames = substring(names(MEs), 3)
  #print(modNames)
  
  MMmax = vector() # List for max MM values
  MMplots = list() # List for MM plots
  for (module in modNames){
    column = match(module, modNames);
    
    # Define dataframe with MM values per module
    df = geneModuleMembership %>%
      filter(moduleColors == module) %>%
      .[column] %>%
      setNames(., 'MM') %>%
      mutate_all(., round, 1) %>%
      abs(.)
    
    if (module != "grey") {
      MMmax <- c(MMmax, max(df$'MM'))
    }
    
    # Make individual plots
    MMplot <- ggplot(df) + 
      geom_bar(mapping = aes(x = MM), fill = module) +
      xlab("Module membership") +
      ylab("Count") +
      ggtitle(paste(module, " n = ", length(df$MM), sep = '')) +
      theme_light() +
      theme(plot.title = element_text(size=13))
    
    # Save plot in list
    MMplots[[module]] <- MMplot
  }
  
  if (plot) {
  rows = ceiling(length(modNames)/5)
  grid <- plot_grid(plotlist = MMplots, ncol = 5, nrow = rows)
  grid
  
  name = paste('Figures/04_ModuleMemberships_',dataType,'_',minMod,'.png', sep='')
  print(name)
  
  ggsave(name, height = rows*2, width = 12)
  }
  
  return(MMmax)
}

plotAvgMM <- function(numMod, MMavgAll, dataType){
  minModSize <- seq(5,20,by=1)
  
  # Make matrix for plotting
  MMmat0 <- cbind(numMod, MMavgAll)
  MMmat <- cbind(minModSize,MMmat0)
  
  if (dataType == "CS_met") {
    title = "Metabolomics" }
  if (dataType == "CS_lip") {
    title = "Lipidomics" }
  
  ggplot(data.frame(MMmat), aes(x = numMod, y = MMavgAll)) +
    geom_point() +
    xlab('Number of modules') +
    ylab('Average max module membership') +
    #title(title) +
    title(sprintf('%s module membership',title)) +
    theme_light() +
    scale_x_continuous(breaks = seq(0, 55, 5))
  
  # Save plot
  name = paste('Figures/04_AverageModuleMembership_',dataType,'.png', sep='')
  print(name)
  ggsave(name, height = 4, width = 6)
  
  return(MMmat)
}

#------------------------------------------------------------------------------
# Annotation
#------------------------------------------------------------------------------

readAnnotationFile <- function(annotFile, datExpr){
  # Analyte IDs and corresponding super- and subpathways
  annot = suppressMessages(read_tsv(annotFile)) %>%
    mutate(SUPER_PATHWAY = replace_na(SUPER_PATHWAY, "Unknown")) %>%
    mutate(SUB_PATHWAY = replace_na(SUB_PATHWAY, "Unknown"))
  IDs = names(datExpr)
  IDs2annot = match(IDs, annot$BIOCHEMICAL)
  return(list(annot, IDs2annot))
}

makeSuperPathwayDataframe <- function(annot, IDs2annot, MEs, dataType){
  # Extract super-pathways
  SUP_PATH = annot$SUPER_PATHWAY[IDs2annot]
  length(unique(SUP_PATH))
  
  # Make dataframe of amount of analytes in each pathway
  if (dataType == "CS_met") {
    SUP_100 <- data.frame(table(SUP_PATH))[c(1:8,10,9),] %>%
      unite("Label", c("SUP_PATH", "Freq"), sep = " n=", remove = FALSE) # Make new column
  }
  if (dataType == "CS_lip") {
    SUP_100 <- data.frame(table(SUP_PATH))[c(1:11),] %>%
      unite("Label", c("SUP_PATH", "Freq"), sep = " n=", remove = FALSE) # Make new column
  }
  names(SUP_100) = c('Label','SUP_PATH','All') # Rename columns

  # Make dataframe of pathway counts per module
  df <- data.frame()
  for (module in substring(names(MEs), 3)){
    # Filter for analytes in module
    modGenes = (moduleColors==module)
    
    # Get super pathways for analytes in module
    modSUP = SUP_PATH[modGenes]
    
    # Make temporary dataframe for pathways and counts of all their analytes
    t = SUP_100
    rownames(t) = SUP_100$SUP_PATH
    
    # Make temporay dataframe with counts for each pathway
    temp = data.frame(table(modSUP))
    rownames(temp) = temp$modSUP
    
    # Merge the two temporary files by pathway to get all pathways for each module
    if (dataType == "CS_met") {
    temp2 = merge(t, temp, by=0, all=T) %>% .[-c(1,4)] %>% .[c(1:8,10,9),] %>% # [c(1:8,10,11,9),] for lipids
      # Add percentage column for each pathway
      mutate(Percent = round(Freq/t$All*100,0)) %>%
      # Add which module the data comes from
      mutate(Module = module)
    }
    if (dataType == "CS_lip") {
      temp2 = merge(t, temp, by=0, all=T) %>% .[-c(1,4)] %>% .[c(1:11),] %>%
        # Add percentage column for each pathway
        mutate(Percent = round(Freq/t$All*100,0)) %>%
        # Add which module the data comes from
        mutate(Module = module)
    }
    # Replace NA with 0
    suppressWarnings(temp2[is.na(temp2)] <- 0) # Makes warnings, but still works
    
    # Row bind to form one big dataframe
    df = rbind(df, temp2)
    
  }
  return(list(df, SUP_PATH))
}

makeSubPathwayDataframe <- function(annot, IDs2annot){
  # Extract super-pathways
  SUB_PATH = annot$SUB_PATHWAY[IDs2annot]
  length(unique(SUB_PATH))
  
  # Make dataframe of amount of analytes in each pathway
  SUB_100 <- data.frame(table(SUB_PATH))[c(1:93,95:98,94),] %>% # Place Unkown last
    unite("Label", c("SUB_PATH", "Freq"), sep = " n=", remove = FALSE) # Make new column
  names(SUB_100) = c('Label','SUB_PATH','All') # Rename colums
  
  # Make dataframe of pathway counts per module
  df <- data.frame()
  for (module in unique(moduleColors)){
    # Filter for analytes in module
    modGenes = (moduleColors==module)
    
    # Get sub pathways for analytes in module
    modSUB = SUB_PATH[modGenes]
    
    # Make temporary dataframe for pathways and counts of all their analytes
    t = SUB_100
    rownames(t) = SUB_100$SUB_PATH
    
    # Make temporay dataframe with counts for each pathway
    temp = data.frame(table(modSUB))
    rownames(temp) = temp$modSUB
    
    # Merge the two temporary files by pathway to get all pathways for each module
    temp2 = merge(t, temp, by=0, all=T) %>% .[-c(1,4)] %>% .[c(1:93,95:98,94),] %>%
      # Add percentage column for each pathway
      mutate(Percent = round(Freq/t$All*100,0)) %>%
      # Add which module the data comes from
      mutate(Module = module)
    
    # Replace NA with 0
    suppressWarnings(temp2[is.na(temp2)] <- 0) # Makes warnings, but still works
    
    # Row bind to form one big dataframe
    df = rbind(df, temp2)
    
  }
  return(list(df, SUB_PATH))
}

#------------------------------------------------------------------------------
# Pathway plots
#------------------------------------------------------------------------------

modsInPathwaysPlot <- function(df, SUP_PATH, minMod, dataType){
  # Make plots of module distributions in pathways
  plots <- list()
  for (pathway in unique(df$SUP_PATH)){
    # Subset dataframe to only include one pathway
    datasub = subset(df, SUP_PATH == pathway)

    # Make plot
    p <- ggplot(data=datasub, aes(x=Module, y=Percent, fill=Module)) +
      geom_bar(stat="identity", fill = datasub$Module) +
      #geom_text(aes(label=Percent), vjust=-0.3, size=3.5) +
      ggtitle(pathway) +
      xlab('') +
      theme_classic() +
      #scale_y_continuous(expand = c(0,0)) + # Remove space between axis and bars
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
            legend.position = 'bottom', plot.title = element_text(size=9)) +
      scale_y_continuous("%", limits = c(0,max(datasub$Percent+5)))
    
    # Save plot to list
    plots[[pathway]] = p
  }
  
  # Plot all pathway barplots in grid (9 pathways)
  grid <- plot_grid(plotlist = plots, ncol = 5)
  
  # Make title for collected plot
  title <- ggdraw() + 
    draw_label(paste(dataType, "minModuleSize = ", minMod), fontface = 'bold', x = 0, hjust = 0, size = 10) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  plot_grid(grid, ncol=1, rel_heights = c(0.1, 1))
  
  name = paste('Figures/04_ModulesInPathways_',dataType,'_',minMod,'.png', sep='')
  print(name)
  
  ggsave(name, height = ceiling(length(plots)/5)*2.5, width = 10)
}

pathwaysInModsPlot <- function(df, moduleColors, minMod, dataType){
  # Make list of plots of pathway distributions in modules
  plots <- list()
  for (module in unique(df$Module)){
    # Subset dataframe to only include one pathway
    datasub = subset(df, Module == module)

    # Make plot
    p <- ggplot(data=datasub, aes(x=factor(Label, unique(Label)), y=Percent, fill=factor(Label, unique(Label)))) +
      geom_bar(stat="identity") +
      scale_fill_brewer(palette="Paired") +
      ggtitle(paste(module,' n=',length(subset(moduleColors, moduleColors==module)), sep='')) + xlab('') +
      theme_light() +
      labs(fill = '') +
      theme(axis.text.x = element_blank(), 
            axis.ticks = element_blank(), 
            plot.title = element_text(size=15),
            legend.position = 'none') +
      scale_y_continuous("%", limits = c(0,max(datasub$Percent+5)))

    # Save plot to list
    plots[[module]] = p
  }
  
  # Wide
  col = 5
  legend_p <- get_legend(p + theme(legend.position="bottom", legend.text=element_text(size=15),
                                   legend.spacing.x = unit(0.4, 'cm')) + 
                           guides(fill=guide_legend(ncol=col)))
  rows = ceiling(length(unique(moduleColors))/col)
  grid_wide <- plot_grid(plotlist = plots, ncol = col)
  plot_grid(grid_wide, legend_p, ncol = 1, rel_heights = c(1, .2))
  name = paste('Figures/04_PathwaysInModules_',dataType,'_',minMod,'_wide.png', sep='')
  print(name)
  ggsave(name, height = rows*2+1, width = col*3.2)
  
  # Long
  col = 2
  legend_p <- get_legend(p + theme(legend.position="bottom", legend.text=element_text(size=13),
                                   legend.spacing.x = unit(0.4, 'cm')) + 
                           guides(fill=guide_legend(ncol=col)))
  rows = ceiling(length(unique(moduleColors))/col)
  grid_long <- plot_grid(plotlist = plots, ncol = col)
  plot_grid(grid_long, legend_p, ncol = 1, rel_heights = c(1, .2))
  name = paste('Figures/04_PathwaysInModules_',dataType,'_',minMod,'_long.png', sep='')
  print(name)
  ggsave(name, height = rows*2+1, width = col*3.5)
}

pathwaysIn1ModPlot <- function(df, moduleColors, minMod, dataType, module){
  # Select data for selected module
  datasub = subset(df, Module == module)
  print(datasub)
  
  p <- ggplot(data=datasub, aes(x=factor(Label, unique(Label)), y=Percent, fill=factor(Label, unique(Label)))) +
    geom_bar(stat="identity") +
    #scale_fill_brewer(palette="Paired") +
    #geom_text(aes(label=Freq), vjust=-0.3, size=3) +
    ggtitle(paste(module,'n =',length(subset(moduleColors, moduleColors==module)))) + xlab('') +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(fill = '') +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size=10)) +
    scale_y_continuous("%", limits = c(0,max(datasub$Percent+5)), expand = c(0,0))
  p
  
  legend_p <- get_legend(p + theme(legend.position="right"))
  
  # Make title for plot
  title <- ggdraw() + 
    draw_label(paste(dataType, "subpathways, minModuleSize =", minMod), fontface = 'bold', x = 0, hjust = 0, size = 10) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  p1 <- plot_grid(title, p, ncol=1, rel_heights = c(0.1, 1))
  p1
  
  plot_grid(p1, legend_p, ncol = 2, rel_widths = c(1,4))
  
  name = paste('Figures/04_PathwaysIn1Module_',dataType,'_',minMod,'.png', sep='')
  print(name)
  
  ggsave(name, height = 7, width = 25)

  
  plot_grid(p1, legend_p, ncol = 2, rel_heights = c(1, 1), rel_widths = c(1,))
  
  name = paste('Figures/04_PathwaysIn1Module_',dataType,'_',minMod,'.png', sep='')
  print(name)
  
  ggsave(name, height = 2+1, width = 10)
}

#------------------------------------------------------------------------------
# Cytoscape
#------------------------------------------------------------------------------

exportCytoscape <- function(TOM, power, minMod, dataType, threshold = 0.2){
  name1 = paste('Networks/04_',dataType,'_cytoscape_edge_',minMod,'.txt', sep='')
  name2 = paste('Networks/04_',dataType,'_cytoscape_node_',minMod,'.txt', sep='')
  print(name1)
  print(name2)
  
  # Export TOMs to Cytoscape
  (exportNetworkToCytoscape(TOM, 
                           edgeFile = name1,
                           nodeFile = name2,
                           nodeAttr = moduleColors, altNodeNames = SUP_PATH, weighted = T,
                           #altNodeNames = moduleColors, 
                           threshold = threshold, 
                           nodeNames = names(datExpr)))
}

#------------------------------------------------------------------------------
# Make lists of IDs per module
#------------------------------------------------------------------------------

exportModuleIDs <- function(annot, IDs2annot, minMod, dataType){
  # Get the corresponding HMDB IDs
  BIOCHEM = annot$BIOCHEMICAL[IDs2annot]
  SUP     = annot$SUPER_PATHWAY[IDs2annot]
  SUB     = annot$SUB_PATHWAY[IDs2annot]
  HMDBIDs = annot$HMDB_ID[IDs2annot]
  KEGGIDs = annot$KEGG_ID[IDs2annot]
  
  HMDB_list = list()
  for (module in unique(moduleColors)){
    modGenes = (moduleColors==module)
    modBIOCHEM = BIOCHEM[modGenes]
    modSUP     = SUP[modGenes]
    modSUB     = SUB[modGenes]
    modHMDBIDs = HMDBIDs[modGenes]
    modKEGGIDs = KEGGIDs[modGenes]
    HMDB_list[[module]] = data.frame("BIOCHEMICAL" = modBIOCHEM,
                                     "SUPER_PATHWAY" = modSUP,
                                     "SUB_PATHWAY" = modSUB,
                                     "HMDB_ID" = modHMDBIDs,
                                     "KEGG_ID" = modKEGGIDs)
  }
  
  # Write excel file with each module in it's own sheet
  name = paste('Data/04_Modules_',dataType,'_',minMod,'.xlsx', sep='')
  print(name)
  openxlsx::write.xlsx(HMDB_list, file = name)
}
