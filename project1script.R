source("prj1functions.R")

genes <- read.table(gzfile('gencode.v29.annotation.gff3.gz'))  
cancer_genes <- read.table('NCG6_cancergenes.tsv', stringsAsFactors = FALSE, fill = TRUE, nrows =-1, skip = 1)
cytoBands <- read.table('cytoBand.txt', stringsAsFactors = FALSE, fill = TRUE) 

# makes a smaller dataframe of just genes
formattedIndeces <- which(as.character(genes$V3)=="gene")

lengthVector <- vector()
nameVector <- vector()
chromosomeEndIndeces <- vector()
chromosomeStartIndeces <- vector()
# index represents chromosome
geneNumber <- vector()
geneStartEndIndeces <- list()
nameListVector <- list()
# global variable for use inside loop
assign("pattern", "gene_name=([^;]*)", envir = .GlobalEnv)
# manual run prior to loop
chromosomeLabels <- vector()
# initial value to overwrite
chromosomeLabel <- "blank"
# chromosome index
i = 0
chromosomalGenesFileIndeces <- list()

for (rowNum in 1:length(formattedIndeces)){
  # when there is an original chromosome read (i.e., 25x)
  if (as.character(genes[formattedIndeces[rowNum],1]) != chromosomeLabel){
    # increments chromosome index
    i = (i + 1)
    print(paste("Loop 1/2: chromosome",i,"/",25))
    # adds the name to a vector for later use
    chromosomeLabel <- as.character(genes[formattedIndeces[rowNum],1])
    chromosomeLabels <- c(chromosomeLabels, chromosomeLabel)
                   
    # gets ALL file indeces for the chromosome whose iteration the loop is on, not just genes
    chromosomeFileIndeces <- which(as.character(genes[,1]) == tail(chromosomeLabels, n=1))
    
    chromosomeI <- genes[chromosomeFileIndeces,4:5]
    # former index 4 is now 1, 5 is now 2, 9 is 6
    chromosomeIStartIndex <- min(as.numeric(chromosomeI[,1]))
    chromosomeIEndIndex <- max(as.numeric(chromosomeI[,2]))

    chromosomeEndIndeces <- c(chromosomeEndIndeces, chromosomeIEndIndex)
    chromosomeStartIndeces <- c(chromosomeStartIndeces, chromosomeIStartIndex)

    # gets the indeces of gene file that are actually genes AND the same chromosome
    chromosomalGenesFileIndeces[[i]] <- (intersect(formattedIndeces, chromosomeFileIndeces))
    # adds the number of genes in the vector index associated w/ the chromosome
    geneNumber <- c(geneNumber, length(chromosomalGenesFileIndeces[[i]]))  
    
    geneStartIndeces <- as.numeric(genes[chromosomalGenesFileIndeces[[i]], 4])
    geneEndIndeces <- as.numeric(genes[chromosomalGenesFileIndeces[[i]], 5])
    # adds all the chromosome's gene starts and ends together as a single element of a list
    geneStartEndIndeces[[i]] = sort(append(geneStartIndeces,geneEndIndeces))
    
    m <- regexec(pattern, genes[chromosomalGenesFileIndeces[[i]],9])
    geneNames <- sapply(regmatches(genes[chromosomalGenesFileIndeces[[i]], 9], m), function(e){return(e[2])})
    nameListVector[[i]] <- geneNames
  } # end if
} # end for

cancer_gene_name_vector <- unique(cancer_genes$V2)

ckListVector <- list()
zoneBounds <- list()
chromosomeLengths <- chromosomeEndIndeces - chromosomeStartIndeces
cancerGeneList <- list()
cancerGeneIndeces <- list()
genesInZones <- list()
cancerGenesInZones <- list()
cancerStartIndeces <- list()
cancerStopStartIndeces <- list()

for (i in 1:length(chromosomeLabels)){
  print(paste("Loop 2/2: chromosome",i,"/",25))
#single iteration test
#for (i in 1){  
  # resets each iteration
  k <- vector()
  kmin <- 2
  kmax <- max(20, min(geneNumber[i]/5, chromosomeLengths[i]/1000000))
  k <- c(kmin, kmax)
  # print(k)
  ckResult <- Ckmeans.1d.dp(geneStartEndIndeces[[i]],k)
  ckListVector[[i]] <- (as.numeric(unlist(ckResult[2])))
  
  zoneBounds[[i]] <- vector()
  # runs the function to change centers to intermediate boundary values
  insideBound <- c(zoneBounds[[i]], boundaries(ckListVector[[i]]))
  # adds the first index of the chromosome for starting bound
  zoneBounds[[i]] <- c(chromosomeStartIndeces[i], insideBound)
  # adds last index of the chromosome for the ending bound
  zoneBounds[[i]] <- c(zoneBounds[[i]], chromosomeEndIndeces[i])

  # get indeces of all cancer genes on chromosome
  
  # name of all the cancer genes on each chromosome
  cancerGeneList[[i]] <- intersect(nameListVector[[i]], cancer_gene_name_vector)
  # gets indeces of matches from nameListVector, to apply to chromosomalGenesFileIndeces (elements correspond to those of nameListVector)
  cancerGeneIndeces[[i]] <- match(cancerGeneList[[i]], nameListVector[[i]]) 
  # cancerGeneIndeces pulls a subset from file indeces for each gene to pull from data table
  cancerStartIndeces[[i]] <- genes[chromosomalGenesFileIndeces[[i]][cancerGeneIndeces[[i]]],4]
  cancerStopStartIndeces[[i]] <- c(cancerStartIndeces[[i]], genes[chromosomalGenesFileIndeces[[i]][cancerGeneIndeces[[i]]],5])

  # makes list of vectors for #of genes in zones, where zones are element indeces to chromosome vectors
  genesInZones[[i]] <- fixedZone(geneStartEndIndeces[[i]], zoneBounds[[i]])
  
  # because if there aren't any cancer genes, the function will break (for mitochondrial chromosome, #25)
  if (length(cancerGeneIndeces[[i]] > 0)){
    cancerGenesInZones[[i]] <- fixedZone(sort(cancerStopStartIndeces[[i]]), zoneBounds[[i]])
  } # end if  
  else{
    # adds a zero to each zone (to fill the vector)
    for (zones in 1:(length(zoneBounds[[i]]) - 1)){
      cancerGenesInZones[[i]] <- rep(0,length(zoneBounds[[i]]) - 1)
    } # end for  
  } # end else
} # end for

pVals <- list()
pValsAdj <- list()
# list elements represent chromosomes, elements of each vector the indeces of enriched zones
enrichedZones <- list()

# sum of all zones in genome
theSum <- 0
for (i in 1:length(genesInZones)){
  if (length(genesInZones[[i]]) > 0){
    theSum <- theSum + length(genesInZones[[i]])
  }
}
# gets avg. number of zones per gene, the factor for p-value adjustment
theFactor <- theSum/length(genesInZones)

for (i in 1:length(chromosomeLabels)){
  pVals[[i]] <- mapply(function(cancerGeneZoneCount, geneZoneCount)
    # first argument is total # of cancer genes on genome, last is total # of general genes on genome
    phyper(cancerGeneZoneCount,length(cancer_gene_name_vector),length(formattedIndeces),geneZoneCount, lower.tail = FALSE),
    # vectors whose elements are fed to mapply function (iterates in parallel)
    cancerGenesInZones[[i]], genesInZones[[i]])
  # multiple testing adjustment per chromosome defaults to divisor of chromosome's zone count 
  # pValsAdj[[i]] <- p.adjust(pVals[[i]], method = "bonferroni")
  # custom multiple testing adjustment -> uses average of zone count per chromosome (theFactor)
  pValsAdj[[i]] <- pVals[[i]]*theFactor
  enrichedZones[[i]] <- which(0.05 > pValsAdj[[i]])
}

final <- data.frame(matrix(ncol = 12, nrow = 0))
columnNames <- c("seqid", "cband.start", "cband.end", "zid", "nzones", "start", "end", "total", "observed", "pval", "pval.adj","contains")
colnames(final) <- columnNames

# makes vector from list (unlist splices out 0s). NOTE: only viable because no chromosomes have more than a single enriched zone
enrichedZoneVector <- vector()
for (i in 1:length(enrichedZones)){
  # for non-enriched chromosomes (since concatenation doesn't work)
  if (length(enrichedZones[[i]]) < 1){
    enrichedZoneVector <- c(enrichedZoneVector, 0)
  }
  else{
    enrichedZoneVector <- c(enrichedZoneVector,enrichedZones[[i]])
  }
}

testFrame <- testOutput(enrichedZones, final)
print(testFrame)

value1=testFrame[,6]
value2=testFrame[,7]
plotLengths <- c(chromosomeLengths[c(1,3,6,6,7,11,11,17,20,22,23)])
plotFrame=data.frame(x=testFrame[,1], 
                     value1 = as.numeric(as.character(value1))/1000000, 
                     value2 = as.numeric(as.character(value2))/1000000,
                     value3 = as.numeric(as.character(plotLengths/1000000)),
                     value4 = as.numeric(as.character(rep(0,11))))

plotFrame$x <- factor(plotFrame$x, levels = c('chr1','chr3','chr6','chr7','chr11','chr17','chr20','chr22','chrX'))

labelXCoordinates <- getXlabel(testFrame)
labelYCoordinates <- getYlabel(testFrame, labelXCoordinates)
labelContent <- getContent(testFrame)

ggplot(plotFrame)+
  geom_segment(aes(x=x, xend=x, y=value4, yend=value3), color = "gray")+
  geom_segment(aes(x=x, xend=x, y=value1, yend=value2), color = "red")+
  xlab("chromosome")+
  ylab("index on chromosome (by million)")+
  ggtitle("Cancer Enriched Zones")+
  scale_y_continuous(breaks = seq(20,260, by=20), expand=expand_scale(mult = c(0, .05)))+
  scale_x_discrete(breaks = c('chr1','chr3','chr6','chr7','chr11','chr17','chr20','chr22','chrX'))+
  annotate('text', y = labelYCoordinates, x = labelXCoordinates, label = labelContent, size = 3)
