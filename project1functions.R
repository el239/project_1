library(Ckmeans.1d.dp)
library(ggplot2)
library(dplyr)
#library(tidyverse)

# function to get string for zone labels
getContent <- function(testFrame){
  stringVector <- vector()
  for(i in 1:nrow(testFrame)){
    label <- paste(as.character(testFrame[i,9]),as.character(testFrame[i,8]), sep = "/")
    stringVector <- c(stringVector, label)
  } # end for
  return (stringVector)
}

# function to get x-axis coordinates for labels
getXlabel <- function(dataFrame){
  xLabelVector <- vector()
  previousChr <- "default"
  count <- 0
  for (i in 1:nrow(dataFrame)){
    if (dataFrame[i,1] != previous){
      count <- count + 1
      xLabelVector <- c(xLabelVector,(count + 0.25))
    } # end if
    else{
      # repeat element (same column on graph)
      xLabelVector <- c(xLabelVector,(count +0.25))
    }
    previous <- dataFrame[i,1]
  } # end for
  return (xLabelVector)
} # end function

# function to get y-axis coordinates for labels 
# checks for closely neighboring zones for output adjustment
getYlabel <- function(dataFrame, xCoordinateVector){
  yLabelVector <- vector()
  xPrevious <- 0
  # default value (longer than longest chromosome to guarantee change)
  yPrevious <- 500
  for (i in 1:nrow(dataFrame)){
    # finds zone midpoint
    average <- (as.numeric(dataFrame[i,6])+as.numeric(dataFrame[i,7]))/2
    # scales to graph
    middle <- average/1000000
    # if on same chromosome, adjusts labels if too close by checking previous values for y
    if ((xCoordinateVector[i] == xPrevious) && (abs(middle - yPrevious) < 10)){
      # finds degree of adjustment
      adjust <- (10 - abs(middle - yPrevious))
      # bumps appropriate zone up or down
      if (yPrevious > middle){
        yLabelVector[i-1] <- (yPrevious + (adjust/2))
        yLabelVector <- c(yLabelVector, (middle - (adjust/2)))
      } # end if
      else{
        yLabelVector[i-1] <- (yPrevious - (adjust/2))
        yLabelVector <- c(yLabelVector, (middle + (adjust/2)))
      } # end else
    } # end if  
    else{
      yLabelVector <- c(yLabelVector, middle)
    } # end else
    xPrevious <- xCoordinateVector[i]
    yPrevious <- middle
  } # end for
  return (yLabelVector)
} # end function

finishedOutput <- function(enrichedZoneVector, dataFrame){
  rowNumber <- 0
  for (i in 1:length(enrichedZoneVector)){
    # i.e., for chromosomes with enriched zones:
    if (enrichedZoneVector[i] > 0){
      dataVector <- vector()
      # determines where to add to the dataframe
      rowNumber <- rowNumber + 1
      # seqid
      dataVector <- c(chromosomeLabels[i])
      # cband.start and cband.end, i is chromosome #
      dataVector <- c(dataVector, getBands(i,enrichedZoneVector[i]))
      # enrichedVector[i] is the enriched zone id
      dataVector <- c(dataVector, enrichedZoneVector[i])
      # just to fetch nzones, genes not important here
      dataVector <- c(dataVector, length(genesInZones[[i]]))
      # start (chromosomal start index of zone, found via indexing zone id)
      zoneStart <- round(zoneBounds[[i]][enrichedZoneVector[i]])
      dataVector <- c(dataVector, zoneStart)
      # end (chromosomal end index of zone, found via indexing zone id)
      zoneEnd <- round(zoneBounds[[i]][enrichedZoneVector[i] + 1])
      dataVector <- c(dataVector, zoneEnd)
      # total (genes in zone)
      dataVector <- c(dataVector, genesInZones[[i]][enrichedZoneVector[i]])
      # observed (cancer genes in zone)
      dataVector <- c(dataVector, cancerGenesInZones[[i]][enrichedZoneVector[i]])
      # pval
      dataVector <- c(dataVector, formatC(pVals[[i]][enrichedZoneVector[i]], format = "e", digits = 2))
      # adj.pval
      dataVector <- c(dataVector, formatC(pValsAdj[[i]][enrichedZoneVector[i]], format = "e", digits = 2))
      # contains (cancer gene names)
      dataVector <- c(dataVector, getCancerGeneNames(i, enrichedZoneVector[i], zoneStart, zoneEnd))
    } # end if
    dataFrame[rowNumber,] <- dataVector
  } # end for
  # print(dataVector)
  return (dataFrame)
} # end function

getCancerGeneNames <- function(chromosome, zoneID, zoneStart, zoneEnd){
  sortedCancerIndeces <- sort(cancerStopStartIndeces[[chromosome]])
  cancerGeneNames <- vector()
  # checks every odd index (a.k.a, the start index)
  for (j in seq(1,length(sortedCancerIndeces)-1, by=2)){
    # if either the start or concomitant stop chromosomal cancer index is within the zone:
    if (((sortedCancerIndeces[j] > zoneStart) && (sortedCancerIndeces[j] < zoneEnd)) || ((sortedCancerIndeces[j + 1] > zoneStart) && (sortedCancerIndeces[j + 1] < zoneEnd))){
      # adds the corresponding name to return vector (cancerGeneList list indeces is half the length of sortedCancerIndeces & names are in order of appearance)
      cancerGeneNames <- c(cancerGeneNames,cancerGeneList[[chromosome]][ceiling(j/2)])
    } # end if
  } # end for
  # because the return must be a single vector for data frame formatting, makes a string
  return (paste(cancerGeneNames,collapse="; "))
} # end function

getBands <- function(chromosome, zoneIndex){
  # for right chromosome:
  chrString <- chromosomeLabels[chromosome]
  # only need zone id, so overall file index is not importants
  indecesSubset <- which(chrString == cytoBands$V1)
  # makes new matrix specific to chromosome
  cytoSubset <- cytoBands[indecesSubset,]
  # extracts vector of chromosomal end indeces for cancer genes in chromosome
  geneEndsVector <- cytoSubset[,3]
  # gets chromosomal index of zone start
  zoneStart <- (zoneBounds[[chromosome]][zoneIndex])
  # gets index of cytoSubset for band containing zone start
  subsetStartIndex <- logSearch(zoneStart,geneEndsVector)
  # gets the name using the index
  startBandID <- cytoSubset[subsetStartIndex,4]
  
  zoneEnd <- (zoneBounds[[chromosome]][zoneIndex + 1])
  subsetEndIndex <- logSearch(zoneEnd,geneEndsVector)
  endBandID <- cytoSubset[subsetEndIndex,4]
  return (c(startBandID,endBandID))
}    

boundaries <- function(clusterCentersVector){
  boundaries <- vector()
  for (i in 1:(length(clusterCentersVector) - 1)){
    # gets the average of every consecutive 2 centers
    boundaries <- c(boundaries, (clusterCentersVector[i] + clusterCentersVector[i + 1])/2)
  }
  return (boundaries)
}

logSearch <- function(number, vector){
  leftBorder <- 1
  rightBorder <- length(vector)
  if (number <= min(vector)){
    return (1)
  }
  if (number >= max(vector)){
    return (length(vector))
  }
  while (leftBorder <= rightBorder){
    # rounds halves down, but that's okay
    median <- round((leftBorder + rightBorder)/2)
    if (vector[median] < number){
      # adjusts left margin for new median next iteration
      leftBorder <- (median + 1)
    } # end if
    # i.e., if baseEnds[median] >= Qstart
    else{
      # adjusts right margin for new median next iteration
      rightBorder <- (median - 1)
    } # end else
  } # end while
  return (rightBorder + 1)
} # end function

# unknown malfunction for chr14; replaced by fixedZone function
zoneTally <- function(indecesVector, zoneBoundsVector){
  j <- 1
  startIn <- TRUE
  zoneGeneCount <- vector()
  geneCount <- 0
  for (i in 1:length(indecesVector)){
    # if the index is within the zone
    if ((indecesVector[i] <= zoneBoundsVector[j + 1]) && (indecesVector[i] >= zoneBoundsVector[j])){
      # for odd indeces (start coordinates)
      if (i %% 2 != 0){
        geneCount <- geneCount + 1
        startIn <- TRUE
      } # end if
      
      # for even indeces (end coordinates)
      else{
        # if the start coordinate is inside, the gene has already been counted
        if (startIn == FALSE){
          geneCount <- geneCount + 1
        } # end if
      } # end else
    } # end if  
    
    # if the index is not inside the zone:
    else{
      # increments zone index, needs loop in case zones are skipped 
      #print(paste("j is: ", j))
      #print(paste("i is: ", i))
      while (!((indecesVector[i] <= zoneBoundsVector[j + 1]) && (indecesVector[i] >= zoneBoundsVector[j]))){
        # pushes the old count onto vector to be saved
        zoneGeneCount <- c(zoneGeneCount, geneCount)
        
        # if we're on a start index, the empty zones need to be marked blank:
        if (i %% 2 != 0){
          # resets the count for new zone
          geneCount <- 0
        } # end if
        
        # otherwise, we're on an end index, and the gene is spanning multiple zones, each with count 1:
        else{
          # count for current zone
          geneCount <- 1
        } # end else
        # skips to next zone
        j <- j + 1
      } # end while
      
      # for odd indeces (start coordinates)
      if (i %% 2 != 0){
        # tallies the count for the new zone
        geneCount <- geneCount + 1
        startIn <- TRUE
      } # end if
      # for even indeces (end coordinates)
      else{
        startIn <- FALSE
      } # end else
    } # end else
  } # end for
  # for last index (an end)
  zoneGeneCount <- c(zoneGeneCount,geneCount)
  # fills ending indeces with zero to full length of all zones
  while (length(zoneGeneCount) < (length(zoneBoundsVector) - 1)){
    zoneGeneCount <- c(zoneGeneCount,0)
  } # end while
  return (zoneGeneCount)
} # end function

fixedZone <- function (sortedCancerIndeces, zoneBoundaries){
  cancerTally <- vector()
  # for each zone
  for (i in 1:(length(zoneBoundaries) - 1)){
    cancerCount <- 0
    # for each stop/start index
    for (j in seq(1,length(sortedCancerIndeces)-1, by=2)){
      # if either the start or concomitant stop chromosomal cancer index is within the zone:
      if (((sortedCancerIndeces[j] > zoneBoundaries[i]) && (sortedCancerIndeces[j] < zoneBoundaries[i + 1])) || ((sortedCancerIndeces[j + 1] > zoneBoundaries[i]) && (sortedCancerIndeces[j + 1] < zoneBoundaries[i + 1]))){
        cancerCount = cancerCount + 1
      } # end if
    } # end for
    cancerTally <- c(cancerTally,cancerCount)
  } # end for
  return (cancerTally)
} # end function

testOutput <- function(enrichedZoneList, dataFrame){
  rowNumber <- 0
  for (i in 1:length(enrichedZoneList)){
    # i.e., for chromosomes with enriched zones:
    if (length(enrichedZoneList[[i]]) > 0){
      for (j in 1:length(enrichedZoneList[[i]])){
        dataVector <- vector()
        # determines where to add to the dataframe
        rowNumber <- rowNumber + 1
        # seqid
        dataVector <- c(chromosomeLabels[i])
        # cband.start and cband.end, i is chromosome #
        dataVector <- c(dataVector, getBands(i,enrichedZoneList[[i]][j]))
        # enrichedVector[i] is the enriched zone id
        dataVector <- c(dataVector, enrichedZoneList[[i]][j])
        # just to fetch nzones, genes not important here
        dataVector <- c(dataVector, length(genesInZones[[i]]))
        # start (chromosomal start index of zone, found via indexing zone id)
        zoneStart <- round(zoneBounds[[i]][enrichedZoneList[[i]][j]])
        dataVector <- c(dataVector, zoneStart)
        # end (chromosomal end index of zone, found via indexing zone id)
        zoneEnd <- round(zoneBounds[[i]][enrichedZoneList[[i]][j] + 1])
        dataVector <- c(dataVector, zoneEnd)
        # total (genes in zone)
        dataVector <- c(dataVector, genesInZones[[i]][enrichedZoneList[[i]][j]])
        # observed (cancer genes in zone)
        dataVector <- c(dataVector, cancerGenesInZones[[i]][enrichedZoneList[[i]][j]])
        # pval
        dataVector <- c(dataVector, formatC(pVals[[i]][enrichedZoneList[[i]][j]], format = "e", digits = 2))
        # adj.pval
        dataVector <- c(dataVector, formatC(pValsAdj[[i]][enrichedZoneList[[i]][j]], format = "e", digits = 2))
        # contains (cancer gene names)
        dataVector <- c(dataVector, getCancerGeneNames(i, enrichedZoneList[[i]][j], zoneStart, zoneEnd))
        dataFrame[rowNumber,] <- dataVector
      } # end for
    } # end if
  } # end for
  # print(dataVector)
  return (dataFrame)
} # end function

'
testZones <- c(2,7,19,21,25)
testIndeces <- c(3,4,5,9,10,11,23,25)
testResult <- zoneTally(testIndeces, testZones)
print(testResult)
# 2,2,0,1

testZones2 <- c(1,7,14,17,29)
testIndeces2 <- c(1,2,6,15)
testResult2 <- zoneTally(testIndeces2, testZones2)
print(testResult2)
# 2,1,1,0

testZones3 <- c(1,3,4,10,11,20)
testIndeces3 <- c(4,6,6,7,12,14)
testResult3 <- zoneTally(testIndeces3, testZones3)
print(testResult3)
#0,1,2,0,1


# test: returns index of first element larger than input
logTest <- logSearch(2800,c(1000,1500,3000,6000))
'