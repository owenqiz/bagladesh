#Arsenic data download from http://www.bgs.ac.uk/home.html
#http://www.bgs.ac.uk/downloads/start.cfm?id=2214
library(readr)
require(foreach)

L <- 10

#read csv file for location and arsenic data
#NationalSurveyData <- read.csv(file.choose(),skip=5)[,c(4,5,9,10,11,15)]


#rename variable from extracting file
colnames(NationalSurveyData) <- c("lat","long","division","district","thana","arsenic")

# Find the index place of arsenic in data frame
# Default to set last variable to be arsenic level
arsenicIndex <- length(NationalSurveyData)


#data cleaning and coding replace "<.5" by ".5" and replace "<6" by "6"
x <- NationalSurveyData$arsenic
x <- replace(x, x == "< 0.5", "0.5")
x <- replace(x, x == "< 6", "6")
x <- as.numeric(levels(x))[x]
NationalSurveyData$arsenic <- x

#drop unused data vector
remove(x)

#Grouped data by their sub-division name, options: division/district/thana
# for division
#groupedData <- split(NationalSurveyData, NationalSurveyData$division)
# for district
# groupedData <- split(NationalSurveyData, NationalSurveyData$district)
# for thana
groupedData <- split(NationalSurveyData, NationalSurveyData$thana)

# Only the thana Chittagong Port has one observation 
#  groupedData$`Chittagong Port` <- NULL

# choose single division
# divisionSurveyData <- subset(NationalSurveyData, subset = division == "Dhaka")
# groupedData <- split(divisionSurveyData, divisionSurveyData$thana, drop = TRUE)


#test Vector
yVector <- groupedData$Lakshmipur$arsenic

#Function to convert an observation into csi score using sample yVector above
csiScore <-function(x,yVector){
# x is argument, i.e. oberservation to be convert
# yVector is a vector of ALL observations
M <- sum(yVector <= L)
n <- length(yVector)
ySort <- sort(yVector)
yTruncate <- ySort[-(1:M)]
LY <- max(0,1-L/x)
coeff <- 2/(n-M-1)
ySum <- sum(pmin(x,yTruncate)/(x + yTruncate))
val<- LY^(coeff*ySum)
return(val)
}

csiStep <- function(x,yVector){
  ySort <- sort(yVector)
  yMax <- max(yVector)
  csiScoreVector <- vapply( ySort, csiScore, yVector = ySort, FUN.VALUE = 0)
  n <- length(ySort) - 1
  M <- sum(yVector <= L)
  val <- 0 # to be discussed
  
  if (x > yMax) return(1)
  
  # first index greater or equal to existed oberservation
  index <- which.max(x <= ySort)
  if (index - 1 == 0) return (val) #WITHOUT THIS LINE MAY CAUSE BUG i.e. EMPTY VALUE
  val <- csiScoreVector[index-1]
  return(val)
}

csiSmooth<-function(y,yVector){
  #smoothing is done for y>L
  if (y <= L) return(0)

  ySort <-sort(yVector)
  n<-length(ySort)
  lam <- n/max(ySort)
    
  #N is the maximum number of points in the sum 
  N <- ceiling(lam * max(ySort))
  mu <- lam * y
  
  k<- 0:N
  kSequence <- L+k/lam
  FTerm <- vapply( kSequence, csiStep, yVector = yVector, FUN.VALUE = 0)
  val <- sum(FTerm*dpois(k,mu)) + 1 - ppois(N,mu)
  return(val)
  }

poissonSmoothDensity <-function(y,yVector){
  n<-length(yVector)
  lam <- n/max(yVector)
  
  #N is the maximum number of points in the sum 
  N <- ceiling(lam * max(yVector))
  mu <- lam * y
  
  Fn <-ecdf(yVector)
  k <- 0:N
  ks <- k/lam 
  ks1<-(k+1)/lam
  
  fTerm <- lam*sum((Fn(ks1)-Fn(ks))*dpois(k,mu))
  return(fTerm)
}

csiIntegration <- function(x, yVector){
  scalarIntegration <-function(x,yVector){
    csiSmooth(x,yVector)*poissonSmoothDensity(x,yVector)
  }
  vapply( x, scalarIntegration, yVector = yVector, FUN.VALUE = 0)
}

# Calculate smooth CSI base on Poisson Smoothing
scsi <- function(yVector){
n <- length(yVector)
if (n < 3) return(mean(pmax(0,1-L/yVector)))
val <- integrate(csiIntegration, lower = L, upper = Inf, yVector = yVector)$value
return(val)
}

#Function of calculating csi for a region
csi <- function(yVector){
  # yVector is a vector of ALL observations
  M <- sum(yVector <= L)
  n <- length(yVector)

  if (n < 3) return(mean(pmax(0,1-L/yVector)))

  ySort <- sort(yVector)
  yTruncate <- ySort[-(1:M)]
  LY <- pmax(0,1-L/yTruncate)
  
  coeff <- 2/(n-M-1)
  yExponent <- rep(NA,length(yTruncate))
  for(i in 1:length(yTruncate)){
    x <- yTruncate[i]
    yExponent[i] <- sum(pmin(x, yTruncate)/(x+yTruncate))  
  }
  yExponent <- coeff*yExponent
  csiScoreVector <- LY^yExponent
  val <- sum(csiScoreVector)/n
  return(val)
}

# Iterator through list of division/district/thana
csiDataFrame <- function(dataList){
  n <- length(dataList)
  regionName <- names(dataList)
  regionCsi <- rep(NA,n)
  for (i in 1:n)
  {
    yVector <- dataList[[i]][[arsenicIndex]]
    regionCsi[i] <- csi(yVector)
  }
  csiDF <- data.frame(regionName, regionCsi)
}

# csiDF <- csiDataFrame(groupedData)


scsiDataFrame <- function(dataList){
  n <- length(dataList)
  regionName <- names(dataList)
  regionCsi <- rep(NA,n)
  # try foreach for parallel computing
  for (i in 1:n)
  {
    yVector <- dataList[[i]][[arsenicIndex]]
    regionCsi[i] <- scsi(yVector)
  }
  scsiDF <- data.frame(regionName, regionCsi)
}

# Time consuming, work in principle, faster for sub-division
# scsiDF <- scsiDataFrame(groupedData)
# Note that the computation is time consuming, we may use other computer to do that
# And export the result so that we can use it
# write.csv(scsiDF, file = "E://University/2017.05.Summer/CUSRA/R/CUSRA Qi/scsiDF.csv")

# Load calculated Smooth CSI from exported file 
# scsiDF <- read.csv(file = "E://University/2017.05.Summer/CUSRA/R/CUSRA Qi/scsiDFDistrict.csv", header = TRUE)[,c(2,3)]
# scsiDF <- read.csv(file.choose(), header = TRUE)[,c(2,3)]
 
 #divisionCompare <- data.frame(regionName = csiDF$regionName, csi = csiDF$regionCsi, scsi = scsiDF$regionCsi)
 #districtCompare <- data.frame(regionName = csiDF$regionName, csi = csiDF$regionCsi, scsi = scsiDF$regionCsi)
 #thanaCompare <- data.frame(regionName = csiDF$regionName, csi = csiDF$regionCsi, scsi = scsiDF$regionCsi)

##################
#Output for Division
##################
#> ptm <- proc.time()
#> scsiDF <- scsiDataFrame(groupedData)
#> proc.time() - ptm
#  user  system elapsed 
#  91560.20    0.03 91560.23 

##################
#Output for District
##################
#> ptm <- proc.time()
#> scsiDF <- scsiDataFrame(groupedData)
#> proc.time() - ptm
# user  system elapsed 
# 3068.27    0.20 3068.93 

##################
#Output for Thana
##################
#> ptm <- proc.time()
#> scsiDF <- scsiDataFrame(groupedData)
#> proc.time() - ptm
#  user  system elapsed 
#  335.32   0.03 335.37 