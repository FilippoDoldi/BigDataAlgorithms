rm(list = ls())
library(tidyverse)
library(car)
Dir.Work <- "./" 
setwd(Dir.Work)

# Test set or final DataSet
debugMode <- FALSE # for debugging / change to false for final

if(debugMode){
  path_Data <- "dataInTiny" 
  real_counts <- "counts_tiny.txt.gz"
} else {
  path_Data <- "dataIn" # for debugging / production
  real_counts <- "counts.txt.gz"
}

# number of files (max) that will be analyzed
nFileMax <- 5000

allFiles <- dir(path = path_Data, pattern = "./*.gz",full.names = "T")
nFiles <- min(nFileMax,length(allFiles))

# parameters used in this code
K <- 5
J <- 1e4

# parameters for hash functions (you can use them by hand)
hash_params <- matrix(data = c(
  3,1561,
  17,277,
  38,394,
  61,13,
  78,246),byrow = T,nrow = 5,ncol = 2)

H <- matrix(data = 0L, nrow = K, ncol = J) # matrix H


# counter for number of total object seen (t in our notation)
nTotWords <- 0L
for(j1 in seq_len(nFiles)){
  cat("number ",j1," out of ",nFiles,": ")
  cat("processing file",allFiles[[j1]],"\n")
  
  # here we read the data (we then analyze it as a stream data)
  streamData <- read_delim(file = allFiles[[j1]],
                           col_names = F,
                           col_types = "i",delim = " ")
  
  # objected arrived in this last file
  objects_stream <- streamData$X1
  
  # counter of total objects read so far
  nTotWords <- nTotWords + length(objects_stream)
  
  # here we apply the k - hash functions 
  for (k in seq_len(K)){
    # reading hash parameters
    a <- hash_params[k,1]
    b <- hash_params[k,2]
    
    # computing the hash function h using the fact that 
    # we have the objects_stream vector
    hashStream <- ((a*objects_stream + b)%%123457)%%J
    
    # idea to speed up the code: group together
    #all the object mapped in the same place
    hashSummary <- tibble(ai = hashStream) %>% group_by(ai) %>% count()
    
    # updating H
    H[k,hashSummary$ai+1L] <- H[k,hashSummary$ai+1L] + hashSummary$n 
    # add 1L since modulus in [0,...n-1]

  }
  
}

# load the true counters
counts_words <- read_delim(real_counts,
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE,col_types = "ii") 

### extract object names 
different_objects <- counts_words$X1 # i = 1, ... , n
### extract frequencies in the whole stream
Fi <- counts_words$X2 # F[i]

# compute hatF for any object: we start with Inf  
hatFi <- rep(x = Inf,length.out = length(different_objects)) # \hat{F}[i]

# and then we take the minimum with the hash
for (k in seq_len(K)){
  # reading hash parameters
  a <- hash_params[k,1]
  b <- hash_params[k,2]
  
  # computing hatFi by taking the minimum
  hatFi <- pmin(hatFi, H[k,((a*different_objects + b)%%123457)%%J+1L] )
  
}

## final plot
ggplot(data = tibble(x = Fi, y = (hatFi-Fi)/Fi)) +
  geom_point(mapping = aes(x=x,y=y), colour = "blue", size = .125) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Object frequency") +
  ylab("Relative error in estimates") +
  annotation_logticks() +
  theme_bw() 

###################
# Now we need to do some regression 
# We start by fitting a line in all of the data linear model
  modellolm<-lm(log((hatFi-Fi)/Fi)~log(Fi))
  summary(modellolm)
  
  #scatterplot of residuals vs estimated values
  plot(modellolm$fitted.values, modellolm$residuals, 
    main='Residui vs valori stimati', lwd=2, xlab='Y stimati', ylab='Residui')
  abline(h=0, lwd=2)
  
  #residuals' QQ Plot
  qqPlot(modellolm$residuals,distribution = "norm",main='QQP dei residui')

# Now we fit a line only to the points where Fi>1000
  aux<-which(Fi>1000)
  Fi2<-Fi[aux]
  hatFi2<-hatFi[aux]
  
  modellolm<-lm(log((hatFi2-Fi2)/Fi2)~log(Fi2))
  summary(modellolm)
  
  #scatterplot of residuals vs estimated values
  plot(modellolm$fitted.values, modellolm$residuals, 
    main='Residui vs valori stimati', lwd=2, xlab='Y stimati', ylab='Residui')
  abline(h=0, lwd=2)
  
  #residuals' QQ Plot
  qqPlot(modellolm$residuals,distribution = "norm",main='QQP dei residui')

# and now where Fi>10000
  aux<-which(Fi>10000)
  Fi3<-Fi[aux]
  hatFi3<-hatFi[aux]
  
  modellolm<-lm(log((hatFi3-Fi3)/Fi3)~log(Fi3))
  summary(modellolm)
  
  #scatterplot of residuals vs estimated values
  plot(modellolm$fitted.values, modellolm$residuals, 
    main='Residui vs valori stimati', lwd=2, xlab='Y stimati', ylab='Residui')
  abline(h=0, lwd=2)
  
  #residuals' QQ Plot
  qqPlot(modellolm$residuals,distribution = "norm",main='QQP dei residui')
  
  