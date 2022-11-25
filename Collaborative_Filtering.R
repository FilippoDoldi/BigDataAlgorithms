# SETUP AND LIBRARIES

rm(list = ls())
library(igraph)
library(Matrix)
library(tidyverse) 
Dir.Work <- "./" 
setwd(Dir.Work)

load("UserShows.RDATA")

# Adjacency matrix
R = Matrix(UserShows, sparse = TRUE) 

# P_ii is the number of objects linked with user i
P <- Diagonal(x = rowSums(R))

# Q_ii is the number of users linked with object i
Q <- Diagonal(x = colSums(R))

Past <- Diagonal(x = 1/sqrt(rowSums(R)))
Qast <- Diagonal(x = 1/sqrt(colSums(R)))

# User/User collaborative system
# We want to advise user 500 on the first 100 shows
vettoreRank563<-(Past%*%tcrossprod(R)%*%Past%*%R)[500,1:100]

# we order user in decreasing order.
suggerClaudiaUtenti <- sort(vettoreRank563, 
                            decreasing = TRUE, 
                            index.return = T)$ix 
# Result
cat(paste(as.vector(Shows$name[suggerClaudiaUtenti[1:5]]),collapse = "\n"))

# Object/object collaborative system
#Sistema collaborativo oggetto-oggetto
# We want to advise user 500 on the first 100 shows
vettoreRank563<-(R%*%Qast%*%crossprod(R)%*%Qast)[500,1:100]

# we order user in decreasing order.
suggerClaudiaOggetti <- sort(vettoreRank563, 
                            decreasing = TRUE, 
                            index.return = T)$ix 
# Result
cat(paste(as.vector(Shows$name[suggerClaudiaOggetti[1:5]]),collapse = "\n"))

#PageRang
# Creates the bipartite graph
# The vertices are in the same vector identified by a binary label
# The first 9985 are the user
g<-graph.incidence(R)

#Logical vector for teleportation on user 500
appoggio <- rep_len(0, length.out = length(V(g)$type))
appoggio[500]=1

#We take only the first 100 shows
vettoreRank563 <- page_rank(graph = g, 
                            personalized = appoggio)$vector[9986:10085]

#we order user in decreasing order.
suggerClaudiaPageRank <- sort(vettoreRank563, 
                             decreasing = TRUE, 
                             index.return = T)$ix 
#Result 
cat(paste(as.vector(Shows$name[suggerClaudiaPageRank[1:5]]),collapse = "\n"))

# Comparison for fraction of corrected suggestions
Performance <- bind_rows(
tibble(y = cumsum(Claudia[suggerClaudiaUtenti])/ seq_len(100)) %>% 
  mutate(x = row_number(), z = "user-user") ,
tibble(y = cumsum(Claudia[suggerClaudiaOggetti])/ seq_len(100)) %>% 
  mutate(x = row_number(), z = "item-item") ,
tibble(y = cumsum(Claudia[suggerClaudiaPageRank])/ seq_len(100)) %>%
  mutate(x = row_number(), z = "page-rank") )
# plot
ggplot(data = Performance %>% filter(x<51)) +
geom_line(aes(y = y, x=x, color = z)) + 
labs(title = "Performances comparison",
                          x = "Number of suggestion",
                          y = "Fraction of corrected suggestion") + 
theme_classic()
