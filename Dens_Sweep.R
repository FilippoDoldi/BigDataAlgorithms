#Setup and libraries

rm(list = ls())
library(igraph)
library(Matrix)
library(tidyverse) # contains ggplot2, dplyr, readr, etc...
library(sparklyr)
Dir.Work <- "./" # Your path
setwd(Dir.Work)
source("./dense_subGraph.R")

sc <- spark_connect(master = "local", version = "2.1")

DataSet <-  spark_read_parquet(sc = sc, 
                    memory = T, # in memory
                    overwrite = T,
                    name = "friends",
      path = paste0(Dir.Work,"/dataHW1") # use "/dataHW1" in final execution 
)

##### set of all vertices: already undirected (no need both V1 and V2)
vertices_tbl <- DataSet %>%
  transmute(id = V1) %>% 
  sdf_drop_duplicates()

##### already undirected V1->V2 and V2->V1 present
edges_tbl <- DataSet %>% 
  transmute(src = V1, dst = V2) %>% 
  sdf_drop_duplicates()

edgR <- edges_tbl %>% collect() # edges in R
vertR <- vertices_tbl %>% collect() # vertices in R
  
# close Spark
spark_disconnect(sc = sc)
rm(edges_tbl,DataSet,sc,vertices_tbl) 

#graph user/friends
g <- graph_from_data_frame(d = edgR, 
                           directed=FALSE,
                           vertices = vertR)

# Finds a dense subgraph (non-directed)
Subgraph <- findDenseSub(vertices = vertR ,edges = edgR) 
gorder(Subgraph) #number of nodes in S
densGr(Subgraph) #density of S
#Comparison
densGr(g) #density of the initial graph g

######################################################
# Sweep on user with id = 1
Components <- components(g)
Components$no # number of connected components
Components$membership[which(vertR$id==1)] # He is in the first component

#Build the graph of the first component
subGr <- induced_subgraph(g,which(Components$membership==1))
Components2 <- components(subGr)
Components2$no

# distribution of restarting: deterministic from LIN
#Position in the subgraph of the node with id = 1
pos<-which(names(Components2$membership)==1)

#Create the vector for Pagerank with restarting from id = 1
S = rep(x = 0,length.out = gorder(subGr))
S[pos] <- 1

PageRank <- page_rank(graph = subGr, damping = 0.85,
                      personalized = S)#restart in pos

# reorder in decreasing score
orderPageRank <- sort(x = PageRank$vector,
                    decreasing = TRUE, 
                    index.return = TRUE)$ix #select the index (not the value)

graphPermuted <- permute(subGr,invPerm(orderPageRank))

# get adiacency matrix. First row <-> LIN, etc...
A1 <- as_adjacency_matrix(graphPermuted,  type = c("both"))

D <- rowSums(A1) # degree 

# check !!!
which((D - degree(graphPermuted)) != 0)

# Conductance computation
# Vol(Ai+1) = Vol(Ai) + di => Vol(Ai+1) - Vol(Ai) = di 
VolAi <- cumsum(D)
# Cut(Ai+1) = Cut(Ai) + di+1 - 2 #{edg from di+1 to Ai}
# Cut(Ai+1) - Cut(Ai) =  di+1 - 2 #{edg from di+1 to Ai}
CutAi <- cumsum(D) - 2*cumsum(rowSums(tril(A1)))

#Plot
ggplot(data = tibble(x = VolAi,y = CutAi) %>%
         mutate(n = row_number()),aes(x=n,y = y/x)) +
  geom_line() +
  labs(x = "Node rank i in decreasing PPR score",
       y = "Conductance") +
  #  ylim(c(0,.5)) +
  theme_classic()




