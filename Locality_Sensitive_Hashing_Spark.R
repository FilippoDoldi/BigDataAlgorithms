#LIBRARIES AND DATA UPLOAD

rm(list = ls())
library(sparklyr)
library(tidyverse) 
Dir.Work <- ""
setwd(Dir.Work)
sc <- spark_connect(master = "local")


typeCols <- paste0("V",1:400,"= 'double'",collapse = ",")
typeCols <- paste("list(",typeCols,")")

# dataBase Images
datasetImages <-  spark_read_csv(sc, name = "datasetImages",
                                 path = "lsh.csv.gz",
                                 columns = eval(parse(text = typeCols)),
                                 header = FALSE, memory = F)
datasetImages %>% head()

# dataSet 20 Images z_i -> TargetImages
# we want to find the six nearest neighbors for each of them

# upload in R
ImagesTestR <- read_csv(file = "lsh2.csv.gz",
                        col_names = FALSE)

names(ImagesTestR) <- paste0("V",1:400)

# upload in Spark
tbfImages <-  copy_to(dest = sc,
                      df = ImagesTestR,
                      name = "tbfImages",
                      memory = F)
tbfImages %>% head()

##############################
#START
##############################
# We transform the datasetImages grouping vectors as one element in the
# dataset and adding an ID to the column for identification

datasetImages_Assembled <- datasetImages %>% 
  ft_vector_assembler(input_cols = paste0("V",1:400),output_col = "vect") %>%
  select(vect) %>% 
  sdf_with_unique_id %>% 
  mutate(id = int(id)) %>%
  sdf_register(name = "datasetImages_Assembled")
# Bring in memory
tbl_cache(sc, "datasetImages_Assembled")
datasetImages_Assembled %>% head()

# We repeat the same procedure for the Target images
datasettbfImages_Assembled <- tbfImages %>% 
  ft_vector_assembler(input_cols = paste0("V",1:400),output_col = "vect") %>%
  select(vect) %>% 
  sdf_with_unique_id %>% 
  mutate(id = int(id)) %>%
  sdf_register(name = "datasettbfImages_Assembled")

tbl_cache(sc, "datasettbfImages_Assembled")
datasettbfImages_Assembled %>% head()

##################
# Building the LSH model
model_LSH <- sc %>%
  ft_bucketed_random_projection_lsh(input_col = "vect",
                                    output_col = "buckets",
                                    bucket_length = 100,
                                    num_hash_tables = 24) %>%
                                    ml_fit(datasetImages_Assembled)

resultID<-matrix(NA,20,7) # ID matrix
distances<-matrix(NA,20,7) # distances matrix
# For each of the 20 images we compute the 7 NN with the function 
# ml_approx_nearest_neighbors. Then we save the results on R.

for(i in 1:20)
{
  NNeigh <- ml_approx_nearest_neighbors(model = model_LSH, 
                            dataset = datasetImages_Assembled, 
                            key = ImagesTestR[i,]%>% unlist() %>% as.vector(), 
                            num_nearest_neighbors = 7 ,
                            dist_col = "distCol")
  
  # Ordering by increasing distance and drop off uninteresting columns
  SevenNNeigh<-NNeigh%>%
    arrange(distCol)%>%
    select(-buckets, -vect)%>%
    collect()
  
  resultID[i,]<-SevenNNeigh$id%>%unlist()%>%as.vector()
  distances[i,]<-SevenNNeigh$distCol%>%unlist()%>%as.vector()
}
# The IDs start from 0, from R the vectors start from 1; we fix this
resultID<-resultID+1

###################################
# RESULTS
###################################

xiast<-matrix(NA,20,400) # 20 Best images
mean_distances<-rep(0,20) # mean distances vector
for(i in 1:20)
{
  z <- read_csv(file = "lsh2.csv.gz",
                col_names = FALSE, show_col_types = FALSE) %>%
    slice(i) %>%
    as.matrix()
  x <- read_csv(file = "lsh.csv.gz",
                col_names = FALSE, show_col_types = FALSE) %>%
    slice(resultID[i,2:7]) %>%
    as.matrix()
  
  mean_distances[i]<-mean(distances[i,])
  xiast[i,] <- x[which(distances[i,] == min(distances[i,])),]
}


z <- read_csv(file = "lsh2.csv.gz",
              col_names = FALSE, show_col_types = FALSE) %>%
  slice(1) %>%
  as.matrix()
image(matrix(z,ncol = 20,byrow = T), axes=FALSE)
image(matrix(xiast[1,],ncol = 20,byrow = T), axes=FALSE)

spark_disconnect(sc)