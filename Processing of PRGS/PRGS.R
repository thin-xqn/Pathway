# Finally, we combined the 15 groups of pathway responsive genes into a pathway response gene set(PRGS).
library(dplyr)
library(magrittr)
path <- "F:/Pathway_Gaussian"
file_list <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
for (file in file_list) {
  file_name <- sub("\\.csv$", "", basename(file))
  assign(file_name, read.csv(file))
}
object_names <- ls()
data_frame_names <- grep("_Gaussian$", object_names, value = TRUE)
print(data_frame_names)
My <- data_frame_names
B <- data.frame()      
for(i in seq_along(My)) {
  data <- get(My[i])  
  B <- rbind(B,data)
}
pathway_response <- B  
setwd("F://Pathway_response")
write.csv(pathway_response, file="pathway_response.csv", row.names = FALSE) 
