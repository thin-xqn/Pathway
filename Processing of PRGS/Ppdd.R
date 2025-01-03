# Ppdd:Per pathway data disposal
# Extracting responsive genes from each pathway separately for further processing.
library(dplyr)
preproccess <- read.csv(file = "F://preproccess.csv")                                              
My=c("Estrogen","H2O2","Hippo","Hypoxia","IL-1","Insulin",
     "JAK-STAT","MAPK+PI3K","Notch","p53","PPAR","TGFb","TNFa","VEGF","Wnt")
for(i in 1:15) {
  b <- My[i]
  c <-preproccess[which(preproccess$p_id==b),]
  assign(My[i], c)
  
}
My <- c("Estrogen", "H2O2", "Hippo", "Hypoxia", "IL-1", "Insulin",
        "JAK-STAT", "MAPK+PI3K", "Notch", "p53", "PPAR", "TGFb",
        "TNFa", "VEGF", "Wnt")
for (i in seq_along(My)) {
  data_frame_name <- My[i]
  data_frame <- get(data_frame_name)
  file_name <- paste(data_frame_name, ".csv", sep="")
  write.csv(data_frame, file = file_name, row.names = FALSE)
}