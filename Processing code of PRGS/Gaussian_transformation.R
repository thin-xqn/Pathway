# We got 15 groups of pathway response genes after the previous processing step(Ppdd).
# There were different number of datasets for each of the pathways.
# For each pathway, gaussian transformation was applied to each dataset.
# Finally, a csv file named "pathway_Gaussian" was obtained.
# As an example, this programme used the Estrogen response gene set for gaussian transformation.
library(dplyr)
library(magrittr)
library(GSVA)
library(GSVAdata)
Estrogen <- read.csv(file = "F://Estrogen.csv") 
E1 <- Estrogen[,c(2,4,7)]
E1 %>% group_by(e_id)%>%mutate(number=n())%T>%str->sum 
sum <- sum %>% distinct(e_id,number,.keep_all = TRUE)  
e_id <- c(sum$e_id)
My=c("GSE12261","GSE20081","GSE22012","GSE23610","GSE25315","GSE26459a","GSE26459b","GSE28006","GSE30931","GSE35287a","GSE35287b","GSE35287c","GSE35287d","GSE3834d","GSE3834e","GSE3834f","GSE28645","GSE1045a","GSE1045b","GSE13458","GSE22593a","GSE22593b","GSE22593c","GSE3834a","GSE3834b","GSE3834c","GSE5258Est","GSE11467","GSE23500","GSE18592","GSE23445","GSE25314a","GSE3834g","GSE3834h","GSE11567a","GSE11791","GSE15717","GSE24592","GSE26834Est","GSE9936a","GSE9936b","GSE9936c","GSE3013b","GSE3013a","GSE3013i","GSE3013j","GSE3013f")
for(i in 1:47) {
  b <- My[i]
  c <-E1[which(E1$e_id ==b),]
  assign(paste0("GSE",i), c)
  
}
My=c("GSE1","GSE2","GSE3","GSE4","GSE5","GSE6","GSE7","GSE8","GSE9","GSE10",
     "GSE11","GSE12","GSE13","GSE14","GSE15","GSE16","GSE17","GSE18","GSE19","GSE20",
     "GSE21","GSE22","GSE23","GSE24","GSE25","GSE26","GSE27","GSE28","GSE29","GSE30",
     "GSE31","GSE32","GSE33","GSE34","GSE35","GSE36","GSE37","GSE38","GSE39","GSE40",
     "GSE41","GSE42","GSE43","GSE44","GSE45","GSE46","GSE47")
for(i in seq_along(My)) {
  b <- get(My[i])
  df <- t(b)
  colnames(df) <- df[1,]  
  dataset_name <- df[2,1]    
  gene_name <- df[1,]     
  df <- df[-c(1:2),]      
  df <- t(df)
  rownames(df) <- dataset_name  
  # Gaussian input
  expr <- df
  sample.idxs <- gene_name
  ## GSVA Gaussian
  n.test.samples <- ncol(expr)
  n.genes <- nrow(expr)
  n.density.samples <- length(sample.idxs)
  rnaseq= FALSE
  gene.density <- NA
  A = .Call("matrix_density_R",
            as.double(t(expr[ ,sample.idxs, drop=FALSE])),
            as.double(t(expr)),
            n.density.samples,
            n.test.samples,
            n.genes,
            as.integer(rnaseq))
  gene.density <- t(matrix(A, n.test.samples, n.genes))
  result <- gene.density
  rownames(result) <- dataset_name   
  colnames(result) <- gene_name     
  assign(paste0("Gaussian",i), result)
  
}
My=c("Gaussian1","Gaussian2","Gaussian3","Gaussian4","Gaussian5","Gaussian6",
     "Gaussian7","Gaussian8","Gaussian9","Gaussian10","Gaussian11","Gaussian12",
     "Gaussian13","Gaussian14","Gaussian15","Gaussian16","Gaussian17","Gaussian18",
     "Gaussian19","Gaussian20","Gaussian21","Gaussian22","Gaussian23","Gaussian24",
     "Gaussian25","Gaussian26","Gaussian27","Gaussian28","Gaussian29","Gaussian30",
     "Gaussian31","Gaussian32","Gaussian33","Gaussian34","Gaussian35","Gaussian36",
     "Gaussian37","Gaussian38","Gaussian39","Gaussian40","Gaussian41","Gaussian42",
     "Gaussian43","Gaussian44","Gaussian45","Gaussian46","Gaussian47")
for(i in seq_along(My)) {
  b <- get(My[i])
  data <- data.frame(t(b))
  assign(paste0("Hist",i), data)
  
}
paste0("Hist", 1:47)
My=c("Hist1","Hist2","Hist3","Hist4","Hist5","Hist6","Hist7","Hist8","Hist9","Hist10",
     "Hist11","Hist12","Hist13","Hist14","Hist15","Hist16","Hist17","Hist18","Hist19",
     "Hist20","Hist21","Hist22","Hist23","Hist24","Hist25","Hist26","Hist27","Hist28",
     "Hist29","Hist30","Hist31","Hist32","Hist33","Hist34","Hist35","Hist36","Hist37",
     "Hist38","Hist39","Hist40","Hist41","Hist42","Hist43","Hist44","Hist45","Hist46",
     "Hist47")
B <- data.frame()       
for(i in seq_along(My)) {
  Hist <- get(My[i])  
  A <- cbind(rownames(Hist),Hist)     
  A$e_id <- c(colnames(Hist))             
  A$Pathway <- c("Estrogen") 
  names(A)[1]="gene"                 
  names(A)[2]="Gaussian"  
  B <- rbind(B,A)
}
B <- B[,c(4,3,1,2)]    
names(B)[1]="p_id"               
Estrogen1 <- B %>%group_by(p_id,gene)%>%mutate(n=n())%T>%str 
Estrogen2 <-Estrogen1 [-which(Estrogen1$n==1),]  
Estrogen3 <- Estrogen2 %>% group_by(p_id, gene)
Estrogen4 <- Estrogen3 %>% mutate(z_mean=mean(Gaussian))
Estrogen5 <- Estrogen4[,-4]
library(dplyr)
Estrogen6 <- Estrogen5 %>% distinct(p_id,gene,n,z_mean,.keep_all = TRUE) 
Estrogen7 <- Estrogen6[,-c(4)]
setwd("F://Pathway_Gaussian")
write.csv(Estrogen7, file="Estrogen_Gaussian.csv", row.names = FALSE) 