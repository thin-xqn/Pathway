# For the pathway response gene sets, we updated all the names of genes to official symbol.
# Then we deleted the non-human genes and got a csv file called PRGS1.
# filter out Wnt activating experiments (unreliable)
# filter out unknown effect
library(dplyr)
signZvalue_vec = function(zscore, effect) {
  return(sapply(1:length(zscore), function(ix)  signZvalue(zscore[ix], effect[ix])))
}
signZvalue = function(zscore, effect) {
  if(effect %in% c("activating", "activation")) {
    return(zscore)
  } else if(effect %in% c("inhibiting", "inhibition", "inhibting")) {
    return(-zscore)
  }
  return(NA)
}
PRGS1 <- read.csv(file = "F://PRGS1.csv")    
to_remove <- c("MAPK", "NFkB", "PI3K", "RAR", "Trail")
PRGS2 <- PRGS1[!PRGS1$p_id %in% to_remove,]
PRGS3 = PRGS2 %>% filter(!(p_id=="Wnt" & effect=="activation"))
PRGS4 = PRGS3%>%
  group_by(p_id,e_id) %>%
  filter(expression > median(expression)) %>%
  ungroup %>%
  filter(effect %in% c("activation", "inhibition")) %>% 
  group_by(p_id, gene, g_id, e_id, effect) %>%
  summarise(zvalue=zvalue[which.max(expression)]) %>%
  ungroup %>%
  mutate(zvalue_signed = signZvalue_vec(zvalue, effect)) %>%
  group_by(p_id, e_id) %>%
  ungroup
write.csv(PRGS4, file="preproccess.csv", row.names = FALSE) 