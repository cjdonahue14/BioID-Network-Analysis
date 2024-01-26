library(ndexr)
ndexcon = ndex_connect() 
network_name <- ndex_find_networks(ndexcon, "HIPPIE") 
View(network_name)
networkId <- network_name[1,"externalId"]
View(networkId)
ndex_network_get_aspect(ndexcon,networkId,"edgeAttributes")->Edgecosts
View(Edgecosts)
library(tidyr)  
pivot_wider(Edgecosts,"propertyOf")->EdgecostW
View(EdgecostW)
library(tidyverse)
EdgecostW%>%mutate(Protein1=str_extract(name,"^[^\\(]+"),Protein2=str_extract(name,"[^\\)]+$"))->EdgecostWS4
View(EdgecostWS4)
cbind(EdgecostWS4$Protein1,EdgecostWS4$Protein2,EdgecostWS4$`Confidence Value`,1)->HIPPIE22
as.data.frame(HIPPIE22)->HIPPIE22DF
transform.data.frame(HIPPIE22DF,V3=as.numeric(V3),V4=as.numeric(V4))->HIPPIE22DFN
View(HIPPIE22DFN)
cbind(HIPPIE22DFN,HIPPIE22DFN[,4]/HIPPIE22DFN[,3])->HIPPIE22Cost
cbind(HIPPIE22Cost$V1,HIPPIE22Cost$V2,HIPPIE22Cost$`HIPPIE22DFN[, 4]/HIPPIE22DFN[, 3]`)->HIP22Final
colnames(HIP22Final)<-c("From","To","Cost")
as.data.frame(HIP22Final)->HIP22FinalDF
View(HIP22FinalDF)
cbind(HIP22FinalDF$From,HIP22FinalDF$To,HIP22FinalDF$Cost)->HIP22Char
View(HIP22Char)
write.csv(HIP22Char,"HIPPIE2022_Conf.csv")