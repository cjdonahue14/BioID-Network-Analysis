#Create hub and spoke network models
library(RCy3)
library(igraph)
library(readxl)
setwd("/path to Extended data folders")
read_excel("Extended data table 2.xlsx")-> BioIDUpAllandGFP
BioIDUpAllandGFP[,1:2]->EdgesBioID2
as.data.frame(EdgesBioID2)->BioID2DF
as.matrix(BioID2DF)->BioIDMx
graph_from_edgelist(BioIDMx)->CytoNet1
plot(CytoNet1)
createNetworkFromIgraph(CytoNet1)
