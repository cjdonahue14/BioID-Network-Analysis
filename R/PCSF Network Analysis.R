library(PCSF)
library(readxl)
#Assemble Interactome
#HIPPIE was extracted from NDexBio using the code in "HIPPIE_pullcode.R". This interactome detailed an erroneous NFX1/NXF1 mix-up interactions (personal communication from BioGRID) which were manually removed from the edge list used to build the final interactome.
#Navigate to the directory containing the edge list for building the HIPPIE 2022 interactome
setwd("/path to directory with edge list for building interactome")
read_excel("HIPPIE2022_ConfR_Corr.xlsx")->HIPX
View(HIPX)
as.data.frame(HIPX)->HIPXDF
construct_interactome(HIPXDF)->HIPNetX
setwd("/path to directory containing extended data table files")
#Load BioID hits from combined VP vs. GFP and other VP datasets
read_excel("Extended data table 2.xlsx")-> BioIDUpAllandGFP
View(BioIDUpAllandGFP)
BioIDUpAllandGFP$`HIPPIE Gene Name`->Genes
BioIDUpAllandGFP$`Saint Score`->Prize
names(Prize)<-Genes
Prize
subnetx<-PCSF(HIPNetX,Prize,w=10,b=2,mu=0.01)
plot.PCSF(subnetx)
#Overlay virus interactions onto network
#Begin overlaying viral protein interactors onto resulting network solution
#Load viral protein interaction table
read_excel("Extended data table 3.xlsx")->BioIDNodeAtts
#Subset dataset so you maintain interactors and interaction scores of each VP with each host protein hit. 
BioIDNodeAtts[,1:3]->AttsOnly
View(AttsOnly)
#Spread dataset so that VPs that do not interact with a host protein are given a score of 0
library(tidyr)
Nodeattssep<- spread(AttsOnly,"VP","Rescore",fill = 0)
View(Nodeattssep)
#Subset so that node attributes are read as a vector
#Subset attributes based on viral protein and load each as a new dataset
cbind(Nodeattssep[,1],Nodeattssep[,2])->LScores
View(LScores)
cbind(Nodeattssep[,1],Nodeattssep[,3])->NPScores
View(NPScores)
cbind(Nodeattssep[,1],Nodeattssep[,4])->VP24Scores
View(VP24Scores)
cbind(Nodeattssep[,1],Nodeattssep[,5])->VP30Scores
View(VP30Scores)
cbind(Nodeattssep[,1],Nodeattssep[,6])->VP35Scores
View(VP35Scores)
cbind(Nodeattssep[,1],Nodeattssep[,7])->VP40Scores
View(VP40Scores)
#Convert individual viral protein interaction scores to numerical data, and map them onto the network as attributes
#Add L Scores as attributes
LScores[,1]<-sapply(strsplit(as.character(LScores[,1]),split = " "),"[[",1)
AttL<-match(V(subnetx)$name,LScores[,1])
V(subnetx)$LInterN<-as.numeric(LScores[,2][AttL])
#Add NP Scores as attributes
NPScores[,1]<-sapply(strsplit(as.character(NPScores[,1]),split = " "),"[[",1)
AttNP<-match(V(subnetx)$name,NPScores[,1])
V(subnetx)$NPInterN<-as.numeric(NPScores[,2][AttNP])
#Add VP35 Scores as attributes
VP35Scores[,1]<-sapply(strsplit(as.character(VP35Scores[,1]),split = " "),"[[",1)
AttVP35<-match(V(subnetx)$name,VP35Scores[,1])
V(subnetx)$VP35InterN<-as.numeric(VP35Scores[,2][AttVP35])
#Add VP40 Scores as attributes
VP40Scores[,1]<-sapply(strsplit(as.character(VP40Scores[,1]),split = " "),"[[",1)
AttVP40<-match(V(subnetx)$name,VP40Scores[,1])
V(subnetx)$VP40InterN<-as.numeric(VP40Scores[,2][AttVP40])
#Add VP24 Scores as attributes
VP24Scores[,1]<-sapply(strsplit(as.character(VP24Scores[,1]),split = " "),"[[",1)
AttVP24<-match(V(subnetx)$name,VP24Scores[,1])
V(subnetx)$VP24InterN<-as.numeric(VP24Scores[,2][AttVP24])
#Add VP30 Scores as attributes  
VP30Scores[,1]<-sapply(strsplit(as.character(VP30Scores[,1]),split = " "),"[[",1)
AttVP30<-match(V(subnetx)$name,VP30Scores[,1])
V(subnetx)$VP30InterN<-as.numeric(VP30Scores[,2][AttVP30])
get.vertex.attribute(subnetx)
#Load resulting network into Cytoscape for aesthetic/ network graphing
library(RCy3)
createNetworkFromIgraph(subnetx,"BioID_PCSF_wHIPPIE22Corr")
#Perform gene set enrichment for network analysis. Note that the below code is adapted from the "enrichment_analysis()" function from PCSF, which is now deprecated 
cluster_edge_betweenness(subnetx)->Clx
data.frame(Membership=Clx$membership, Genes=Clx$names)
x <- data.frame(Membership=Clx$membership, Genes=Clx$names) %>% dplyr::arrange(Membership)
unique(x$Membership)
unique(x$Membership) %>% lapply(function(i) {x[x$Membership==i, "Genes"]})
y <- unique(x$Membership) %>% lapply(function(i) {x[x$Membership==i, "Genes"]})
y[[1]]
library(enrichR)
setEnrichrSite("Enrichr")
Enrdb<-(c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023","KEGG_2021_Human","Reactome_2022"))
encl1<-enrichr(y[[1]],Enrdb)
encl2<-enrichr(y[[2]],Enrdb)
encl3<-enrichr(y[[3]],Enrdb)
encl4<-enrichr(y[[4]],Enrdb)
encl5<-enrichr(y[[5]],Enrdb)
encl6<-enrichr(y[[6]],Enrdb)
encl7<-enrichr(y[[7]],Enrdb)
encl8<-enrichr(y[[8]],Enrdb)
encl9<-enrichr(y[[9]],Enrdb)
encl10<-enrichr(y[[10]],Enrdb)
encl11<-enrichr(y[[11]],Enrdb)
encl12<-enrichr(y[[12]],Enrdb)
encl13<-enrichr(y[[13]],Enrdb)
encl14<-enrichr(y[[14]],Enrdb)
encl15<-enrichr(y[[15]],Enrdb)
encl16<-enrichr(y[[16]],Enrdb)
encl17<-enrichr(y[[17]],Enrdb)
encl18<-enrichr(y[[18]],Enrdb)
encl19<-enrichr(y[[19]],Enrdb)
encl20<-enrichr(y[[20]],Enrdb)
encl21<-enrichr(y[[21]],Enrdb)
Cluster1<-rbind(encl1$GO_Biological_Process_2023,encl1$GO_Cellular_Component_2023,encl1$GO_Molecular_Function_2023,encl1$KEGG_2021_Human,encl1$Reactome_2022)
Cluster2<-rbind(encl2$GO_Biological_Process_2023,encl2$GO_Cellular_Component_2023,encl2$GO_Molecular_Function_2023,encl2$KEGG_2021_Human,encl2$Reactome_2022)
Cluster3<-rbind(encl3$GO_Biological_Process_2023,encl3$GO_Cellular_Component_2023,encl3$GO_Molecular_Function_2023,encl3$KEGG_2021_Human,encl3$Reactome_2022)
Cluster4<-rbind(encl4$GO_Biological_Process_2023,encl4$GO_Cellular_Component_2023,encl4$GO_Molecular_Function_2023,encl4$KEGG_2021_Human,encl4$Reactome_2022)
Cluster5<-rbind(encl5$GO_Biological_Process_2023,encl5$GO_Cellular_Component_2023,encl5$GO_Molecular_Function_2023,encl5$KEGG_2021_Human,encl5$Reactome_2022)
Cluster6<-rbind(encl6$GO_Biological_Process_2023,encl6$GO_Cellular_Component_2023,encl6$GO_Molecular_Function_2023,encl6$KEGG_2021_Human,encl6$Reactome_2022)
Cluster7<-rbind(encl7$GO_Biological_Process_2023,encl7$GO_Cellular_Component_2023,encl7$GO_Molecular_Function_2023,encl7$KEGG_2021_Human,encl7$Reactome_2022)
Cluster8<-rbind(encl8$GO_Biological_Process_2023,encl8$GO_Cellular_Component_2023,encl8$GO_Molecular_Function_2023,encl8$KEGG_2021_Human,encl8$Reactome_2022)
Cluster9<-rbind(encl9$GO_Biological_Process_2023,encl9$GO_Cellular_Component_2023,encl9$GO_Molecular_Function_2023,encl9$KEGG_2021_Human,encl9$Reactome_2022)
Cluster10<-rbind(encl10$GO_Biological_Process_2023,encl10$GO_Cellular_Component_2023,encl10$GO_Molecular_Function_2023,encl10$KEGG_2021_Human,encl10$Reactome_2022)
Cluster11<-rbind(encl11$GO_Biological_Process_2023,encl11$GO_Cellular_Component_2023,encl11$GO_Molecular_Function_2023,encl11$KEGG_2021_Human,encl11$Reactome_2022)
Cluster12<-rbind(encl12$GO_Biological_Process_2023,encl12$GO_Cellular_Component_2023,encl12$GO_Molecular_Function_2023,encl12$KEGG_2021_Human,encl12$Reactome_2022)
Cluster13<-rbind(encl13$GO_Biological_Process_2023,encl13$GO_Cellular_Component_2023,encl13$GO_Molecular_Function_2023,encl13$KEGG_2021_Human,encl13$Reactome_2022)
Cluster14<-rbind(encl14$GO_Biological_Process_2023,encl14$GO_Cellular_Component_2023,encl14$GO_Molecular_Function_2023,encl14$KEGG_2021_Human,encl14$Reactome_2022)
Cluster15<-rbind(encl15$GO_Biological_Process_2023,encl15$GO_Cellular_Component_2023,encl15$GO_Molecular_Function_2023,encl15$KEGG_2021_Human,encl15$Reactome_2022)
Cluster16<-rbind(encl16$GO_Biological_Process_2023,encl16$GO_Cellular_Component_2023,encl16$GO_Molecular_Function_2023,encl16$KEGG_2021_Human,encl16$Reactome_2022)
Cluster17<-rbind(encl17$GO_Biological_Process_2023,encl17$GO_Cellular_Component_2023,encl17$GO_Molecular_Function_2023,encl17$KEGG_2021_Human,encl17$Reactome_2022)
Cluster18<-rbind(encl18$GO_Biological_Process_2023,encl18$GO_Cellular_Component_2023,encl18$GO_Molecular_Function_2023,encl18$KEGG_2021_Human,encl18$Reactome_2022)
Cluster19<-rbind(encl19$GO_Biological_Process_2023,encl19$GO_Cellular_Component_2023,encl19$GO_Molecular_Function_2023,encl19$KEGG_2021_Human,encl19$Reactome_2022)
Cluster20<-rbind(encl20$GO_Biological_Process_2023,encl20$GO_Cellular_Component_2023,encl20$GO_Molecular_Function_2023,encl20$KEGG_2021_Human,encl20$Reactome_2022)
Cluster21<-rbind(encl21$GO_Biological_Process_2023,encl21$GO_Cellular_Component_2023,encl21$GO_Molecular_Function_2023,encl21$KEGG_2021_Human,encl21$Reactome_2022)
#Merge into excel doc
library(xlsx)
write.xlsx(Cluster1,sheetName = c("Cluster1"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE)
write.xlsx(Cluster2,sheetName = c("Cluster2"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster3,sheetName = c("Cluster3"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster4,sheetName = c("Cluster4"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster5,sheetName = c("Cluster5"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster6,sheetName = c("Cluster6"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster7,sheetName = c("Cluster7"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster8,sheetName = c("Cluster8"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster9,sheetName = c("Cluster9"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster10,sheetName = c("Cluster10"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster11,sheetName = c("Cluster11"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster12,sheetName = c("Cluster12"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster13,sheetName = c("Cluster13"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster14,sheetName = c("Cluster14"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster15,sheetName = c("Cluster15"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster16,sheetName = c("Cluster16"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster17,sheetName = c("Cluster17"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster18,sheetName = c("Cluster18"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster19,sheetName = c("Cluster19"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster20,sheetName = c("Cluster20"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)
write.xlsx(Cluster21,sheetName = c("Cluster21"),file = "BioID HIPPIE22 EnrichmentI.xlsx",row.names = TRUE,append = TRUE)