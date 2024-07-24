library(RCy3)
library(igraph)
library(readxl)
library(tidyr)
#Downloaded Metascape cluster Cytoscape file, converted to edge list, and uplaoded edgelist into R
setwd("C:/Users/Comma/Desktop/LaCount_Collabs/Decapping IF/Publication/Review/Post Review")
#Upload edgelist 
read_excel("MetascapeEdges_UpdatedPrey.xlsx", sheet=2)->NewEdges
View(NewEdges)
as.matrix(NewEdges)->NewEdgesM
View(NewEdgesM)
#Create igraph network from edge list
graph_from_edgelist(NewEdgesM)->Metascape
plot(Metascape)
#Begin attributing network
setwd("C:/Users/Comma/Desktop/LaCount_Collabs/Decapping IF/Publication/Review")
read_excel("Extended data table 1 updated 07-05-24 for publication_PCSFEdit.xlsx", sheet = 6)->BioIDNodeAtts
#Subset dataset so you maintain interactors and interaction scores of each VP with each host protein hit. 
BioIDNodeAtts[,1:3]->AttsOnly
View(AttsOnly)
#Spread dataset so that VPs that do not interact with a host protein are given a score of 0
library(tidyr)
Nodeattssep<- spread(AttsOnly,"Bait","Rescore",fill = 0)
View(Nodeattssep)
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
LScores[,1]<-sapply(strsplit(as.character(LScores[,1]),split = " "),"[[",1)
AttL<-match(V(Metascape)$name,LScores[,1])
V(Metascape)$LInterN<-as.numeric(LScores[,2][AttL])
#Add NP Scores as attributes
NPScores[,1]<-sapply(strsplit(as.character(NPScores[,1]),split = " "),"[[",1)
AttNP<-match(V(Metascape)$name,NPScores[,1])
V(Metascape)$NPInterN<-as.numeric(NPScores[,2][AttNP])
#Add VP35 Scores as attributes
VP35Scores[,1]<-sapply(strsplit(as.character(VP35Scores[,1]),split = " "),"[[",1)
AttVP35<-match(V(Metascape)$name,VP35Scores[,1])
V(Metascape)$VP35InterN<-as.numeric(VP35Scores[,2][AttVP35])
#Add VP40 Scores as attributes
VP40Scores[,1]<-sapply(strsplit(as.character(VP40Scores[,1]),split = " "),"[[",1)
AttVP40<-match(V(Metascape)$name,VP40Scores[,1])
V(Metascape)$VP40InterN<-as.numeric(VP40Scores[,2][AttVP40])
#Add VP24 Scores as attributes
VP24Scores[,1]<-sapply(strsplit(as.character(VP24Scores[,1]),split = " "),"[[",1)
AttVP24<-match(V(Metascape)$name,VP24Scores[,1])
V(Metascape)$VP24InterN<-as.numeric(VP24Scores[,2][AttVP24])
#Add VP30 Scores as attributes  
VP30Scores[,1]<-sapply(strsplit(as.character(VP30Scores[,1]),split = " "),"[[",1)
AttVP30<-match(V(Metascape)$name,VP30Scores[,1])
V(Metascape)$VP30InterN<-as.numeric(VP30Scores[,2][AttVP30])
get.vertex.attribute(Metascape)
createNetworkFromIgraph(Metascape,"Annotated Metascape")
