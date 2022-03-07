#Analyses Metabolic rates
#R version 3.5.1

library(ape)
library(geiger)
library(phytools)
library(castor)
library(paleotree)
library(ggtree)
library(phytools)
library(strap)


# Time scaling phylogeny ----

#Read the topology
phy <- read.nexus("Final_tree.nex")


#Read the stratigraphic ranges of species (FAD and LAD)
taxonTimes <- read.table("TaxonTimes.txt",header=TRUE,row.names=1)
#Read the minimum and maximum ages of each time interval
intervalTimes <- read.table("IntervalTimes.txt",header=TRUE,row.names=1)

timeList <-list(intervalTimes,taxonTimes)



#in addition to the tips of the tree, we added calibration ages to certain nodes: 


#- Synapsida vs diapsida (Edaphosaurus and mammals on one hand and all the rest) 
  #should be constrained between 312.3 and 330.4 Mya --> NODE 56
#- Archosauromorpha vs lepidosauromorpha (crocodiles+dinos and birds, and lizards) 
  #should be between 259.7 and 299.8--> NODE 76
#- Birds and crocodiles (birds+dinosa and crocs) should be between 
  #235 and 250.4 --> NODE 86
#- Crown birds (Crypturellus and the rest) should be between 66 and 86.5 --> NODE 101
#- Crown mammals (all mammals excluding Edaphosaurus) should be 
  #between 162.5 and 191.1 --> NODE 58
#- Marsupials plus placentals 124.6 and 138.4 --> NODE 60
#- Placentals should be between 95.3 and 113 --> NODE 63
#- pterosaurs (rhamphorinchoid + Pteranodon) and dinosaurs to 251.2-247.2 million years --> NODE 88
#- dinosaurs to 233.23-231.4 million years --> NODE 90


#Based on estimated age in timetree.org
# divergence of Mustela and Puma around 54 Ma --> NODE 68
# divergence of Cetacea and Bos around 56 Ma --> NODE 70
# divergence of Ochotona and Aotus 90 Ma --> NODE 71
# divergence of Ochotona and Rattus 82 Ma --> NODE 72
# divergence of Chironectes and Marmosa 30 Ma --> NODE 62
# divergence of Tachyglossus and Platypus 46 Ma --> NODE 59
# divergence of Iguana and Agama 157 Ma --> NODE 78
# divergence of Furcifer and Agama 101 Ma --> NODE 79
# divergence of Ara and Paser 82 Ma --> NODE 109
# divergence of Oceanodroma and Spheniscus 72.1 Ma --> NODE 108




pdf("nodelabels.pdf")
plot(phy, cex=0.5)
nodelabels(col= "red", cex=0.5, frame="none")
dev.off()

set.seed(15)

nodemins <- rep("NA", phy$Nnode)
nodemins[56-length(phy$tip.label)] <- sample(c(312.3:330.4),1)
nodemins[76-length(phy$tip.label)] <- sample(c(259.7:299.8),1)
nodemins[86-length(phy$tip.label)] <- sample(c(235:250.4),1)
nodemins[101-length(phy$tip.label)] <- sample(c(66.3:86.5),1)
nodemins[58-length(phy$tip.label)] <- sample(c(162.5:191.1),1)
nodemins[60-length(phy$tip.label)] <- sample(c(124.6:138.4),1)
nodemins[63-length(phy$tip.label)] <- sample(c(95.3:113),1)
nodemins[88-length(phy$tip.label)] <- sample(c(247.2:251.2),1)
nodemins[90-length(phy$tip.label)] <- sample(c(233.23:231.4),1)

nodemins[68-length(phy$tip.label)] <- 54
nodemins[70-length(phy$tip.label)] <- 56
nodemins[71-length(phy$tip.label)] <- 90
nodemins[72-length(phy$tip.label)] <- 82
nodemins[62-length(phy$tip.label)] <- 30
nodemins[59-length(phy$tip.label)] <- 46
nodemins[78-length(phy$tip.label)] <- 157
nodemins[79-length(phy$tip.label)] <- 101
nodemins[109-length(phy$tip.label)] <- 82
nodemins[108-length(phy$tip.label)] <- 72.1

nodemins <- as.numeric(nodemins)


#Time scaling of 1000 trees
bintrees <- bin_timePaleoPhy(phy, timeList=timeList, type = "mbl", 
                             vartime = 2, ntrees = 1000, nonstoch.bin = FALSE, randres = T, 
                             timeres = F,sites = NULL, point.occur = FALSE, add.term = T, 
                             inc.term.adj = T, dateTreatment = "firstLast", node.mins = nodemins, 
                             noisyDrop = TRUE, plot = F)


#Consensus edges

consensus_phy <- consensus.edges(bintrees,method="mean.edge", if.absent="ignore")
consensus_phy$root.time <- bintrees[[1]]$root.time

Consensus_phy_lad <- ladderize(consensus_phy, right = FALSE)


pdf("scaled_phylo_final.pdf")
geoscalePhylo(Consensus_phy_lad,units=c("Period", "Epoch", "Age"), cex.tip=0.2, width= 1, show.tip.label=T)
dev.off()


writeNexus(Consensus_phy_lad, "consensus_Paleotree_final.tre")





# +Ancestral reconstruction  ----

dat_MR <- read.table("TaxonTimes.txt",header=TRUE,row.names=1)

MetRate <- dat_MR$Calculated.MRs..mL.O2...1.h...1.g.
names(MetRate) <- dat_MR$Taxon


phydata <- treedata(Consensus_phy_lad, MetRate, sort=TRUE)

MetRate_phy <- phydata$data
colnames(MetRate_phy) <- "MetRate"



fit <-fastAnc(phydata$phy,MetRate_phy,vars=TRUE,CI=TRUE)
fit$ace

ramp_rate <- colorRampPalette(c("#02b2ce","#ffd004",  "#e52920"), bias=1.5) 

contmap_MR <- contMap(phydata$phy, MetRate_phy[,1], method = "user", 
                      anc.states = fit$ace)
contmap_MR$cols[] <- ramp_rate(1001)



pdf("contmap_phy.pdf")
plot(contmap_MR, lwd = 2, fsize=0.5, outline=F)
dev.off()



td <- data.frame(node = nodeid(phydata$phy, names(MetRate_phy[,1])),
                 trait = MetRate_phy[,1])
nd <- data.frame(node = as.numeric(names(fit$ace)),
                 trait = fit$ace)
d <- rbind(td, nd)


write.table(d, "ASR_MetRate.csv", sep=";", dec=".")







######

