# Code to plot phylogenetic tree of Eo and Ac
# load libraries
library(ape)
library(phytools)
library(ggtree)
library(ggplot2)

# clear environment
rm(list = ls())
# close plots
graphics.off()

load("~/obis_tree_F1pt8_log10Ac_3emin2pmax.RData")
#load("~/obis_tree_F1pt8_Eo_3emin2pmax.RData")

# label phyla nodes on the tree?
makelab=0

# plot
#dev.new(width=4, height=20, unit="in")
if (trait == "E") {
  p2 =
    ggtree(fish.tree, ,ladderize = FALSE, size=.6) + 
    geom_tree(aes(color=trait), continuous=T, size=.25) + 
    #ggtree(fish.tree, layout='circular',ladderize = FALSE, size=.5) + 
    #geom_tree(aes(color=trait), continuous=T, size=.3) + 
    #scale_color_gradientn(colours=c("blue","yellow","red")) +
    scale_color_gradientn(colours=inlmisc::GetColors(length(seq(0,20,5)),scheme='turbo'),limits=c(-.5,1))  + 
    #scale_color_gradientn(colours=inlmisc::GetColors(length(seq(0,20,5)),scheme='BuRd'),limits=c(-1,1))  + 
    #scale_color_gradientn(low="yellow",high="red",limits=c(-.5,1))  + 
  #  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  #scale_color_gradient2(limits=c(-2,2))  + 
    xlim(0, 1.2) +
    theme(legend.position = c(0,0), plot.margin = margin(2,2,2,2,"cm"))
  
  print(p2)
} else {
  p =
    ggtree(fish.tree,ladderize = FALSE, size=.6) + 
    geom_tree(aes(color=trait), continuous=T, size=.25) + 
    #ggtree(fish.tree, layout='circular',ladderize = FALSE, size=.5) + 
    #geom_tree(aes(color=trait), continuous=T, size=.3) + 
    #scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
    scale_color_gradientn(colours=inlmisc::GetColors(length(seq(0,20,5)),scheme='turbo'),limits=c(0.5,1.5))  + 
    xlim(0, 1.2) +
    theme(legend.position = c(0,0), plot.margin = margin(2,2,2,2,"cm"))
  
  print(p)
}

if (makelab == 1){
# label phyla
phytools::plotTree(my_tree,fsize=0.1,ftype="i",show.tip.label=TRUE)
  
  idx<-which(my_tree$node.label=="Demospongiae")
  nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
  
idx<-which(my_tree$node.label=="Vertebrata (subphylum in Deuterostomia)")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Cnidaria")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Arthropoda")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Mollusca")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Bryozoa")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Tunicata")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Echinodermata")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
idx<-which(my_tree$node.label=="Deuterostomia")
nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5,frame='none')# label all nodes on tree
}
