# Code to plot phylogenetic tree of Eo and Ac
# load libraries
library(brms)
library(ape)
library(dplyr)
library(phytools)
library(tidybayes)
library(bayestestR)
library(ggtree)
library(treeio)
library(ggplot2)
library(corrplot)
library(visreg)

# clear environment
rm(list = ls())
# close plots
graphics.off()


# tick tock
ptm <- proc.time()

# parameters
brl_power = 0.4 #grafen's rho, 0.2-1
F1_min <- 0.8 # F1 fit filter
nobs_min <- 0 # nobs filter
trait<-"A"  # A for Ac or E for Eo
Amin <- 10000 # Ac threshold (no effect)
rmout <- 0 # saturate outliers
lgAc<- 1 # log10(1) or raw (0) Ac?

# read species list from obis data (filtered)
data <- read.csv("obis_Eoeff_z10.csv")
print(data)

# remove outlier Ac?
Amax<-which(data[,3]>Amin)
data[Amax,2]<-NaN
data[Amax,3]<-NaN

# filter poor F1 fits
data[data[,5]<F1_min,2:3]<-NaN # F1>F1_min
data[data[,6]<nobs_min,2:3]<-NaN # nobs>nobs_min

# species with trait data
speciesq <- data[is.nan(data[,2])==FALSE,1]

# remove duplicated species 
#speciesq <- speciesq[!duplicated(speciesq)]

# extract species OTTs from tree of life
resolved_names <- rotl::tnrs_match_names(names = speciesq)
resolved_names

# remove species with approximate matches 
my_ott_ids <-resolved_names[,4]
my_matches <-resolved_names[,3]
my_matches[is.na(my_matches)]=TRUE
my_species <-resolved_names[,2]
my_ott_ids <-my_ott_ids[my_matches==FALSE]
my_species_pruned <-my_species[my_matches==FALSE]

# remove duplicates
#my_ott_ids <- my_ott_ids[!duplicated(my_species_pruned)]
#my_species_pruned <- my_species_pruned[!duplicated(my_species_pruned)]

# ott of species in obis but not found in induced subtree
probchild_ott<-c(3594036,1088122,151938,165111,191186,238110,238117,2828488,2848924,2848934,2850353,2850442,2850448,2851044,2851338,2851341,2851367,2852284, 2852393,2852465,2852674,2852749,2853292,2854977,2854991,2855001,2855060,2855447,2855657,2855872,2855974,2855990,2856007,2856052,2856139,2857609,2859772,2859990,2860070,2862316,2866162,2874044,2874096,2874129,2874142,2874147,2874162,2874176,2874187,2874252,2874255,2874276,2874346,2874377,2874406,2874409,2874414,2874424,2874451,2874479,2874481,2874495,2874497,2874498,2874545,2874561,2874625,2874666,2876686,2876688,2876700,2876701,2876702,2876721,2876743,2876747,2876767,2876818,2876870,2876931,2876936,2876994,2877025,2877043,2877206,2877212,2877239,2877258,2877293,2877403,2877549,2877556,2877583,2877638,2877657,2877674,2877872,2878166,2878221,2878296,2878297,2878740,2878766, 2878768,2878900,2878938,2883462,2886582,2887165,2887251,2894511,2894527,2894544,2894607,2894748,2894751,2894755,2894776,2894778,2894842,2894850,2894859,2894862,2894864,2894928,2894929,2894933,2894934,2894962,2895000,2895038,2895048,2895050,2895056,2895160,2895244,2895259,2895332,2895434,2895436,2900781,2901763,2902802,2902893,2904401,2909241,2909730,2912660,2913799,2915325,2915593,2915716,2917521,2917774,2918238,2918478,2920547,2920594,2920897,2922204,2922581,2935069,2935672,29472,29493,2958645,2958819,2961353,2961393,2961461,2963030,2963133,2963250,2963301,2963424,2964057,2965354,2965423,2965713,2965972,2966117,2966411,2966955,2983168,2989673,3001299,354630,3592081,3594036,3594073,3633170,3633670,3633813,3633831,3634612,3634622,3634685,3635438,3635440,3635473,3635479,3635490,3635498,3635508,3635524,3635528,3635531,3635534,3635539,3635580,3636496,3638528,3638533,3638536,3642351,3642353,3642357,3647943,3676445,3676809,375424,406416,4146663,4147074,4150760,4150763,4150842,4150884,4150938,4150998,4151000,4151136,4154865,4157095,4157530,420469,432262,445485,460905,476444,4956670,4958553,4960257,4967537,513965,5254156,5501286,5501537,5776365,635612,6383234,651758,691870,7493273,751007,780814,819585,864448,927296,965149,978677,987433,4157095,2828488,2895436,4157095,5254156,6383234,3633157,3634951,4162784,7072584,4146661,2877621,2878444,2915842,2966871,3635578,3647821,4146661,5501093,7064037,2963293,2874237,2874237,2918249,2964633,2966422,3635501,3636499,4958549,4147010,3633673,3636502,3594050,255562,2854749,2863603,2874149,2874511,2876981,2877104,2892580,2894580,2894747,2914896,2916891,2916892,2922635,2934888,2961411,2962281,2963199,2963207,2963272,2963621,2966118,2966769,3594050,3633629,3633633,3633650,3633737,3633789,3633824,3634608,3634610,3634689,3635494,3635521,3635536,3635576,3636490,3636494,3639841,3642352,3681815,4148765,88225,36423523,882253,3642352)
idxtoremove=probchild_ott*0
nn<-1
for (xx in 1:length(probchild_ott)) {
  
  # index of problem child
  probchild_idx1 <- which(my_ott_ids==probchild_ott[xx]) # pruned species list (resolved names w exact match)
  # name of problem child
  probchild<-my_species_pruned[probchild_idx1]
  
  # find index in full data set to remove them 
  probchild_idx3<-which(data[,1]==probchild)# full species list (all resolved names)
  # save index
  if (length(probchild_idx3)>0) {
    idxtoremove[nn]<-probchild_idx3
    nn<-nn+1
  } else {
    
    print(xx)
    print(probchild_ott[xx])
  }
  
}

# remove these species (small fraction of final tree)
for (xx in 1:length(idxtoremove)) {
  if (idxtoremove[xx]==0) {
  } else {
    # remove species
    data[idxtoremove[xx],2:3] <- NaN
  }
}

# other species not recognized
data[which(data[,1]=="Ervilia castanea"),2:3]=NaN
data[which(data[,1]=="Hesperisternia multangulus"),2:3]=NaN
data[which(data[,1]=="Cirsotrema magellanicum"),2:3]=NaN
data[which(data[,1]=="Strigamia maritima"),2:3]=NaN
data[which(data[,1]=="Ervilia nitens"),2:3]=NaN
data[which(data[,1]=="Cheilopogon pinnatibarbatus"),2:3]=NaN
data[which(data[,1]=="Lycodes diapterus"),2:3]=NaN
data[which(data[,1]=="Notoscopelus elongatus"),2:3]=NaN
data[which(data[,1]=="Adeonellopsis subsulcata"),2:3]=NaN
data[which(data[,1]=="Osmerus mordax"),2:3]=NaN

# mischaracterized 
data[which(data[,1]=="Fritillaria tenella"),2:3]=NaN
data[which(data[,1]=="Pontogeneia inermis"),2:3]=NaN
data[which(data[,1]=="Porina gracilis"),2:3]=NaN
data[which(data[,1]=="Ophidion rochei"),2:3]=NaN
data[which(data[,1]=="Notopogon lilliei"),2:3]=NaN


#try again after filtering out otts not in induced subtree
# species with trait data
speciesq <- data[is.nan(data[,2])==FALSE,1]
Eo <- data[is.nan(data[,2])==FALSE,2]
Ac <- data[is.nan(data[,2])==FALSE,3]
dEdT <- data[is.nan(data[,2])==FALSE,4]
F1 <- data[is.nan(data[,2])==FALSE,5]
nobs <- data[is.nan(data[,2])==FALSE,6]
phyla <- data[is.nan(data[,2])==FALSE,7]

# remove duplicated species 
#Ac <- Ac[!duplicated(speciesq)]
#Eo <- Eo[!duplicated(speciesq)]
#F1 <- F1[!duplicated(speciesq)]
#speciesq <- speciesq[!duplicated(speciesq)]

# extract species OTTs from tree of life
resolved_names <- rotl::tnrs_match_names(names = speciesq)
resolved_names

# Eo for species with resolved names
resolved_names_Eo <- Eo

# Ao for species with resolved names
resolved_names_Ac <- Ac

# remove species with approximate matches 
my_ott_ids <-resolved_names[,4]
my_matches <-resolved_names[,3]
my_matches[is.na(my_matches)]=TRUE
my_species <-resolved_names[,2]
my_ott_ids <-my_ott_ids[my_matches==FALSE]
my_Eo <- resolved_names_Eo[my_matches==FALSE]
my_Ao <- resolved_names_Ac[my_matches==FALSE]
my_species_pruned <-my_species[my_matches==FALSE]

# remove duplicates
#my_Ao <- my_Ao[!duplicated(my_species_pruned)]
#my_Eo <- my_Eo[!duplicated(my_species_pruned)]
#my_ott_ids <- my_ott_ids[!duplicated(my_species_pruned)]
#my_species_pruned <- my_species_pruned[!duplicated(my_species_pruned)]

# get subtree with queried taxa
my_tree <- rotl::tol_induced_subtree(ott_ids = my_ott_ids,label_format = "name")

# plot subtree
ape::plot.phylo(my_tree, cex =.05,show.node.label=TRUE) # or just plot(my_tree, cex = 2)

# Wilco's circle tree plot
phytools::plotTree(my_tree,fsize=0.1,ftype="i",type="fan",show.tip.label=TRUE)

# label phyla 
#idx<-which(my_tree$node.label=="Vertebrata (subphylum in Deuterostomia)")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Cnidaria")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Arthropoda")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Mollusca")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Bryozoa")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Tunicata")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Echinodermata")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Deuterostomia")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree
#idx<-which(my_tree$node.label=="Ctenophora")
#nodelabels(text=my_tree$node.label[idx],node=idx+Ntip(my_tree),cex=0.5)# label all nodes on tree


# get names in format
my_species_pruned_nospace <- my_species_pruned
for (xx in 1:length(my_species_pruned)) {
  # remove spaces
  my_species_pruned_nospace[xx] <- gsub(" ","_",my_species_pruned[xx])
}

# choose trait
if (trait == "E") {
  treetiptrait<-unique(data.frame(label=my_species_pruned_nospace,traitvalue=(my_Eo)))
} else {
  if (lgAc==1){
    treetiptrait<-unique(data.frame(label=my_species_pruned_nospace,traitvalue=log10(my_Ao)))
  } else {
    treetiptrait<-unique(data.frame(label=my_species_pruned_nospace,traitvalue=(my_Ao)))
  }
}


# extract species traits
traitdata<-treetiptrait$traitvalue

# extract species names
names(traitdata)<-treetiptrait$label

# estimate trait data for higher taxa
fit <- phytools::fastAnc(compute.brlen(my_tree, power=brl_power),traitdata,vars=TRUE,CI=TRUE)

# get node information for each trait
td <- data.frame(node = nodeid(my_tree, names(traitdata)),trait = traitdata)

# higher taxa nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# combine nodes with and w/o data
dd <- rbind(td, nd)

# extract node info as numeric
dd$node <- as.numeric(dd$node)

# saturate outliers for visualization
if (rmout==1){
  if (trait=="A") {
    dd[dd[,2]>1.5,2]<-1.5 # set outliers to upper 95th percentile for ease of visualization
    dd[dd[,2]<0.5,2]<-0.5 # set outliers to upper 95th percentile for ease of visualization
    #dd[dd[,2]>20,2]<-20 # set outliers to upper 95th percentile for ease of visualization
    
  }
  if (trait=="E") {
    dd[dd[,2]>1,2]<-1 # set outliers to upper 95th percentile for ease of visualization
    dd[dd[,2]<(-.5),2]<- (-.5) # set outliers to upper 95th percentile for ease of visualization
  }
}
# construct tree 
fish.tree <- full_join(compute.brlen(my_tree, power=brl_power), dd, by = 'node')

# lamda statistics
psL<-phylosig(compute.brlen(my_tree, power=brl_power),traitdata,method="lambda",test=TRUE)
psL_1<-phylosig(compute.brlen(my_tree, power=0.2),traitdata,method="lambda",test=TRUE)
psL_2<-phylosig(compute.brlen(my_tree, power=1),traitdata,method="lambda",test=TRUE)
print(psL)

print(proc.time() - ptm)

# save output 
save.image(file='obis_tree_F1pt8_log10Ac_3emin2pmax.RData')
