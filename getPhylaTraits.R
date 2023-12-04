rm(list = ls())

# load package
library(kites)

# hypoxia tolerance
data <- read.csv("Fig5B.csv")
dsum<-colSums(data)
dsumtot<-max(dsum)
plot_kite(data, 40,xlabtext='Active hypoxia tolerance, Aeco (1/atm)',xticks= seq(0,20,5))


# temperature sensitivity
data <- read.csv("Fig5D.csv")
dsum<-colSums(data)
dsumtot<-max(dsum)
plot_kite(data, 40,xlabtext='Temperature sensitivity, Eeco (eV)',xticks= seq(-2,2,0.5))
