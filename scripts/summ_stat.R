
library(PopGenome)

Hcan <- readData("data")
Hcan
Hcan@n.sites
Hcan@region.names

get.sum.data(Hcan)

#get.individuals(Hcan)
#Hcan <- set.populations(Hcan,list(c("CcanDQ062140"), 
#                                  c("Ccan11", "Ccan17", "Ccan18", "Ccan01", "Ccan27", "Ccan30", "Ccan55", "Ccan48",
#                                    "Ccan08", "Ccan15", "Ccan26", "Ccan29", "Ccan19", "Ccan06", "Ccan50", "Ccan24", 
#                                    "Ccan03", "Ccan31", "Ccan49", "Ccan53", "Ccan25", "Ccan21", "Ccan10", "Ccan23",
#                                    "Ccan05", "Ccan37", "Ccan04", "Ccan16", "Ccan12")))


#Hcan  <- detail.stats(Hcan)
#get.detail(Hcan)[[1]]

Hcan  <- diversity.stats(Hcan, pi = T)
get.diversity(Hcan)[[1]]

Hcan  <- neutrality.stats(Hcan, do.R2=T)
get.neutrality(Hcan, theta=T)[[1]]


ms <- MS(Hcan,detail=T,thetaID="Tajima",neutrality=T)
ms@prob.equal[[1]]
ms@prob.less[[1]]
MS_getStats(ms)





library(ape)

dist_Hcan <- read.dna("data/cytBtotal.fas", format = "fasta")
dist_Hcan
dist.dna(dist_Hcan)
dist_HcanPY<-dist.dna(dist_Hcan, model = "F84", as.matrix = T)[-1,1]
mean(dist_HcanPY)
dist_HcanBR<-dist.dna(dist_Hcan, model = "F84", as.matrix = T)[-1,-1]
diag(dist_HcanBR)<-NA
mean(dist_HcanBR, na.rm = T)





library(pegas)
haps_Hcan <- haplotype(dist_Hcan)
haps_Hcan
net_Hacan <- haploNet(haps_Hcan)
ind.hap<-with(stack(setNames(attr(haps_Hcan, "index"), rownames(haps_Hcan))),
              table(hap=ind, individuals=rownames(dist_Hcan)[values]))

mydata <- as.data.frame(ind.hap)
good <- mydata[mydata$Freq == 1,]
locations <- c(rep("PNE-BR", 29), "Paraguay")
new.hap <- table(good$hap, locations)

png(filename = "Hap.png", width = 20, height = 20, units = "cm", bg = "transparent", res = 600)
plot(net_Hacan, size = attr(net_Hacan, "freq"), scale.ratio = 2, cex = 0.8, pie=new.hap, label=F, bg=c("black", "gray"))
legend(x=-30, y=9, c("Paraguay", "PNE-BR"), pch = 16, pt.cex = 2, cex = 1.1, col = c("black", "gray"), bty="n")
dev.off()





library(rhierbaps)
baps_Hcan <- load_fasta(dist_Hcan)
resbaps_Hcan <- hierBAPS(baps_Hcan, max.depth = 2, n.pops = 2, n.extra.rounds = 100, assignment.probs = T)
resbaps_Hcan$partition.df






