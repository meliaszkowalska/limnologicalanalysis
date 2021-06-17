####---------------squere root transf

save.image("C:/Users/Monika/Desktop/najnowsze analizy do 2 art/the ultimate terror/rrrr.RData")

if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}


library("car")


library(MASS)
library(rcompanion)
opis
gatunki5
# run CA using vegan command
my.ca2 <- cca(gatunki5)

plot(my.ca2)

my.ca2$CA$eig/my.ca$tot.chi

# run CA using vegan command
my.cca2 <- cca(gatunki5 ~ SFI + pora.roku + przewodnosc, data=opis)

plot(my.cca2)

my.cca2$CCA$eig/my.cca$tot.chi
a <- my.ca2$CA$eig/my.ca$tot.chi
b <- my.cca2$CCA$eig/my.cca$tot.chi

# Test what fraction of CA1 variance is explained in CCA1
b[1]/a[1]

# get Euclidean distance
d.euc5 <- dist(gatunki5)

# get Bray-Curtis distances (default for Vegan)
d.bray5 <- vegdist(gatunki5)

d.bray5
# get Chi-square distances using vegan command
# we will extract chi-square distances from correspondence analysis
my.ca5 <- cca(gatunki5)
d.chisq5 <- as.matrix(dist(my.ca5$CA$u[,1:2]))
# get the indices of the 10 dominant OTUs
otu.ix5 <- order(colMeans(gatunki5),decreasing=TRUE)[1:10]

# PCA on the subset including 10 OTUs
pca.otus5 <- princomp(gatunki5[,otu.ix5])$scores

# For comparison, do PCoA on Euclidean distances
pc.euc5 <- cmdscale(dist(gatunki5[,otu.ix5]))

# plot the PC1 scores for PCA and PCoA
# note: we might have to flip one axis because the directionality is arbitrary
# these are perfectly correlated
plot(pc.euc5[,1], pca.otus5[,1])
# NMDS using Vegan package
####proby przynajmniej 1% w 2 probach
pr12pr<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
pr12pr
mds.bray <- metaMDS(pr12pr)$points
gatunki
mds.bray5=metaMDS(pr12pr,k=3,trymax=100000000000)$points
# makes a gradient from red to blue
my.colors5 <- colorRampPalette(c('red','blue'))(10)

layer5 <- opis[,'pH']
plot(mds.bray5[,1], mds.bray5[,2], mds.bray5[,3], col=my.colors5[layer5], cex=3, pch=16)
library(xlsx)
write.xlsx(mds.bray5,file = "myworkbook10.xlsx",
           sheetName = "mds", append = FALSE)
library(plot3D)
scatter3D(mds.bray5[,1], mds.bray5[,2], mds.bray5[,3], theta = 15, phi = 20)
scatter3D(mds.bray5[,1], mds.bray5[,2], mds.bray5[,3], phi = 0, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed")
text3D(mds.bray5[,1], mds.bray5[,2], mds.bray5[,3],  labels = rownames(mds.bray5),
       add = TRUE, colkey = FALSE, cex = 0.5)

####matrixdissimilarity
#przygotowanie do diss

mds.bray5

de




# Run PCoA (not PCA)
pc.bray <- cmdscale(d.bray,k=2)
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)



# use matrix multiplication to calculated weighted average of each taxon along each axis
wa <- t(gatunki3) %*% pc.bray

# plot the PCoA sample scores
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)

# get the indices of the 10 dominant OTUs
otu.ix <- order(colMeans(gatunki3),decreasing=TRUE)[1:10]
#####poprawic
# PCA on the subset including 10 OTUs
pca.otus <- princomp(gatunki3[,otu.ix])$scores

# For comparison, do PCoA on Euclidean distances
pc.euc <- cmdscale(dist(gatunki3[,otu.ix]))

# plot the PC1 scores for PCA and PCoA
# note: we might have to flip one axis because the directionality is arbitrary
# these are perfectly correlated
plot(pc.euc[,1], pca.otus[,1])
# First get a normalized version of the OTU table
# where each OTU's relative abundances sum to 1
otus.norm <- sweep(gatunki3,2,colSums(gatunki3),'/')

# use matrix multiplication to calculated weighted average of each taxon along each axis
wa <- t(otus.norm) %*% pc.bray

# plot the PCoA sample scores
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)

# add points and labels for the weighted average positions of the top 10 dominant OTUs
points(wa[otu.ix,1],wa[otu.ix,2],cex=2,pch=1)
text(wa[otu.ix,1],wa[otu.ix,2],colnames(gatunki)[otu.ix],pos=1, cex=.5)

# Keep only species present in at least 5% of samples
# This is fairly aggressive but will reduce the clutter in biplots
gatunki4 <- gatunkia[colMeans(gatunki2>0)>.05,]
gatunki4
scaleformclust <- scale(gatunki4)
scaleformclust2 <- scale(gatunki3)
library(mclust)

#NMDS 2020 ciag dalszy
#NMDS using Vegan package
gatunki3


#Use the function pairwise.adonis() with following arguments
#
#x = community table
#
#factors = a column or vector with all factors to be tested pairwise
#
#sim.function = which function to calculate the similarity matrix. eg 'daisy' or 'vegdist' default is 'vegdist'
# NOTE that if you wnat to use daisy, you need to install library 'cluster'
#
#sim.method = similarity method from daisy or vegdist: default is 'bray' 
#
#p.adjust.m = the p.value correction method, one of the methods supported by p.adjust(); default is 'bonferroni'
#
#The function will return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.value
#
# load in runnig R session with: 
# source('pairwise.adonis.txt')
#

# example:
# data(iris)
# pairwise.adonis(iris[,1:4],iris$Species)
#
#[1] "Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1"
#                    pairs    F.Model        R2 p.value p.adjusted sig
#1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *

# similarity euclidean from vegdist and holm correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')

# similarity manhattan from daisy and bonferroni correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',sim.method='manhattan',p.adjust.m='bonferroni')
##############################

 ###dopublikacji nmds! 
  
  example_NMDS=metaMDS(pr12pr,k=3,trymax=100)$points
library(vegan)
  stressplot(example_NMDS)
  
  scatter3D(example_NMDS)
  
  ordiplot(example_NMDS,type="n")
  orditorp(example_NMDS,display="species",col="red",air=0.01)
  orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
  
  ordiplot(example_NMDS,type="n")
  ordihull(example_NMDS,groups=opis$SFI,draw="polygon",col="grey90",
           label=FALSE)
  
  example_NMDS
  
  
  
  
  library(vegan)
  dune
  
  data(dune)
  data(dune.env)
  
  # default test by terms
  
  pr12pr
  
  adonis2(pr12pr ~ pora.roku*SFI*stanowisko, data = opis)
  
  # overall tests
  
  adonis2(pr12pr ~ pora.roku*SFI*stanowisko, data = opis, by = NULL)
  
  
  
  
  
  
  
  
  
  
  addonis <- adonis(formula = pr12pr ~ pora.roku*SFI*stanowisko, data = opis, permutations = 999) 

addonis
co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()


for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1")
return(pairw.res)

} 


#### adonis proba na nmds matrix from 3d calculation
NMDSmatrix3d <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
NMDSmatrix3d
dis<-as.dist(NMDSmatrix3d)
dis
opis
NMDSmatrix3d<-as.dist(NMDSmatrix3d)
adonis <- adonis(NMDSmatrix3d ~ pora.roku+jezioro, data = opis)
adonis
#czy stanowiska siê ró¿ni¹ pomiêdzy sob¹ jeœli zrobimy poprawkê na jezioro
adonisstanowiska<- adonis2(NMDSmatrix3d ~ jezioro/pora.roku, data = opis)
adonisstanowiska
#same stanowiska bez jeziora
adonisstanowiska2<- adonis2(NMDSmatrix3d ~ stanowisko, data = opis)
adonisstanowiska2
##SFI z poprawk¹ na jeziora
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ pora.roku+stanowisko/jezioro, strata = opis$jezioro, data = opis)
adonisstanowiska3
#dyspersja

dis
mod <- betadisper(dis, opis$jezioro)
mod
plot(mod)
mod2 <- betadisper(dis, opis$pora.roku)
plot(mod2)
opis
opis
betadisper(dis, opis$jezioro)
######z pubblikacji przyklad
## Perform test
anova(mod)

## Permutation test for F
permutest(mod2, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod2.HSD <- TukeyHSD(mod2))
plot(mod2.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod2)

## can also specify which axes to plot, ordering respected
plot(mod, axes = c(3,1))

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)

## `scores` and `eigenvals` also work
scrs <- scores(mod)
str(scrs)
head(scores(mod, 1:4, display = "sites"))
# group centroids/medians 
scores(mod, 1:4, display = "centroids")
# eigenvalues from the underlying principal coordinates analysis
eigenvals(mod) 

## try out bias correction; compare with mod3
(mod3B <- betadisper(dis, opis$stanowiskoxjezioro, type = "median", bias.adjust=TRUE))

## should always work for a single group
group <- factor(rep("stanowiskoxjezioro", NROW(opis)))
(tmp <- betadisper(dis, group, type = "median"))
(tmp <- betadisper(dis, group, type = "centroid"))


## Using group centroids
mod3 <- betadisper(dis, opis$stanowiskoxjezioro, type = "centroid")
mod3B
permutest(mod3, permutations = 999)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3))


######---- po permanovie jakie testy
library(ggplot2)
library(vegan)
library(grid)
opis

#######stanowiska siê od siebie ró¿ni¹ bior¹c pod uwagê wszystkie pory roku
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ stanowisko, strata = NULL, data = opis)
adonisstanowiska3
###jeziora siê od siebie ró¿ni¹
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ jezioro, strata = NULL, data = opis)
adonisstanowiska3
####pory roku siê od siebie ró¿ni¹ bez wzglêdu na jezioro
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ pora.roku, strata = NULL, data = opis)
##adonisstanowiska3
##model bior¹cy pod uwagê stanowiska, jeziora, pora roku<- stanowisko nic nie wnosi do modelu, za to jezioro i pora roku tak
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ jezioro * pora.roku, strata = NULL, data = opis)
adonisstanowiska3
##czy stanowiska na jeziorze w tej samej porze roku siê od siebie ró¿ni¹? - ze strat¹ i bez - nie ma ró¿nic
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ stanowisko/jezioro/pora.roku, strata = jezioro/pora.roku, data = opis)
adonisstanowiska3
##czy jeziora z poprawk¹ na porê roku siê ró¿ni¹? - ze strata i bez siê ró¿ni¹
adonisstanowiska3<- adonis2(NMDSmatrix3d ~ jezioro/pora.roku, strata = NULL, data = opis)
adonisstanowiska3


####z palca funkcja
pairwise.adonis3 <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis2(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

jeziora<-opis$jezioro


rm(jezioro)
########################
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
NMDSmatrix3d
adonisstanowiska3pairwise<- pairwise.adonis2(NMDSmatrix3d ~ jezioro, strata = NULL, data=opis)
adonisstanowiska3pairwise
write.xlsx(adonisstanowiska3pairwise,file = "adonisstanowiska3pairwise.xlsx",
           sheetName = "adonisstanowiska3pairwise", append = FALSE)
adonisstanowiska3pairwise2<- pairwise.adonis2(NMDSmatrix3d ~ jezioro, strata = 'jezioro', data=opis)
adonisstanowiska3pairwise2

write.xlsx(adonisstanowiska3pairwise2,file = "adonisstanowiska3pairwise2.xlsx",
           sheetName = "adonisstanowiska3pairwise2", append = FALSE)

matrixgat <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
opis

adonisstanowiska3pairwise<- pairwise.adonis2(matrixgat ~ stanowisko, strata = NULL, data=opis2)
adonisstanowiska3pairwise
write.xlsx(adonisstanowiska3pairwise,file = "adonisstanowiska3pairwise.xlsx",
           sheetName = "adonisstanowiska3pairwise", append = FALSE)
adonisstanowiska3pairwise2<- pairwise.adonis2(matrixgat ~ stanowiskoxjezioro, strata = 'jezioro', data=opis)
adonisstanowiska3pairwise2

##########chemistry pe³ne
##envdata <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
##ioj2009 <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
##ioj2019 <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
##ocena2009 <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
##ocena2019 <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
envdata<-normalize(envdata, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
ioj2009<-normalize(ioj2009, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
ioj2019<-normalize(ioj2019, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
ocena2009<-normalize(ocena2009, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
ocena2019<-normalize(ocena2019, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")


envdist<-dist(envdata, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ioj2009dist<-dist(ioj2009, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ioj2019dist<-dist(ioj2019, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ocena2009dist<-dist(ocena2009, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ocena2019dist<-dist(ocena2019, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

adonisstanowiska3pairwise<- pairwise.adonis2(envdist ~ stanowisko, strata = 'jezioro', data=opis2)
adonisstanowiska3pairwise
opis2<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

####test dla par ioj2009x2019

library(tidyverse)
library(ggpubr)
library(rstatix)
iojxioj<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
ocenaxocena<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
iojxioj2<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
ocenaxocena2<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
ggboxplot(iojxioj, x = "IOJ", y = "war", add = "jitter")
ggboxplot(ocenaxocena, x = "OCENA", y = "war", add = "jitter")

res.fried <- iojxioj %>% friedman_test(war ~ IOJ |id)
res.fried
res.fried <- ocenaxocena %>% friedman_test(war ~ OCENA |id)
res.fried

iojxioj %>%
  group_by(IOJ) %>%
  identify_outliers(war)
iojxioj %>%
  group_by(IOJ) %>%
  shapiro_test(war)
ggqqplot(iojxioj, "war", facet.by = "IOJ")
res.aov <- anova_test(data = iojxioj, dv = war, wid = id, within = IOJ)
get_anova_table(res.aov)
pwc <- iojxioj %>%
  pairwise_t_test(
    war ~ IOJ, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
# Subset weight data before treatment
iojxioj2<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
iojxioj
before <- iojxioj2$IOJ.2009
before
# subset weight data after treatment
after <- iojxioj2$IOJ.2019
after
# Plot paired data
library(PairedData)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()
res <- wilcox.test(before, after, paired = TRUE)
res

wilcox.test(war ~ IOJ, data = iojxioj, paired = TRUE,
            alternative = "greater")

###dla oceny
before <- ocenaxocena2$OCENA.2009
before
# subset weight data after treatment
after <- ocenaxocena2$OCENA.2019
after
# Plot paired data
library(PairedData)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()
res <- wilcox.test(before, after, paired = TRUE)
res
library('BiocManager')
#Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
library(ape)

#This package will also help us more easily manipulate our data
library(dplyr)
library(ggplot2)

#This package is used to calculate and plot Venn diagrams as well as heatmaps
library(gplots)
#Linear mixed-effects models like repeated measures analysis
library(lme4)
## Loading required package: Matrix
#used to read in mothur-formatted files
library(phangorn)

#The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(phyloseq)

#A package to create interactive web graphics of use in 3D plots
library(plotly)
library(tidyr)
library(vegan)
simper(matrixgat, opis$pora.roku, permutations=999)
kruskal.test(matrixgat$cymbella0affiniformis ~ opis$pora.roku)
matrixgat
library(clustsig)

######ró¿nice w porach roku ioj poszczególne stanowiska
# Subset weight data before treatment
#iojxioj3<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
iojxioj
before <- iojxioj3$jesien
before
# subset weight data after treatment
after <- iojxioj3$lato
after
# Plot paired data
library(PairedData)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()
res <- wilcox.test(before, after, paired = TRUE)
res

wilcox.test(war ~ pora.roku, data = iojxioj3, paired = TRUE,
            alternative = "greater")
iojxioj3
ggboxplot(iojxioj3, x = 'stanowisko', y = 'war', 
          color = "stanowisko", 
          
          ylab = "Weight", xlab = "Treatment")
####kruskal i odpowiednik anovy z post hoc
kruskal.test(war ~ jezioro, data = iojxioj3)
pairwise.wilcox.test(iojxioj3$war, iojxioj3$stanowisko,
                     p.adjust.method = "BH")
library(ggplot2)
library(multcomp)

iojxioj3$stanowisko <- as.factor(iojxioj3$stanowisko)
glht(model, mcp(stanowisko="Tukey"))

model2 = glm(war ~ stanowisko,
             data = iojxioj3, family='poisson')

model2
model = lm(war ~ stanowisko,
           data = iojxioj3)
iojxioj3
war

model
model3<-glm(formula = war ~ jezioro, family = poisson, data = opis)
model3

###############################RDA

envdata
write.xlsx(envdata,file = "myworkbook10.xlsx",
           sheetName = "mds", append = FALSE)
##ograniczenie zmiennych
envdataogr<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
#pca

pca12 <- PCA(envdataogr)
summary(pca12)
print(pca12)

fviz_eig(pca12)
fviz_pca_var(pca12, col.var = "black")

fviz_pca_biplot(pca12, repel = T, col.var = "blue", col.ind = "black")

fviz_pca_biplot(pca12, col.ind = "black", palette = "BrBG", col.var = "black", repel = T, habillage = 1)

pcawig
#po normalizacji
library(BBmisc)
IOJ2009<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
IOJ2019<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
envdataogr<-normalize(envdataogr, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
IOJ2019<-normalize(IOJ2019, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
IOJ2009<-normalize(IOJ2009, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")

envdataogrt<-as.data.frame(envdataogrt)
x <- rda(envdataogr~ IOJ2019$IOJ.2019)
plot(x)

plot(ord)

y<- prcomp(envdataogr, scale=FALSE)
biplot(y)
ord <- rda(envdataogr ~ IOJ2019$IOJ.2019)
plot(ord, las=1, bty="l", col="red", pch=19)
help("plot")
ord
ord2009 <- rda(envdataogr ~ IOJ2009$IOJ.2009)
plot(ord2009)
ord2009
anova(ord)
anova(ord2009)
anova(ord, by="axis", permutations=499)
anova(ord2009, by="axis", permutations=499)
 #################dbRDA
ord <- capscale(envdataogr ~ IOJ2019$IOJ.2019)
plot(ord)
ord
ord2009 <- capscale(envdataogr ~ IOJ2009$IOJ.2009)
plot(ord2009)
ord2009
anova(ord) ## overall test of the significance of the analysis
anova(ord, by="axis", perm.max=500) ## test axes for significance
anova(ord, by="terms", permu=200) ## test for sig. environ. variables

anova(ord2009) ## overall test of the significance of the analysis
anova(ord2009, by="axis", perm.max=500) ## test axes for significance
anova(ord2009, by="terms", permu=200) ## test for sig. environ. variables


###########obliczenia PCA wygenerowane
#wstawienie tablicy
#jezioro1
opis2
opis2 <- read.table(file="clipboard", sep = "\t", dec = ",", head = F, row.names = 1)
jezioro1<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
jezioro1 = sqrt(jezioro1)
library(vegan)
d.bray.jezioro1 <- vegdist(jezioro1)


adonis <- adonis(jezioro1 ~ opis2$saplehistory, data = opis2)
adonis
pc.jezioro1 <- cmdscale(d.bray.jezioro1,k=2)

library(geosphere)
###odleglosc od srodka
poczatek<-pc.jezioro1[1:4,1:2]
wszystkie<-pc.jezioro1[,1:2]
koniec<-centroid(pc.jezioro1[,1:2])
newtable<-rbind(poczatek,wszystkie,koniec)
newtable
library(xlsx)
write.xlsx(newtable,file = "centroids.xlsx",
           sheetName = "centroids1", append = FALSE)
my.colors <- c("#00AFBB")

plot(pc.jezioro1[,1], pc.jezioro1[,2], col=my.colors, cex=3, pch=16)
points(pc.jezioro1[c(1:4, 1001:1004, 2001:2004, 3001:3004, 4001:4004, 5001:5004, 6001:6004, 7001:7004, 8001:8004, 9001:9004,10001:10004, 11001:11004),1], 
       pc.jezioro1[c(1:4, 1001:1004, 2001:2004, 3001:3004, 4001:4004, 5001:5004, 6001:6004, 7001:7004, 8001:8004, 9001:9004,10001:10004, 11001:11004),2], col="#FC4E07", cex=3, pch=16)

text(pc.jezioro1[c(1:4, 1001:1004, 2001:2004, 3001:3004, 4001:4004, 5001:5004, 6001:6004, 7001:7004, 8001:8004, 9001:9004,10001:10004, 11001:11004),1],
     pc.jezioro1[c(1:4, 1001:1004, 2001:2004, 3001:3004, 4001:4004, 5001:5004, 6001:6004, 7001:7004, 8001:8004, 9001:9004,10001:10004, 11001:11004),2], labels=rownames(opis2), cex=0.9, font=20)

#################################

glmofcentroid<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
glmofcentroid$jezioro<-as.factor(glmofcentroid$jezioro)
glmcentr <- lm(glmofcentroid$cent ~ glmofcentroid$jezioro,
                       data = glmofcentroid)
ggplot(glmofcentroid, aes(x = jezioro, y = cent)) + geom_boxplot()
summary(glmcentr)
pairwise.t.test(glmofcentroid$cent ,glmofcentroid$jezioro, p.adj = "bonf")
library(car)
muma<-Anova(glmcentr, type="III")
res.aov111 <- aov(glmofcentroid$cent ~ glmofcentroid$jezioro, data = glmofcentroid)
TukeyHSD(res.aov111)


summary(lm(glmofcentroid$cent~ glmofcentroid$pora.roku,
           data = glmofcentroid))
options(scipen=4)
attributes(glmcentr)
plot(glmcentr)
# print default contrasts
contrasts(glmofcentroid$pora.roku)
library(multcomp)

#################################

glmofsrednia<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

glmsrodek <- lm(glmofsrednia$odsredniej ~ glmofsrednia$rok +glmofsrednia$pora,
               data = glmofsrednia)
ggplot(glmofsrednia, aes(x = rok, y = odsredniej)) + geom_boxplot()
summary(glmsrodek)
Anova(model, type="III")
summary(lm(glmofsrednia$odsredniej ~ glmofsrednia$rok,
           data = glmofsrednia))
options(scipen=4)
attributes(glmsrodek)
plot(glmsrodek)



##########################################################

#histogramy

histogramjezioro1<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
histogramrzecz<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

hist(histogramjezioro1$IOJ, breaks=20, main="With breaks=30")

stripchart(histogramrzecz,
           method = "jitter",
           pch = 23,
           bg = "red",
           add = TRUE)
########################################################
library(pairwiseAdonis)
write.xlsx(iojxioj3,file = "iojxioj3.xlsx",
           sheetName = "mds", append = FALSE)
iojxioj3<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
iojpairwise<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
iojdist<-dist(iojxioj3$warioj, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
iojpairwise$IOJ2019
opis$IOJ
adonissimple<-adonis(formula = iojdist ~ iojxioj3$rok, data = opis, permutations = 999)
adonissimple
adonisiojpairwise<- pairwise.adonis(iojdist, opis$IOJ, sim.method = "bray",
                                    p.adjust.m = "none")

adonisiojpairwise

ggboxplot(iojxioj3, x = 'jezioro', y = 'war', 
          color = "jezioro", 
          
          ylab = "IOJ2019", xlab = "Lake")
ggboxplot(iojxioj3, x = 'pora.roku', y = 'war', 
          color = "pora.roku", 
          
          ylab = "IOJ2019", xlab = "season")
opis
dis
