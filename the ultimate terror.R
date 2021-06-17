###----------------------------------------------------
save.image("C:/Users/Monika/Desktop/najnowsze analizy do 2 art/the ultimate terror/rrrr.RData")
##zapis pliku

##----jackard sorensen i bray curtis i czego tam jeszcze matka nie miala

library(vegan)
#przygotowanie do jaccarda

utjacc <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
utjacc

jacc = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  
  return ((M.11) / (M.11 + M.10 + M.01))
}

input.variables = utjacc

m = matrix(data = NA, nrow = length(input.variables), ncol = length(input.variables))
for (r in 1:length(input.variables)) {
  for (c in 1:length(input.variables)) {
    if (c == r) {
      m[r,c] = 1
    } else if (c > r) {
      m[r,c] = jacc(input.variables[,r], input.variables[,c])
    }
  }
}

variable.names = sapply(input.variables, attr, "label")
colnames(m) = variable.names
rownames(m) = variable.names   

jaccards = m


####------------------sorensen
Sorensen = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  
  return ((M.10+M.01) / (M.11+ M.11 + M.10 + M.01))
}

input.variables = utjacc

a = matrix(data = NA, nrow = length(input.variables), ncol = length(input.variables))
for (r in 1:length(input.variables)) {
  for (c in 1:length(input.variables)) {
    if (c == r) {
      a[r,c] = 1
    } else if (c > r) {
      a[r,c] = Sorensen(input.variables[,r], input.variables[,c])
    }
  }
}

variable.names = sapply(input.variables, attr, "label")
colnames(a) = variable.names
rownames(a) = variable.names   

sorensens = a
library(CommEcol)
dist<-dis.nness(DF, m=1, ness=TRUE)
dist2<-vegdist(DF, method="morisita")

library(xlsx)
# Write the first data set in a new workbook
write.xlsx(mds.bray, file = "myworkbook2.xlsx",
           sheetName = "sorensens", append = FALSE)
DF<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

braycurtis<-vegdist(DF, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
bcmatrix<-as.matrix(braycurtis)
bcmatrix
bcdf<-as.data.frame(bcmatrix)
bcdf
write.xlsx(bcdf, file = "myworkbook3.xlsx",
           sheetName = "braycurtis", append = FALSE)
##turnover and nestedness
DF<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

#####----------------
buba = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  if (M.10<M.01){
    M.2=M.10
  }else {
    M.2=M.01
  }
  
  return (M.2/ (M.11 + M.2))
}

input.variables = utjacc

d = matrix(data = NA, nrow = length(input.variables), ncol = length(input.variables))
for (r in 1:length(input.variables)) {
  for (c in 1:length(input.variables)) {
    if (c == r) {
      d[r,c] = 1
    } else if (c > r) {
      d[r,c] = buba(input.variables[,r], input.variables[,c])
    }
  }
}

variable.names = sapply(input.variables, attr, "label")
colnames(d) = variable.names
rownames(d) = variable.names   

simpson = d


write.xlsx(simpson, file = "myworkbook4.xlsx",
           sheetName = "simpsonbeta", append = FALSE)
####---------------------obliczenia ze stezeniami
stezenia <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
stezenia

# IMPORT ------------------------------------------------------------------
wig <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

wig$Azotyny <- as.numeric(wig$Azotyny) 
wig$SFI <- as.factor(wig$SFI)
wig$pora.roku <- as.factor(wig$pora.roku)

install.packages("ggplot2")

library(ggplot2)
library(FactoMineR)
library(factoextra)

summary(wig)
# PCA ---------------------------------------------------------------------

pca1 <- PCA(wig, quali.sup = c(1,16,17))
summary(pca1)

fviz_eig(pca1)
fviz_pca_var(pca1, col.var = "black")

fviz_pca_biplot(pca1, repel = T, col.var = "blue", col.ind = "black")

fviz_pca_biplot(pca1, col.ind = "black", palette = "BrBG", col.var = "black", repel = T, habillage = 1)



# GLM ---------------------------------------------------------------------

betaana <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
betaana
plot(betaana)
hist(wig$Azotany..V.)

m1 <- glm(wig$Azotany..V. ~ wig$SFI)
summary(m1)

ggplot(wig, aes(x = SFI, y = Azotany..V.)) + geom_boxplot()


m2 <- glm(wig$Chlorki ~wig$SFI)
summary(m2)

ggplot(wig, aes(x = SFI, y = Chlorki)) + geom_boxplot()


##----------------obliczenia do IOJ

pr12pr <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
danekosmiczne
plot(danekosmiczne)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(factoextra)
library(FactoMineR)

pca2 <- PCA(danekosmiczne, quali.sup = c(1,16,17))
summary(pca1)

fviz_eig(pca2)
fviz_pca_var(pca2, col.var = "black")

fviz_pca_biplot(pca2, repel = T, col.var = "blue", col.ind = "black")


fviz_pca_biplot(pca2, col.ind = "black", palette = "BrBG", col.var = "black", repel = T, habillage = 1)

ultimateterror
library(ggplot2)
# Boxplot iq versus condition
boxplot(ultimateterror$DI ~ ultimateterror$group, main = "Boxplot",
        xlab = "group", ylab = "DI")

# Plot the data
ggplot(data = ultimateterror, aes(x = group, y = DI)) +
  geom_line() + geom_point()
##zabawy z PCoA i NMDS

install.packages(c('biom','vegan'),repo='http://cran.wustl.edu')
library('biom')
library('vegan')
# load mapping file
gatunki <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
pr12pr2 <- sweep(pr12pr,2,colSums(pr12pr),'/')
# load mapping file
opis <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)


# find the overlapping samples
common.ids <- intersect(rownames(opis), rownames(gatunki2))
# get just the overlapping samples
gatunki2 <- gatunki2[common.ids,]
opis <- opis[common.ids,]

library(xlsx)
write.xlsx(gatunki3, file = "myworkbook2.xlsx",
           sheetName = "sorensens", append = FALSE)
# Keep only species present in at least 5% of samples
# This is fairly aggressive but will reduce the clutter in biplots
gatunki5 <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)

gatunki3 <- gatunki2[,colMeans(gatunki2>0)>.05]
gatunki4<- gatunki[,colMeans(gatunki2>0)>.5]
# run CA using vegan command
my.ca <- cca(pr12pr2)

plot(my.ca)

my.ca$CA$eig/my.ca$tot.chi

# run CA using vegan command
my.cca <- cca(pr12pr2 ~ pora.roku + stanowisko + SFI, data=opis)

plot(my.cca)

my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi

# Test what fraction of CA1 variance is explained in CCA1
b[1]/a[1]

# get Euclidean distance
d.euc <- dist(pr12pr2)

# get Bray-Curtis distances (default for Vegan)
d.bray <- vegdist(pr12pr2)

d.bray
# get Chi-square distances using vegan command
# we will extract chi-square distances from correspondence analysis
my.ca <- cca(pr12pr2)
d.chisq <- as.matrix(dist(my.ca$CA$u[,1:2]))
# get the indices of the 10 dominant OTUs
otu.ix <- order(colMeans(pr12pr2),decreasing=TRUE)[1:10]

# PCA on the subset including 10 OTUs
pca.otus <- princomp(pr12pr2[,otu.ix])$scores

# For comparison, do PCoA on Euclidean distances
pc.euc <- cmdscale(dist(pr12pr2[,otu.ix]))

# plot the PC1 scores for PCA and PCoA
# note: we might have to flip one axis because the directionality is arbitrary
# these are perfectly correlated
plot(pc.euc[,1], pca.otus[,1])
# NMDS using Vegan package
mds.bray <- metaMDS(pr12pr)$points
# makes a gradient from red to blue
my.colors <- colorRampPalette(c('red','blue'))(10)

layer <- opis[,'pH']
plot(mds.bray[,1], mds.bray[,2], col=my.colors[layer], cex=3, pch=16)
# Run PCoA (not PCA)
pc.bray <- cmdscale(d.bray,k=2)
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)


# use matrix multiplication to calculated weighted average of each taxon along each axis
wa <- t(pr12pr) %*% pc.bray

# plot the PCoA sample scores
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)

# get the indices of the 10 dominant OTUs
otu.ix <- order(colMeans(pr12pr),decreasing=TRUE)[1:10]
#####poprawic
# PCA on the subset including 10 OTUs
pca.otus <- princomp(pr12pr[,otu.ix])$scores

# For comparison, do PCoA on Euclidean distances
pc.euc <- cmdscale(dist(pr12pr[,otu.ix]))

# plot the PC1 scores for PCA and PCoA
# note: we might have to flip one axis because the directionality is arbitrary
# these are perfectly correlated
plot(pc.euc[,1], pca.otus[,1])
# First get a normalized version of the OTU table
# where each OTU's relative abundances sum to 1
otus.norm <- sweep(pr12pr,2,colSums(gatunki3),'/')

# use matrix multiplication to calculated weighted average of each taxon along each axis
wa <- t(otus.norm) %*% pc.bray

# plot the PCoA sample scores
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)

# add points and labels for the weighted average positions of the top 10 dominant OTUs
points(wa[otu.ix,1],wa[otu.ix,2],cex=2,pch=1)
text(wa[otu.ix,1],wa[otu.ix,2],colnames(pr12pr)[otu.ix],pos=1, cex=.5)

# Keep only species present in at least 5% of samples
# This is fairly aggressive but will reduce the clutter in biplots
gatunki4 <- gatunkia[colMeans(gatunki2>0)>.05,]
gatunki3
scaleformclust <- scale(pr12pr)
scaleformclust2 <- scale(gatunki3)
Mclust(scaleformclust, G = NULL, modelNames = NULL, 
       prior = NULL, 
       control = emControl(), 
       initialization = NULL, 
       warn = mclust.options("warn"), 
       x =  NULL 
      )
mc <- Mclust(scaleformclust)
mc2 <- Mclust(scaleformclust2) 
summary(mc) 
summary(mc2) 
mc$modelName
mc2$modelName# Optimal selected model ==> "VVV"
mc$G                        # Optimal number of cluster => 3
head(mc$z, 30)              # Probality to belong to a given cluster
head(mc$classification, 30) # Cluster assignement of each observation
library(factoextra)
# BIC values used for choosing the number of clusters
fviz_mclust(mc, "BIC", palette = "jco")
# Classification: plot showing the clustering
fviz_mclust(mc, "classification", geom = "point", 
            pointsize = 1.5, palette = "jco")
# Classification uncertainty
fviz_mclust(mc, "uncertainty", palette = "jco")

#NMDS 2020 ciag dalszy
#NMDS using Vegan package
gatunki3
##------transformacja


mds.bray2 <- metaMDS(...)$points
# makes a gradient from red to blue
my.colors <- colorRampPalette(c('red','blue'))(10)

layer <- opis[,'pH']
plot(mds.bray[,1], mds.bray[,2], col=my.colors[layer], cex=3, pch=16)
# Run PCoA (not PCA)
pc.bray <- cmdscale(d.bray,k=2)
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)


# use matrix multiplication to calculated weighted average of each taxon along each axis
wa <- t(gatunki3) %*% pc.bray

# plot the PCoA sample scores
plot(pc.bray[,1], pc.bray[,2], col=my.colors[layer], cex=3, pch=16)



####PERMANOVA with Adonis net
x = data
adonis(x~factors)
##Now, pairwise.adonis() for my factor location
factors=factors$Location
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
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
PW.Adonis=pairwise.adonis(x,factors,sim.method="bray",p.adjust.m = "bonferroni")
write.table(PW.Adonis,"Adonis-Results.csv",sep=",")

