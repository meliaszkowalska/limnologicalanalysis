###analysis of beta and alpha diversity changes with chemical and physical variables of water
save.image("C:/Users/Monika/Desktop/uzupe³nienie artyku³u 2 marzec 2021/r lion fight/lionfightenvironment2.RData")

##species analysis to calculate beta indices
##read table fromplik dla wojtka/dane kosmiczne - tylko gatunki 1%w2próbkach species which where at least 1%@2samples
#species<- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
#species

#install.packages("labdsv")
library(labdsv)
##hellinger transformation for species - for other methods done in those methods
species_h<-hellinger(species)
species_h_data_frame<-as.data.frame(species_h)
library(betapart)
getwd()

install.packages("writexl")
library("writexl")
#presence absence data import
#speciespresence_absence <-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1) 
###the idea of this is to make big dataframe with beta, alpha and chemistry distaces are variables and pairs of samples 
#are cases

#sorensen
betapresence_absence_sorensen<-beta.pair(speciespresence_absence, index.family="sorensen")

m.sim <- as.matrix(betapresence_absence_sorensen$beta.sim)
m.sne <- as.matrix(betapresence_absence_sorensen$beta.sne)
m.sor <- as.matrix(betapresence_absence_sorensen$beta.sor)
m.sim

######flatt your corellation matrix so useful, why on earth I didn't find it earlier ^^'''''
#it goes like ths -> ->
#                 -> -> _> and so on
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

m.sim<-flattenCorrMatrix(m.sim)
m.sim
m.sne<-flattenCorrMatrix(m.sne)
m.sor<-flattenCorrMatrix(m.sor)

write.csv(m.sim, "sim.csv")
write.csv(m.sne, "sne.csv")
write.csv(m.sor, "sor.csv")

m.simdf <- as.data.frame(m.sim)
m.snedf <- as.data.frame(m.sne)
m.sordf <-as.data.frame(m.sor)

write_xlsx(m.simdf, "sim.xlsx")
write_xlsx(m.snedf, "sne.xlsx")
write_xlsx(m.sordf, "sor.xlsx")

#jaccard
betapresence_absence_jaccard<-beta.pair(speciespresence_absence, index.family="jaccard")
betapresence_absence_jaccard


m.simj <- as.matrix(betapresence_absence_jaccard$beta.jtu)
m.snej <- as.matrix(betapresence_absence_jaccard$beta.jne)
m.sorj <- as.matrix(betapresence_absence_jaccard$beta.jac)

write.csv(m.simj, "simj.csv")
write.csv(m.snej, "snej.csv")
write.csv(m.sorj, "sorj.csv")

m.simdfj <- as.data.frame(m.simj)
m.snedfj <- as.data.frame(m.snej)
m.sordfj <-as.data.frame(m.sorj)

write_xlsx(m.simdfj, "simj.xlsx")
write_xlsx(m.snedfj, "snej.xlsx")
write_xlsx(m.sordfj, "sorj.xlsx")

#####calculeted beta for presence absence



######################################################################
##############calculation for LCBD
#install.packages("adespatial")
library("adespatial")
species

localbetaindices<-beta.div(species, method = "ab.sorensen", sqrt.D = FALSE, samp = TRUE,
         nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
localbetaindiceshellinger<-beta.div(species, method = "hellinger", sqrt.D = TRUE, samp = TRUE,
                           nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
localbetaindices

localbetaindiceshellinger



if(require("vegan", quietly = TRUE) & require("adegraphics", quietly = TRUE)){
  res = beta.div(species, "hellinger", nperm=999)
  # # Plot a map of the LCDB indices using the Cartesian coordinates
  # data(mite.xy)
  # s.value(mite.xy, res$LCBD, symbol = "circle", col = c("white", "brown"), main="Map of mite LCBD")
  ### Example using the mite abundance data and the percentage difference dissimilarity
  res = beta.div(species, "percentdiff", nperm=999, clock=TRUE)
  # Plot a map of the LCDB indices
  signif = which(res$p.LCBD <= 0.05) # Which are the significant LCDB indices?
  nonsignif = which(res$p.LCBD > 0.05) # Which are the non-significant LCDB indices?
  # g1 <- s.value(mite.xy[signif,], res$LCBD[signif], ppoint.alpha = 0.5, plegend.drawKey = FALSE,
                # symbol = "circle", col = c("white", "red"), main="Map of mite LCBD (red = significant indices)")
  # g2 <- s.value(mite.xy[nonsignif,], res$LCBD[nonsignif], ppoint.alpha = 0.5,
                # symbol = "circle", col = c("white", "blue"))
  # g2+g1
}
signif
nonsignif
betadivcompquanttrBS<-beta.div.comp(species, coef = "BS", quant = TRUE, save.abc = FALSE)
BSm<-as.matrix(betadivcompquanttrBS$repl)
BS<-as.matrix(betadivcompquanttrBS$repl)
BSm
BS<-flattenCorrMatrix(BS)

BSdf<-as.data.frame(BS) 
# 
write_xlsx(BSdf, "BSdf.xlsx")
###rich
BSrich<-as.matrix(betadivcompquanttrBS$rich)
BS<-as.matrix(betadivcompquanttrBS$rich)
BSrichdif<-BS
BS<-flattenCorrMatrix(BS)
BSdfrich<-as.data.frame(BS) 
# 
write_xlsx(BSdfrich, "BSdfrich.xlsx")

betadivcompquanttrBS
#D
BSD<-as.matrix(betadivcompquanttrBS$D)
BS<-as.matrix(betadivcompquanttrBS$D)
BS<-flattenCorrMatrix(BS)
BS

BSdD<-as.data.frame(BS) 
# 
write_xlsx(BSdD, "BSdfD.xlsx")
########################alfa matrixes-distance wyp³aszczenie i osi¹gniêcie pe³nej tablicy do analiz
#alpha_matrix <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
#last_table <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
alpha_matrix

TaxaS<-dist(alpha_matrix$Taxa_S, method = "euclidean",diag = FALSE, upper = FALSE)
TaxaS
TaxaSMatrix<-as.matrix(TaxaS)
flattentaxaS<-flattenCorrMatrix(TaxaSMatrix)
last_table$TaxaS<-flattentaxaS$cor
last_table
####
Dominance_D<-dist(alpha_matrix$Dominance_D, method = "euclidean",diag = FALSE, upper = FALSE)
Dominance_D
Dominance_DMatrix<-as.matrix(Dominance_D)
flattenDominance_D<-flattenCorrMatrix(Dominance_DMatrix)
last_table$Dominance_D<-flattenDominance_D$cor
######
Simpson_1.D<-dist(alpha_matrix$Simpson_1.D, method = "euclidean",diag = FALSE, upper = FALSE)
Simpson_1.DMatrix<-as.matrix(Simpson_1.D)
flattenSimpson_1.D<-flattenCorrMatrix(Simpson_1.DMatrix)
last_table$Simpson_1.D<-flattenSimpson_1.D$cor
last_table

########Shannon_H
Shannon_H<-dist(alpha_matrix$Shannon_H, method = "euclidean",diag = FALSE, upper = FALSE)
Shannon_HMatrix<-as.matrix(Shannon_H)
flattenShannon_H<-flattenCorrMatrix(Shannon_HMatrix)
last_table$Shannon_H<-flattenShannon_H$cor
last_table
#######Evenness_e.H.S
Evenness_e.H.S<-dist(alpha_matrix$Evenness_e.H.S, method = "euclidean",diag = FALSE, upper = FALSE)
Evenness_e.H.SMatrix<-as.matrix(Evenness_e.H.S)
flattenEvenness_e.H.S<-flattenCorrMatrix(Evenness_e.H.SMatrix)
last_table$Evenness_e.H.S<-flattenEvenness_e.H.S$cor
last_table
#######Brillouin
Brillouin<-dist(alpha_matrix$Brillouin, method = "euclidean",diag = FALSE, upper = FALSE)
BrillouinMatrix<-as.matrix(Brillouin)
flattenBrillouin<-flattenCorrMatrix(BrillouinMatrix)
last_table$Brillouin<-flattenBrillouin$cor
last_table
#######Menhinick
Menhinick<-dist(alpha_matrix$Menhinick, method = "euclidean",diag = FALSE, upper = FALSE)
MenhinickMatrix<-as.matrix(Menhinick)
flattenMenhinick<-flattenCorrMatrix(MenhinickMatrix)
last_table$Menhinick<-flattenMenhinick$cor
last_table
#######Margalef
Margalef<-dist(alpha_matrix$Margalef, method = "euclidean",diag = FALSE, upper = FALSE)
MargalefMatrix<-as.matrix(Margalef)
flattenMargalef<-flattenCorrMatrix(SMargalefMatrix)
last_table$Margalef<-flattenMargalef$cor
last_table
######Equitability_J
Equitability_J<-dist(alpha_matrix$Equitability_J, method = "euclidean",diag = FALSE, upper = FALSE)
Equitability_JMatrix<-as.matrix(Equitability_J)
flattenEquitability_J<-flattenCorrMatrix(Equitability_JMatrix)
last_table$Equitability_J<-flattenEquitability_J$cor
last_table
######Fisher_alpha
Fisher_alpha<-dist(alpha_matrix$Fisher_alpha, method = "euclidean",diag = FALSE, upper = FALSE)
Fisher_alphaMatrix<-as.matrix(Fisher_alpha)
flattenFisher_alpha<-flattenCorrMatrix(Fisher_alphaMatrix)
last_table$Fisher_alpha<-flattenFisher_alpha$cor
last_table
######Chao.1
Chao.1<-dist(alpha_matrix$Chao.1, method = "euclidean",diag = FALSE, upper = FALSE)
Chao.1Matrix<-as.matrix(Chao.1)
flattenChao.1<-flattenCorrMatrix(Chao.1Matrix)
last_table$Chao.1<-flattenChao.1$cor
last_table
namesforrows <- paste(BS$row, BS$column)
namesforrows <-as.matrix(namesforrows)
namesforrows
rownames(last_table) = namesforrows[,1]
last_table
####chemia do dist
#wigchemia_SFI_pory <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
str(wigchemia_SFI_pory)
wigchemia_SFI_pory[,c(1:13)]
#oznaczenie które kategoryczne
wigchemia_SFI_pory$SFI <- as.factor(wigchemia_SFI_pory$SFI)
wigchemia_SFI_pory$season <- as.factor(wigchemia_SFI_pory$season)

#standarization of environmental and biotic data
wigchemia_SFI_pory[,c(1:13)]<-scale(wigchemia_SFI_pory[,c(1:13)], center = TRUE, scale = TRUE)
wigchemia_SFI_pory[,c(1:13)]
wigchemia_SFI_pory
### doloz pc1ipc2 i chemie
########chlorki
chlorki<-dist(wigchemia_SFI_pory$chlorides, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(chlorki)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$chlorki<-flattenchlorki$cor
last_table
#######weglany carbonate tymczasowe zmienne
carbonate<-dist(wigchemia_SFI_pory$carbonate, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(carbonate)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$weglany<-flattenchlorki$cor
last_table
#######siarczany sulphate
sulphate<-dist(wigchemia_SFI_pory$sulphate, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(sulphate)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$siarczany<-flattenchlorki$cor
last_table
#######azotany nitrates
nitrates<-dist(wigchemia_SFI_pory$nitrates, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(nitrates)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$azotany<-flattenchlorki$cor
last_table
#######fosforany phosphate
phosphate<-dist(wigchemia_SFI_pory$phosphate, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(phosphate)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$fosforany<-flattenchlorki$cor
last_table
######lit Lithium
Lithium<-dist(wigchemia_SFI_pory$Lithium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(Lithium)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$lit<-flattenchlorki$cor
last_table
######sod Sodium
Sodium<-dist(wigchemia_SFI_pory$Sodium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(Sodium)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$sod<-flattenchlorki$cor
last_table
######amon Ammonium
Ammonium<-dist(wigchemia_SFI_pory$Ammonium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(Ammonium)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$amon<-flattenchlorki$cor
last_table
#####potas Potassium
chlorki<-dist(wigchemia_SFI_pory$Potassium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(chlorki)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$potas<-flattenchlorki$cor
last_table
#####magnez Magnesium
chlorki<-dist(wigchemia_SFI_pory$Magnesium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(chlorki)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$magnez<-flattenchlorki$cor
last_table
####wapn Calcium
Calcium<-dist(wigchemia_SFI_pory$Calcium, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(Calcium)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$wapn<-flattenchlorki$cor
last_table
####ph pH
chlorki<-dist(wigchemia_SFI_pory$pH, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(chlorki)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$ph<-flattenchlorki$cor
last_table

####przewodno conductivity
chlorki<-dist(wigchemia_SFI_pory$conductivity, method = "euclidean",diag = FALSE, upper = FALSE)
chlorkiMatrix<-as.matrix(chlorki)
flattenchlorki<-flattenCorrMatrix(chlorkiMatrix)
last_table$przewodnosc<-flattenchlorki$cor
last_table
last_table$namesforrows <- paste(BS$row, BS$column)
write_xlsx(last_table, "last_table.xlsx")
#
#
#########################################################################################
#### and this is ready to go table

distance_dataframe = subset(last_table, select = -c(namesforrows) )

#########################################################################################
#second table - sample/lake groups
#the beta indices values where grouped in lakes/seasons, mean value for every point where calculated
###
#dataframe_singlevalues <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
dataframe_singlevalues

##chemistry
##wigchemia_SFI_pory <- read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
str(dataframe_singlevalues)
#oznaczenie które kategoryczne
dataframe_singlevalues$SFI <- as.factor(dataframe_singlevalues$SFI)
dataframe_singlevalues$pora.roku <- as.factor(dataframe_singlevalues$pora.roku)
#standarization of environmental and biotic data
sub<-dataframe_singlevalues[,0:30]
sub
sub<-scale(sub, center = TRUE, scale = TRUE)
sub
sub2<-dataframe_singlevalues[,31:32]

dataframe_singlevalues_scaled<- cbind(sub,sub2)

dataframe_singlevalues_scaled
#one hot encoding
library("mltools")
library(data.table)


chemistryonehot <- one_hot(as.data.table(dataframe_singlevalues_scaled[,31:32]))
chemistryonehot
sub
dataframe_singlevalues_onehot<- cbind(sub,chemistryonehot)
#pearson correlation between variables - to delete >0.7
cormat<-cor(sub[,c(18:30)], method = "pearson", use = "complete.obs")
cormat
chemistrywithoutcorr<-sub[,c(18:30)]
cormat<-cor(chemistrywithoutcorr[,c(3,4,5,8,11)], method = "pearson", use = "complete.obs")
cormat
chemistryonehot
###bin >0.8
chemistrywithoutcorr08<-chemistrywithoutcorr[,c(1,3,4,5,8,11,12)]
###bin >0.7
chemistrywithoutcorr07<-chemistrywithoutcorr[,c(3,4,5,8,11)]
indscaled<-dataframe_singlevalues_onehot[,c(1:17)]
Chemistry08<-cbind(chemistrywithoutcorr08, indscaled,chemistryonehot)
Chemistry07<-cbind(chemistrywithoutcorr07, indscaled,chemistryonehot)

#######two tables with all single point data
###0.8 corellation border
Chemistry08
#####0.7 corellation border
Chemistry07



library(corrplot)
corrplot(cormat, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

###PCA 
#####do³ó¿ dane z PCA!
names<-read.table(file="clipboard", sep = "\t", dec = ",", head = T, row.names = 1)
names
nam<-names[,c(1)]
row.names(Chemistry07)<-nam
Chemistry07
chem<-cbind(names,Chemistry07)
Ch07PCA<-chem[,c(2:6)]
Ch08PCA<-Chemistry08[,c(1:7)]

PCA07<-prcomp(Ch07PCA, scale = FALSE)
PCA08<-prcomp(Ch08PCA, scale = FALSE)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library(factoextra)

fviz_eig(PCA07)
fviz_eig(PCA08)

fviz_pca_ind(PCA07,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(PCA07,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(PCA08,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
PCA07

fviz_pca_biplot(PCA07, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(PCA07)
eig.val

# Results for Variables
res.var <- get_pca_var(PCA07)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(PCA07)
###first two is chemistry combined
coordin<-res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 
####RDA for chemistry and wskaŸniki
#install.packages("ggplot2")

library(ggplot2)
library(FactoMineR)
library(factoextra)

summary(Chemistry07)
onlychemicaldata<-chem[,2:6]
onlyindices<-chem[,7:23]
library(vegan)
library(ggplot2)
library(ggrepel)
###############################RDA
#LCBDALL

x <- rda(onlychemicaldata~ onlyindices$LCBDall)
plot(x)

plot(ord)

y<- prcomp(onlychemicaldata, scale=FALSE)
biplot(y)
ord <- rda(onlychemicaldata~ onlyindices$LCBDall)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################LCBDrepl

x <- rda(onlychemicaldata~ onlyindices$LCBDrepl)
plot(x)

plot(ord)

ord <- rda(onlychemicaldata~ onlyindices$LCBDrepl)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################LCBDnest

ord <- rda(onlychemicaldata~ onlyindices$LCBDnest)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################Beta.sorensen

ord <- rda(onlychemicaldata~ onlyindices$Beta.sorensen)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)

#####################################Beta.sorensen.turn

ord <- rda(onlychemicaldata~ onlyindices$Beta.sorensen.turn)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################alpha taxa

ord <- rda(onlychemicaldata~ onlyindices$Taxa)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################dominance

ord <- rda(onlychemicaldata~ onlyindices$Dominance)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################alpha simson

ord <- rda(onlychemicaldata~ onlyindices$Simpson)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
#####################################alpha shannon

ord <- rda(onlychemicaldata~ onlyindices$Shannon)
plot(ord, las=1, bty="l", col="red", pch=19)
ord

anova(ord)
anova(ord, by="axis", permutations=499)
####idziemy z mrm testem, s¹ ju¿ wszystkie matrixy
#wzi¹æ poszczególne matrixy ze sob¹ i przeprowadziæ testy korelacji

####MRM test
# variables in the non-redundant abiotic and biotic variable sets (Martiny et al. 2011).
# Then, we removed the non-significant variables from this initial MRM test and re-ran
# the test. The significance of the partial regression was tested 999 times by a matrix
# permutation

# data(graze)
# # Abundance of this grass is related to forest cover but not location
# MRM(dist(LOAR10) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=10)
# # Abundance of this legume is related to location but not forest cover
# MRM(dist(TRRE3) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=10)
# # Compare to presence/absence of grass LOAR10 using logistic regression
# LOAR10.presence <- ifelse(graze$LOAR10 > 0, 1, 0)
# MRM(dist(LOAR10.presence) ~ dist(sitelocation) + dist(forestpct),
#     data=graze, nperm=10, method="logistic")
library(ecodist)
dist_chem<-dist(onlychemicaldata, method = "euclidean",diag = FALSE, upper = FALSE)
sulphate
nitrates
phosphate
Ammonium
Calcium
onlyindices
m.simj <- as.matrix(betapresence_absence_jaccard$beta.jtu)
m.snej <- as.matrix(betapresence_absence_jaccard$beta.jne)
m.sorj <- as.matrix(betapresence_absence_jaccard$beta.jac)
MRM(betapresence_absence_jaccard$beta.jtu ~ sulphate + nitrates , nperm=1000)
MRM(betapresence_absence_jaccard$beta.jne ~ sulphate   +Calcium, nperm=1000)
MRM(betapresence_absence_jaccard$beta.jac ~  nitrates  +Calcium, nperm=1000)
MRM(betadivcompquanttrBS$repl ~  nitrates + phosphate  , nperm=1000)
MRM(betadivcompquanttrBS$rich ~  nitrates, nperm=1000)
MRM(betadivcompquanttrBS$D ~  nitrates + phosphate, nperm=1000)
MRM(BSrepl) ~ (dist_chem + nitrates,  nperm=1000)

####LBCD
##replacement matrix
BSrepl<-as.matrix(betadivcompquanttrBS$repl)
###nestedness
BSrich<-as.matrix(betadivcompquanttrBS$rich)
####dissimilarity matrix
BSD<-as.matrix(betadivcompquanttrBS$D)

MRM(dist(LOAR10) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=10)


###random forest for indexes and variables

distance_dataframe










####single values

dataframe_singlevalues



#single values for lakes seasons and sfi glm!
dataframe_singlevalues


#Then, an optimal number of 2000
# trees was produced using cross-validation (Elith et al. 2008). The importance of each
# predictor variable was determined by its frequency of selection (for splitting)
# weighted by a measure of improvement of the model given each split and averaged
# across all the trees (contributions were scaled to sum to 100). To reduce the effect of
# spurious relationships between variables, we first ran the RF test with all the selected
# variables. Then, we removed the variable with the lowest contribution and re-ran the
# test until the lowest contribution of each variable was greater than 5%.







#####varpart for seasons, chemistry, sfi
#We performed variation partitioning analyses (Anderson and Cribble 1998) to
# reveal the effects of the spatial, environmental and biotic variables on beta diversity,
# the LCBD and their components. All the significant environmental and biotic
# variables were selected by forward selection against the biological characteristics data
# with 9999 permutations for all three taxonomic groups.


# These above analyses were performed in the R environment using the following
# packages, such as ‘betapart’ V1.5.1 (Baselga et al. 2018), ‘randomForestSRC’ V2.9.0
# (Liaw and Wiener 2002), ‘vegan’ V2.5-4 (Borcard and Legendre 2002), ‘ecodist’
# V2.0.1 (Goslee and Urban 2007), and ‘SpatialEpi’ V1.2.3 package (Kim and
#                                                                  Wakefield 2010).



####INDVAL
