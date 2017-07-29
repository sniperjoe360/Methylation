########
# Run the RandomForest Classifier
# for Bone Sarcomas. July 29, 2017
# 
# Developer: S. Peter Wu, M.D. 
# peter.wu10@nyumc.org
#
# for those who wish to process their own beta-values, 
# the (very large!) .idat files, please see GSE
# 
# for those who wish to use the beta-values within this file 
# (ie. beta_tCGA, beta_Tirado, beta_TARGET_OS)
# please contact me at above email address 
# 
# the csv files for the beta-values of the test samples (1,2,3)
# are found within this branch of the edit. 
#
# this is the final version of this code and will not be maintained further, 
# the next will be published on CRAN. 
# 
########

# read the gene csv with the B-values for the methylation profile
#bset = read.csv("R:\\karajannislabspace\\Projects\\Project 07 Sarcoma\\Results_Gold_7_11_2016_4\\tracks_and_tables_data\\csv\\betas_1.csv")
# find the 400 most variable genes

beta<-read.csv('H:\\beta_UltimateTest3.csv') # Fill in the directory where this file is stored
bset_test3<-data.frame(beta)
bset_test3$ID<-bset_test3$X
bset_test3$X <- NULL

beta<-read.csv('H:\\beta_UltimateTest2.csv')
bset_test2<-data.frame(beta)
bset_test2$ID<-bset_test2$X
bset_test2$X <- NULL

beta<-read.csv('H:\\beta_ultimatetest1.csv')
bset_test1<-data.frame(beta)
bset_test1$ID<-bset_test1$X
bset_test1$X <- NULL

load('H:\\bset_tcga')
bset2<-beta_tCGA
load('H:\\bset_target_OS.R')
bset3<-beta_tagetOS
load('H:\\bset_tirado_EWS.R')
bset4<-beta_tirado_EWS

data.dir <-"R:\\"

Grouping <- read.csv(file.path(data.dir, "ClassifierGrouping_GoldStandard_Final_edited.csv"))
# Fill in the directory where this file is stored

###### Variance 3 ########


# # merge all data into one set
# bset_all<-merge(bset,bset2,by="ID")
# bset_all<-merge(bset_all,bset3,by="ID")
# bset_all<-merge(bset_all,bset4,by="ID")


# bset_all<-bset_all[,c(1,6:44,49:58,63:234,239:253)]


# # order by variance
# # different version of variance 
# # select the top 400 most variable probes
# variance = apply(bset_all[1:length(bset_all)], 1, var)
# bset_all$variance = variance
# M<-bset_all[order(variance,decreasing=TRUE),]
# n=400
# B<-M[1:n,]

# B$variance = NULL

# bset = B
# # # merge all data into one set
# bset<-merge(B[,1:39],bset2,by="ID")
# bset<-merge(bset,bset3,by="ID")
# bset<-merge(bset,bset4,by="ID")

# bset<-bset[,c(6:44,49:58,63:234,239:253)]


########### Variance 1 ########

# select probes that have maximum variance
# across all values however have minimum variance
# within each set

bset_training_SS  <-bset[,c(1,6,10,13,28:35)]
bset_training_OS  <-bset[,c(1,14,15,24:27,36:44)]
bset_training_EWS <-bset[,c(1,8,11,16:23)]

mydata = data.frame(cbind(bset_training_SS,bset_training_OS,bset_training_EWS))
row.names(mydata)=mydata$ID
mydata <- mydata[,c(2:12,14:28,30:39)]

# calculate the one way ANOVA
library(limma)

Group = c('SS','SS','SS','SS','SS','SS','SS','SS','SS','SS','SS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','OS','EWS','EWS','EWS','EWS','EWS','EWS','EWS','EWS','EWS','EWS')

# transform into M-values
mydata_M <- log2(mydata/(1-mydata))

dmat<-model.matrix(~ Group)
fit<-lmFit(mydata_M,dmat)
eBfit <- eBayes(fit)
topTable(eBfit) # top 10 hits

# significant
nrow(topTable(eBfit, number = Inf, p.value = 0.01/dim(beta)[1]))

# top  hits
rows<-row.names(topTable(eBfit, number = Inf, p.value = 0.01/dim(beta)[1]))

# # write a fun

# bset[,c(1,6,10,13,28:35),'Group'] = 'Synovial'
# bset[,c(1,14,15,24:27,36:44),'Group'] = 'Osteosarcoma'
# bset[,c(1,8,11,16:23),'Group'] = 'Ewings'



# # order by variance
# # different version of variance 
# # select the top 400 most variable probes
# variance_SS = apply(bset_training_SS[2:length(bset_training_SS)], 1, var)
# variance_OS = apply(bset_training_OS[2:length(bset_training_OS)], 1, var)
# variance_EWS = apply(bset_training_EWS[2:length(bset_training_EWS)], 1, var)
# variance_all = apply(bset[6:length(bset)], 1, var)

# # create now maximum variance across values
# variance <- variance_all - variance_SS - variance_OS - variance_EWS


# instead of variance use the p-val from one-way ANOVA
M<-bset[bset$ID %in% rows,]
n=nrow(M)
B<-M[1:n,]


#annotate the genes
annogene<-read.csv("H:\\Snuderl\\450KManifest.csv")

# # merge all data into one set
bset<-merge(B,bset2,by="ID")
bset<-merge(bset,bset3,by="ID")
bset<-merge(bset,bset4,by="ID")
bset<-merge(bset,bset_test,by="ID")
bset<-merge(bset,bset_test1,by="ID")
bset<-merge(bset,bset_test2,by="ID")
bset<-merge(bset,bset_test3,by="ID")


#bset_merge<-bset[,c(1,6:44,49:58,63:234,239:253)]
write.csv(merge(B,annogene,by='ID'),"H:\\Snuderl\\Annogene_5508.csv")

write.csv(H:\\Snuderl\\Annogene_5508.csv")


bset<-bset[,c(6:44,49:58,63:234,239:length(bset))]

########### Variance 2 ########

# # order by variance
# # different version of variance 
# # select the top 400 most variable probes
# variance = apply(bset[6:length(bset)], 1, var)
# bset$variance = variance
# M<-bset[order(variance,decreasing=TRUE),]
# n=400
# B<-M[1:n,]

# # # merge all data into one set
# bset<-merge(B,bset2,by="ID")
# bset<-merge(bset,bset3,by="ID")
# bset<-merge(bset,bset4,by="ID")

# bset<-bset[,c(6:44,49:58,63:234,239:253)]


####### return to algorithm


# remove all columns except the values
# create column based on groups

# load the groupings
groups <- read.csv('H:\\Snuderl\\ClassifierGrouping_GoldStandard_Final_edited.csv')[,c(1,2,3)]
# define a new datatable

B2 <- t(bset)# transpose

rownames(groups) <- groups[,1]
B3<-merge(B2,groups,by="row.names")

# remove the last several columns, retain grouping
B4 <- B3[,c(2:nrow(bset),nrow(bset)+3)]
#B4 <- B3[,c(3:n,n+4)] # only for variation 3


# remove NA (mixed group)
B5<-B4[complete.cases(B4),]
colnames(B5)[ncol(B5)] <- 'Grouping'
# get the random forest running


# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite(c("plyr"))
# biocLite(c("minfi","randomForest","gridExtra","ggplot2", "glmnet","IlluminaHumanMethylation450kmanifest","IlluminaHumanMethylationEPICmanifest","SummarizedExperiment","conumee","CopyNumber450kData","limma"))
# biocLite(c("limma"))
# biocLite(c("openssl"))
# biocLite(c("preprocessCore"))

# install.packages('RColorBrewer')
# install.packages('randomForest')
# install.packages('gplots')
library(gplots)

library(RColorBrewer)
library(randomForest)

save(B5,file="H:\\Snuderl\\B6_Final_Data_input_1_31_2016.Rdata")
load("H:\\Snuderl\\B6_Final_Data_input_1_31_2016.Rdata")

B5_SR <- B5[B5$Grouping %in% c("EWS","Osteosarcoma","Synovial sarcoma","TARGET_OS","Unknown3","Unknown2","Unknown1"),]

B6<-droplevels(B5_SR[c(1:36),])
B6<-droplevels(B5_SR[c(1:36,57,123:125),])
B6<-droplevels(B5_SR[c(1:45,157),])
B6<-droplevels(B5_SR[c(1:36,seq(62,233,2)),])

B6<-droplevels(B5[c(1:94,seq(94, 265, 2)),])

colnames(B6)[ncol(B6)] <- 'Grouping'

# 55, 68 are the problematic ones

B_rf <- randomForest(Grouping~.,data=B6,mtry=28,proximity=TRUE,ntree=1000)
#ROC curve
bestmtry <- tuneRF(B6,B6$Grouping, ntreeTry=42, 
     stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)

# generate random forest tree
table(predict(B_rf),B6$Grouping)
plot(margin(B_rf,B6$Grouping))
legend(20,0.6, legend=levels(B6$Grouping),pch=15,col=1:5,cex=0.7)

levels(B6$Grouping) <- c("EWS (N=10)","Osteosarcoma (N=15)","Synovial Sarcoma (N=11)","TARGET-40-PASEFS","Unknown #1","Unknown #2","Unknown #3")

#B6$Grouping[37] <- "TARGET-40-PASEFS"
darkcols <- brewer.pal(12, "Set2")


MDSplot(B_rf,B6$Grouping,pch=15,bg=as.numeric(B6$Grouping)+1,cex=1.2, main = "Multidimensional Scaling Plot \n(Full Model)",xlab="", ylab = "",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,palette=darkcols)
legend(-0.45,0.5, legend=levels(B6$Grouping),pch=15,col=darkcols,cex=0.7)

# run the randomforest on the remaining samples
# for single sample
predict(B_rf,B5_SR[c(57,123:125),],'prob')

# for all the non-training samples
predict(B_rf,B5_SR[c(37:61,seq(62, 233, 2),233),],'prob')


# All test sequences w/ uterine
table(predict(B_rf,B5[c(28:59,69:94,seq(94, 265, 2)),]),B5[c(28:59,69:94,seq(94, 265, 2)),]$Grouping)

# All test sequences
table(predict(B_rf,B5[c(37:61,seq(62, 233, 2),233),]),B5[c(37:61,seq(62, 233, 2),233),]$Grouping)

# original
heatmap(as.matrix(t(B5[c(1:36),1:nrow(bset)-1])),col = bluered(128),scale='column')

# Rwings
heatmap(as.matrix(t(B5[c(1:36,47:61),1:nrow(bset)-1])),col = bluered(128),scale='column')

png("H:\\heatmap.png",width=100,height=80)

par(mar=c(20,5,20,5))
heatmap(as.matrix(t(B5[c(1:61,seq(62, 233, 2),233),1:nrow(bset)-1])),col = bluered(90),Colv=NA)

dev.off()

heatmap(as.matrix(t(B5[c(1:36,47:61),1:nrow(bset)-1])),col = redblue(128),scale='column')

heatmap(as.matrix(t(B5[c(1:61,seq(62, 233, 2),233),1:nrow(bset)-1])),col = redblue(128))
save.image()

#install.packages(ROCR)
library(ROCR)

B_rf_pr = predict(B_rf,newdata=B6,'prob')

# find importance
importance(B_rf)
varImpPlot(B_rf)
