# DAPC
library(adegenet)
library(ggplot2)
library(reshape2)
library(dplyr)
#setwd("/Users/lucianocosme/Dropbox/Albopictus/manuscript_chip/data/no_autogenous/analysis/adegenet")
#for windows
#setwd("C:/Users/degop/Dropbox/Albopictus/manuscript_chip/data/no_autogenous/analysis/adegenet")

# get the right format plink --bfile albo --recodeA --out albo

#import the data
####################################################################################
albo<-read.PLINK("./adegenet/albo.raw", quiet = FALSE, chunkSize = 1000, # map.file = "file1.map", (after plink file) 
                   parallel = require("parallel"), n.cores = 6)

#for windows
#albo<-read.PLINK("albo.raw", quiet = FALSE, chunkSize = 1000, parallel=FALSE)


# We specify that we want to evaluate up to k = 20 
grp <- find.clusters(albo, max.n.clust=15, n.pca = 200 )

#
names(grp)
head(grp$Kstat, 10)
grp$stat
head(grp$grp, 10)
grp$size

# number of mosquitoes per group
table(pop(albo),grp$grp)

groups2<-table(pop(albo),grp$grp)
write.csv(groups2, file = "groups.csv")



table.value(table(pop(albo), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("ori", 1:10))

# We run the analysis on the previous toy dataset, using the inferred groups stored in grp$grp:
dapc1 <- dapc(albo, grp$grp)

dapc1

?scatter.dapc
#scatter plot

scatter(dapc1)

#different options of the plots
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)


myCol <- c("darkblue","purple","green","orange","red","blue")
# #eb5fd1 pink south
# #915951 alboown east
# #ffb235 orange north-central
# #1eb81e green  north
# #4eb5ff blue  central-southeast
myCol <- c("#e88ad2",
           "#9e4725",
           "#ffb235",
           "#4eb5ff",
           "#e35235",
           "#c2db43",
           "#610ffa",
           "#5cd168",
           "#b2aae0",
           "#524d50",
           "#648a58",
           "#d6a080")


scatter(dapc1, posi.da="bottomright", bg="white",
        pch=15:19, cstar=0, col=myCol, scree.pca=FALSE, posi.pca="topright", legend = FALSE)
# legend(1,1,legend=c("North", "North-Central", "Central-Southeast", "South", "East"),
#        fill=c("#1eb81e","#ffb235","#4eb5ff","#eb5fd1","#915951"), cex=1,box.lwd = 0,box.col = "transparent",bg = "transparent")

dev.off()


# Restore default clipping rect
par(mar=c(5, 4, 4, 2) + 0.1)



jpeg("DAPC_pca2.jpeg", width = 6, height = 5, units='in', res = 600)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol,solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:10))
add.scatter.eig(dapc1$eig,15,1,2, posi="bottomright", inset=.02)
dev.off()


##
# scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0,
#         col=myCol, solid=.4, cex=3, clab=0, mstree=TRUE, scree.da=FALSE, 
#         posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:10))
# par(xpd=TRUE)
# points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,cex=3, lwd=8, col="black")
# points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,cex=3, lwd=2, col=myCol)
# 
# myInset <- function(){
#   temp <- dapc1$pca.eig
#   temp <- 100* cumsum(temp)/sum(temp)
#   plot(temp, col=rep(c("black","lightgrey"),c(dapc1$n.pca,1000)), ylim=c(0,100), 
#        xlab="PCA axis", ylab="Cumulated variance (%)",
#        cex=1, pch=20, type="h", lwd=2)}
# add.scatter(myInset(), posi="bottomright", inset=c(-0.03,-0.01), ratio=.28,
#             bg=transp("white"))
# 

# figure paper
jpeg("DF1.jpeg", width = 6, height = 4, units='in', res = 600)

scatter(dapc1,1,1, col=myCol, bg="white", scree.da=FALSE, legend=FALSE, solid=.4)
legend("top", legend=c(1:10),
       fill=c(myCol), cex=0.8,box.lwd = 0,box.col = "transparent",bg = "transparent", ncol=5)
#
dev.off()

jpeg("1.MACr.jpeg", width = 20, height = 10, units='in', res = 600)
scatter(dapc1,2,2, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc1,3,3, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc1,4,4, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

# myPal <- colorRampPalette(c("blue","gold","red"))
# scatter(dapc1, col=transp(myPal(10)), scree.da=FALSE,cell=1.5, cex=2, bg="white",cstar=0)

#paper
######################################################
# Loading plot
######################################################
set.seed(4)
options(scipen=999)
contrib <- loadingplot(dapc1$var.contr, axis=2,
                      lab.jitter=.2,thres=0.0003)

# Get the loci 

# chech the name of the loci 
albo$loc.names

# to get the list of loci above the threshold
contrib$var.idx

# get the ones above the threshold in the plot, example of one locus
albo$loc.names[3305]
# PS. it always show the reference allele on the output

# get the ones above the threshold in the plot, all loci
library(stringr)
albo$loc.names[c(contrib$var.idx)]
lociTop<-as.data.frame(albo$loc.names[c(contrib$var.idx)])
head(lociTop)

colnames(lociTop)<- c("SNP")

loadingSNps<- strsplit(lociTop$SNP, "_", fixed=TRUE)
head(loadingSNps)

dfSNPs<- do.call(rbind.data.frame, loadingSNps)
colnames(dfSNPs)<- c("SNP", "refA")
dfSNPs2<- cbind(dfSNPs,contrib$var.idx)
colnames(dfSNPs2)<- c("SNP", "refA", 'dapcID')
head(dfSNPs2)




######################################################
# Get the bim file
######################################################
bimFile<-read.table("albo.bim", header=F, row.names=NULL, sep="\t", stringsAsFactors = FALSE, colClasses = "character")
head(bimFile)
colnames(bimFile)<-c("chrom", "SNP", "V3", "bp", "REF", "ALT")
# remove column V3
bimFile2<- subset(bimFile, 
               select=c(-V3))

# merge with DAPC data
head(bimFile2)
head(head(dfSNPs))
bimFile3 <- merge(bimFile2, dfSNPs2, by= "SNP", all = FALSE)
head(bimFile3)
# get columns we want
bimFile4<- subset(bimFile3, 
                  select=c("dapcID", "chrom", "SNP", "bp", 'REF', "ALT"))
head(bimFile4)

# convert genlight to genind
######################################################

# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DART")

library(dartR)
albo2<-gl2gi(albo, probar = TRUE, verbose = NULL)


# change the order of the factors, following fastStructure
levels(albo2$pop)
albo2$pop <- factor(albo2$pop, levels = c('GELc',
                                          'BENc',
                                          'KUNc',
                                          'JAFc',
                                          'CHAc',
                                          'LAMc',
                                          'SUUc',
                                          'SUFc',
                                          'KATc',
                                          'KAGc',
                                          'LIBc',
                                          'AWKc',
                                          'TROc',
                                          'MORc',
                                          'MADc',
                                          'BARc',
                                          'TIRc',
                                          'ROMc',
                                          'IMPc',
                                          'KRAc',
                                          'NEGc',
                                          'PORc',
                                          'GRAc',
                                          'MAUc',
                                          'NOVc',
                                          'TUCc',
                                          'RECc',
                                          'SAIc',
                                          'PALc',
                                          'MACc',
                                          'HOUc',
                                          'BEAc',
                                          'PEOc',
                                          'COLc',
                                          'BERc',
                                          'SPRc',
                                          'NEWc',
                                          'MANc',
                                          'LOSc'))
levels(albo2$pop)                 

# 
# las
# numeric in {0,1,2,3}; the style of axis labels.
# 0: always parallel to the axis [default],
# 1: always horizontal,
# 2: always perpendicular to the axis,
# 3: always vertical.

head(bimFile4)
# check margins of plot
par("mar")

# get the frequencies
freq_3305<- tab(genind2genpop(albo2[loc=c("AX-583301796_T")]),freq=TRUE)
freq_27472<- tab(genind2genpop(albo2[loc=c("AX-582523725_C")]),freq=TRUE)
freq_57371<- tab(genind2genpop(albo2[loc=c("AX-581354204_G")]),freq=TRUE)
freq_52296<- tab(genind2genpop(albo2[loc=c("AX-580721505_C")]),freq=TRUE)
freq_36169<- tab(genind2genpop(albo2[loc=c("AX-579475215_G")]),freq=TRUE)

#######################################################
# trying to make it with ggplot
######################################################
# rbind the dataframes
tomelt<-cbind(freq_3305, freq_27472, freq_57371, freq_52296, freq_36169)

# rename
head(tomelt)
colnames(tomelt)<-c("AX-583301796_T", "AX-583301796_C",
                     "AX-582523725_C", "AX-582523725_T",
                     "AX-581354204_G", "AX-581354204_A",
                     "AX-580721505_C", "AX-580721505_T",
                     "AX-579475215_G", "AX-579475215_A")
head(tomelt)

library(reshape2)
melted_dapc <- melt(tomelt,id.var="id")
head(melted_dapc)

# split the melted file to get each allele
meltedAlleles<- strsplit(as.character(melted_dapc$Var2), "_", fixed=TRUE)
head(meltedAlleles)

# make it a dataframe and rename the column names
meltedAlleles2<- do.call(rbind.data.frame, meltedAlleles)
head(meltedAlleles2)
colnames(meltedAlleles2)<- c("SNP", "Allele")
head(meltedAlleles2)

# merge with the melted_dapc
meltedAlleles3 <- cbind(melted_dapc, meltedAlleles2)
head (meltedAlleles3)
colnames(meltedAlleles3)<- c("Population", "Var2", "Frequency", "SNP", "Allele")
head (meltedAlleles3)

# keep the columns needed
meltedAlleles4<- subset(meltedAlleles3, 
                  select=c("Population", "SNP", "Allele", 'Frequency'))
head(meltedAlleles4)


#plot
# ggplot(melted_dapc, aes(x=Var1,y=value,group=Var1,colour=Var2)) +
#   geom_point()+
#   geom_line(aes(color=Var2))
# 
# ggplot(melted_dapc, aes(x=Var1, y=value, group=Var2, colour=Var2)) +
#   geom_point() +
#   geom_line(aes(lty=Var2)) +
#   theme_classic()

head(meltedAlleles4)
base <- ggplot(meltedAlleles4, aes(x=Population, y=Frequency, group=Allele, color = Allele, fill = Allele)) + #, group=Allele, colour=SNP
  geom_blank() + 
  xlab(NULL) + 
  ylab(NULL)

p <- base + facet_wrap(~SNP, ncol = 1,  strip.position="right") +  #geom_point() +
  #geom_label(aes(label = Allele), size = 4, alpha = 1, label.size = NA, fill = "lightgray")+
  #geom_label(aes(label = Allele), size = 4, alpha = .1, label.size = NA)+
  geom_line(aes(linetype = Allele), size = .6, show.legend = FALSE) +
  geom_label(aes(label=Allele), colour= "white", label.padding=unit(0.05,"lines"), size=5, alpha = 1) +
  geom_text(aes(label=Allele), size=4, color =  "black") +
  scale_fill_manual(values = c("A" = "#fccee9", "G" = "#a0d3ed","C" = "#ffcc74", "T" = "#b5f283")) +
  scale_linetype_manual(values = c("A" = "solid", "T" = "dashed", 
                                     "C" = "dotted", "G" = "twodash")) +
  scale_color_manual(values = c("A" = "#f29ce6", "G" = "#2626f0","C" = "#f09a19", "T" = "#397511")) +
  theme(legend.position="none") +  
  theme_bw() +
  theme( 
    strip.background = element_rect(fill="#e6e2d6"),
    strip.text = element_text(face="bold")
  ) + theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45)) 

# theme(strip.text.y = element_text(angle = 0)
p
#
ggsave(filename = "loading_loci.pdf", p, width = 14, height = 8, units = "in")

#######################################################
#                   diapause SNPs
#######################################################
# import the file with the SNPs on diapause genes
diapause<-read.table("diapause.txt", header=T, row.names=NULL, sep="\t", stringsAsFactors = FALSE, colClasses = "character")
head(diapause)

diapause2 <- merge(bimFile2, diapause, by= "SNP", all = FALSE)
head(diapause2)

# list of SNPs that I used for dapc, some of the SNPs were filtered out
# snps with data
head(diapause2)
# paste(df$n,df$s,sep="-")
diapause2$allele <- paste(diapause2$SNP,diapause2$REF,sep="_")
diapause2$alleleA <- paste(diapause2$SNP,diapause2$ALT,sep="_")

head(diapause2$allele)
head(diapause2$alleleA)

# # doing it for just a few SNPs
# freq01<- tab(genind2genpop(albo2[loc=c("AX-580336192_A")]),freq=TRUE)
# freq02<- tab(genind2genpop(albo2[loc=c("AX-580336357_A")]),freq=TRUE)
# freq03<- tab(genind2genpop(albo2[loc=c("AX-580336395_A")]),freq=TRUE)
# freq04<- tab(genind2genpop(albo2[loc=c("AX-580336499_T")]),freq=TRUE)
# freq05<- tab(genind2genpop(albo2[loc=c("AX-580336526_C")]),freq=TRUE)
# get all frequencies
freq888<- tab(genind2genpop(albo2[loc=c(diapause2$allele)]),freq=TRUE)
head(freq888)

# rbind the dataframes
tomelt2<-cbind(freq888)

# rename, first check the columns
head(tomelt2)
head(diapause2$allele)
head(diapause2$alleleA)


# rename a few alleles
head(tomelt2)
# colnames(tomelt2)<-c("AX-580336192_A", "AX-580336192_T",
#                     "AX-580336357_A", "AX-580336357_T",
#                     "AX-580336395_A", "AX-580336395_T",
#                     "AX-580336499_T", "AX-580336499_A",
#                     "AX-580336526_C", "AX-580336526_T")

# all alleles
# make new column with the loci names (R and A)
diapause2$lociNames<-lapply(1:nrow(diapause2), function(i) c(diapause2$allele[i], diapause2$alleleA[i]))
head(diapause2$lociNames)

head(tomelt2)

# melt the column with the loci names
diapauseNames<-melt(diapause2$lociNames)
head(diapauseNames)

# change the names
colnames(tomelt2)<-diapauseNames$value

# check it
head(tomelt2)

# melt if with reshape
melted_dapc2 <- melt(tomelt2,id.var="id")
head(melted_dapc2)

# split the melted file to get each allele
meltedAlleles2<- strsplit(as.character(melted_dapc2$Var2), "_", fixed=TRUE)
head(meltedAlleles2)

# make it a dataframe and rename the column names
meltedAlleles3<- do.call(rbind.data.frame, meltedAlleles2)
head(meltedAlleles3)
colnames(meltedAlleles3)<- c("SNP", "Allele")
head(meltedAlleles3)

# merge with the melted_dapc
meltedAlleles4 <- cbind(melted_dapc2, meltedAlleles3)
head (meltedAlleles4)
colnames(meltedAlleles4)<- c("Population", "Var2", "Frequency", "SNP", "Allele")
head (meltedAlleles4)

# keep the columns needed
meltedAlleles5<- subset(meltedAlleles4, 
                        select=c("Population", "SNP", "Allele", 'Frequency'))
head(meltedAlleles5)



head(meltedAlleles5)
base <- ggplot(meltedAlleles5, aes(x=Population, y=Frequency, group=Allele, color = Allele, fill = Allele)) + #, group=Allele, colour=SNP
  geom_blank() + 
  xlab(NULL) + 
  ylab(NULL)

p2 <- base + facet_wrap(~SNP, ncol = 3,  strip.position="right") +  #geom_point() +
  #geom_label(aes(label = Allele), size = 4, alpha = 1, label.size = NA, fill = "lightgray")+
  #geom_label(aes(label = Allele), size = 4, alpha = .1, label.size = NA)+
  geom_line(aes(linetype = Allele), size = .6, show.legend = FALSE) +
  geom_label(aes(label=Allele), colour= "white", label.padding=unit(0.05,"lines"), size=5, alpha = 1) +
  geom_text(aes(label=Allele), size=4, color =  "black") +
  scale_fill_manual(values = c("A" = "#fccee9", "G" = "#a0d3ed","C" = "#ffcc74", "T" = "#b5f283")) +
  scale_linetype_manual(values = c("A" = "solid", "T" = "dashed", 
                                   "C" = "dotted", "G" = "twodash")) +
  scale_color_manual(values = c("A" = "#f29ce6", "G" = "#2626f0","C" = "#f09a19", "T" = "#397511")) +
  theme(legend.position="none") +  
  theme_bw() +
  theme( 
    strip.background = element_rect(fill="#e6e2d6"),
    strip.text = element_text(face="bold")
  ) + theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45)) 

# theme(strip.text.y = element_text(angle = 0)
p2
#
ggsave(filename = "diapause_loci.pdf", p2, width = 25, height = 20, units = "in")


#######################################################
#                   immunity SNPs
#######################################################
# import the file with the SNPs on diapause genes
immunity<-read.table("immunity_genes_palatini.txt", header=T, row.names=NULL, sep="\t", stringsAsFactors = FALSE, colClasses = "character")
head(immunity)

immunity2 <- merge(bimFile2, immunity, by= "SNP", all = FALSE)
head(immunity2)

# list of SNPs that I used for dapc, some of the SNPs were filtered out
# snps with data
head(immunity2)

immunity2$allele <- paste(immunity2$SNP,immunity2$REF,sep="_")
immunity2$alleleA <- paste(immunity2$SNP,immunity2$ALT,sep="_")

head(immunity2$allele)
head(immunity2$alleleA)

# get all frequencies
freq999<- tab(genind2genpop(albo2[loc=c(immunity2$allele)]),freq=TRUE)
head(freq999)

# rbind the dataframes
tomelt3<-cbind(freq999)

# rename
head(tomelt3)
head(immunity2$allele)
head(immunity2$alleleA)

# make new column with the loci names (R and A)
immunity2$lociNames<-lapply(1:nrow(immunity2), function(i) c(immunity2$allele[i], immunity2$alleleA[i]))
head(immunity2$lociNames)

# melt the column with the loci names
immunityNames<-melt(immunity2$lociNames)

colnames(tomelt3)<-immunityNames$value

head(tomelt3)


melted_dapc3 <- melt(tomelt3,id.var="id")
head(melted_dapc3)

# split the melted file to get each allele
meltedAlleles3<- strsplit(as.character(melted_dapc3$Var2), "_", fixed=TRUE)
head(meltedAlleles3)

# make it a dataframe and rename the column names
meltedAlleles4<- do.call(rbind.data.frame, meltedAlleles3)
head(meltedAlleles4)
colnames(meltedAlleles4)<- c("SNP", "Allele")
head(meltedAlleles4)

# merge with the melted_dapc
meltedAlleles5 <- cbind(melted_dapc3, meltedAlleles4)
head (meltedAlleles5)
colnames(meltedAlleles5)<- c("Population", "Var2", "Frequency", "SNP", "Allele")
head (meltedAlleles5)

# keep the columns needed
meltedAlleles6<- subset(meltedAlleles5, 
                        select=c("Population", "SNP", "Allele", 'Frequency'))
head(meltedAlleles6)



head(meltedAlleles6)
base <- ggplot(meltedAlleles6, aes(x=Population, y=Frequency, group=Allele, color = Allele, fill = Allele)) + #, group=Allele, colour=SNP
  geom_blank() + 
  xlab(NULL) + 
  ylab(NULL)

p3 <- base + facet_wrap(~SNP, ncol = 3,  strip.position="right") +  #geom_point() +
  #geom_label(aes(label = Allele), size = 4, alpha = 1, label.size = NA, fill = "lightgray")+
  #geom_label(aes(label = Allele), size = 4, alpha = .1, label.size = NA)+
  geom_line(aes(linetype = Allele), size = .6, show.legend = FALSE) +
  geom_label(aes(label=Allele), colour= "white", label.padding=unit(0.05,"lines"), size=5, alpha = 1) +
  geom_text(aes(label=Allele), size=4, color =  "black") +
  scale_fill_manual(values = c("A" = "#fccee9", "G" = "#a0d3ed","C" = "#ffcc74", "T" = "#b5f283")) +
  scale_linetype_manual(values = c("A" = "solid", "T" = "dashed", 
                                   "C" = "dotted", "G" = "twodash")) +
  scale_color_manual(values = c("A" = "#f29ce6", "G" = "#2626f0","C" = "#f09a19", "T" = "#397511")) +
  theme(legend.position="none") +  
  theme_bw() +
  theme( 
    strip.background = element_rect(fill="#e6e2d6"),
    strip.text = element_text(face="bold")
  ) + theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45)) 

# theme(strip.text.y = element_text(angle = 0)
p3
#
ggsave(filename = "immunity_loci.pdf", p3, width = 25, height = 20, units = "in")
#



  
#######################################################
# ggplot plot looks better but going to leave this code using the matplot
######################################################
par(mfrow=c(5,1), mar=c(1,2,3,1),las=3)  # from mar=c(5.1,4.1,4.1,.1) to par(mar=c(1,1,1,1))

#freq_3305
matplot(freq_3305, pch=c("T","C"), type="b",
        xlab="Population",ylab="allele frequency", xaxt="n",
        main="SNP AX-583301796")
# legend(x = c(0.567, 1.432), y = c(1.542015, 1.432), legend="", xpd=T, bg="grey")
# mtext("SNP AX-583301796", side=1, line=-1)
# freq_27472
matplot(freq_27472, pch=c("C","T"), type="b",
        xlab="Population",ylab="allele frequency", xaxt="n",
        main="SNP AX-582523725")
# freq_57371
matplot(freq_57371, pch=c("G","A"), type="b",
        xlab="Population",ylab="allele frequency", xaxt="n",
        main="SNP AX-581354204")
# freq_52296
matplot(freq_52296, pch=c("C","T"), type="b",
        xlab="Population",ylab="allele frequency", xaxt="n",
        main="SNP AX-580721505")
# freq_36169
matplot(freq_36169, pch=c("G","A"), type="b",
        xlab="Population",ylab="allele frequency", xaxt="n",
        main="SNP AX-579475215")

axis(side=1, at=1:39, lab=levels(albo2$pop))



# DAPC commands not usefull
######################################################

## label individuals at the periphery
scatter(dapc1, label.inds = list(air = 2, pch = NA))

## only ellipse, custom labels, use insets
scatter(dapc1, cell=2, pch="", cstar=0, posi.pca="topright", posi.da="bottomright", scree.pca=TRUE,
        inset.pca=c(.01,.3), label=paste(c("Central-Southeast","South", "North-Central", "Northwest","East")), axesel=FALSE, col=terrain.colors(10))

scatter(dapc1, col=transp(myPal(6)), scree.da=FALSE, cell=1.5, cex=2, bg="white",cstar=0,
        legend=TRUE)


scatter(dapc1, col=transp(myPal(6)), scree.da=T, cell=1.5, cex=2, bg="white",cstar=0, scree.pca = TRUE,
        legend=TRUE)
######################################################

# get some statistics
######################################################
#Interpreting group memberships
class(dapc1$posterior)
## [1] "matrix"
dim(dapc1$posterior)
## [1] 600   6
round(head(dapc1$posterior),3)

summary(dapc1)


# The slot assign.per.pop indicates the proportions of successful reassignment
assignplot(dapc1, subset=1:50)


# structure like plot 
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:10), lab="",
          ncol=1, xlab="individuals", col=funky(10))


#We can also have a closer look at a subset of individuals; for instance, for the first 50 individuals
compoplot(dapc1, subset=1:50, posi="bottomright", txt.leg=paste("Cluster", 1:5), lab="",
          ncol=2, xlab="individuals", col=funky(5))


# Obviously, we can use the power of R to lead our investigation further.
# For instance, which are the most ’admixed’ individuals? Let us consider as
# admixed individuals having no more than 90% of probability of membership in a single cluster:

# not working
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

compoplot(dapc1, subset=temp, posi="bottomright", txt.leg=paste("Cluster", 1:5),
          ncol=2, col=funky(5))


# When and why group memberships can be unreliable   5PCs
albo2
temp <- summary(dapc(albo2, n.da=100, n.pca=5))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual population",horiz=TRUE, las=1)


png("reassignment.png", width = 6, height = 10, units='in',res = 600)
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual population",horiz=TRUE, las=1)
dev.off()
##################################################################################

# Using the a-score

#dapc2 <- dapc(albo,n.pca=200,n.da=15)
dapc2 <- dapc(albo2, n.da=100, n.pca=200) 
temp <- a.score(dapc2)
names(temp)

temp$tab[1:5,1:5]
temp$pop.score
temp$mean

# The number of retained PCs can be chosen so as to optimize the a-score; this is achived by optim.a.score:
dapc2 <- dapc(albo2, n.da=100, n.pca=300)
temp <- optim.a.score(dapc2)


pdf("ideal_number_PCs.pdf", width = 9, height = 5)
temp <- optim.a.score(dapc2)
dev.off()


#  We perform the analysis with 20 PCs retained, and then map the membership probabilities as before:
dapc3 <- dapc(albo2, n.da=100, n.pca=51)
myCol <- rainbow(15)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc3, lab="", posi=list(x=30,y=-.01), cleg=.4)

# does not work with our dataset because low admixture
############################################################################################
# And as before, we can further investigate admixed individuals, 
# which we arbitrarily define as those having no more than 0.5 probability of membership to any group:

temp <- which(apply(dapc3$posterior,1, function(e) all(e<0.1)))
temp

lab <- pop(albo2)
par(mar=c(8,4,5,1), xpd=TRUE)
compoplot(dapc3, subset=temp, cleg=.6, posi=list(x=0,y=1.2), lab=lab)


# Using supplementary individuals
set.seed(2)
kept.id <- unlist(tapply(1:nInd(albo2), pop(albo2),
                         function(e) sample(e, 2, replace=FALSE)))
x.sup <- albo2[-kept.id]
x <- albo2[kept.id]
nInd(x)

# We perform the DAPC of x, and use predict to predict results for the supplementary individuals:
dapc4 <- dapc(x,n.pca=20,n.da=15)
pred.sup <- predict.dapc(dapc4, newdata=x.sup)
names(pred.sup)

head(pred.sup$assign)
pred.sup$ind.scores[1:5,1:3]
round(pred.sup$posterior[1:5, 1:5],3)


col <- rainbow(length(levels(pop(x))))
col.points <- transp(col[as.integer(pop(x))],.2)
scatter(dapc4, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, xlim=c(-10,10), legend=F)

par(xpd=TRUE)
points(dapc4$ind.coord[,1], dapc4$ind.coord[,2], pch=20, col=col.points, cex=5)
col.sup <- col[as.integer(pop(x.sup))]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=15,
       col=transp(col.sup,.7), cex=2)
add.scatter.eig(dapc4$eig,15,1,2, posi="topleft", inset=.02)
# Add legend to top right, outside plot region
legend("topleft", inset=c(3.3,0), legend=levels(pop(x)), pch=20, title="Population", col = col, cex=1)


mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))


table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))

#############################################################################################
# repeat the first steps but now with the ideal number of PCs
# We specify that we want to evaluate up to k = 20 
grp3 <- find.clusters(albo2, max.n.clust=15, n.pca = 51)


names(grp3)
head(grp3$Kstat, 10)
grp3$stat
head(grp3$grp, 10)
grp3$size

# number of mosquitoes per group
table(pop(albo),grp3$grp)

groups4<-table(pop(albo2),grp3$grp)
write.csv(groups4, file = "groups.csv")

head(groups4)

table.value(table(pop(albo), grp3$grp),
            col.lab=paste("Cluster", 1:10),
            clegend = 0,
            row.lab=paste(row.names(groups4)))#row.lab=paste("pop", 1:39))

# We run the analysis on the previous toy dataset, using the inferred groups stored in grp$grp:
dapc4 <- dapc(albo2, grp3$grp, n.pca=51, n.da = 5 )

dapc4

myCol <- c("#e88ad2",
           "#9e4725",
           "#ffb235",
           "#4eb5ff",
           "#e35235",
           "#c2db43",
           "#610ffa",
           "#5cd168",
           "#b2aae0",
           "#524d50",
           "#648a58",
           "#d6a080")

scatter(dapc4, posi.da="bottomleft",  bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topright")

# DAs plot
scatter(dapc4,1,1, col=myCol, bg="white", scree.da=FALSE, legend=FALSE, solid=.4)
legend("top", legend=c(1:10),
       fill=c(myCol), cex=0.8,box.lwd = 0,box.col = "transparent",bg = "transparent", ncol=5)


scatter(dapc4,2,2, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc4,3,3, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc4,4,4, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc4,5,5, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)



scatter(dapc4, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid =.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:10))

# 
# pdf("DAPC_final.pdf", width = 8, height = 6)
scatter(dapc4, col=myCol,
        bg="white", scree.pca=TRUE, posi.da="bottomleft",
        posi.pca="topright",
        cell=1.5, pch=20, cex=2,cstar=0)#, leg=TRUE, txt.leg=paste("Cluster",1:10))
# dev.off()

# pdf("DAPC_final.pdf", width = 8, height = 6)
scatter(dapc4, scree.pca=TRUE, posi.da="bottomleft",posi.pca="topleft",
        bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid =.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:10))
# dev.off()


# make it ggplot
# https://grunwaldlab.github.io/Population_Genetics_in_R/clustering_plot.html
my_df <- as.data.frame(dapc4[[length(dapc4)]]$ind.coord)
my_df <- as.data.frame(dapc4$ind.coord)
head(my_df)
my_df$Group<- dapc4$assign
head(my_df)

# sort by row idex
xx<- my_df[order(as.numeric(rownames(my_df))),,drop=FALSE]
head(xx)

# check the name of each individual on albo2$pop
xxx <- as.data.frame(albo2$pop)
head(xxx)

# merge by row index
xy <- cbind(xx,xxx)
head(xy)
colnames(xy)<-c("LD1", "LD2", "LD3", "LD4", "LD5", "Group", "ID1")
head(xy)

# now import matadata to make shape by continet or range
yy2<- read.table("metadata.txt", header = TRUE, sep = "\t")
head(yy2)

# merge with xy by ID1
xy2<- merge(xy, yy2, by = "ID1")
head(xy2)

# make plot with ggplot
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_classic()
p2 <- p2 + scale_color_manual(values=c(myCol))
p2 <- p2 + scale_fill_manual(values=c(paste(myCol, "66", sep = "")))
p2
#
ggsave(filename = "PCA_scatter_ggplot.pdf", p2, width = 8, height = 6, units = "in")



#to see all shapes -> plot shapes - para escolher os simbolos
N = 100; M = 1000
good.shapes = c(1:25,33:127)
foo = data.frame( x = rnorm(M), y = rnorm(M), s = factor( sample(1:N, M, replace = TRUE) ) )
ggplot(aes(x,y,shape=s ), data=foo ) +
  scale_shape_manual(values=good.shapes[1:N]) +
  geom_point()

# make plot reflecting the continents using different shapes
head(xy2)
# make plot with ggplot, change the linear discriminant (LD) to see the clustering changes
p2 <- ggplot(xy2, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = Continent ))
#p2 <- ggplot(xy2, aes(x = LD2, y = LD3, color = Group, fill = Group, shape = Continent ))
p2 <- p2 + geom_point(size = 3)
p2 <- p2 + scale_shape_manual(values = c(21:24))
p2 <- p2 + theme_classic()
p2 <- p2 + scale_color_manual(values=c(myCol))
p2 <- p2 + scale_fill_manual(values=c(paste(myCol, "66", sep = "")))
p2

ggsave(filename = "PCA_scatter_ggplot_by_continent_LD1_LD2.pdf", p2, width = 8, height = 6, units = "in")
ggsave(filename = "PCA_scatter_ggplot_by_continent_LD2_LD3.pdf", p2, width = 8, height = 6, units = "in")

# I will drop the 3 individuals from NEG that are very different, might help visualize the points better,
# maybe I have to do it since before running DAPC
head(xy2)
xy3<- xy2 %>%
  filter(LD1 < 60)

# plot it again to see how it changes, change the linear discriminant (LD) to see the clustering changes
p2 <- ggplot(xy3, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = Continent ))
#p2 <- ggplot(xy2, aes(x = LD2, y = LD3, color = Group, fill = Group, shape = Continent ))
p2 <- p2 + geom_point(size = 3)
p2 <- p2 + scale_shape_manual(values = c(21:24))
p2 <- p2 + theme_classic()
p2 <- p2 + scale_color_manual(values=c(myCol))
p2 <- p2 + scale_fill_manual(values=c(paste(myCol, "66", sep = "")))
p2

ggsave(filename = "PCA_scatter_ggplot_by_continent_LD1_LD2_dropping_outliers.pdf", p2, width = 8, height = 6, units = "in")


# make it 3D
#install.packages("threejs")
library(threejs)
scatterplot3js(as.matrix(my_df[,1:3]),col=myCol[my_df$Group],size=0.3)

library(plotly)
p <- plot_ly(my_df, x=~LD1, y=~LD2, 
             z=~LD3, color=~Group) %>%
  add_markers(size=1.5)
print(p)

#https://stackoverflow.com/questions/45052188/how-to-plot-3d-scatter-diagram-using-ggplot
#devtools::install_github("AckerDWM/gg3D")
library("gg3D")

## An empty plot with 3 axes
qplot(x=0, y=0, z=0, geom="blank") + 
  theme_void() +
  axes_3D()

## Axes can be populated with points using the function stat_3D.
ggplot(my_df, aes(x=LD1, y=LD2, z=LD3, color=Group)) +
  scale_color_manual(values=c(myCol)) +
  # scale_fill_manual(values=c(myCol))+
  theme_void() +
  axes_3D() +
  stat_3D()
#
ggsave(filename = "3D_PCA.pdf", p, width = 8, height = 5, units = "in")


