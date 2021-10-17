library(ggplot2)
library(ggcorrplot)
#install.packages("dplyr")
#library(dplyr) ## use group.by function
#library(data.table)
library(scales) ##needed for percentin barplot
#install.packages('GGally')
library(GGally)
#install.packages('stringi')
library(moments)
library(nortest)
## setwd

setwd("~/MSc_project_R/Webserver_1-69_results_compiled")

## Read in all files
fnames <- list.files(path = getwd(), pattern= "*_results.csv", full.name = TRUE, recursive = FALSE)
fnames

## All csv into one dataframe
df_comb <- NULL
for (aa in fnames)
{df <- read.csv(file=aa,header=TRUE)
  #df<- df$Mutant
  pdb <- strsplit(aa, "_comp_results.csv") ##file name up to underscore
  pdb <- substr(pdb, nchar(pdb) - 3, nchar(pdb)) ## last 4 char of pdb string
  pdb_list = replicate(nrow(df), pdb) ## replicate pdb according to rows in csv
  df$pdb <- pdb_list
  if (is.null(df_comb)) {
      df_comb <- df
      
  } else
  { df_comb <- rbind(df_comb, df)}
  
  
}

#View(df_comb)
#View(df_prop)
nrow(df_comb)

## A function to create data frames to map properties to AAs (attach properties to Mutant AAs)
create_df_prop <- function(groups, names, header) {
  for (i in 1:length(names)) {
    prop_list <- replicate(nchar(groups[i]), names[i])## create vector of repeated prop names
    temp_df  <- cbind(unlist(strsplit(groups[i],"")), prop_list) ## join to vector of AAs in temp_df
    colnames(temp_df) <- c("AA", header)
    if (i==1) { ## joing all props into one df
      return_df <- temp_df
    } else {
      return_df <- rbind(return_df, temp_df)
    }
  }
  #prop_list
  #return_df
  return(return_df) ## returned into df_prop
}



##function to create charts
create_plots <- function(df_data, df_groups, title){
  ## Join the two dataframes on Mutant(df_comb) and AA(df_prop)
  df_result <- merge(x=df_data, y = df_groups, by.x = "Mutant", by.y = "AA")
  df_result <- df_result[order(df_result$pdb),]

  
  ## Barplot based on properties.
  jpeg(paste(title, "_count.jpg"))
  ## Barplot with y as counts
  plot <- ggplot(df_result,aes(pdb)) + geom_bar(aes_string(fill = title))+scale_fill_brewer(palette="Set1") + theme(legend.title = element_blank())
  print(plot)
  dev.off()
  jpeg(paste(title, "_percent.jpg"))
  ## Barplot with y as percentage
  plot<- ggplot(df_result,aes(pdb)) + geom_bar(aes_string(fill = title),position="fill")+scale_fill_brewer(palette="Set1")+scale_y_continuous(labels = percent,"Percentage") + theme(legend.title = element_blank())
  print(plot)
  dev.off()
  return (df_result) #returns the dataframe with the mapped properties
}

###Part 1 - features of mutant amino acids

##a) Chemical property

## Chemical property groups
## Assigning chemical property to amino acid

aliphatic <-"AGILPV"
aromatic <-"FWY"
sulphur <-"CM"
hydroxyl <-"ST"
basic <-"RHK"
acidic <-"DE"
amide <-"NQ"
head = "chem_props"
groups <- c(aliphatic, aromatic, sulphur, hydroxyl, basic, acidic, amide)
group_names <- c("aliphatic", "aromatic", "sulphur", "hydroxyl", "basic", "acidic", "amide")
df_prop<- create_df_prop(groups, group_names, head)
#View(df_prop)
df_all <- df_comb
df_all <- create_plots(df_all, df_prop, head)
#View(df_all)

##b) Size of aa
small <-"GASCDPNT"
medium <-"QEHV"
large <-"RILKMFWY"
head = "size"
groups <- c(small, medium, large)
group_names <- c("small", "medium", "large")
df_prop<- create_df_prop(groups, group_names, head)
df_all <- create_plots(df_all, df_prop, head )


##c) Polarity
polar <-"RNDQEHKSTY"
non_polar <-"ACGILMFPWV"
head = "polarity"
groups <- c(polar, non_polar)
group_names <- c("polar", "non-polar")
df_prop<- create_df_prop(groups, group_names, head)
df_all <- create_plots(df_all, df_prop, head )

##d) Hydrogen bond
donor <- "RKW"
acceptor <- "DE"
donor_and_acceptor <-"NQHSTY"
none <- "ACGILMFPV"
head = "h_bond"
groups <- c(donor, acceptor, donor_and_acceptor,none)
group_names <- c("donor", "acceptor","donor_and_acceptor","none")
df_prop<- create_df_prop(groups, group_names, head)
df_all <- create_plots(df_all, df_prop, head )

##e) Hydrophobicity
hydrophobic <-"ACILMFWV"
neutral <-"GHPSTY"
hydrophilic <-"RNDQEK"
head = "hydrophobicity"
groups <- c(hydrophobic, neutral, hydrophilic)
group_names <- c("hydrophobic", "neutral", "hydrophilic")
df_prop<- create_df_prop(groups, group_names, head)
df_all <- create_plots(df_all, df_prop, head )
View(df_all)
#nrow(df_all) #check this tallies with df_prop #134
#write.csv(df_all, "df_all.csv")
  
###########################################################################################################
### Part 2) - Distance to interface v change in Gibbs free energy


### DDG from all webservers pooled. Not using GGplot as DDG awkward to label.
x <- df_all$mcsm.ppi2.prediction.distance.to.interface
y1 <- df_all$mcsm.ppi2.prediction.ddG..kCal.mol.##red
y2 <- -(df_all$MutaBind2.DDG) ##Blue ## Reverse sign
y3 <- -(df_all$SAAMBE3D..ddG..kCal.mol.) ##Green ## Reverse sign
data_ggp <-data.frame(x=x,y=c(y1,y2,y3))
par(mfrow=c(1,1))
plot(data_ggp,xlab=expression("Minimum distance to the interface"(ring(A))),ylab="\u0394\u0394G",col = c(y1="orange",y2="blue",y3="green"), ylim=c(-2,2),pch=16)
legend("topright", c(y1="mCSM-PPI2",y2="MutaBind2",y3="SAAMBE3D"), col = c(y1="orange",y2="blue",y3="green"),pch=16,lty = 1,bty = "n",cex=0.7)
abline(h=mean(data_ggp$y),col="red",lwd=3)

### DDG for each of the webservers
par(mfrow=c(2,2),font=2,ps=10,family="sans")
##mCSM-PPI2
plot(x=df_all$mcsm.ppi2.prediction.distance.to.interface,y=df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,xlab=expression("Minimum distance to the interface"(ring(A))),ylab="\u0394\u0394G",col="orange",main="mCSM-PPI2",ylim=c(-2,2),pch=16,cex.lab=1)
abline(h=mean(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.),col="red",lwd=2)
##MutaBind2 - using distance to interface data from mCSM-PPI2 output (reverse sign)
plot(x=df_all$mcsm.ppi2.prediction.distance.to.interface,y=-(df_all$MutaBind2.DDG),xlab=expression("Minimum distance to the interface"(ring(A))),ylab="\u0394\u0394G",col="blue",main="MutaBind2",ylim=c(-2,2),pch=16,cex.lab=1)
abline(h=mean(-df_all$MutaBind2.DDG),col="red",lwd=2)
#SAAMBE3D - using distance to interface data from mCSM-PPI2 output (reverse sign)
plot(x=df_all$mcsm.ppi2.prediction.distance.to.interface,y=-(df_all$SAAMBE3D..ddG..kCal.mol.),xlab=expression("Minimum distance to the interface"(ring(A))),ylab="\u0394\u0394G",col="green",main="SAAMBE3D",ylim=c(-2,2),pch=16,cex.lab=1)
abline(h=mean(-df_all$SAAMBE3D..ddG..kCal.mol.),col="red",lwd=2)


##################################
## Part 3) Correlation plots: Gibbs free energy: comparing output from the 3 webservers

### Is there a correlation between DDG output by the 3 webservers?

## reverse sign for mutabind2 and saambe3d
par(mfrow=c(2,1),font=2,ps=10,family="sans")

## GGPLOT needs df
DDG_corr_comp <- data.frame(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,-(df_all$MutaBind2.DDG),-(df_all$SAAMBE3D..ddG..kCal.mol.))
#View(DDG_corr_comp)
pairs(DDG_corr_comp,labels = c("mCSM_PPI2","MutaBind2","SAAMBE3D"),pch=16)
require(GGally)
ggpairs(DDG_corr_comp,columnLabels = c("mCSM-PPI2","MutaBind2","SAAMBE3D"),upper = list(continuous = wrap("cor", method = "spearman")))

### Ranked correlation plot
DDG_corr_comp_rank <- data.frame(rank(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.),rank(-(df_all$MutaBind2.DDG)),rank(-(df_all$SAAMBE3D..ddG..kCal.mol.)))
pairs(DDG_corr_comp_rank,labels = c("mCSM_PPI2","MutaBind2","SAAMBE3D"),pch=16)
ggpairs(DDG_corr_comp_rank,columnLabels = c("mCSM-PPI2","MutaBind2","SAAMBE3D"),upper = list(continuous = wrap("cor", method = "spearman")))



#####A ggcorplot
corr <-round(cor(DDG_corr_comp,method ="spearman"),2)
colnames(corr) <- c("mCSM-PPI2","MutaBind2","SAAMBE3D")
rownames(corr) <- c("mCSM-PPI2","MutaBind2","SAAMBE3D")
corr
p.mat <-cor_pmat(DDG_corr_comp)
ggcorrplot(corr,method="square",lab=TRUE,lab_size = 3)


##############
##Correlation tests

cor.test(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,-df_all$MutaBind2.DDG,method=c('kendall'))
cor.test(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,-df_all$SAAMBE3D..ddG..kCal.mol.,method=c('kendall'))
cor.test(-df_all$MutaBind2.DDG,-df_all$SAAMBE3D..ddG..kCal.mol.,method=c('kendall'))

cor.test(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,-df_all$MutaBind2.DDG)
cor.test(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.,-df_all$SAAMBE3D..ddG..kCal.mol.)
cor.test(-df_all$MutaBind2.DDG,-df_all$SAAMBE3D..ddG..kCal.mol.)
### Display results as a table underneath




###Plotting distribution of change in free energy

par(mfrow=c(2,2),font=2,ps=10,family="sans")
my.breaks <- seq(-2, 2, by=0.5)

hist(df_all$mcsm.ppi2.prediction.ddG..kCal.mol., breaks=my.breaks,xlab = "\u0394\u0394G",main = "mCSMPPI2",col="orange")
## MutaBind2 and SAAMBE3D - decreasing affinity >0, so reverse sign
hist(-(df_all$MutaBind2.DDG), breaks=my.breaks,xlab = "\u0394\u0394G",main = "MutaBind2",col="blue")
hist(-(df_all$SAAMBE3D..ddG..kCal.mol.), breaks=my.breaks, xlab = "\u0394\u0394G",main = "SAAMBE3D",col="green")

## empirical cumulative dist function
plot(ecdf(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.),main = "mCSMPPI2",verticals=TRUE)##heavy tail
plot(ecdf(df_all$MutaBind2.DDG), main = "MutaBind2",verticals=TRUE) ## blocky
plot(ecdf(df_all$SAAMBE3D..ddG..kCal.mol.),main = "SAAMBE3D",verticals=TRUE)

##May need different threshold for each
sd(df_all$mcsm.ppi2.prediction.ddG..kCal.mol.) ##[1] 0.4244935
sd(df_all$MutaBind2.DDG) #[1] 0.666546
sd(df_all$SAAMBE3D..ddG..kCal.mol.)##[1] 0.4480055





##Tests on Gibbs Free Enery Distributions
#  THe servers are giving different answers. Check for statistical differences in output.

#1.  test for normality
mcsm <- df_all$mcsm.ppi2.prediction.ddG..kCal.mol.
mutabind <- -df_all$MutaBind2.DDG
saambe <- -df_all$SAAMBE3D..ddG..kCal.mol.
shapiro.test(mcsm) #not normal - probably outliers
shapiro.test(mutabind) #fat not normal at 5% level
shapiro.test(saambe) #p-value = 0.01594,not normal, negative skew


##non-parametric test ##nortest and moments
ad.test(mcsm) #not normal
ad.test(mutabind) #not normal
ad.test(saambe) ##p=0.06528. Just about,accept null, data from normal distribution
#p-values are bigger than shapiro - less powerful test

#calculate variances
var(mcsm)
var(mutabind)
var(saambe)
#looks like mutabind variance is different.  Dependent samples and not normal so
#no obvious formal test.  But most reasonable to conclude different variances.

#means and medians
mean(mcsm)
mean(mutabind)
mean(saambe)
median(mcsm)
median(mutabind)
median(saambe)

#Visually quite different mean from Saambe but can run tests from the others
#paired t test requires difference to be normal
shapiro.test(mcsm - saambe) #p-value = 0.05586
shapiro.test(mutabind - saambe) #p-value = 0.122
shapiro.test(mutabind - mcsm) #p-value = 0.2591
#none great so care with t test. Use Wilcoxon signed ranked test instead.
#but both done here
t.test(mcsm, mutabind, paired = TRUE)
wilcox.test(mcsm, mutabind, paired = TRUE)
#both tests show different means for MCSM and mutabind
t.test(mcsm, saambe, paired = TRUE)
wilcox.test(mcsm, saambe, paired = TRUE)
t.test(saambe, mutabind, paired = TRUE)
wilcox.test(saambe, mutabind, paired = TRUE)
#obvious from histograms that saambe significantly different

##Look at correlations in more detail - 
##do the results get more in line with each other for the more extreme changes in Gibbs Free Energy? 

df <- df_all[, c("pdb","Position","Wild.Type","Mutant","mcsm.ppi2.prediction.distance.to.interface","Interface.","chem_props","size","polarity", "h_bond","hydrophobicity","mcsm.ppi2.prediction.ddG..kCal.mol.", "MutaBind2.DDG", "SAAMBE3D..ddG..kCal.mol.")]
names(df)[names(df)=="mcsm.ppi2.prediction.distance.to.interface"] <- "Distance_to_interface"
names(df)[names(df)=="Interface."] <- "MutaBind_Interface_Ind"
names(df)[names(df)=="mcsm.ppi2.prediction.ddG..kCal.mol."] <- "MCSM"
names(df)[names(df)=="MutaBind2.DDG"] <- "MutaBind"
names(df)[names(df)=="SAAMBE3D..ddG..kCal.mol."] <- "SAAMBE"
df$SAAMBE <- -df$SAAMBE
df$MutaBind <- -df$MutaBind
#n = nrow(df)
df$MCSM_rank <- rank(df$MCSM) 
df$MutaBind_rank <- rank(df$MutaBind) 
df$SAAMBE_rank <- rank(df$SAAMBE) 
#mean(df$MutaBind_PC[df$MCSM_PC>.75])
View(df) 
#write.csv(df, file="final_list_of_interest.csv") ### Final table in excel -thesis- total=119 mutations 


conditional_mean_lt <- function(rank, vecAve, vecCond) {  
  #vecAve - the vector of DDGs that you are taking average of
  #vecCond = the vector of DDGs that decides which mutations to average over
  #rank = this is the threshold i.e. the function takes the average of all the ranks in vecAve, where the rank in vecCond is lower than this "rank" 
  ans <- mean(vecAve[vecCond <= rank])
  return(ans)
}

#conditional_mean_gt <- function(rank, vecAve, vecCond) {
# same as conditional_mean_lt but with "greater than" to look at increasing affinity - positive DDG
#  ans <- mean(vecAve[vecCond >= rank])
#  return(ans)
#}


#df$Muta_MCSM_gt <- lapply(df$MCSM_rank,conditional_mean_gt, vecAve = df$MutaBind_rank  , vecCond = df$MCSM_rank)
#df$SAAMBE_MCSM_gt <- lapply(df$MCSM_rank,conditional_mean_gt, vecAve = df$SAAMBE_rank  , vecCond = df$MCSM_rank)
#df$MCSM_Muta_gt <- lapply(df$MutaBind_rank,conditional_mean_gt, vecAve = df$MCSM_rank  , vecCond = df$MutaBind_rank)
#plot(df$MCSM_rank, df$Muta_MCSM_gt, ylim = c(10,67), col = 4)
#points(df$MCSM_rank, df$SAAMBE_MCSM_gt, col = 2)
#legend("topright", c("MutaBind", "SAAMBE"), lty = 1, col = c(4,2))

#Look at the negative changes
par(mfrow=c(1,1),font=2,ps=10,family="sans")
df$Muta_MCSM_lt <- lapply(df$MCSM_rank,conditional_mean_lt, vecAve = df$MutaBind_rank  , vecCond = df$MCSM_rank)
df$SAAMBE_MCSM_lt <- lapply(df$MCSM_rank,conditional_mean_lt, vecAve = df$SAAMBE_rank  , vecCond = df$MCSM_rank)
df$MCSM_Muta_lt <- lapply(df$MutaBind_rank,conditional_mean_lt, vecAve = df$MCSM_rank  , vecCond = df$MutaBind_rank)
df$SAAMBE_Muta_lt <- lapply(df$MutaBind_rank,conditional_mean_lt, vecAve = df$SAAMBE_rank  , vecCond = df$MutaBind_rank)
df$MCSM_SAAMBE_lt <- lapply(df$SAAMBE_rank,conditional_mean_lt, vecAve = df$MCSM_rank  , vecCond = df$SAAMBE_rank)
df$Muta_SAAMBE_lt <- lapply(df$SAAMBE_rank,conditional_mean_lt, vecAve = df$MutaBind_rank  , vecCond = df$SAAMBE_rank)
plot(df$MCSM_rank, df$Muta_MCSM_lt, xlab="mCSM-PPI2 ranked",ylab="Conditional average rank",ylim= c(0,134), col = 4)
points(df$MCSM_rank, df$SAAMBE_MCSM_lt, col = 2)
legend("topright", c("MutaBind", "SAAMBE"), col = c(4,2),pch = 1,bty = "n")

#plot(df$MutaBind_rank, df$MCSM_Muta_gt, ylim = c(0,67), col = 4)
#legend("topright", c("MCSM"), lty = 1, col = 4)

plot(df$MutaBind_rank, df$MCSM_Muta_lt, xlab="MutaBind2 ranked",ylab="Conditional average rank",ylim = c(0,134), col = 4)
points(df$MutaBind_rank, df$SAAMBE_Muta_lt, col = 2)
legend("topright", c("MCSM", "SAAMBE"), pch = 1, col = c(4,2),bty = "n")


plot(df$SAAMBE_rank, df$MCSM_SAAMBE_lt, xlab="SAAMBE ranked",ylab="Conditional average rank",ylim = c(0,134), col = 4)
points(df$SAAMBE_rank, df$Muta_SAAMBE_lt, col = 2)
legend("topright", c("MCSM", "Muta"), pch = 1, col = c(4,2),bty = "n")

