library("rJava")
library(rcdk)
library(readxl)
library(caret)
library(dplyr)
library(ggplot2)
library(latex2exp)


setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD")


data <- read.csv('DATABASES/SolCbio3Database_merged.csv')
descs <- read.csv('DATABASES/descriptors_SolCbio3Database_merged.csv')
descs <- descs[,-1]




p <- ggplot(data) + 
  geom_histogram(aes(x=diff_logS),colour='black',fill='blue4',binwidth = 0.1) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) + 
  scale_x_continuous(breaks = round(seq(min(data$diff_logS),max(data$diff_logS),by=1),1)) + 
  xlab(TeX("$\\Delta \\log{S}$"))
ggsave(p,file='PLOTS/threshold_diff/diff_his.pdf',width=270,height = 170,units = 'mm')
ggsave(p,file='PLOTS/threshold_diff/diff_his.png',width=270,height = 170,units = 'mm')


##CALCULAR CONDICIONALES CON VALORES DE diff ENTRE 0 y 1
for (i in 1:20) {
new <- ifelse(data$diff_logS > i/10,0,1)  
data[,ncol(data)+1] <- new
colnames(data)[ncol(data)] <- paste0('cond', i/10)
}




#reduction step#
descs <- descs[, !apply(descs, 2, function(x) any(is.na(x)) )] #remove NAs
descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )] #remove constant columns



#remove correlated descriptors
r2 <- which(cor(descs)^2 > .6, arr.ind=TRUE)
r2 <- r2[ r2[,1] > r2[,2] , ]
descs <- descs[, -unique(r2[,2])]



for (i in 0:19) {
  new <- data[,13 + i]
  descs[,ncol(descs)+1] <- new
  colnames(descs)[ncol(descs)] <- paste0('cond',i/10)
}






pval <- data.frame(matrix(ncol=20,nrow=225))
cols <- c()
for(n in 1:20){
  new <- c()
  for(i in 1:225){
    new <- c(new,t.test(ifelse(descs[,225+n]==1,descs[,i],NA),ifelse(descs[,225+n]==0,descs[,i],NA))[[3]])
  }
  pval[,n] <- new
  cols <- c(cols,paste0('pvalue',n/10))
}
colnames(pval) <- cols



pvalues <- data.frame(
  'diff' = c(1:20)/10
)


mols <- c()
nmols <- c()

for(i in 1:20){
  mols <- c(mols,nrow(pval[pval[,i]<0.05,]))
  nmols <- c(nmols,nrow(descs[descs[,225+i] == 0,]))
}
pvalues$descs <- mols
pvalues$numMols <- nmols
pvalues$typeN <- 'Number of entries'





p <- ggplot(pvalues,aes(x=diff)) + 
  geom_col(aes(y=numMols), fill = 'blue4',colour='black',alpha=0.5,width = 0.05) + 
  geom_line(aes(y=descs*15),linetype=2,size=1) + 
  geom_point(aes(y=descs*15),size=6,shape=21, fill = 'gray50') + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = c(0.8,0.4)) +
  scale_y_continuous(name='Number of entries',sec.axis=sec_axis(trans = ~./15,name='Descriptors with p < 0.05')) + 
  scale_x_continuous(name = TeX("$\\Delta \\log{S}$ threshold"), breaks = c(0:20)/5)

ggsave(p,file='PLOTS/threshold_diff/threshold.pdf',width=200,height = 130,units = 'mm')
ggsave(p,file='PLOTS/threshold_diff/threshold.png',width=200,height = 130,units = 'mm')






