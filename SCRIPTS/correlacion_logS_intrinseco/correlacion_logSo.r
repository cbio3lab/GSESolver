library(readxl)
library(dplyr)
library(latex2exp) #enable LaTeX notation 


setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD")


data <- read.csv('DATABASES/SolCbio3Database_merged.csv')



data$dev <- ifelse(data$diff_logS>1,'YES','NO')

data$num <- NA
for(i in 1:nrow(data)){
  if(data$dev[i]=='YES'){
    data$num[i] <- as.character(rownames(data)[i])
  }
}

#resultados de la GSE
library(Metrics)
library(ggplot2)


linear <- lm(data$logS~data$GSE)


rmse <- round(rmse(data$logS, data$GSE),2)
mae <- round(mae(data$logS,data$GSE),2)
mse<-round(sum((data$logS-data$GSE)/length(data$logS)),2)
r2 <- round(cor(data$logS,data$GSE)^2,2)

#correlation plot
p <- ggplot(data) +
  #geom_point(aes(y=logS-0.05,x=GSE-0.05),colour="darkgray",size=4,alpha=0.7) + 
  geom_point(aes(y=logS,x=GSE,colour=dev),size=2,alpha=0.7) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = c(0.7,0.1)) + 
  scale_color_manual(values=c('blue4','gray50'),name='|logS (exp) - logS (GSE)| > 1?') +
  geom_abline(slope = coef(linear)[2],intercept = coef(linear)[1]) +
  geom_abline(slope = 1,intercept = 0,linetype='dashed') +
  ylab('logS Experimental') +
  xlab('logS GSE') +
  coord_cartesian(xlim = c(-15,2), ylim = c(-15,2)) +
  annotate('text',x=-14,y=-5,label=paste("R^2 == ",r2),parse=TRUE) + 
  annotate('text',x=-14,y=-5.8,label=paste('RMSE =',rmse)) + 
  annotate('text',x=-14,y=-6.6,label=paste('MAE =',mae)) + 
  annotate('text',x=-14,y=-7.4,label=paste('MSE =',mse))

ggsave(p, filename = 'PLOTS/correlacion_logSo/logSovsGSE.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/logSovsGSE.pdf', width = 7, height = 6, units = 'in')



dev_count <- data %>% count(dev)
p <- ggplot(data) +
  geom_bar(aes(y = dev, fill = dev), color = 'black') +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = 'none') + 
  annotate('text',x=dev_count[2,2]/2,y=2,label=paste("n =",dev_count[2,2]), size = 5, color = 'white', fontface = 'bold') + 
  annotate('text',x=dev_count[1,2]/2,y=1,label=paste("n =",dev_count[1,2]), size = 5 ,color = 'white', fontface = 'bold') + 
   scale_fill_manual(values=c('blue4','gray50')) +
  ylab(TeX("$\\Delta \\log{S} > 1$?")) +
  xlab('')

ggsave(p, filename = 'PLOTS/correlacion_logSo/barplot_logS.png', width = 3, height = 2, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/barplot_logS.pdf', width = 3, height = 2, units = 'in')


#residuals plot de los outliers
outliers <- data %>% filter(dev == 'YES')

outliers$res <- outliers$GSE - outliers$logS


p <- ggplot(outliers) + 
  geom_segment(aes(x = GSE, y = res, xend = GSE, yend = 0), color = 'gray50', linewidth = 0.1) +
  geom_point(aes(x = GSE, y = res), size = 2, color = 'gray50', alpha = 0.8) + 
  geom_hline(yintercept = 0, linewidth = 1) +
  #geom_text(aes(x = GSE, y = res, label = num), nudge_y = 0.25, size = 3, alpha = 0.5) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = 'logS GSE', y = 'Residuals (logS units)') + 
  coord_cartesian(xlim = c(-12.27,2.6), ylim = c(-4.4,7.3))

ggsave(p, filename = 'PLOTS/correlacion_logSo/residuals_outliers.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/residuals_outliers.pdf', width = 7, height = 6, units = 'in')






########### CHEMICAL SPACE ##################


descs <- read.csv('DATABASES/descriptors_SolCbio3Database_merged.csv')
descs$dev <- data$dev
descs$logS <- data$logS

p1 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=MW,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=MW,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y='MW (g/mol)')

p2 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=HBA2_Obabel,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=HBA2_Obabel,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y='HBA count')


p3 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=HBD_Obabel,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=HBD_Obabel,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y='HBD count')

p4 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=ALogP,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=ALogP,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y='clogP')


p5 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=TPSA_rdkit,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=TPSA_rdkit,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y=bquote('TPSA'~(ring(A)^2)),parse=TRUE)

p6 <- ggplot(descs) + 
  geom_violin(aes(x=dev,y=logS,fill=dev),alpha=0.7) + 
  geom_boxplot(aes(x=dev,y=logS,fill=dev),alpha=0.5,size=1, width=0.1) + 
  scale_fill_manual(values=c("blue4", "gray50")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=TeX("Is $\\Delta \\log{S} > 1$?"),y='logS',parse=TRUE)


violins <- list(p1,p2,p3,p4,p5,p6)
ncol <- 6
nrow <- 1

library(gridExtra) #create a grid plot

#grid plot creation
grid_violin <- grid.arrange(do.call(arrangeGrob, c(violins,ncol = ncol, nrow = nrow)))
ggsave('PLOTS/correlacion_logSo/chemspace_violins.pdf', plot = grid_violin, width = 11, height = 3.5, units = 'in')
ggsave('PLOTS/correlacion_logSo/chemspace_violins.png', plot = grid_violin, width = 11, height = 3.5, units = 'in')





#### PCA OF THE CHEMICAL SPACE #######

library(ggrepel)


####### PCA with Lipinski descriptors
pca_descs <- data.frame('MW' = descs$MW,
                        'HBD' = descs$HBD_Obabel,
                        'HBA' = descs$HBA2_Obabel,
                        'TPSA' = descs$TPSA_rdkit,
                        'logS' = descs$logS)

pca_descs$dev <- data$dev



  # Perform PCA
  pca_result <- prcomp(pca_descs[,-ncol(pca_descs)], scale = TRUE)
  
  # Extract scores and loadings
  scores <- as.data.frame(pca_result$x)
  loadings <- as.data.frame(pca_result$rotation)
  
  # Calculate variance explained
  variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)
  
  # Add grouping variable if provided

    scores$Group <- data$dev
  
  

  
  # Add biplot vectors 
    # Scale loadings for better visualization
    loadings_scaled <- loadings
    loadings_scaled$PC1 <- loadings_scaled$PC1 * 0.8 * max(abs(scores$PC1))
    loadings_scaled$PC2 <- loadings_scaled$PC2 * 0.8 * max(abs(scores$PC2))
    loadings_scaled$Variable <- rownames(loadings_scaled)
    

   #create plot 
    p <- ggplot() +
      theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 16),
            legend.title = element_text(size=16),
            legend.text = element_text(size=14),
            legend.position = c(0.2,0.8)) +
      labs(x = paste0("PC1 (", variance_explained[1], "%)"),
           y = paste0("PC2 (", variance_explained[2], "%)")) + 
      geom_point(data = scores, 
                   aes(x = PC1, y = PC2, color = Group),
                   alpha = 0.7, size = 4) + 
      scale_color_manual(TeX("Is $\\Delta \\log{S} > 1$?"),values = c('blue4', 'gray50')) + 
      geom_segment(data = loadings_scaled,
                   aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "black", size = 1) +
      geom_text_repel(data = loadings_scaled,
                      aes(x = PC1, y = PC2, label = Variable),
                      color = "black", size = 4,
                      max.overlaps = Inf)
    
    
    
    
ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_lipinski.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_lipinski.pdf', width = 7, height = 6, units = 'in') 



          ####PCA with all descriptors    
descs <- read.csv('DATABASES/descriptors_SolCbio3Database_merged.csv') #reload the descs dataframe
descs <- descs[,-1]


#reduction step#
descs <- descs[, !apply(descs, 2, function(x) any(is.na(x)) )] #remove NAs
descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )] #remove constant columns



#remove correlated descriptors
r2 <- which(cor(descs)^2 > .6, arr.ind=TRUE)
r2 <- r2[ r2[,1] > r2[,2] , ]
descs <- descs[, -unique(r2[,2])]

#add deviation variable
descs$dev <- data$dev

# Perform PCA
pca_result <- prcomp(descs[,-ncol(descs)], scale = TRUE)

# Extract scores and loadings
scores <- as.data.frame(pca_result$x)
loadings <- as.data.frame(pca_result$rotation)

# Calculate variance explained
variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Add grouping variable if provided

scores$Group <- data$dev




# Add top contributors biplot vectors 
loadings$contribution <- loadings$PC1^2 + loadings$PC2^2
top_vars <- loadings %>%
  arrange(desc(contribution)) %>%
  head(5)

# Scale loadings for better visualization
top_vars_scaled <- top_vars
top_vars_scaled$PC1 <- top_vars_scaled$PC1 * 2 * max(abs(scores$PC1))
top_vars_scaled$PC2 <- top_vars_scaled$PC2 * 2 * max(abs(scores$PC2))
top_vars_scaled$Variable <- rownames(top_vars_scaled)


#create plot 
p <- ggplot() +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = c(0.2,0.8)) +
  labs(x = paste0("PC1 (", variance_explained[1], "%)"),
       y = paste0("PC2 (", variance_explained[2], "%)")) + 
  geom_point(data = scores, 
             aes(x = PC1, y = PC2, color = Group),
             alpha = 0.7, size = 4) + 
  scale_color_manual(TeX("Is $\\Delta \\log{S} > 1$?"),values = c('blue4', 'gray50')) + 
  geom_segment(data = top_vars_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", size = 1) +
  geom_text_repel(data = top_vars_scaled,
                  aes(x = PC1, y = PC2, label = Variable),
                  color = "black", size = 4,
                  max.overlaps = Inf)




ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_all.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_all.pdf', width = 7, height = 6, units = 'in')



#### PCA WITH ONLY SIGNIFICANT DESCRIPTORS

#Dividir los dfs con los descriptores entre los que tienen cond de 0 o 1
descsYES <- descs %>% filter(dev=='YES') %>% select(-dev)
descsNO <- descs %>% filter(dev=='NO') %>% select(-dev)

means <- data.frame(
  'descs' = c(colnames(descsYES))
)

#calculate means and standard deviation of every descriptor
for (i in 1:length(means[,1])){
  means$meanYES[i] = mean(descsYES[,i])
  means$meanNO[i] = mean(descsNO[,i])
  means$diff_m[i] = abs(means$meanYES[i]-means$meanNO[i])
  means$rsdYES[i] = sd(descsYES[,i])/sqrt(length(descsYES[,i]))
  means$rsdNO[i] = sd(descsNO[,i])/sqrt(length(descsNO[,i]))
  means$diff_rsd[i] = sqrt(means$rsdYES[i]^2+means$rsdNO[i]^2)
}

means <- means %>% filter(diff_rsd<diff_m) #remove descriptors in which delta_sd > delta_mean



#create dataframe with filtered descriptors
descs_fil <- data.frame(matrix(NA,
                             nrow = nrow(descs),
                             ncol = ncol(descsYES)))

for (i in 1:(ncol(descs)-1)){
  for (n in 1:nrow(means)){
    if (colnames(descs)[i] == means$descs[n]) {
      descs_fil[,i] <- descs[,i]
    }
  }
}
descs_fil <- descs_fil %>% select_if(~ !any(is.na(.)))
colnames(descs_fil) <- means$descs
descs_fil$dev <- descs$dev



#_______________________________WELCH'S T-TEST_________________________________
pvalue <- c()
for (i in 1:(ncol(descs_fil)-1)){
  pvalue <- c(pvalue,t.test(ifelse(descs_fil[,ncol(descs_fil)]=='NO',descs_fil[,i],NA),ifelse(descs_fil[,ncol(descs_fil)]=='YES',descs_fil[,i],NA))[[3]])
}
means$welchs_p <- pvalue
means <- means %>% filter(welchs_p<0.05)

descs_fil <- data.frame(matrix(NA,
                               nrow = nrow(descs),
                               ncol = ncol(descsYES)))

for (i in 1:(ncol(descs)-1)){
  for (n in 1:nrow(means)){
    if (colnames(descs)[i] == means$descs[n]) {
      descs_fil[,i] <- descs[,i]
    }
  }
}
descs_fil <- descs_fil %>% select_if(~ !any(is.na(.)))
colnames(descs_fil) <- means$descs



#Correlation plot for the descriptors after UFS
library(ggcorrplot)

corr <- round(cor(descs_fil),2)
p <- ggcorrplot(corr,show.diag = F,type='lower', lab = TRUE, lab_size=1, tl.cex = 6) + 
  scale_fill_gradient2(low='blue4',high='red4',breaks=c(-1,-0.5,0,0.5,1),limit=c(-1,1),name='Correlation')

ggsave(p,file='PLOTS/correlacion_logSo/correlation_plot.png',width = 4000, height = 4000,units = 'px')
ggsave(p,file='PLOTS/correlacion_logSo/correlation_plot.pdf',width = 4000, height = 4000,units = 'px')



####

descs_fil$dev <- descs$dev

write.csv(means, file = "DATABASES/descs_list.csv")




# Perform PCA
pca_result <- prcomp(descs_fil[,-ncol(descs_fil)], scale = TRUE)

# Extract scores and loadings
scores <- as.data.frame(pca_result$x)
loadings <- as.data.frame(pca_result$rotation)

# Calculate variance explained
variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Add grouping variable if provided

scores$Group <- data$dev




# Add top contributors biplot vectors 
loadings$contribution <- loadings$PC1^2 + loadings$PC2^2
top_vars <- loadings %>%
  arrange(desc(contribution)) %>%
  head(5)

# Scale loadings for better visualization
top_vars_scaled <- top_vars
top_vars_scaled$PC1 <- top_vars_scaled$PC1 * 2 * max(abs(scores$PC1))
top_vars_scaled$PC2 <- top_vars_scaled$PC2 * 2 * max(abs(scores$PC2))
top_vars_scaled$Variable <- rownames(top_vars_scaled)


#create plot 
p <- ggplot() +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = c(0.8,0.8)) +
  labs(x = paste0("PC1 (", variance_explained[1], "%)"),
       y = paste0("PC2 (", variance_explained[2], "%)")) + 
  geom_point(data = scores, 
             aes(x = PC1, y = PC2, color = Group),
             alpha = 0.7, size = 4) + 
  scale_color_manual(TeX("Is $\\Delta \\log{S} > 1$?"),values = c('blue4', 'gray50')) + 
  geom_segment(data = top_vars_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", size = 1) +
  geom_text_repel(data = top_vars_scaled,
                  aes(x = PC1, y = PC2, label = Variable),
                  color = "black", size = 4,
                  max.overlaps = Inf)




ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_significant.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/PCA_significant.pdf', width = 7, height = 6, units = 'in')


descs_fil$smiles <- data$smiles
write.csv(descs_fil,file = 'DATABASES/filtered_descriptors_SolCbio3Database_merged.csv')





#### t-SNE dimension reduction
library(Rtsne)
X <- descs_fil[,-c(ncol(descs_fil)-1,ncol(descs_fil))]
y <- descs_fil$dev

set.seed(1234) #set seed for a consistent t-SNE

tsne_results <- Rtsne(X, dims = 2, perplexity = 25, verbose = TRUE, max_iter = 1500, check_duplicates = FALSE)

tsne_df <- data.frame(
  X = tsne_results$Y[, 1],
  Y = tsne_results$Y[, 2],
  digit = y
)


p <- ggplot(tsne_df) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = c(0.87,0.13)) +
  labs(x = paste0("t-SNE Dimension 1"),
       y = paste0("t-SNE Dimension 2")) + 
  geom_point(aes(x = X, y = Y, color = factor(digit)),
             alpha = 0.6, size = 3) + 
  scale_color_manual(TeX("Is $\\Delta \\log{S} > 1$?"),values = c('blue4','gray50'))

ggsave(p, filename = 'PLOTS/correlacion_logSo/tSNE.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/tSNE.pdf', width = 7, height = 6, units = 'in')
