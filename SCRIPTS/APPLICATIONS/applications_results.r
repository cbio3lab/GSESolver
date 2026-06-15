library(readxl)
library(dplyr)
library(latex2exp) #enable LaTeX notation 


setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD")


data <- read_excel('SCRIPTS/APPLICATIONS/LogS agrochemical.xlsx')


data$dev <- ifelse(data$DlogS>1, 'YES','NO')
data$dev_binary <- ifelse(data$dev=='NO', 0,1)
data$GSESolver_binary <- ifelse(data$`GSESolver prediction`=='Use GSE', 0, 1)

#tabulate confusion matrices values into vectors for each ML model

#indexes: 1. true positive, 2. false negative, 3. false positive, 4. true negative
confusion <- c(table(data$dev_binary,data$GSESolver_binary)[1],
          table(data$dev_binary,data$GSESolver_binary)[2],
          table(data$dev_binary,data$GSESolver_binary)[3],
          table(data$dev_binary,data$GSESolver_binary)[4])

accuracy <- round((confusion[1]+confusion[4])/sum(confusion),2)
sensitivity <- round((confusion[4])/(confusion[3]+confusion[4]),2)
specificity <- round((confusion[1])/(confusion[1]+confusion[2]),2)


library(Metrics)
library(ggplot2)


linear <- lm(data$logS~data$logS_GSE)


rmse <- round(rmse(data$logS, data$logS_GSE),2)
mae <- round(mae(data$logS,data$logS_GSE),2)
mse<-round(sum((data$logS-data$logS_GSE)/length(data$logS)),2)
r2 <- round(cor(data$logS,data$logS_GSE)^2,2)

#correlation plot
p <- ggplot(data) +
  geom_point(aes(y=logS,x=logS_GSE,colour=Type, shape = `GSESolver prediction`),size=4,alpha=0.7) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = 'none') + 
  scale_color_manual(values=c('red4','green4', 'magenta4', 'orange2'),name='Type') +
  scale_shape_manual(values = c(17,19)) + 
  geom_abline(slope = coef(linear)[2],intercept = coef(linear)[1]) +
  geom_abline(slope = 1,intercept = 0,linetype='dashed') +
  ylab('logS Experimental') +
  xlab('logS GSE') +
  coord_cartesian(xlim = c(-8,2.5), ylim = c(-8,2.5)) +
  annotate('text',x=-6.5,y=2.5,label=paste("R^2 == ",r2),parse=TRUE) + 
  annotate('text',x=-6.5,y=2,label=paste('RMSE =',rmse)) + 
  annotate('text',x=-6.5,y=1.5,label=paste('MAE =',mae)) + 
  annotate('text',x=-6.5,y=1,label=paste('MSE =',mse)) + 
  annotate('text',x=-6.5,y=0.5,label=paste('Accuracy =',accuracy)) + 
  annotate('text',x=-6.5,y=0,label=paste('Sensitivity =',sensitivity)) + 
  annotate('text',x=-6.5,y=-0.5,label=paste('Specificity =',specificity))
  

ggsave(p, filename = 'PLOTS/APPLICATIONS/corr_agrochemicals.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/APPLICATIONS/corr_agrochemicals.pdf', width = 7, height = 6, units = 'in')






#residuals plot 

data$res <- data$logS_GSE - data$logS

p <- ggplot(data = data ) + 
  # Add the filled rectangle first so it appears behind everything else
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = 1,
           alpha = 0.2, fill = "blue4") +
  geom_segment(aes(x = logS_GSE, y = res, xend = logS_GSE, yend = 0), 
               color = 'gray50', linewidth = 0.1) +
  geom_point(aes(x = logS_GSE, y = res, color = Type, shape = `GSESolver prediction`), 
             size = 4, alpha = 0.8) + 
  geom_hline(yintercept = 0, linewidth = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  labs(x = 'logS GSE', y = 'Residuals (logS units)') + 
  coord_cartesian(xlim = c(-8,2.6), ylim = c(-2.5,5)) + 
  scale_color_manual(values=c('red4','green4', 'magenta4', 'orange2'), name='Type') +
  scale_shape_manual(values = c(17,19)) + 
  scale_y_continuous(position='right')
  
  

ggsave(p, filename = 'PLOTS/APPLICATIONS/res_agrochemicals.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/APPLICATIONS/res_agrochemicals.pdf', width = 7, height = 6, units = 'in')






########### CHEMICAL SPACE ##################



#### t-SNE dimension reduction
library(Rtsne)
app_descs <- read.csv('SCRIPTS/APPLICATIONS/APPLICATIONS_descs.csv')
data_descs <- read.csv('DATABASES/filtered_descriptors_SolCbio3Database_merged.csv')

data_descs$GSESolver.prediction <- ifelse(data_descs$output==0,'Use GSE','Do not use GSE')
data_descs$Type <- "SolCbio3Database"
data_descs <- data_descs %>% select(-c(output,dev))

app_descs <- app_descs %>% select(colnames(data_descs))

tsne_data <- rbind(data_descs,app_descs)

X <- tsne_data[,- c((ncol(tsne_data)-2):ncol(tsne_data))]
y <- tsne_data$GSESolver.prediction

set.seed(1234) #set seed for a consistent t-SNE

tsne_results <- Rtsne(X, dims = 2, perplexity = 25, verbose = TRUE, max_iter = 1500, check_duplicates = FALSE)

tsne_df <- data.frame(
  X = tsne_results$Y[, 1],
  Y = tsne_results$Y[, 2],
  digit = y,
  Type = tsne_data$Type
)

solcbio3 <- tsne_df %>% filter(Type=='SolCbio3Database')
apps <-  tsne_df %>% filter(Type !='SolCbio3Database')


p <- ggplot() +
  theme(
    panel.background = element_blank(), 
    panel.border = element_rect(fill = NA, colour = "black", size = 1),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = 'right'
  ) +
  labs(
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) + 
  geom_point(
    data = solcbio3,
    aes(x = X, y = Y, color = factor(digit)),
    alpha = 0.2, size = 3
  ) + 
  scale_color_manual("SolCbio3Database", values = c("blue4", "gray50")) + 
  geom_point(
    data = apps,
    aes(x = X, y = Y, shape = factor(digit), fill = Type),
    alpha = 0.9, size = 5, color = 'black'
  ) + 
  scale_fill_manual(
    "Agrochemicals", 
    values = c('red4', 'green4', 'black', 'magenta4', 'orange2')
  ) + 
  scale_shape_manual("GSESolver prediction", values = c(24, 21)) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,  # Use a shape that shows fill (like 21 = filled circle)
        color = NA,  # Remove the black border in the legend
        size = 3
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "gray70",  # Add a fill to shape legend for visibility
        size = 3
      )
    )
  )

ggsave(p, filename = 'PLOTS/APPLICATIONS/tSNE.png', width = 10, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/APPLICATIONS/tSNE.pdf', width = 10, height = 6, units = 'in')




##violin plots
app_descs <- read.csv('SCRIPTS/APPLICATIONS/APPLICATIONS_descs.csv')
data_descs <- read.csv('DATABASES/descriptors_SolCbio3Database_merged.csv')
apps_full <- read_excel('SCRIPTS/APPLICATIONS/APPLICATIONS_FULL.xlsx')
solcbio3 <- read.csv('DATABASES/SolCbio3Database_merged.csv')

app_descs$logS <- apps_full$logS
data_descs$logS <- solcbio3$logS

p1 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=MW),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = MW, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=MW),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y='MW (g/mol)')

p2 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=HBA2_Obabel),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = HBA2_Obabel, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=HBA2_Obabel),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y='HBA')


p3 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=HBD_Obabel),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = HBD_Obabel, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=HBD_Obabel),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y='HBD') 

p4 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=ALogP),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = ALogP, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=ALogP),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y='clogP')

p5 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=TPSA_rdkit),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = TPSA_rdkit, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=TPSA_rdkit),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y=bquote('TPSA'~(ring(A)^2)),parse=TRUE) 

p6 <- ggplot() + 
  geom_violin(data = data_descs, aes(x = factor(1), y=logS),alpha=0.7, fill = 'blue4') + 
  geom_jitter(data = app_descs, alpha = 0.7, aes(x = factor(1), y = logS, color = Type), size = 2) + 
  geom_boxplot(data = data_descs, aes(x = factor(1), y=logS),alpha=0.5, size = 1, width = 0.1, fill = 'blue4') + 
  scale_color_manual("",values = c('red4','green4', 'black', 'magenta4', 'orange2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x="",y='logS')





violins <- list(p1,p2,p3,p4,p5,p6)
ncol <- 6
nrow <- 1

library(gridExtra) #create a grid plot

#grid plot creation
grid_violin <- grid.arrange(do.call(arrangeGrob, c(violins,ncol = ncol, nrow = nrow)))
ggsave('PLOTS/APPLICATIONS/chemspace_violins.pdf', plot = grid_violin, width = 11, height = 3.5, units = 'in')
ggsave('PLOTS/APPLICATIONS/chemspace_violins.png', plot = grid_violin, width = 11, height = 3.5, units = 'in')



