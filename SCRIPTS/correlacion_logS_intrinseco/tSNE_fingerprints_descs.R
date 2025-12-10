library(readxl)
library(tidyverse)
library(latex2exp) #enable LaTeX notation 


setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD")


data <- read.csv('DATABASES/SolCbio3_descs_fps.csv')



#### t-SNE dimension reduction
library(Rtsne)
X <- data[,-c(1,103)]
y <- data$dev

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
  scale_color_manual(TeX("Is $\\Delta \\log{S} > 1$?"),values = c('gray50','blue4'))

ggsave(p, filename = 'PLOTS/correlacion_logSo/tSNE_fps_descs.png', width = 7, height = 6, units = 'in')
ggsave(p, filename = 'PLOTS/correlacion_logSo/tSNE_fps_descs.pdf', width = 7, height = 6, units = 'in')
