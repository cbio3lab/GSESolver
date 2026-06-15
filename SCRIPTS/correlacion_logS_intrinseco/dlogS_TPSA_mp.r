library(readxl)
library(tidyverse)
library(latex2exp) #enable LaTeX notation 


setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD")


data <- read.csv('DATABASES/SolCbio3Database_merged.csv')



data$dev <- ifelse(data$diff_logS>1,'YES','NO')




descs <- read.csv('DATABASES/descriptors_SolCbio3Database_merged.csv')
descs$diff_logS <- data$diff_logS
descs$dev <- data$dev
descs$logS <- data$logS
descs$mp <- data$mp


p1 <- ggplot(descs) + 
  geom_point(aes(x = TPSA_Obabel, y = diff_logS, color = dev), alpha = 0.4, size = 3) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'none') + 
  labs(x=bquote('TPSA'~(ring(A)^2)),y=TeX("$\\Delta \\log{S}$"), parse = TRUE) + 
  scale_color_manual(name = TeX("Is $\\Delta \\log{S} > 1$?"), values=c("blue4", "gray50"))

p2 <- ggplot(descs) + 
  geom_point(aes(x = mp, y = diff_logS, color = dev),alpha = 0.4, size = 3) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA,colour="black",size=1),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = c(0.8,0.8)) + 
  labs(x= 'mp (°C)',y=TeX("$\\Delta \\log{S}$")) + 
  scale_color_manual(name = TeX("Is $\\Delta \\log{S} > 1$?"), values=c("blue4", "gray50"))


violins <- list(p1,p2)
ncol <- 2
nrow <- 1

library(gridExtra) #create a grid plot

#grid plot creation
grid_violin <- grid.arrange(do.call(arrangeGrob, c(violins,ncol = ncol, nrow = nrow)))
ggsave('PLOTS/correlacion_logSo/mp_TPSA_vslogS.pdf', plot = grid_violin, width = 11, height = 5, units = 'in')
ggsave('PLOTS/correlacion_logSo/mp_TPSA_vslogS.png', plot = grid_violin, width = 11, height = 5, units = 'in')
