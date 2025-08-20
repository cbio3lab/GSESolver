library(tidyverse)

sol <- read.csv('SolCbio3Database_merged.csv')

sol_sinlogPNA <- sol %>% filter(!is.na(logPN))
sol_sinmpNA <- sol_sinlogPNA %>% filter(!is.na(mp))

write.csv(sol_sinmpNA, file = "SolCbio3Database_merged_noNA.csv")
