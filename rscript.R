library(tidyverse)
library(reshape2)
library(gprofiler2)

out <- read_csv("output.csv")

go <- gprofiler2::gost(out$gene, evcodes = TRUE)
glutamate_genes <- go[["result"]][14,]
glutamate_genes <- glutamate_genes$intersection

out_2 <- out[9:228]

#313-314 776 1650 3086-3088 grik5
#856-857 1210-1211 1810-1811 3370 -3378 slc12a5
#1302 2621-2623 SORCS1
# 2347-2350 ZC3H12A
#4275-4277 SLC1A1

mout <- melt(out_2)

resumen <- mout %>% 
  group_by(variable) %>% 
  summarize(avg_value = mean(value)) %>% 
  arrange(desc(avg_value)) %>% 

resumen$variable <-  as.factor(resumen$variable)
resumen$variable <- fct_reorder(resumen$variable, desc(avg_value))

pdf("test.pdf", height = 10, width = 45)
ggplot(resumen, aes(reorder(variable, avg_value), avg_value))+
  theme_minimal() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))
dev.off()
  
out_2 <- out_2 %>% select(-strand)
mout <- melt(out_2)

mout_2 <- mout %>% 
  group_by(gene) %>% 
  mutate(abs_value = sum(abs(value))) %>% 
  mutate(dir_value = sum(value)) %>% 
  select(-value, -variable)

mout_2 <- unique(mout_2)

mout_2$abs_value <- as.numeric(mout_2$abs_value)
mout_2$dir_value <- as.numeric(mout_2$dir_value)

ggplot(mout_2, aes(abs_value)) %>% 
  theme_minimal() +
  geom_point()
