library(tidyverse)
library(reshape2)
library(gprofiler2)

out <- read_csv("output.csv")

go <- gprofiler2::gost(out$gene, evcodes = TRUE)
glutamate_genes <- go[["result"]][14,]
glutamate_genes <- glutamate_genes$intersection
glutamate_genes <- strsplit(glutamate_genes, split = ",")

#Divide output in glutamate signaling variants (and those that arent)
glut_res <- function(exp_output, glut){
  exp_output <- exp_output[9:228]
  if ( glut == TRUE){
    out <- exp_output %>% 
      filter(gene %in% glutamate_genes[[1]])
  } else if (glut == FALSE) {
    out <- exp_output %>% 
      filter(!gene %in% glutamate_genes[[1]])
  }
  mout <- melt(out)
  return(mout)
}

avg_tissue <- function(df){
  out <- df %>% 
    group_by(variable) %>% 
    summarize(avg_value = mean(value))
  return(out)
}

avg_gene <- function(df){
  out <- df %>% 
    group_by(gene, variable) %>% 
    summarize(avg_value = mean(value))
  return(out)
}

#Test for glut vs the rest, per tissue
glut <- glut_res(out, TRUE)
non_glut <- glut_res(out, FALSE)
gluttiss <-avg_tissue(glut)
NONgluttiss <-avg_tissue(non_glut)

#Is taking into account 0s so probably skewed
test2 <- cbind(gluttiss, NONgluttiss)
colnames(test2) <- c("gluttis", "avgglut", "nongluttis", "avgnonglut")
t.test(test2$avgnonglut, test2$avgglut, 
       data=test2,
       paired=TRUE,
       conf.level=0.95)

#Generates dataframe for plot
boxplot <-  NULL
boxplot$values <- c(test2$avgglut, test2$avgnonglut) 
boxplot$ids <- NA
boxplot$ids <- rep("Gutamatergic signalling", 218)
boxplot$ids[219:436] <- rep("Control", 218)
boxplot <- as.data.frame(boxplot)

#Plot averages
ggas("fig1.png")
ggplot(boxplot, aes(ids, values))+ 
  theme_minimal() +
  geom_violin() +
  geom_jitter(aes(color = "red"), width = 0.2, alpha = 0.6) +
  labs(title="Average per tissue",
       x ="Category", y = "Average predicted value") +
  theme(legend.position = "none")
dev.off()

# Including everything
glut <- glut_res(out, TRUE)
non_glut <- glut_res(out, FALSE)
boxplot_everything <- NULL
boxplot_everything$values <- c(glut$value, non_glut$value)
boxplot_everything$ids <- NA
boxplot_everything$ids[1:length(glut$value)] <- rep("Glutamate", length(glut$value))
boxplot_everything$ids[length(glut$value)+1:length(non_glut$value)] <- rep("Non glutamate", length(non_glut$value))
boxplot_everything <- as.data.frame(boxplot_everything)

#Plot 
pdf("boxplot.pdf")
aaa <- ggplot(boxplot_everything, aes(ids, values))+ 
  theme_minimal() +
  geom_boxplot() +
  labs(title="Average per tissue",
       x ="Category", y = "Average predicted value") +
  theme(legend.position = "none")
dev.off()


#Test per gene
glutgene <-avg_gene(glut)
NONglutgene<-avg_gene(non_glut)

test <- NULL
test <- rbind(glutgene, NONglutgene)
test$id <- NA
test$id[1:length(glutgene$avg_value)] <- "glut"
test$id[length(glutgene$avg_value)+1:length(NONglutgene$avg_value)] <- "NONglut"

model <- lm(avg_value ~ id,
           data = test)
summary(model)

test <- test %>% 
  filter(avg_value < -0.01)
  
a <- ggplot(test, aes(avg_value, gene))+ 
  theme_minimal() +
  geom_violin() +
  geom_jitter() +
  labs(title="Average per tissue",
       x ="Category", y = "Average predicted value") +
  theme(legend.position = "none")

a <- a + facet_wrap(~ id)
a

# w
test <- NULL
test <- rbind(glut, non_glut)
test$id <-NA
test$id[1:length(glut$value)] <- rep("glut", length(glut$value))
test$id[length(glut$value)+1:length(non_glut$value)] <- rep("nonglut", length(non_glut$value))

# Generates dataframe for cone plot
genecone <- test %>% 
  group_by(gene, id, variable) %>% 
  mutate(sum = sum(value)) %>% 
  mutate(abs = sum(abs(value))) %>% 
  select(gene, variable, id, sum, abs) %>% 
  unique()
  #mutate(highlight, ifelse(abs > 10 && , TRUE, FALSE))


# Extract brain tissues
braintissues <- colnames(out[11:228])
indx <- c(3,8,9, 10, 11, 12, 
       13, 14, 15, 16,17,
       18,20, 42, 61, 60, 70, 81, 106, 107, 136, 
       137, 145, 161, 186, 187, 
       189, 209)
braintissues <- braintissues[indx]



pdf("test5.pdf")
a1 <- ggplot(genecone, aes(sum, abs)) +
  theme_minimal() +
  coord_flip() +  
  geom_point(colour = 'grey', size = 0.1 ) 

a1 <- a1 + geom_point(colour = (ifelse((genecone$id == "glut" & !(genecone$variable %in% braintissues)), 
                      '#440154FF', NA)), size = 0.3) +
  geom_point(colour = (ifelse((genecone$id == "glut" & (genecone$variable %in% braintissues)), 
                                 '#20A387FF', NA)), size = 0.3) +
  labs(title="", x ="Sum of absolute expression", y = "Absolut sum of absolute expression")

dev.off()

nrs <- c(-0.2282983)

h <- genecone %>%
  arrange(sum)

h$gene <- stringr::str_replace_all(h$gene, "ENSG00000163618", "CADPS")
#test <- gprofiler2::gconvert(genecone$gene)
h$gene <- stringr::str_replace_all(h$gene, "ENSG00000178235", "SLITRK1")



a1 + 
  geom_label_repel(data = h[3,], aes(x=sum, y=abs, label=paste(gene, "-", variable)), box.padding = 0.5, max.overlaps = Inf) +
  geom_label_repel(data = h[90408,], aes(x=sum, y=abs, label=paste(gene, "-", variable)), box.padding = 0.5, max.overlaps = Inf) 
  
#313-314 776 1650 3086-3088 grik5
#856-857 1210-1211 1810-1811 3370 -3378 slc12a5
#1302 2621-2623 SORCS1
# 2347-2350 ZC3H12A
#4275-4277 SLC1A1
all <- rbind(glut, non_glut)
all2 <- all %>% 
  group_by(gene) %>% 
  summarize(value = sum(value)) %>% 
  filter(value != 0)

all2  <- all2 %>% 
  mutate(z_scores = value-mean(value)/sd(value))

fdr.out <- fdrtool(all2[[3]], statistic="normal")

h <- all2 %>%
  filter(glutamateornot == TRUE)

DE <- all2 %>% 
  filter(fdr < 0.01)

pdf("DEplot.pdf")
ggplot(data=all2, aes(x=value, y=-1*log10(fdr))) +
  theme_minimal() +
  geom_point() +
  geom_point(data = DE, aes(x=z_scores, y=-1*log10(fdr)), color = "red") +
  geom_label_repel(data = h, aes(x=z_scores, y=-1*log10(fdr), label=name), box.padding = 1, max.overlaps = Inf) +
  labs(title="Differentially expressed genes", x ="Aggregated expression value", y ="1-log10(FDR)")
dev.off()


