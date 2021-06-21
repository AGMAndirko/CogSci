library(tidyverse)
library(ggrepel)
library(reshape2)
library(gprofiler2)
library(fdrtool)

#Reads ExPecto output, including all genes under positive selection according to Peyregne et al. 2017
out <- read_csv("output.csv")

#Does GO enrichment analysis of genes on those genes
go <- gprofiler2::gost(out$gene, evcodes = TRUE)
glutamate_genes <- go[["result"]][14,]
glutamate_genes <- glutamate_genes$intersection
glutamate_genes <- strsplit(glutamate_genes, split = ",")

# Divides ExPecto output into glutamate signaling variants (argument 2 = TRUE) and non-glutamatergic ones (arg2 = FALSE)
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

# Averages values by tissue, regardless of gene
avg_tissue <- function(df){
  out <- df %>% 
    group_by(variable) %>% 
    summarize(avg_value = mean(value))
  return(out)
}

# Averages values by gene and tissue
avg_gene <- function(df){
  out <- df %>% 
    group_by(gene, variable) %>% 
    summarize(avg_value = mean(value))
  return(out)
}

#Test for glut vs the rest, per tissue
# Creates two dfs, one for glutamatergic genes and one that isn't
glut <- glut_res(out, TRUE)
non_glut <- glut_res(out, FALSE)

# Average by tissue
gluttiss <-avg_tissue(glut)
NONgluttiss <-avg_tissue(non_glut)

#Test2 is a dataframe to test average value variance between the two
test2 <- cbind(gluttiss, NONgluttiss)
colnames(test2) <- c("gluttis", "avgglut", "nongluttis", "avgnonglut")
t.test(test2$avgnonglut, test2$avgglut, 
       data=test2,
       paired=TRUE,
       conf.level=0.95)
# Significant

# Generates ex profeso dataframe for boxplot
boxplot <-  NULL
boxplot$values <- c(test2$avgglut, test2$avgnonglut) 
boxplot$ids <- NA
boxplot$ids <- rep("Gutamatergic signalling", 218)
boxplot$ids[219:436] <- rep("Control", 218)
boxplot <- as.data.frame(boxplot)

# Violin plot with averages per tissue (red dots) in glutamatergic and non glutamatergic signalling genes
pdf("fig1.pdf")
ggplot(boxplot, aes(ids, values))+ 
  theme_minimal() +
  geom_violin() +
  geom_jitter(aes(color = "red"), width = 0.2, alpha = 0.6) +
  labs(title="Average per tissue",
       x ="Category", y = "Average predicted value") +
  theme(legend.position = "none")
dev.off()

# Generates ex profeso dataframe for an extra boxplot with every variant
glut <- glut_res(out, TRUE)
non_glut <- glut_res(out, FALSE)
boxplot_everything <- NULL
boxplot_everything$values <- c(glut$value, non_glut$value)
boxplot_everything$ids <- NA
# column "ids" denotes the group to which the variant belongs (glutamate, not glutamate)
boxplot_everything$ids[1:length(glut$value)] <- rep("Glutamate", length(glut$value))
boxplot_everything$ids[length(glut$value)+1:length(non_glut$value)] <- rep("Non glutamate", length(non_glut$value))
boxplot_everything <- as.data.frame(boxplot_everything)

# Figure 2: boxplot with all variants
# Note: Might be heavy or slow due to the number of variants
pdf("boxplot.pdf")
aaa <- ggplot(boxplot_everything, aes(ids, values))+ 
  theme_minimal() +
  geom_boxplot() +
  labs(title="Average per tissue",
       x ="Category", y = "Average predicted value") +
  theme(legend.position = "none")
dev.off()

# Testing significance per gene
# Average contrast groups by tissue and gene
glutgene <-avg_gene(glut)
NONglutgene<-avg_gene(non_glut)

# Test = ex profeso dataframe
test <- NULL
test <- rbind(glutgene, NONglutgene)
test$id <- NA
test$id[1:length(glutgene$avg_value)] <- "glut"
test$id[length(glutgene$avg_value)+1:length(NONglutgene$avg_value)] <- "NONglut"

# Linear model
model <- lm(avg_value ~ id,
           data = test)
summary(model)
# Significant

# Plot those values that are lower than -0.01 in both groups
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

# Creates another ex profeso dataframe, resets "test" as null to avoid confusion
test <- NULL
test <- rbind(glut, non_glut)
test$id <-NA
test$id[1:length(glut$value)] <- rep("glut", length(glut$value))
test$id[length(glut$value)+1:length(non_glut$value)] <- rep("nonglut", length(non_glut$value))

# Generates dataframe for cone plot
genecone <- test %>% 
  group_by(gene, id, variable) %>% 
  mutate(sum = sum(value)) %>% # directional value per gene, tissue
  mutate(abs = sum(abs(value))) %>%  # absolut value per gene, tissue
  select(gene, variable, id, sum, abs) %>% 
  unique()
  #mutate(highlight, ifelse(abs > 10 && , TRUE, FALSE)) # Creates a cuttof point to highlight

# Labels & highlights
# Extracts brain tissues of interest to highlight in plot
braintissues <- colnames(out[11:228])
indx <- c(3,8,9, 10, 11, 12, 
       13, 14, 15, 16,17,
       18,20, 42, 61, 60, 70, 81, 106, 107, 136, 
       137, 145, 161, 186, 187, 
       189, 209)
braintissues <- braintissues[indx]

# Arrange dataframe by value, extract names from gene ids of top and bottom values for plot:
lab <- genecone %>%
  arrange(sum)
lab$gene <- stringr::str_replace_all(lab$gene, "ENSG00000163618", "CADPS")
lab$gene <- stringr::str_replace_all(lab$gene, "ENSG00000178235", "SLITRK1")

# Coneplot, highlighting brain tissues though ifelse conditionals + labels for extreme values
pdf("coneplot.pdf")
a1 <- ggplot(genecone, aes(sum, abs)) +
  theme_minimal() +
  coord_flip() +  
  geom_point(colour = 'grey', size = 0.1 ) + 
  geom_point(colour = (ifelse((genecone$id == "glut" & !(genecone$variable %in% braintissues)), 
                      '#440154FF', NA)), size = 0.3) +
  geom_point(colour = (ifelse((genecone$id == "glut" & (genecone$variable %in% braintissues)), 
                                 '#20A387FF', NA)), size = 0.3) +
  labs(title="", x ="Sum of absolute expression", y = "Absolut sum of absolute expression") + 
  geom_label_repel(data = lab[3,], aes(x=sum, y=abs, label=paste(gene, "-", variable)), box.padding = 0.5, max.overlaps = Inf) +
  geom_label_repel(data = lab[90408,], aes(x=sum, y=abs, label=paste(gene, "-", variable)), box.padding = 0.5, max.overlaps = Inf) 
dev.off()


# join both glutamate and non-glutamate dataframes (should be the same as "out")
all <- rbind(glut, non_glut)

# Take out 0's because they mess up the z score
all <- all %>% 
  group_by(gene) %>% 
  summarize(value = sum(value)) %>% 
  filter(value != 0) %>% 
  mutate(glutamateornot = gene %in% glutamate_genes[[1]])


# Calculates z scores and local FDR
all  <- all %>% 
  mutate(z_scores = value-mean(value)/sd(value)) 
fdr.out <- fdrtool(all[[4]], statistic="normal")
all$fdr <- fdr.out$lfdr

# Creates df for highlights
highlight <- all %>%
  filter(glutamateornot == TRUE)

namedf <- gprofiler2::gconvert(highlight$gene)
#This should give all TRUE:
table(namedf$input == highlight$gene)
highlight$name <- namedf$name
highlight <- highlight %>% 
  filter(fdr < 0.01)


#Filters genes by false discovery rate cutoff
DE <- all %>% 
  filter(fdr < 0.01)

# Plot genes with higher than expected expression, plus label for the glutamate genes
pdf("DEplot.pdf")
ggplot(data=all, aes(x=value, y=-1*log10(fdr))) +
  theme_minimal() +
  geom_point() +
  geom_point(data = DE, aes(x=z_scores, y=-1*log10(fdr)), color = "red") +
  geom_label_repel(data = highlight, aes(x=z_scores, y=-1*log10(fdr), label=name), box.padding = 1, max.overlaps = Inf) +
  labs(title="Differentially expressed genes", x ="Aggregated expression value", y ="1-log10(FDR)")
dev.off()


