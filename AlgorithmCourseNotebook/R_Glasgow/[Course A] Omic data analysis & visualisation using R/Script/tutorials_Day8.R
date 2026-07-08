#Input library & Data
library(ggplot2)
library(reshape2)
File_path = "C:/Users/96328/Downloads/data/data/"
em_symbols = read.table("C:/Users/96328/Downloads/data/data/em_symbols.csv", header=TRUE, sep="\t")
em_scaled = read.table(paste0(File_path,"em_scaled.csv"),header = TRUE, sep = "\t")
sample_sheet = read.table("C:/Users/96328/Downloads/data/data/sample_sheet.csv", header=TRUE, sep="\t")
master = read.table(paste0(File_path,"master.csv"), header=TRUE, sep="\t")

#Data clean
gene = "Cdc45"
gene_data = em_symbols[gene,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = sample_sheet$SAMPLE_GROUP
names(gene_data) = c("expression","sample_group")

gene_data$sample_group = factor(gene_data$sample_group, levels=c("gut","duct","node"))

vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression))+
  geom_boxplot()
vilolin_ggp

vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression))+
  geom_boxplot()+ 
  geom_violin()

vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression))+
  geom_boxplot()+ 
  geom_jitter()

vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression))+
  geom_boxplot()+ 
  geom_jitter()

colour = "red"
fill = "black"

vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression,colour = sample_group,fill = sample_group))+
  geom_boxplot()+ 
  geom_jitter(width = 0.1, colour = "red")
vilolin_ggp

#Task2
sorted_order = order(master[,"p"], decreasing=FALSE)
master = master[sorted_order,]
gene1 = row.names(master)[1]
gene2 = row.names(master)[2]
gene3 = row.names(master)[3]

gene_data = em_symbols[gene2,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = sample_sheet$SAMPLE_GROUP
names(gene_data) = c("expression","sample_group")
vilolin_ggp = ggplot(gene_data,aes(x=sample_group,y=expression,colour = sample_group,fill = sample_group))+
  geom_boxplot()+ 
  geom_jitter(width = 0.1, colour = "red")
vilolin_ggp

#Task3
candidate_genes = row.names(master)[1:10]
gene_data = em_scaled[candidate_genes,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = sample_sheet$SAMPLE_GROUP
gene_data = gene_data[,-11]
#gene_data.m = melt(as.matrix(gene_data), id_vars = sample_group, candidate_genes)

gene_data.m = melt(as.matrix(gene_data))

#try and create a vector that contains sample groups 10x
sample_groups_vector = sample_sheet$SAMPLE_GROUP
sample_groups_vector = rep(sample_groups_vector, times = 10)
sample_groups_vector
gene_data.m$sample_group = sample_groups_vector

vilolin_ggp = ggplot(gene_data.m,aes(x = Var2, y=value,colour = sample_group,fill = sample_group))+
  geom_boxplot()+ 
  geom_jitter(width = 0.1, colour = "red")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
vilolin_ggp

#Task4
vilolin_ggp = ggplot(gene_data.m,aes(x = sample_group, y=value,colour = sample_group,fill = sample_group))+
  geom_boxplot()+ 
  geom_jitter(width = 0.1, colour = "red")+
  facet_wrap(~Var2, ncol=5)
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
vilolin_ggp
