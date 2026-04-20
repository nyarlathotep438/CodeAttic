#Package & Data input####
library(ggplot2)
library(reshape2)
library(amap)
Path_File = "C:/Users/96328/Downloads/data/data/"
em_symbols_sig = read.table(paste0(Path_File,"em_symbols_sig.csv"),header = TRUE,sep = "\t")
em_scaled_sig = read.table(paste0(Path_File,"em_scaled_sig.csv"),header = TRUE,sep = "\t")
sample_sheet = read.table(paste0(Path_File,"sample_sheet.csv"), header=TRUE, sep="\t")

#Data Clean####
em_sig = em_symbols_sig[1:100,]
hm.matrix = as.matrix(em_sig)
hm.matrix = melt(hm.matrix)

heatmap_ggp = ggplot(hm.matrix,aes(x = Var1,y = Var2,fill = value))+
  geom_tile()
heatmap_ggp

em_sig = em_scaled_sig[1:100,]
hm.matrix = as.matrix(em_sig)
hm.matrix = melt(hm.matrix)

heatmap_ggp = ggplot(hm.matrix,aes(x = Var2,y = Var1,fill = value))+
  geom_tile()
heatmap_ggp

colours = c("red", "orange", "yellow","green","blue")
colorRampPalette(colours)(100)

heatmap_ggp = ggplot(hm.matrix,aes(x = Var2,y = Var1,fill = value))+
  geom_tile()+
  scale_fill_gradientn(colours = colorRampPalette(colours)(100))
heatmap_ggp

hm.matrix = as.matrix(em_sig)
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)
hm.matrix_clustered = hm.matrix[y.order,]
hm.matrix_clustered = melt(hm.matrix_clustered)
heatmap_ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = colorRampPalette(colours)(100))
heatmap_ggp

x.dist = Dist(t(hm.matrix), method = "spearman")
x.cluster = hclust(x.dist, method = "average")
x.dd = as.dendrogram(x.cluster)
x.order = order.dendrogram(x.dd)
hm.matrix_clustered = hm.matrix[y.order, x.order]
hm.matrix_clustered = melt(hm.matrix_clustered)
heatmap_ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = colorRampPalette(colours)(100)) +
  ylab("") +
  xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'))
heatmap_ggp

#Task2####
# rug for discrete variable
groups_data = as.matrix(as.numeric(as.factor(sample_sheet$SAMPLE_GROUP)))
groups_data = melt(groups_data)

#rug
rug_colours = c("red","green","blue")
ggp = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = rug_colours) +
  theme(plot.margin=unit(c(0,1,1,1), "cm"), axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
ggp

# Function that gets the default ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# uses the function to get the colors
rug_colours = gg_color_hue(3)