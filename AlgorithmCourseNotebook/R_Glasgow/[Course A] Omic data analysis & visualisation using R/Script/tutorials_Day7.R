###Package & Data Input###
library(ggplot2)
ss = read.table("C:/Users/96328/Downloads/data/data/sample_sheet.csv", header=TRUE, sep="\t")
em_symbols = read.table("C:/Users/96328/Downloads/data/data/em_symbols.csv", header=TRUE, sep="\t")
em_scaled = read.table("C:/Users/96328/Downloads/data/data/em_scaled.csv", header=TRUE, sep="\t")

###Data clean1###
em_symbols = em_symbols[ , ss$SAMPLE]

em_symbols = as.matrix(sapply(em_symbols,as.numeric))
pca = prcomp(t(em_symbols))
pca_coordinates = data.frame(pca$x)

pca_coordinates$color_group = ss$SAMPLE_GROUP

vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100
x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
                     
pca_ggp = ggplot(pca_coordinates,aes(x=PC1, y=PC2,colour = color_group)) + 
  geom_point() +
  scale_colour_manual(values = c("black","red", "blue"), labels=c("duct","gut","node")) +
  labs(x= x_axis_label ,y= y_axis_label)
pca_ggp

###Data clean2###
em_scaled = em_scaled[ , ss$SAMPLE]

em_scaled = as.matrix(sapply(em_scaled,as.numeric))
pca = prcomp(t(em_scaled))
pca_coordinates = data.frame(pca$x)

pca_coordinates$color_group = ss$SAMPLE_GROUP

vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100
x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")

pca_ggp = ggplot(pca_coordinates,aes(x=PC1, y=PC2,colour = color_group)) + 
  geom_point() +
  scale_colour_manual(values = c("black","red", "blue"), labels=c("duct","gut","node")) +
  labs(x= x_axis_label ,y= y_axis_label)
pca_ggp

#expression density plot
ed_plot_gut1 = ggplot(em_symbols,aes(x = log10(gut_r1))) +
  xlim(c(0,6))+
  geom_density(fill = "red",size = 0,alpha = 0.5)
ed_plot_gut1

install.packages("reshape2")
library(reshape2)

em_symbols.m = melt(em_symbols)
ed_plot = ggplot(em_symbols.m,aes(x = log10(value))) +
  geom_density(fill = "red",size = 0,alpha = 0.5)+
  facet_wrap(~Var2, ncol=9)+
  theme(panel.spacing = unit(1, "lines"))
ed_plot


