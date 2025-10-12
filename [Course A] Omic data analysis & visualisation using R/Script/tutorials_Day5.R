#Task1
install.packages("ggplot2")
library(ggplot2)

#Task2
ggp = ggplot(master, aes(x=master$log2fold,y= master$mlog10p))
ggp = ggplot(master, aes(x=master$log2fold,y= master$mlog10p)) + geom_point()

ggp = ggplot(master, aes(x=master$log2fold,y= master$mlog10p)) + geom_point(colour="red")

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + geom_point(colour="black") + geom_point(colour="red")
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + geom_point(colour="red") + geom_point(colour="black")
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + geom_point(colour="black") + geom_point(data= master_sig,colour="red")

master_sig_up = subset(master, log2fold > 0 & sig == TRUE)
master_sig_down = subset(master, log2fold < 0 & sig == TRUE)
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + geom_point(colour="black") + geom_point(data= master_sig_up,colour="red") + geom_point(data = master_sig_down, colour = "blue")

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue")

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p")

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  theme_dark()

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5)

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  xlim(c(-20, 20)) + 
  ylim(c(0, 50))

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  xlim(c(-20, 20)) + 
  ylim(c(0, 50)) +
  geom_point(size = 1, alpha = 0.5, shape = 3)

master_sig_up_top5 = master_sig_up[1:5,]
master_sig_down_top5 = master_sig_down[1:5,]

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  xlim(c(-20, 20)) + 
  ylim(c(0, 50)) +
  geom_point(size = 1, alpha = 0.5, shape = 3) +
  geom_text(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")

install.packages("ggrepel")
library(ggrepel)
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  xlim(c(-20, 20)) + 
  ylim(c(0, 50)) +
  geom_point(size = 1, alpha = 0.5, shape = 3) +
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")

png("C:/Users/96328/Downloads/data/data/master_plot_volcano.png", height = 400, width = 400)
print(ggp)
dev.off()

ggp2 = ggplot(master, aes(x= log10(mean_em), y=log2fold)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "The change of gene expression",x="log10mean_em",y="log2fold") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  xlim(c(-20, 20)) + 
  ylim(c(-50, 50)) +
  geom_point(size = 1, alpha = 0.5, shape = 3) +
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")

png("C:/Users/96328/Downloads/data/data/master_plot_MA.png", height = 400, width = 400)
print(ggp2)
dev.off()
