library(ggplot2)
library(ggrepel)

#Task1
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(colour="black") + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "Volcano",x="log2fold",y="-log10p") +
  geom_vline(xintercept= -1, linetype="dashed", color = "grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(aes(colour="black")) + 
  geom_point(data= master_sig_up,colour="red") + 
  geom_point(data = master_sig_down, colour = "blue") + 
  labs(title = "Volcano",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept= -1, linetype="dashed", color = "grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", size=0.5) +
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5), colour = "red")) +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5), colour = "blue"))

ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  geom_point(aes(colour="a")) + 
  geom_point(data= master_sig_up, aes(colour="b")) + 
  geom_point(data = master_sig_down, aes( colour = "c")) + 
  labs(title = "The change of gene expression",x="log2fold",y="-log10p") +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  geom_point(size = 1, alpha = 0.5, shape = 1) +
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5), colour = "b")) +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5), colour = "c")) +
  scale_colour_manual(values = c("black", "red", "blue"), labels=c("No significate","up","down"))



my_theme = theme(
  plot.title = element_text(size=20),
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title.x = element_text(size=16),
  axis.title.y = element_text(size=16)
)

#Task2
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + geom_point()
ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = sig)) + geom_point()

ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = sig)) + 
  geom_point()+ 
  scale_colour_manual(values = c("red", "blue"), labels=c("True","False"))

master$sig = factor(master$direction, levels = c("a", "b","c"))
levels(master$direction)

master_non_sig = subset(master, sig==FALSE)
master_non_sig$direction = "a"
master_sig_up$direction = "b"
master_sig_down$direction = "c"
master = rbind(master_non_sig, master_sig_up, master_sig_down)

ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = direction)) + 
  geom_point()+ 
  scale_colour_manual(values = c("black","red", "blue"), labels=c("No significate","up","down"))

ggp2 = ggplot(master, aes(x= log10(mean_em), y=log2fold,colour = direction)) + 
  geom_point()+
  scale_colour_manual(values = c("black","red", "blue"), labels=c("No significate","up","down"))

ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = direction)) + 
  geom_point()+ 
  scale_colour_manual(values = c("black","red", "blue"), labels=c("No significate","up","down"))+
  my_theme +
  geom_vline(xintercept= -1, linetype="dashed", color = "grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  labs(title = "Volcano of gene expression",x="log2fold",y="-log10p") +
  xlim(c(-30, 30)) + 
  ylim(c(0, 50))+
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")

ggp2 = ggplot(master, aes(x= log10(mean_em), y=log2fold,colour = direction)) + 
  geom_point()+
  scale_colour_manual(values = c("black","red", "blue"), labels=c("No significate","up","down"))+
  my_theme +
  geom_vline(xintercept= -1, linetype="dashed", color = "grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) +
  labs(title = "MA of gene expression",x="log10mean_em",y="log2fold") +
  xlim(c(-20, 20)) + 
  ylim(c(-50, 50))+
  geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)), colour = "red") +
  geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)), colour = "blue")
