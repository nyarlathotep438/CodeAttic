#Task1
master = na.omit(master)

#Task2
sorted_order = order(master[,"p"], decreasing=FALSE)
master = master[sorted_order,]

#Task3
master$mean_em = rowMeans(master[,2:10])

#Task4
master$mlog10p = -log10(master$p)

#Task5
master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0)

#Task6
em_scaled = t(scale(t(em_symbols)))
em_scaled = na.omit(em_scaled)
em_scaled = data.frame(em_scaled)

#Task7
master_sig = subset(master,p.adj < 0.05 & abs(log2fold) > 1.0)
sig_gene = rownames(master_sig)

#Task8
em_symbols_sig = em_symbols[sig_gene,]
em_scaled_sig = em_scaled[sig_gene,]

#Task9
write.table(master,file="C:/Users/96328/Downloads/data/data/master.csv", sep="\t")
write.table(master_sig,file="C:/Users/96328/Downloads/data/data/master_sig.csv", sep="\t")
write.table(em_symbols,file="C:/Users/96328/Downloads/data/data/em_symbols.csv", sep="\t")
write.table(em_symbols_sig,file="C:/Users/96328/Downloads/data/data/em_symbols_sig.csv", sep="\t")
write.table(em_scaled,file="C:/Users/96328/Downloads/data/data/em_scaled.csv", sep="\t")
write.table(em_scaled_sig,file="C:/Users/96328/Downloads/data/data/em_scaled_sig.csv", sep="\t")