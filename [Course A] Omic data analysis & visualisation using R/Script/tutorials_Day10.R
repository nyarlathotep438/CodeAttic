#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)
Path_File = "C:/Users/96328/Downloads/data/data/"
annotations = read.table(paste0(Path_File,"annotations.csv"),header = TRUE,sep = "\t")
em_symbols_sig = read.table(paste0(Path_File,"em_symbols_sig.csv"),header = TRUE,sep = "\t")
master = read.table(paste0(Path_File,"master.csv"), header=TRUE, sep="\t")
sig_genes = row.names(em_symbols_sig)

sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

ggp = barplot(ora_results, showCategory=10)
ggp

ggp = dotplot(ora_results, showCategory=10)
ggp

ggp = goplot(ora_results, showCategory = 10)
ggp

ggp = cnetplot(ora_results, categorySize="pvalue")
ggp

ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont = "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.10)


gene_sig_up = row.names(subset(master,sig==TRUE & log2fold > 0))
gene_sig_down = row.names(subset(master,sig==TRUE & log2fold < 0))
em_symbols_sig_up = em_symbols_sig[gene_sig_up,]
em_symbols_sig_down = em_symbols_sig[gene_sig_down,]
sig_genes_entrez = bitr(gene_sig_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

# we want the log2 fold change
gsea_input = master$log2fold

# add gene names the vector
names(gsea_input) = row.names(master)

# omit any NA values
gsea_input = na.omit(gsea_input)

# sort the list in decreasing order (required for clusterProfiler)
gsea_input = sort(gsea_input, decreasing = TRUE)

gse_results = gseGO(geneList=gsea_input,
                    ont ="CC",
                    keyType = "SYMBOL",
                    nPerm = 10000,
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = org.Mm.eg.db,
                    pAdjustMethod = "none")
ggp = ridgeplot(gse_results)
ggp
