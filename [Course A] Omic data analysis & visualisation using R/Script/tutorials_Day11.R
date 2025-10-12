#Input data
gene_sets = ora_results$geneID
description = ora_results $Description
p.adj = ora_results$p.adjust

#Data clean
ora_results_table = data.frame(
  GeneSet = gene_sets,
  AdjustedP = p.adj,
  row.names = description
)

enriched_gene_set = as.character(ora_results_table [1,1])
candidate_genes = unlist(strsplit(enriched_gene_set, "/"))

indices = match(candidate_genes, rownames(master))
candidate_genes_table = master[indices, , drop = FALSE]
gsea_input = candidate_genes_table$log2fold
names(gsea_input) = row.names(candidate_genes_table)

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

most_enriched_index = which.max(gse_results$NES)
most_enriched_description = gse_results$Description[most_enriched_index]
candidate_genes_string <- gse_results$core_enrichment[most_enriched_index]
candidate_genes <- unlist(strsplit(candidate_genes_string, "/"))
expression_values <- em_symbols_sig[candidate_genes, ]

most_significant_go_index <- which.max(gse_results$NES) 
candidate_genes_string <- gse_results$core_enrichment[most_significant_go_index]
candidate_genes <- unlist(strsplit(candidate_genes_string, "/"))

#BiocManager::install("STRINGdb")
library(STRINGdb)
candidate_genes_table = data.frame(candidate_genes)

names(candidate_genes_table) = "gene"
string_db = STRINGdb$new( version="11.5",
                          species=9606, 
                          score_threshold=200, 
                          network_type="full", 
                          input_directory="")
string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
string_db$plot_network(string_mapped)
