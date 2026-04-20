#Set_Work_Dictionary####
Work_Dictionary = "C:/Users/96328/Documents/assessment_data"
setwd(Work_Dictionary) 

#Manual Image Save & Load####
save.image("image_assessment_2.RData")
load("image_assessment_2.RData")

#Library####
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)

#Function####

##Gene Expression Differential Analysis
#Make new Differential Expression table with significant gene mark for volcano plot
Mark_Significant_Gene = function(DE)
{
  #Set significant range:p.adj < 0.05 and |log2fold|>1
  DE$sig = as.factor(DE$p.adj < 0.05 & abs(DE$log2fold) > 1)
  
  #Set up three subsets: significantly up-regulated, significantly down-regulated, and no significant change
  DE_no_sig = subset(DE,DE$sig == FALSE)
  DE_sig_up = subset(DE,DE$sig == TRUE & DE$log2fold > 1)
  DE_sig_down = subset(DE,DE$sig == TRUE & DE$log2fold < -1)
  
  #Mark three group with a,b,c
  DE_no_sig$direction = "a"
  DE_sig_up$direction = "b"
  DE_sig_down$direction = "c"
  
  #Re-bind the DE table
  DE_direction = rbind(DE_no_sig, DE_sig_up, DE_sig_down)
  return(DE_direction)
}

#Get Candidate Gene Name List for analysis
Get_Candidate_Gene = function(DE_direction)
{
  #Record the list of significant genes and select the five top genes
  Sig_up_gene = row.names(subset(DE_direction,direction == "b"))
  Sig_up_gene_top5 = Sig_up_gene[1:5]
  Sig_down_gene = row.names(subset(DE_direction,direction == "c"))
  Sig_down_gene_top5 = Sig_down_gene[1:5]
  
  #Sig_gene: All significant genes Candidate_gene:Genes used in volcano plots and box plots
  Sig_gene = c(Sig_up_gene,Sig_down_gene)
  Candidate_gene = c(Sig_up_gene_top5,Sig_down_gene_top5)
  
  #Set an empty list to load all variables for output
  Sig_gene_name = list("Sig_up_gene" = list(),
                       "Sig_down_gene" = list(),
                       "Sig_up_gene_top5" = list(),
                       "Sig_down_gene_top5" = list(),
                       "Sig_gene" = list(),
                       "Candidante_gene" = list())
  
  #Load result into list
  Sig_gene_name$Sig_up_gene = Sig_up_gene
  Sig_gene_name$Sig_down_gene = Sig_down_gene
  Sig_gene_name$Sig_up_gene_top5 = Sig_up_gene_top5
  Sig_gene_name$Sig_down_gene_top5 = Sig_down_gene_top5
  Sig_gene_name$Sig_gene = Sig_gene
  Sig_gene_name$Candidante_gene = Candidate_gene
  
  return(Sig_gene_name)
}

#Volcano Plot for the Different Expression Table after direction
Make_Volcano_Plot = function(DE_direction,DE_sig_up_top5,DE_sig_down_top5)
{
  #Calculate -log10p and record in DE table
  DE_direction$mlog10p = -log10(DE_direction$p.adj)
  
  #Select five significant genes and set up subset
  DE_sig_up_top5 = subset(DE_direction,DE_direction$direction == "b")[1:5,]
  DE_sig_down_top5 = subset(DE_direction,DE_direction$direction == "c")[1:5,]
  
  #Input necessary packages
  library(ggplot2)
  library(ggrepel)
  
  #Set the image style
  Self_theme = theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)
  )
  
  #Plotting a volcano
  ggp_volcano = ggplot(DE_direction, aes(x=log2fold, y=mlog10p,colour = direction)) + 
    Self_theme+
    geom_point() + 
    scale_colour_manual(name = "",
                        values = c("grey","red", "blue"),
                        labels=c("No significate","UP","DOWN")) +
    
    # X-axis & Y-axis adjusted as needed
    xlim(c(-5, 5)) + 
    #ylim(c(0, 13))+
    
    #Set the horizontal and vertical axis labels
    labs(x="log2fold",y="-log10p") +
    
    #Set significant range mark lines
    geom_vline(xintercept= -1, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", linewidth =0.5) +
    geom_vline(xintercept=1, linetype="dashed", color = "grey", linewidth =0.5) +
    
    #Mark five significantly up-regulated or down-regulated genes
    geom_text_repel(data=DE_sig_up_top5, aes(label=row.names(DE_sig_up_top5)),
                    colour = "darkred",
                    size = 4,
                    nudge_x = 1,
                    hjust = 0,
                    direction = "y",
                    segment.size = 0.2, 
                    force = 10) +
    geom_text_repel(data=DE_sig_down_top5, aes(label=row.names(DE_sig_down_top5)),
                    colour = "darkblue",
                    size = 4,
                    nudge_x = -1,
                    hjust = 1, 
                    direction = "y",
                    segment.size = 0.2,
                    force = 10)
  
  return(ggp_volcano)
}

#For Volcano Plot & Candidate Gene List(Call above three functions)
Show_Sig_Gene = function(DE)
{
  #Call the above functions to obtain the volcano plot and the list of significant genes
  DE_direction = Mark_Significant_Gene(DE)
  Candidate_gene_list = Get_Candidate_Gene(DE_direction)
  Volcano_Plot = Make_Volcano_Plot(DE_direction)
  
  #Create the list "Sig_Gene_result" and add variables to it
  Sig_Gene_result = list("sig_gene_name" = list(), "sig_gene_plot" = list())
  Sig_Gene_result$sig_gene_name = Candidate_gene_list
  Sig_Gene_result$sig_gene_plot = Volcano_Plot
  return(Sig_Gene_result)
}

#Boxplot for Candidate Genes
Make_Boxplot = function(em_scaled,candidate_gene)
{
  #Input necessary packages
  library(reshape2)
  library(ggplot2)
  
  #Transpose the scaled expression matrix table and add the sample group column
  gene_data = em_scaled[candidate_gene,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = OD$ss$sample_group
  
  #Melt the table with sample group as id
  gene_data.m = melt(gene_data, 
                     id.vars = "sample_group",
                     variable.name = "Gene",
                     value.name = "Expression")
  
  #Draw the violin box plot
  ggp_boxplot = ggplot(gene_data.m, aes(x = sample_group, y = Expression, fill = sample_group)) +
    geom_violin(trim = FALSE, scale = "width", width = 0.7) +  
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +  
    scale_fill_manual(values = c("green", "purple", "orange")) +  
    facet_wrap(~ Gene, ncol = 5, scales = "free_y") +  
    labs(x = "Sample Group", y = "Scaled Expression") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none") 
  
  return(ggp_boxplot)
}

#Heatmap for the Candidate Genes
Make_Heatmap = function(em_scaled, candidate_gene, sample_sheet) {
  # Load required packages
  library(amap)
  library(ggplot2)
  library(reshape2)
  
  # Subset candidate genes
  gene_data = em_scaled[candidate_gene, , drop = FALSE]
  hm_matrix = as.matrix(gene_data)
  
  # Cluster rows (genes)
  row_dist = Dist(hm_matrix, method = "spearman")
  row_clust = hclust(row_dist, method = "average")
  row_order = row_clust$order
  
  # Sort columns by sample group from sample_sheet
  sample_groups = sample_sheet[colnames(hm_matrix), "sample_group"]
  col_order = order(sample_groups)
  grouped_matrix = hm_matrix[row_order, col_order]
  
  # Prepare data for ggplot
  melted_data = melt(grouped_matrix)
  colnames(melted_data) = c("Gene", "Sample", "Expression")
  
  # Create rug data for group boundaries
  group_factor = factor(sample_groups[col_order], levels = unique(sample_groups[col_order]))
  rug_data = data.frame(
    xstart = c(1, cumsum(table(group_factor))[-length(table(group_factor))] + 0.5),
    xend = cumsum(table(group_factor)) + 0.5,
    group = levels(group_factor)
  )
  
  # Create color palette
  color_palette = colorRampPalette(c("blue", "gray", "red"))(100)
  
  # Create heatmap with rug marks
  heatmap_plot = ggplot(melted_data, aes(x = Sample, y = Gene)) +
    geom_tile(aes(fill = Expression)) +
    # Add rug marks at bottom
    geom_segment(data = rug_data, 
                 aes(x = xstart, xend = xend, y = -0.5, yend = -0.5, color = group),
                 inherit.aes = FALSE, linewidth = 2) +
    scale_fill_gradientn(colors = color_palette) +
    scale_color_manual(values = c("GB_1" = "green", "GB_2" = "purple", "HC" = "orange")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "italic", hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm"), 
      plot.margin = unit(c(1,1,2,1), "lines")
    ) +
    guides(
      fill = guide_colorbar(title.position = "top", 
                            barwidth = unit(0.5, "cm"), 
                            barheight = unit(4, "cm")), 
      color = guide_legend(title.position = "top") 
    ) +
    labs(x = "Samples", y = "Genes", fill = "Expression", color = "Sample Group")
  
  return(heatmap_plot)
}

##Pathway Analysis
#Over-Representation Analysis
Make_ORA = function(candidate_gene, ont = "BP")
{
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  #Convert symbols to entrez for analysis
  candidate_genes_entrez = bitr(candidate_gene,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
  
  #ORA Analysis
  ora_results = enrichGO(gene = candidate_genes_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         readable = T,
                         ont = ont,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.10)
  
  #ORA Key Information Table
  gene_sets = ora_results$geneID
  description = ora_results$Description
  p.adj = ora_results$p.adjust
  ora_results_table = data.frame(cbind(gene_sets, description, p.adj))
  
  #ORA plots
  ggp.barplot = barplot(ora_results, showCategory=10)
  ggp.dotplot = dotplot(ora_results, showCategory=10)
  
  #Creating list for tables and plots
  Total_ORA_results = list("tables" = list(), "plots" = list())
  
  # store tables
  Total_ORA_results$tables$ora_results = ora_results
  Total_ORA_results$tables$ora_results_table = ora_results_table
  
  # store plots
  Total_ORA_results$plots$ggp.barplot = ggp.barplot
  Total_ORA_results$plots$ggp.dotplot = ggp.dotplot
  
  return(Total_ORA_results)
}

##General Function
#General Figure save function
Save_Plot= function(plot, plot.path, plot.height, plot.width)
{
  tiff(plot.path, width = plot.width, height = plot.height)
  print(plot)
  dev.off()
}

#Input all Original Data in work direction
read_csv_to_list = function(folder_path, sep = "\t") 
{
  #Find all CSV Files
  file_list = list.files(path = folder_path, 
                          pattern = "\\.csv$", 
                          full.names = TRUE, 
                          recursive = TRUE)
  Original_data = list()
  
  #Read all CSV files as data frame
  for (file_path in file_list) {
    file_name = tools::file_path_sans_ext(basename(file_path))
    df = read.table(file = file_path,
                     header = TRUE,
                     row.names = 1,
                     sep = sep,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
    
    Original_data[[file_name]] = df
    message(file_name, " Input complete")
  }
  
  return(Original_data)
}


#Data Clean####
#Input Original Data
OD = read_csv_to_list(folder_path = Work_Dictionary,sep = "\t")

#Confirm and display candidate genes (obtain genes lists and volcano plots)
DE_GB_1_vs_HC_result = Show_Sig_Gene(OD$DE_GB_1_vs_HC)
DE_GB_2_vs_HC_result = Show_Sig_Gene(OD$DE_GB_2_vs_HC)
DE_GB_2_vs_GB_1_result = Show_Sig_Gene(OD$DE_GB_2_vs_GB_1)

#Scale Gene Expression Matrix
em_scaled = na.omit(data.frame(t(scale(t(OD$em)))))

#Boxplot for candidate gene
Boxplot_1vH = Make_Boxplot(em_scaled, DE_GB_1_vs_HC_result$sig_gene_name$Candidante_gene)
Boxplot_2vH = Make_Boxplot(em_scaled, DE_GB_2_vs_HC_result$sig_gene_name$Candidante_gene)
Boxplot_2v1 = Make_Boxplot(em_scaled, DE_GB_2_vs_GB_1_result$sig_gene_name$Candidante_gene)

#Heatmap for candidate gene
Heatmap_all = Make_Heatmap(em_scaled,row.names(em_scaled),OD$ss)
Heatmap_1vH = Make_Heatmap(em_scaled, DE_GB_1_vs_HC_result$sig_gene_name$Sig_gene,OD$ss)
Heatmap_2vH = Make_Heatmap(em_scaled, DE_GB_2_vs_HC_result$sig_gene_name$Sig_gene,OD$ss)
Heatmap_2v1 = Make_Heatmap(em_scaled, DE_GB_2_vs_GB_1_result$sig_gene_name$Sig_gene,OD$ss)

##Pathway analysis
#ORA Analysis
ORA_1vH_sig_up = Make_ORA(DE_GB_1_vs_HC_result$sig_gene_name$Sig_up_gene, "BP")
ORA_1vH_sig_down = Make_ORA(DE_GB_1_vs_HC_result$sig_gene_name$Sig_down_gene, "BP")

ORA_2vH_sig_up = Make_ORA(DE_GB_2_vs_HC_result$sig_gene_name$Sig_up_gene, "BP")
ORA_2vH_sig_down = Make_ORA(DE_GB_2_vs_HC_result$sig_gene_name$Sig_down_gene, "BP")

ORA_2v1_sig_up = Make_ORA(DE_GB_2_vs_GB_1_result$sig_gene_name$Sig_up_gene, "BP")
ORA_2v1_sig_down = Make_ORA(DE_GB_2_vs_GB_1_result$sig_gene_name$Sig_down_gene, "BP")

#Save Plots to be displayed in the report
#List all useful plots
Plots_useful = list(
  "1vH_volcano" = DE_GB_1_vs_HC_result$sig_gene_plot,
  "2vH_volcano" = DE_GB_2_vs_HC_result$sig_gene_plot,
  "2v1_volcano" = DE_GB_2_vs_GB_1_result$sig_gene_plot,
  "Boxplot_1vH" = Boxplot_1vH,
  "Boxplot_2vH" = Boxplot_2vH,
  "Boxplot_2v1" = Boxplot_2v1,
  "1vH_sig_up" = ORA_1vH_sig_up$plots$ggp.barplot,
  "1vH_sig_down" = ORA_1vH_sig_down$plots$ggp.barplot,
  "2vH_sig_up" = ORA_2vH_sig_up$plots$ggp.barplot,
  "2vH_sig_down" = ORA_2vH_sig_down$plots$ggp.barplot,
  "2v1_sig_up" = ORA_2v1_sig_up$plots$ggp.barplot,
  "2v1_sig_down" = ORA_2v1_sig_down$plots$ggp.barplot
)

#Save all useful plots
for (plot_name in names(Plots_useful)) {
  Save_Plot(
    plot = Plots_useful[[plot_name]],
    plot.path = paste0(plot_name, ".tiff"), 
    plot.height = 549,
    plot.width = 537
  )
}
