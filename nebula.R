

# functions ---------------------------------------------------------------

nebula_dotplot = ggplot(dotplot_df[dotplot_df$signif, ]) +
  aes(x = name, y = feature, colour = logFC_capped, size = logFDR) +
  #viridis::scale_colour_viridis() + 
  #scale_colour_gradient2(low = "blue", mid = "grey", high = "red", limits = c(-1.01, 3), midpoint = 0,
  #                       name = "Norm.\nenrichment\nscore", na.value = "black") +
  scale_colour_gradientn(values = scales::rescale(c(-1, 0, 1, max(dotplot_df$logFC_capped))), colours =  c("darkblue", "white", "red", "darkred")) +
  scale_size_area(na.value = 3) +
  geom_point(stroke =2) +
  guides(colour = guide_colourbar(title = 'logFC', order = 1),
         size = guide_legend(title = "-log10(q)", order = 2)) +
  xlab("Cell type") +
  ylab("TE") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))


# mouse atlas -------------------------------------------------------------

library(ArchR)
library(tidyverse)
library(jj)
library(nebula)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
proj$cluster_annotation_level0 = from_to(vec = proj$Clusters_0.5, 
                                         old_new_map_vec = c(
                                           C1='Plasma cell', 
                                           C2='undefined', 
                                           C3='Plasma cell', 
                                           C4='undefined', 
                                           C5='B cell', C6='B cell', 
                                           C7='T cell', 
                                           C8='ILC', 
                                           C9='ILC', 
                                           C10='NK cell', 
                                           C11='T cell', 
                                           C12='T cell', 
                                           C13='T cell', 
                                           C14='T cell', 
                                           C15='T cell', 
                                           C16='Mast cell/Eosinophil/Basophil',
                                           C17='Neutrophil', 
                                           C18='Macrophage', 
                                           C19='Macrophage',
                                           C20='Macrophage', 
                                           C21='Macrophage',
                                           C22='Monocyte', 
                                           C23='Macrophage',
                                           C24='DC', 
                                           C25='DC', 
                                           C26='DC',
                                           C27='DC',
                                           C28='DC' 
                                         ))
proj = proj[!proj$cluster_annotation_level0 == 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')
count_df = prepare_comb_norm_df(pconfig, proj, normalize = F)
# # count_df[, sample(seq(ncol(count_df)), 100)] %>%
# #   as.data.frame %>% 
# #   # Reshape
# #   gather(key = name, value = val) %>%
# #   # Basic chart
# #   ggplot(aes(x = val)) +
# #   geom_histogram(colour = "darkgreen", fill = "gray") +
# #   facet_wrap(~name, nrow = 10, scales = 'free') +
# #   ## Theme and looks 
# #   ggtitle("Histograms") +
# #   theme(strip.background = element_rect(fill = "gray80", colour = "black",
# #                                         size = 0.5, linetype = "solid"),
# #         strip.text = element_text(face = "bold"))
# 
# #https://divingintogeneticsandgenomics.rbind.io/post/negative-bionomial-distribution-in-single-cell-rnaseq/
# ##negative binomial distribution?
# ##zero-inflated?
# ##poisson distribution?
# ##what else
# # library(sparseMatrixStats)
# # count_df_raw = base::t(prepare_comb_norm_df(pconfig, proj, normalize = F))
# # TE_means<- sparseMatrixStats::rowMeans2(count_df_raw)
# # TE_vars<- sparseMatrixStats::rowVars(count_df_raw)
# # df<- bind_cols(TE_means = TE_means, TE_vars = TE_vars)
# # df %>% ggplot(aes(x = log10(TE_means), y = log10(TE_vars))) +
# #   geom_point() +
# #   theme_classic(base_size = 14) +
# #   ggtitle("TE mean and var")
# # 
# # model<- lm(TE_vars ~  1* TE_means + I(TE_means^2) + 0, data =df )
# # summary(model)
# # predicted_df<- data.frame(mean = df$TE_means, var_predict = 
# #                             df$TE_means + coef(model) * (df$TE_means)^2 )
# # 
# # df %>%  ggplot(aes(x = log10(TE_means), y = log10(TE_vars))) +
# #   geom_point() +
# #   geom_line(color = "red", data = predicted_df, aes(x = log10(TE_means), y =log10(var_predict))) + 
# #   theme_classic(base_size = 14) +
# #   ggtitle("TE mean and var with negative binomial prediction")
# # 
# # phi <- 1/coef(model)
# # zeros_nb<- (phi/(TE_means + phi))^phi
# # zeros_observed<- apply(count_df_raw, 1, function(x) mean(x ==0))
# # 
# # data.frame(zeros_nb = zeros_nb, zeros_observed = zeros_observed, 
# #            TE_means = TE_means) %>%
# #   ggplot(aes(x =log10(TE_means), y = zeros_observed)) +
# #   geom_point() +
# #   geom_line(aes(x = log10(TE_means), y = zeros_nb), color = "red") +
# #   theme_classic(base_size = 14) +
# #   ggtitle("Observed and predicted zero fraction vs TE means")
# 
# 
# 
# 
# # data(sample_data)
# # df = model.matrix(~X1+X2+cc, data=sample_data$pred)
# # re = nebula(sample_data$count,sample_data$sid,pred=df)
# 
# #aim: find TEs differentially accessible in peripheral tissues versus the spleen accounting for
# #Sample, celltype, and nFrags
# 
# #important: first level of the factor is used as reference for the logFC (shown in intercept)
# dr_df$peripheral_tissue = factor(ifelse(dr_df$Tissue == 'Spleen', 'spleen',  'peripheral'), levels = c('spleen', 'peripheral'))
# dr_df$cluster_annotation_level0 = factor(dr_df$cluster_annotation_level0, levels = c('B cell','T cell','DC','ILC','Macrophage','Mast cell/Eosinophil/Basophil','Monocyte','Neutrophil','NK cell','Plasma cell'))
# nebula_count_df = Matrix::t(count_df) #[, c('AmnSINE1','B1_Mm', 'RLTR48A', 'RMER6C')])
# df = model.matrix(~peripheral_tissue+cluster_annotation_level0, data = dr_df)
# 
# te_dr_df = cbind(jj_get_reduction_coords(proj, 'UMAP'), prepare_comb_norm_df(pconfig, proj, normalize = T))
# jj_plot_features(reduction = te_dr_df, meta_features = 'RLTR48A', cap_top = 'auto', cont_or_disc = 'c')
# jj_plot_numeric_by_group(te_dr_df, feature_column = 'RLTR48A', group_column = 'cluster_annotation_level3', order = T,flip_coordinates = T)
# 
# #data_g = group_cell(count=nebula_count_df,id=dr_df$Sample,pred=df)
# #re = nebula(data_g$count,data_g$id,pred=data_g$pred)
# #do not use normalized data since nebula models the raw counts. Instead specify offset for scaling.
# re = nebula(nebula_count_df, id = dr_df$Sample, pred=df, offset = dr_df$nFrags)#, model = 'PMM') #for scATAC, one could also try the PMM (Poisson distribution)
# #write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-04-13-nebula_te/nebula_peripheral_vs_spleen.RDS')
# #write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-01-peripheral_markers/nebula_peripheral_vs_spleen_nbmm.RDS')
# 
# re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-04-13-nebula_te/nebula_peripheral_vs_spleen.RDS')
# 
# # dr_df$cell_of_interest = factor(ifelse(dr_df$cluster_annotation_level0 == 'T cell', 'T cell',  'Other'), levels = c('Other', 'T cell'))
# # df = model.matrix(~Tissue+cell_of_interest+nFrags, data = dr_df)
# # re = nebula(nebula_count_df, id = dr_df$Sample, pred=df, verbose = T, model = 'PMM') #for scATAC, one could also try the PMM (Poisson distribution)
# 
# res_df = as.data.frame(t(re$summary))
# colnames(res_df) = re$summary$gene
# pval_df = res_df[grepl('p_', rownames(res_df)), ]
# rownames(pval_df) = gsub('p_', '', rownames(pval_df))
# logFC_df = res_df[grepl('logFC_', rownames(res_df)), ] #%>% dplyr::arrange(desc(Ccr8))
# rownames(logFC_df) = gsub('logFC_', '', rownames(logFC_df))
# fdr_df = as.data.frame(matrix(p.adjust(as.numeric(unlist(as.list(pval_df))), method = 'fdr'), ncol = ncol(pval_df), dimnames = list(rownames(pval_df), colnames(pval_df))))
# 
# pval_thres = 0.05#1e-20 #1e-50
# logFC_thres = 0
# pval_df['peripheral_tissueperipheral', ] < pval_thres & logFC_df['peripheral_tissueperipheral', ] > logFC_thres
# #logical_df = pval_df[grepl('Tissue', rownames(pval_df)), ] < pval_thres & logFC_df[grepl('Tissue', rownames(logFC_df)), ] > 0.2
# 
# 
# res_df = re$summary
# res_df[res_df$logFC_peripheral_tissueperipheral > logFC_thres & res_df$p_peripheral_tissueperipheral < pval_thres, ] %>% 
#   dplyr::arrange(desc(logFC_peripheral_tissueperipheral)) %>% 
#   dplyr::select(gene, logFC_peripheral_tissueperipheral, p_peripheral_tissueperipheral )
# 
# res_df[res_df$`logFC_cell_of_interestT cell` > logFC_thres & res_df$`logFC_cell_of_interestT cell` < pval_thres, ] %>% 
#   dplyr::arrange(desc(`logFC_cell_of_interestT cell`)) %>% 
#   dplyr::select(gene, `logFC_cell_of_interestT cell`, `p_cell_of_interestT cell` )
# 
# fdr_tdf = as.data.frame(t(fdr_df))
# res_sig = res_df[res_df$logFC_peripheral_tissueperipheral > logFC_thres & fdr_tdf$peripheral_tissueperipheral < pval_thres, ] %>% 
#   dplyr::arrange(p_peripheral_tissueperipheral ) %>% 
#   dplyr::select(gene, logFC_peripheral_tissueperipheral , p_peripheral_tissueperipheral )
# jj_volcano_plot(res_df, 
#                 logfc_column = 'logFC_peripheral_tissueperipheral', 
#                 pval_column = 'p_peripheral_tissueperipheral', 
#                 symbol_column = 'gene', marker_thres = c(0.5, 5))
# wilcoxon_tes = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-01-peripheral_markers/wilcoxon_peripheral_tes.txt')
# jj_compare_vectors(wilcoxon_tes, res_sig$gene)
# catp(intersect(wilcoxon_tes, res_sig$gene))
# 
# logFC_df %>% dplyr::filter(RLTR48A > 0)
# fdr_df %>% dplyr::filter(RLTR48A < 1e-50)
# fdr_df$RLTR48A < 1e-50 & logFC_df$RLTR48A > 0.2
# #In each cluster, we defined marker genes as those with CPC > 1, 
# # the logarithm of fold change (log(FC)) > 0.2, and p < 1e-50 in the differential expression analysis. 
# # These criteria ensured that these marker genes had abundant expression and showed significantly 
# # higher expression in that cluster than in the remaining clusters. 
# # We then ranked them according to their log(FC) and used the top 40 gene
# logical_df = pval_df < pval_thres & logFC_df > logFC_thres
# 
# 
# library(presto)
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# dr_df$cluster_annotation_use = from_to(dr_df$cluster_annotation_level0, c(Macrophage='Monocyte/Macrophage', Monocyte='Monocyte/Macrophage'))
# cells_compare = c('B cell', 'Monocyte/Macrophage', 'T cell', 'NK cell')
# res_list = list()
# for(i in seq_along(cells_compare)){
#   cell_i = cells_compare[i]
#   count_df_use = count_df[dr_df$cluster_annotation_use == cell_i, ]
# 
#   group_vec = dr_df$Tissue[dr_df$cluster_annotation_use== cell_i]
#   group_vec = factor(ifelse(group_vec == 'Spleen', 'Spleen', 'Peripheral'), levels = c('Peripheral', 'Spleen'))
#   res = wilcoxauc(Matrix::t(count_df_use), y = group_vec)
#   res_list[[cell_i]] = res %>% dplyr::filter(group == 'Peripheral' & padj < 0.05 & logFC > 0) %>% dplyr::arrange(padj)
# }
# res_list_tes = lapply(res_list, '[[', 'feature')
# jj_plot_upsetr(res_list_tes)
# common_tes = Reduce('intersect', res_list_tes)
# writeLines(common_tes, paste0(storeFigPath, 'wilcoxon_peripheral_tes.txt'))


# nebula on TEs in the T/NK/ILC subset ---------------------------------------------------------

library(ArchR)
library(tidyverse)
library(jj)
library(nebula)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

proj = proj[!proj$cluster_annotation_fine == 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')
count_df = t(prepare_comb_norm_df(pconfig, proj, normalize = F))

dr_df$cluster_annotation_fine = factor(dr_df$cluster_annotation_fine, levels = c('CD4_Tnaive', unique(dr_df$cluster_annotation_fine)[!unique(dr_df$cluster_annotation_fine) == 'CD4_Tnaive']))
dr_df$Tissue = factor(dr_df$Tissue, levels = c('Spleen', 'Colon', 'VAT', 'Skin'))

df = model.matrix(~Tissue+cluster_annotation_fine, data = dr_df)

#cpc (default = 0.005)
#A non-negative threshold for filtering low-expressed genes. Genes with counts per cell smaller than the specified value will not be analyzed.
# Remove  118996  genes having low expression.
# Analyzing  135549  genes with  8  subjects and  10676  cells.
re = nebula(count_df, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005)
#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_tnkilc_20230605.RDS')


# te_olap_df = jj_make_overlap_df(da_te_list[c('Colon','Skin', 'VAT')])
# motifs_plot = unique(unlist(da_tf_list))
# motifs_plot = na.omit(unique(c(unlist(motif_olap_df$overlaps[,!(colnames(motif_olap_df$overlaps) %in% c('a','b','c','d','e'))]))))

# heatmap_binding_domain = rep('Other', ncol(plot_mat))
# for(i in binding_domains){
#   heatmap_binding_domain[grepl(i, colnames(plot_mat))] = i
# }
# ha = columnAnnotation('Binding Domain'= heatmap_binding_domain, col = list('Binding Domain' = jj_get_colours(heatmap_binding_domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))
# 


dotplot_df = res_df %>% rownames_to_column('feature') %>% dplyr::filter(feature %in% tes_plot)
dotplot_df = dotplot_df %>% pivot_longer(cols = -1, names_to = 'name', values_to = 'value') %>% dplyr::mutate(info_type = ifelse(grepl('_FDR', name), 'FDR', 'logFC')) %>%
  dplyr::mutate(name = gsub('_FDR|_logFC', '', name)) %>% 
  spread(key = info_type, value = value)
dotplot_df$logFDR <- (-1) * log10(dotplot_df$FDR)
dotplot_df$signif = dotplot_df$FDR < 0.05
dotplot_df$logFC_capped = jj_cap_vals(dotplot_df$logFC, cap_bottom = -1, cap_top = NULL)
anyNA(dotplot_df$logFDR)
anyNA(dotplot_df$logFC_capped)

# nebula_dotplot = ggplot(dotplot_df) +
#   aes(x = name, y = feature, colour = logFC_capped, size = logFDR, shape = signif) +
#   #viridis::scale_colour_viridis() + 
#   #scale_colour_gradient2(low = "blue", mid = "grey", high = "red", limits = c(-1.01, 3), midpoint = 0,
#   #                       name = "Norm.\nenrichment\nscore", na.value = "black") +
#   scale_colour_gradientn(values = scales::rescale(c(-1, 0, 1, max(dotplot_df$logFC_capped))), colours =  c("darkblue", "white", "red", "darkred")) +
#   scale_shape_manual(breaks = c(T, F), values = c(16, 3)) +
#   scale_size_area(na.value = 3) +
#   geom_point(stroke =2) +
#   guides(colour = guide_colourbar(title = 'logFC', order = 1),
#          size = guide_legend(title = "-log10(q)", order = 2),
#          shape = guide_legend(title = 'FDR < 0.05')) +
#   xlab("Cell type") +
#   ylab("TE") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 1, vjust = 0.5),
#         axis.text.y = element_text(colour = "black"),
#         axis.ticks = element_line(colour = "black"))

pdf(paste0(storeFigPath, 'mouse_tnkilc_te_celltype_nebula_dotplot_tisTreg_sig_tes.pdf'), width =8, height=8)
nebula_dotplot
dev.off()

for(k in marker_te_groups){
  message(k)
  tes_plot =  head(da_te_list[[k]], 5)
  te_umaps = jj_plot_features(te_dr_df, features= tes_plot, cap_top = 'auto', return_gg_object = T)
  pdf(paste0(storeFigPath,  sprintf('mouse_tnkilc_top5_te_umaps_%s_tes.pdf', k)),width = 8, height = 6)
  print(te_umaps)
  dev.off()
}



# plot_mat = res_df[rownames(res_df) %in% tes_plot, grepl('_logFC', colnames(res_df))]
# colnames(plot_mat) = gsub('_logFC', '', colnames(plot_mat))
# plot_mat = plot_mat[, colnames(plot_mat) %in% c('Colon', 'Skin', 'VAT')]
# 
# library(circlize)
# Heatmap(t(plot_mat)) #define range based on automatic colorscale
# col_fun = colorRamp2(c(0, 0.1, 0.25, 0.5), c("darkblue", "white", "red", "darkred"))
# h1 <- Heatmap(t(plot_mat), col = col_fun,
#               name = 'logFC vs spleen',
#               column_names_gp = gpar(fontsize = 8), 
#               row_names_gp = gpar(fontsize = 8))#,
#               #top_annotation = ha)
# 
# pdf(paste0(storeFigPath, 'mouse_peripheral_te_heatmap.pdf'), width =10, height=3)
# png(paste0(storeFigPath, 'mouse_peripheral_te_heatmap.pdf'), width = 8, height = 3, res = 400, units = 'in')
# h1
# dev.off()

# te_dr_df$tissue_celltype = paste(te_dr_df$Tissue, te_dr_df$cluster_annotation_fine, sep ='_')
# te_dr_df$tissue_celltype = replace_if(te_dr_df$tissue_celltype, count_below = 30, 'undefined')
# jj_plot_numeric_by_group(te_dr_df[te_dr_df$tissue_celltype != 'undefined',],
#                          feature_column = 'L1Md_F', group_column = 'tissue_celltype', order = T, flip_coordinates = T)
# 
# res_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0.5)
# da_te_list = lapply(res_list, '[[', 'feature')
# mean_df = jj_summarize_sparse_mat(Matrix::t(prepare_comb_norm_df(pconfig, proj, normalize = T)), summarize_by_vec = proj$Tissue)
# mean_df = as.data.frame(mean_df)
# mean_thres = 1
# mean_df = mean_df[with(mean_df, Colon > mean_thres | Skin > mean_thres | VAT > mean_thres), ]
# da_te_list = lapply(da_te_list, function(x) x[x %in% rownames(mean_df)])
# #significantly associated tes *fdr < 0.05 with logFC > 0.5 and norm mean(expression) > 1 in any of the peripheral tissues
# jj_plot_upsetr(da_te_list, nintersects = 20) #reduce so that only intersects > 1 are shown
# Reduce(intersect, da_te_list[c('Colon','Skin', 'VAT')])


### save the t/nk/ilc nebula results for all tes with logFC > 0 and fdr < 0.05 as xlsx file
res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_fine', '', colnames(res_df))
colnames(res_df) = gsub('Tissue', '', colnames(res_df))
sig_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0)
#variables = unique(gsub('_FDR', '', gsub('_logFC', '', colnames(res_df))))
# sig_list = list()
# for(i in seq_along(variables)){
#   peaks_keep = res_df[res_df[, paste0(variables[i], '_FDR') ] < 0.05 & res_df[, paste0(variables[i], '_logFC') ] > 0, ]
#   peaks_keep$feature = rownames(peaks_keep)
#   sig_list[[ variables[i] ]] = peaks_keep[, c('feature', paste0(variables[i], '_FDR'), paste0(variables[i], '_logFC'))]
#   colnames(sig_list[[ variables[i] ]]) = c('feature', 'FDR', 'logFC')
# }
for(i in seq_along(sig_list)){
  other_variables_da_peaks = unlist(sapply(sig_list[-i], '[[', 'feature'))
  sig_list[[i]]$exclusive = !sig_list[[i]]$feature %in% other_variables_da_peaks
}
celltype_da_peaks = unique(unlist(sapply(sig_list[!names(sig_list) %in% c('Colon','Skin','VAT')], '[[', 'feature')))
for(i in c("Colon", "VAT",  "Skin")){
  sig_list[[i]]$tissue_exclusive = !sig_list[[i]]$feature %in% celltype_da_peaks
}
sig_list = lapply(sig_list, function(x) x[order(x$FDR), ])
#jj_save_excel(sig_list, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc0_fdr05.xlsx')


# nebula on peaks in the complete atlas ---------------------------------------------------------


library(ArchR)
library(tidyverse)
library(jj)
library(nebula)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
proj$cluster_annotation_level0 = from_to(vec = proj$Clusters_0.5, 
                                         old_new_map_vec = c(
                                           C1='Plasma cell', 
                                           C2='undefined', 
                                           C3='Plasma cell', 
                                           C4='undefined', 
                                           C5='B cell', C6='B cell', 
                                           C7='T cell', 
                                           C8='ILC', 
                                           C9='ILC', 
                                           C10='NK cell', 
                                           C11='T cell', 
                                           C12='T cell', 
                                           C13='T cell', 
                                           C14='T cell', 
                                           C15='T cell', 
                                           C16='Mast cell/Eosinophil/Basophil',
                                           C17='Neutrophil', 
                                           C18='Macrophage', 
                                           C19='Macrophage',
                                           C20='Macrophage', 
                                           C21='Macrophage',
                                           C22='Monocyte', 
                                           C23='Macrophage',
                                           C24='DC', 
                                           C25='DC', 
                                           C26='DC',
                                           C27='DC',
                                           C28='DC' 
                                         ))
proj = proj[!proj$cluster_annotation_level0 == 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')

#saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_filtered_no_doublets_corrected/ArchR_proj_nonbinary_peak_mat', load = T, dropCells = T)

###
pgr_use = getPeakSet(proj)
proj = addPeakSet(proj, peakSet = pgr_use, force=T)
#check current peak set with proj@peakSet
proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchR_proj_nonbinary_peak_mat')  
proj <- addPeakMatrix(proj,binarize = FALSE, ceiling = 50)
#the variance is 2.5 * mean
proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchR_proj_nonbinary_peak_mat')
pmat = get_peak_mat(proj, binarize = FALSE)
pmat = ceiling(pmat / 2) #fixes problem of equal counts being more abundant, but still var = 1.25*mean
max(pmat) #25

###
dr_df = jj_get_reduction_coords(proj)
#dr_df$peripheral_tissue = factor(ifelse(dr_df$Tissue == 'Spleen', 'spleen',  'peripheral'), levels = c('spleen', 'peripheral'))
ref_celltype = 'B cell'
dr_df$cluster_annotation_level0 = factor(dr_df$cluster_annotation_level0, levels = c(ref_celltype, unique(dr_df$cluster_annotation_level0)[!unique(dr_df$cluster_annotation_level0) == ref_celltype]))
dr_df$Tissue = factor(dr_df$Tissue, levels = c('Spleen', 'Colon', 'VAT', 'Skin'))
df = model.matrix(~Tissue+cluster_annotation_level0, data = dr_df) 

# Note that the Poisson mixed model (PMM) should not be used to test a cell-level predictor because it only estimates the subject-level overdispersion
# also the fragment counts should be used instead of read counts or cut sites: https://twitter.com/TuXinming/status/1522293071215415296
# the data is better modeled by a binomial distribution than poisson, since variance is greater than the mean
#pmat_bin = pmat_use
#pmat_bin[pmat_bin > 1] = 1

#cpc (default = 0.005)
#A non-negative threshold for filtering low-expressed genes. Genes with counts per cell smaller than the specified value will not be analyzed.
#run on the cluster with
re = nebula(pmat, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005, ncore =1, verbose = T)  
#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_mouse_atlas_nonbinary_20230623.RDS') 

# pmat = get_peak_mat(proj)
# rsums = Matrix::rowSums(pmat)
# #peaks accessible in at least 1% of cells
# peak_acc_1pct = rsums >= (ncol(pmat) / 100)
# pmat_use = pmat[peak_acc_1pct, ]
# dr_df$cluster_annotation_level0 = factor(dr_df$cluster_annotation_level0, levels = c('B cell','T cell','DC','ILC','Macrophage','Mast cell/Eosinophil/Basophil','Monocyte','Neutrophil','NK cell','Plasma cell'))
# df = model.matrix(~peripheral_tissue+cluster_annotation_level0, data = dr_df)
# re = nebula(pmat_use, id = dr_df$Sample, pred=df, offset = dr_df$nFrags)#, model = 'PMM') #for scATAC, one could also try the PMM (Poisson distribution)
## write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/nebula_peaks_peripheral_vs_spleen.RDS')

re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_mouse_atlas_nonbinary_20230623.RDS')
table(re$convergence) #convergence worse when binary matrix is used
re$summary

res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_level0', '', colnames(res_df))
colnames(res_df) = gsub('Tissue', '', colnames(res_df))

variables = unique(gsub('_FDR', '', gsub('_logFC', '', colnames(res_df))))
sig_list = list()
for(i in seq_along(variables)){
  sig_list[[ variables[i] ]] = rownames(res_df)[res_df[, paste0(variables[i], '_FDR') ] < 0.05 & res_df[, paste0(variables[i], '_logFC') ] > 0]
}
jj_plot_upsetr(sig_list)

##add gene annotation
library(genomic_region_tools)
sig_list = list()
for(i in seq_along(variables)){
  peaks_keep = res_df[res_df[, paste0(variables[i], '_FDR') ] < 0.05 & res_df[, paste0(variables[i], '_logFC') ] > 0, ]
  sig_list[[ variables[i] ]] = cbind(rownames(peaks_keep), peaks_keep[, c(paste0(variables[i], '_FDR'), paste0(variables[i], '_logFC'))])
  colnames(sig_list[[ variables[i] ]]) = c('feature','FDR', 'logFC')
  rownames(sig_list[[ variables[i] ]]) = NULL
  
}

for(i in seq_along(sig_list)){
  other_variables_da_peaks = unlist(sapply(sig_list[-i], '[[', 'feature'))
  sig_list[[i]]$exclusive = !sig_list[[i]]$feature %in% other_variables_da_peaks
  sig_list[[i]]$query_region = NULL
}

tissues = c('Colon','Skin','VAT')
celltype_da_peaks = unique(unlist(sapply(sig_list[!names(sig_list) %in% tissues], '[[', 'feature')))
for(i in tissues){
  sig_list[[i]]$tissue_exclusive = !sig_list[[i]]$feature %in% celltype_da_peaks
}
sig_list = lapply(sig_list, function(x) x[order(x$FDR), ])

tissue_peaks = lapply(sig_list[tissues], function(x) x[x$tissue_exclusive & x$logFC > 0.5, ])
sapply(tissue_peaks, nrow)
tissue_peaks_gr_list = lapply(tissue_peaks, function(x) convert_granges(x$feature))

proj = archr_add_peak_signatures(proj, signature_list = tissue_peaks_gr_list, signature_name = 'nebula_atlas_tissue')

# nebula on peaks T/NK/ILC subset ---------------------------------------------------------


library(ArchR)
library(tidyverse)
library(jj)
library(nebula)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')


cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/treg_fine_annotation.csv')
cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cluster_annotation_fine %in% c('tisTregST2_prog', 'Treg_naive', 'tisTregST2','pTreg'), ]
#jj_compare_vectors(cluster_annotation_fine_df$cellNames, proj$cellNames[proj$cluster_annotation_level3 == 'Treg'])
cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cellNames %in% proj$cellNames[proj$cluster_annotation_level3 == 'Treg'], ]
proj$cluster_annotation_fine = proj$cluster_annotation_level3
proj$cluster_annotation_fine[match(cluster_annotation_fine_df$cellNames, proj$cellNames)] = cluster_annotation_fine_df$cluster_annotation_fine

#jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level3')
proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'undef'] = 'undefined'
proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'

proj = proj[!proj$cluster_annotation_fine == 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')

pmat = get_peak_mat(proj)
rsums = Matrix::rowSums(pmat)
# #peaks accessible in at least 5% of cells
# peak_acc_pct_filter = rsums >= (ncol(pmat) * 0.05)
peak_acc_pct_filter = T

pgr = getPeakSet(proj)
pgr_use = pgr[peak_acc_pct_filter]

pmat_use = pmat[peak_acc_pct_filter, ]



###
setwd('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/')
proj = addPeakSet(proj, peakSet = pgr_use, force=T)
#check current peak set with proj@peakSet
saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_t_nk_ilc_subset_nonbinary_peakset')
proj <- addPeakMatrix(proj,binarize = FALSE, ceiling = 50)
#the variance is 2.5 * mean
proj = loadArchRProject('mouse_normal/ArchRProject_t_nk_ilc_subset_nonbinary_peakset')
pmat = get_peak_mat(proj, binarize = FALSE)
pmat_use = ceiling(pmat / 2) #fixes problem of equal counts being more abundant, but still var = 1.25*mean

###

dr_df = jj_get_reduction_coords(proj)
#dr_df$peripheral_tissue = factor(ifelse(dr_df$Tissue == 'Spleen', 'spleen',  'peripheral'), levels = c('spleen', 'peripheral'))
dr_df$cluster_annotation_fine = factor(dr_df$cluster_annotation_fine, levels = c('CD4_Tnaive', unique(dr_df$cluster_annotation_fine)[!unique(dr_df$cluster_annotation_fine) == 'CD4_Tnaive']))
dr_df$Tissue = factor(dr_df$Tissue, levels = c('Spleen', 'Colon', 'VAT', 'Skin'))

df = model.matrix(~Tissue+cluster_annotation_fine, data = dr_df)

# Note that the Poisson mixed model (PMM) should not be used to test a cell-level predictor because it only estimates the subject-level overdispersion
# also the fragment counts should be used instead of read counts or cut sites: https://twitter.com/TuXinming/status/1522293071215415296
# the data is better modeled by a binomial distribution than poisson, since variance is greater than the mean
#pmat_bin = pmat_use
#pmat_bin[pmat_bin > 1] = 1

#cpc (default = 0.005)
#A non-negative threshold for filtering low-expressed genes. Genes with counts per cell smaller than the specified value will not be analyzed.
# Remove  118996  genes having low expression.
# Analyzing  135549  genes with  8  subjects and  10676  cells.
re = nebula(pmat_use, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005) #for scATAC, one could also try the PMM (Poisson distribution)
#results based on peak matrix with non-binary fragment counts (not cut sites)
#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/nebula_peaks_tnkilc_full_nonbinary_peakset_20230530.RDS')
#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/nebula_peaks_tnkilc.RDS')
#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/nebula_peaks_binary_tnkilc.RDS')
re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/nebula_peaks_tnkilc_full_nonbinary_peakset_20230530.RDS')

table(re$convergence) #convergence worse when binary matrix is used
re$summary

res_df = as.data.frame(t(re$summary))
colnames(res_df) = re$summary$gene
pval_df = res_df[grepl('p_', rownames(res_df)), ]
rownames(pval_df) = gsub('p_', '', rownames(pval_df))
logFC_df = res_df[grepl('logFC_', rownames(res_df)), ] #%>% dplyr::arrange(desc(Ccr8))
rownames(logFC_df) = gsub('logFC_', '', rownames(logFC_df))
logFC_df = logFC_df[-1, ]
logFC_df = t(logFC_df)
logFC_df = apply(logFC_df, 2, as.numeric)

fdr_df = as.data.frame(matrix(p.adjust(as.numeric(unlist(as.list(pval_df))), method = 'fdr'), ncol = ncol(pval_df), dimnames = list(rownames(pval_df), colnames(pval_df))))
fdr_df = fdr_df[-1, ] #exclude intercept
sum(apply(fdr_df, 2, function(x) all(x < 0.05))) #6
fdr_df = as.data.frame(t(fdr_df))
sig_list = list()
for(i in seq_along(fdr_df)){
  sig_list[[ colnames(fdr_df)[i] ]] = rownames(fdr_df)[fdr_df[, i] < 0.05 & logFC_df[, i] > 0]
}
names(sig_list) = gsub('cluster_annotation_fine', '', names(sig_list))
jj_plot_upsetr(sig_list)
sapply(sig_list, length) %>% sort(., decreasing = T) %>% data.frame(var = names(.), n = .) %>% 
  dplyr::filter(grepl('cluster', var))

# olaps = list()
# for(i in seq_along(sig_list)){
#   olaps[[i]] = jj_compare_vectors(sig_list[[i]], sig_list_nonbin[[i]])
# } 


##add gene annotation
library(genomic_region_tools)
sig_list = list()
for(i in seq_along(fdr_df)){
  peaks_keep = fdr_df[, i] < 0.05 & logFC_df[, i] > 0
  sig_list[[ colnames(fdr_df)[i] ]] = data.frame(feature = rownames(fdr_df[peaks_keep, ]), FDR = fdr_df[peaks_keep, i], logFC = logFC_df[peaks_keep, i])
}

sig_list_anno = lapply(sig_list, function(x) return(cbind(x, closest_feature(convert_granges(x$feature), annotations_gr = geneAnnoMm10$genes))))
for(i in seq_along(sig_list_anno)){
  other_variables_da_peaks = unlist(sapply(sig_list_anno[-i], '[[', 'feature'))
  sig_list_anno[[i]]$exclusive = !sig_list_anno[[i]]$feature %in% other_variables_da_peaks
  sig_list_anno[[i]]$query_region = NULL
}

celltype_da_peaks = unlist(sapply(sig_list_anno[-c(1:3)], '[[', 'feature'))
for(i in c("TissueColon", "TissueVAT",  "TissueSkin")){
  sig_list_anno[[i]]$tissue_exlusive = !sig_list_anno[[i]]$feature %in% celltype_da_peaks
}
names(sig_list_anno) = gsub("Tissue|cluster_annotation_fine", "", names(sig_list_anno))
#jj_save_excel(sig_list_anno, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/tnkilc_nebula_fragment_count_da_peaks_lfc0_fdr05.xlsx')
jj_save_excel(sig_list_anno, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/tnkilc_nebula_fragment_count_da_peaks_lfc0_fdr05_complete_nonbinary_peakset.xlsx')



# other stuff -------------------------------------------------------------


markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix',
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df = archr_get_markers_as_df(markers_se, proj, cutOff = "FDR <= 0.05 & Log2FC >=0")
marker_list = split(marker_df, marker_df$comparison)
marker_list = lapply(marker_list, '[[', 'feature')
sapply(marker_list, length) %>% sort(., decreasing = T) %>% as.data.frame
jj_plot_upsetr(marker_list)

for(i in seq_along(names(marker_list))){
  print(names(marker_list)[i])
  print(jj_compare_vectors(marker_list[[i]], sig_unique_list[[i]]))
} #agreement is bad ~1%

library(genomic_region_tools)
sig_gr_list = lapply(sig_list, convert_granges)
proj = archr_add_peak_signatures(proj, signature_list = sig_gr_list, signature_name = 'nebula')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
colnames(dr_df) = replace_special_chars(colnames(dr_df))
gg = jj_plot_features(reduction = dr_df, meta_features = grep('z_', colnames(dr_df), value = T), cap_top = 'q99',cap_bottom = 'q01', return_gg_object = T)
jj_arrange_ggplots(gg, nplots = 12, cols = 4)

sig_unique_list = list()
for(i in seq_along(sig_list)){
  sig_peaks = sig_list[[i]]
  other_peaks = unlist(sig_list[-i])
  sig_unique_list[[ names(sig_list)[i] ]] = sig_peaks[!sig_peaks %in% other_peaks]
}
sapply(sig_unique_list, length)
sig_unique_gr_list = lapply(sig_unique_list, convert_granges)
proj = archr_add_peak_signatures(proj, signature_list = sig_unique_gr_list, signature_name = 'nebula_unique')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
colnames(dr_df) = replace_special_chars(colnames(dr_df))
gg = jj_plot_features(reduction = dr_df, meta_features = grep('z_', colnames(dr_df), value = T), cap_top = 'q99',cap_bottom = 'q01', return_gg_object = T)
jj_arrange_ggplots(gg, nplots = 12, cols = 4)


tistreg_peaks_fdr_sorted = fdr_df[rownames(fdr_df) %in% sig_unique_list$tisTregST2, ] %>% dplyr::arrange(cluster_annotation_finetisTregST2) %>% rownames
tistreg_sig = convert_granges(sig_unique_list$tisTregST2)
gannot = annotate_surrounding_genes(tistreg_sig, gene_annotation = geneAnnoMm10$genes)
p = plotBrowserTrack(proj, region = convert_granges(tistreg_peaks_fdr_sorted[1])+1000, groupBy = 'cluster_annotation_fine', tileSize = 10)
plot(p)

#transpose = T results in wrong clustering of cell types, nLabel = 0 is not possible...
plotMarkerHeatmap(markers_se, 
                  cutOff = "FDR <= 0.01 & Log2FC >=0.5",  
                  transpose = F, clusterCols = T, binaryClusterRows = T, #clutering not working with transpose = T
                  nLabel = 1) 
# dr_df$cell_of_interest = factor(ifelse(dr_df$cluster_annotation_level0 == 'T cell', 'T cell',  'Other'), levels = c('Other', 'T cell'))
# df = model.matrix(~Tissue+cell_of_interest+nFrags, data = dr_df)
# re = nebula(nebula_count_df, id = dr_df$Sample, pred=df, verbose = T, model = 'NBGMM') 

data_df = as.data.frame(as.matrix(t(pmat_use[1:2, ])))
model = glmer(`chr1-4785512-4786012` ~ 1, data = data_df, family = 'binomial')



### check the distributions (based on fragment count matrix)
# Generate 100 occurrences of flipping 10 coins, each with 30% probability
#rbinom(100,10,.3)

library(sparseMatrixStats)
peak_means<- sparseMatrixStats::rowMeans2(pmat_use)
peak_vars<- sparseMatrixStats::rowVars(pmat_use)
peaksums = rowSums2(pmat_use)
zero_counts = apply(pmat_use, 1, function(x) sum(x==0))
zero_fraction = zero_counts / ncol(pmat_use)
expected_poisson = exp(-peak_means)
expected_binomial = dbinom(x = 0, size = 1, prob = peak_means)
dbinom(0,1,peak_means[1]) #prob of 0 successes in 1 coin flip given the probability from peak_means
1 - mean(rbinom(ncol(pmat_use), 1, peak_means[1])) #same from sampling: 10000 repetitions of 1 coin flip with prob peak_means -> count successes and invert to get failures



df = data.frame(peak_means = peak_means, zero_fraction = zero_fraction, poisson_exp = expected_poisson, bin_expected = expected_binomial)
df %>% ggplot() + geom_point(aes(x = log(peak_means), y = zero_fraction)) + geom_line(aes(x = log(peak_means), y = poisson_exp),color='red')


df<- bind_cols(peak_means = peak_means, peak_vars = peak_vars)
df %>% ggplot(aes(x = peak_means, y = peak_vars)) +
  geom_point() + geom_abline(slope = 1, intercept = 0, colour = 'red') + 
  theme_classic(base_size = 14) +
  ggtitle("Mean fragment number and their variance per peak") + coord_fixed()


#binomial prediction
model<- lm(peak_vars ~  1* peak_means + I(peak_means^2) + 0, data =df )
summary(model)
predicted_df<- data.frame(mean = df$peak_means, var_predict =
                            df$peak_means + coef(model) * (df$peak_means)^2 )
predicted_df$var_predict_poisson = predicted_df$mean

df %>%  ggplot(aes(x = log10(peak_means), y = log10(peak_vars))) +
  geom_point() +
  geom_line(color = "red", data = predicted_df, aes(x = log10(peak_means), y =log10(var_predict)), size = 1) +
  geom_line(color = "blue", data = predicted_df, aes(x = log10(peak_means), y =log10(var_predict_poisson)), size = 1) +
  #geom_abline(slope = 1, intercept = 0, colour = 'blue') + 
  theme_classic(base_size = 14) +
  ggtitle("TE mean and var with negative binomial prediction (red) and poisson prediction (blue)")

phi <- 1/coef(model)
zeros_nb<- (phi/(TE_means + phi))^phi
zeros_observed<- apply(count_df_raw, 1, function(x) mean(x ==0))

data.frame(zeros_nb = zeros_nb, zeros_observed = zeros_observed,
           TE_means = TE_means) %>%
  ggplot(aes(x =log10(TE_means), y = zeros_observed)) +
  geom_point() +
  geom_line(aes(x = log10(TE_means), y = zeros_nb), color = "red") +
  theme_classic(base_size = 14) +
  ggtitle("Observed and predicted zero fraction vs TE means")


# plot peak signatures ----------------------------------------------------

pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
list.dirs(pconfig$ARCHR_DIR, recursive = F)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
proj$cluster_annotation_level0 = from_to(vec = proj$Clusters_0.5, 
                                         old_new_map_vec = c(
                                           C1='Plasma cell', 
                                           C2='undefined', 
                                           C3='Plasma cell', 
                                           C4='undefined', 
                                           C5='B cell', C6='B cell', 
                                           C7='T cell', 
                                           C8='ILC', 
                                           C9='ILC', 
                                           C10='NK cell', 
                                           C11='T cell', 
                                           C12='T cell', 
                                           C13='T cell', 
                                           C14='T cell', 
                                           C15='T cell', 
                                           C16='Mast cell/Eosinophil/Basophil',
                                           C17='Neutrophil', 
                                           C18='Macrophage', 
                                           C19='Macrophage',
                                           C20='Macrophage', 
                                           C21='Macrophage',
                                           C22='Monocyte', 
                                           C23='Macrophage',
                                           C24='DC', 
                                           C25='DC', 
                                           C26='DC',
                                           C27='DC',
                                           C28='DC' 
                                         ))
dr_df = jj_get_reduction_coords(proj, 'UMAP')

#load nebula results
nebula_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/tnkilc_nebula_fragment_count_da_peaks_lfc0_fdr05.xlsx')
nebula_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/tnkilc_nebula_fragment_count_da_peaks_lfc0_fdr05_complete_nonbinary_peakset.xlsx')
nebula_list = lapply(nebula_list, function(x) x[order(x$FDR), ])
nebula_list = lapply(nebula_list[1:3], function(x) x[x$tissue_exlusive & x$logFC > 0.5, ])
nebula_list = lapply(nebula_list[1:3], function(x) x[x$exclusive & x$logFC > 1, ])
nebula_sig_list = lapply(nebula_list, '[[', 'feature')
jj_plot_upsetr(nebula_sig_list)


# nebula on nonbinary peaks (max 25) with full mouse atlas ----------------

re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_mouse_atlas_nonbinary_20230623.RDS') 
