
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_areg_gfp.yaml')

#print config variables
pconfig

set.seed(1)
#addArchRThreads(threads = 1) 
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
getwd()

#view available data
list.dirs(getwd(), recursive = F)


proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks'))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F)
jj_plot_features(reduction=dr_df, meta_features='Clusters_0.5', label = T, pt.size = 0.5)

proj$cluster_annotation = from_to(annot_vec = proj$Clusters_0.5, 
                                  old_new_map_vec = c(
                                    C1='myeloid', 
                                    C2='myeloid',
                                    C3='myeloid', 
                                    C4='CD8',
                                    C5='CD4', C6='CD4', C7='CD4',
                                    C8='ILC',
                                    C9='NK', 
                                    C10='undefined', 
                                    C11='CD4',C12='CD4',C13='CD4',C14='CD4',
                                    C15='B', 
                                    C16='B' 
                                  ))
dr_df = get_reduction_coords(proj)
jj_plot_features(reduction=dr_df, meta_features='cluster_annotation', label = F)

# chromvar ----------------------------------------------------------------

signature_list = read_rds('/abi/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
signature_list = signature_list[-c(4)]
signature_list = lapply(signature_list, makeGRangesFromDataFrame)
proj = archr_add_peak_signatures(proj, signature_list, signature_name = 'tisTreg')
gg = jj_plot_features(reduction = dr_df,
                      pt.size = 0.5,
                      meta_features = grep('z_', colnames(dr_df), value=T),
                      cont_or_disc = 'c',
                      colorScale = 'viridis',
                      cap_top = 'q95', cap_bottom = 'q05',
                      custom_theme = theme_minimal(), return_gg_object = T)

pdf(paste0(storeFigPath, 't_subset_tisTreg_signatures.pdf'), width=12, height=8)
  plot_grobs(gg, grobs = 4, cols = 2)
dev.off()

ex_sig = readLines("/abi/data2/simonma/projects/scATAC/analysis/2021-03-03-homer_mouse_cd8/mouse_cd8_tissue_cd8_a_3_10_11_14_15_16_vs_all_peaks.txt")
ex_sig = ex_sig[!ex_sig=='NA']
ex_sig2 = readLines("/abi/data2/simonma/projects/scATAC/analysis/2020-08-20-satpathy_all_t/cd8_cl0_vs_1_6_7_peaks.txt")
ex_sig_list= list(ex_sig = StringToGRanges(ex_sig), ex_sig2 = StringToGRanges(ex_sig2))
proj = archr_add_peak_signatures(proj, ex_sig_list, signature_name = 'hcc_ex_sig_mouse')
dr_df = get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features=c('z_ex_sig', 'z_ex_sig2'), colorScale = 'viridis',cap_top = 'q95', cap_bottom = 'q05')

res = archr_plot_markers(proj, marker_genes = c('Havcr2','Ifng','Pdcd1','Lag3','Tigit','Gzmb','Tnf','Areg','Sell'), use_magic = T)
plot_grobs(res,grobs = 9, cols=3)


# summarize cell types ----------------------------------------------------

#proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks'))
dr_df = get_reduction_coords(proj, redname='UMAP')
jj_plot_features(reduction=dr_df, meta_features = 'Tissue')
jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation')
summary_df = summarise_fractions(dr_df[, c('Tissue', 'cluster_annotation')], summarise_by = 'Tissue')
total_cells = sum(summary_df$ncells)
summary_df$total_frac = summary_df$ncells / total_cells

ggplot(summary_df, aes(x=Tissue, fill = cluster_annotation, y = fraction)) +
  geom_bar(stat='identity', position="fill") + 
  #facet_wrap(.~clinical) +
  coord_flip()  + theme_minimal() + #NoLegend() +
  scale_fill_manual(values=msPickSampleColors(summary_df$cluster_annotation))

summary_table = summary_df %>% dplyr::select(Tissue, cluster_annotation, fraction) %>% 
  pivot_wider(names_from = 'cluster_annotation', values_from = 'fraction')
summary_table[is.na(summary_table)] = 0
summary_table = summary_table %>% dplyr::mutate_at(vars(-Tissue), round, 3)
summary_table = summary_table[, c(T, apply(summary_table[, 2:ncol(summary_table)],2, max) >= 0.001)]
summary_table %>% View


# plot markers ------------------------------------------------------------

gg = jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F, return_gg_object = T,
                 custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
res = archr_plot_markers(proj, marker_genes = c('Cd3d','Cd3e','Cd4', 'Cd8a', 'Cd8b1', 'Sell'))
table(dr_df$singler_label)
dr_df$prediction = dr_df$singler_label
dr_df$prediction[dr_df$prediction %in% c('Fibroblasts','Mast cells', 'Microglia','Neutrophils','Epithelial cells')] = 'other'
gg2 = jj_plot_features(reduction=dr_df, meta_features='prediction', label = F, return_gg_object = T)


pdf(paste0(storeFigPath, 'mouse_areg_gfp.pdf'), width = 6, height= 5)
print(add_umap_arrows(gg[[1]], theme_use = theme_minimal2) )
print(add_umap_arrows(gg2[[1]], theme_use = theme_minimal2) )
print(plot_grobs(res,grobs = 2, cols=2))
dev.off()

jj_plot_categorical_by_group(dr_df, 'prediction', 'Tissue', flip_coordinates = T)

