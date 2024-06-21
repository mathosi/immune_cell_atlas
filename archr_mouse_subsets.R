library(jj)
library(ArchR)

plot_sigs = function(dr_df, plot_tissue = F, features_plot =c('core tisTregST2 signature' = 'z_core_tisTreg_sig',
                                             'Skin tisTregST2 signature' = 'z_skin_treg_sig',
                                             'VAT tisTregST2 signature' = 'z_fat_treg_sig',
                                             'Colon tisTregST2 signature' = 'z_colon_treg_sig',
                                             'Early progenitor signature' = 'z_early_precursor_sig',
                                             'Late progenitor signature' = 'z_late_precursor_sig')){
  
  gg = list(); for(i in features_plot){
    gg[[i]] = jj_plot_features(reduction=dr_df, meta_features=i,my_title = names(features_plot)[features_plot==i],
                               colorScale = 'viridis', cap_top = 'q95', return_gg_object = T, pt.size=0.5)[[1]] + 
      theme_minimal2 +  theme(legend.position="none")
  }
  if(plot_tissue){
    gg_tissue = jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5, return_gg_object = T,
                                 custom_colors = jj_get_colours(dr_df$Tissue,'/abi/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
    gg = c(gg_tissue, gg)
  }
  return(gg)
}

#mouse subset analysis\
#addArchRThreads(threads = 1) 
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))

#fix tissue annotation
proj$Sample_old = proj$Sample
proj$Sample_new = plyr::mapvalues(proj$Sample, from=c('MD-scATAC_71_Spleen',
                                                      'MD-scATAC_72_Spleen',
                                                      'MD-scATAC_73_Colon',
                                                      'MD-scATAC_74_Colon',
                                                      'MD-scATAC_75_VAT',
                                                      'MD-scATAC_76_VAT',
                                                      'MD-scATAC_77_Skin',
                                                      'MD-scATAC_78_Skin'),
                                  to=c('MD-scATAC_71_Spleen',
                                       'MD-scATAC_72_Spleen',
                                       'MD-scATAC_73_VAT',
                                       'MD-scATAC_74_VAT',
                                       'MD-scATAC_75_Skin',
                                       'MD-scATAC_76_Skin',
                                       'MD-scATAC_77_Colon',
                                       'MD-scATAC_78_Colon'))
proj$Tissue = sapply(strsplit(proj$Sample_new, '_'), '[[', 3)

dr_df = get_reduction_coords(proj)
gg1 = jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F, return_gg_object = T)[[1]]
gg2 = jj_plot_features(reduction=dr_df, meta_features='Clusters_0.5', label = T, pt.size = 0.5, return_gg_object = T)[[1]]
plot_grobs(list(gg1, gg2))

proj$cluster_annotation = from_to(annot_vec = proj$Clusters_0.5, 
                                  old_new_map_vec = c(
                                    C1='Plasma cell',
                                    C2='undefined', 
                                    C3='Plasma cell',
                                    C4='undefined',
                                    C5='B cell', C6='B cell', 
                                    C7='NKT/Tgd cell',
                                    C8='ILC', 
                                    C9='ILC',
                                    C10='NK cell', 
                                    C11='CD4 T cell',
                                    C12='CD4 T cell', 
                                    C13='naive T cell', 
                                    C14='CD8 T cell', 
                                    C15='CD8 T cell',
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


dr_df = jj_get_reduction_coords(proj)
jj_plot_features(reduction=dr_df, meta_features='cluster_annotation', label = F)
jj_plot_features(reduction=dr_df, meta_features='cluster_annotation', label = F, custom_colors = jj_get_colours(dr_df$cluster_annotation, '/abi/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'))

#transfer annotation from t nk ilc subset: 
#proj$cluster_annotation_fine = proj$cluster_annotation
#proj$cluster_annotation_fine[match(dr_df_t$cellNames, proj$cellNames)] = dr_df_t$cluster_annotation_fine

#add fine singler annoation from immgen
singler_pred = read_csv('/omics/groups/OE0436/internal/msimon/scATAC/imm_atlas_mouse_normal_gene_activity_mat_singler_ref_immgendata_label_own_predictions.csv')
singler_pred = singler_pred[match(rownames(dr_df), singler_pred$X1), ]
stopifnot(identical(singler_pred$X1, rownames(dr_df)))
proj$singler_label_fine = singler_pred$labels

dr_df$singler_label_fine = singler_pred$labels
summarise_fractions(dr_df[, c('Clusters_0.5', 'singler_label_fine')], 'Clusters_0.5') %>%
  dplyr::filter(Clusters_0.5 == 'C7')

#add deviations for tisTreg signatures to cellColData
sig_se = getMatrixFromProject(proj, 'tisTregMatrix')
z_score_mat = t(assays(sig_se)[['z']])
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))
proj = add_df_to_cellcoldata(proj, z_score_mat)
dr_df = jj_get_reduction_coords(proj, redname = 'UMAP')
gg = plot_sigs(dr_df)

pdf(paste0(storeFigPath, 'mouse_atlas_tissue.pdf'), width = 6, height= 5)
  print(add_umap_arrows(gg[[1]] + labs(color=''), theme_use = theme_minimal2))
dev.off()
pdf(paste0(storeFigPath, 'mouse_atlas_tisTreg_signatures.pdf'), width = 8, height= 5)
  print(plot_grobs(gg[2:7], grobs=6, cols=3))
dev.off()

#saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets_corrected', load=T, overwrite = T, dropCells = T)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected', pconfig$DONOR))

### make bigwigs für Clusters_0.5
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features='Clusters_0.5', label = F, pt.size = 0.5)
res = getGroupBW(
  ArchRProj = proj,
  groupBy = "Clusters_0.5",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 1000,
  ceiling = 4
)

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tissue.pdf'),  width = 10, height = 8)
jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F, pt.size = 0.5, 
                 custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                 custom_theme = theme_minimal()) 

dev.off()


# update cluster annotation -----------------------------------------------

#update the whole atlas with fine grained annotations from the T/NK/ILC subset
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
dr_df = jj_get_reduction_coords(proj, 'UMAP')

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected', pconfig$DONOR))
dr_df = jj_get_reduction_coords(proj, 'UMAP')

proj$cluster_annotation_level1 = proj$cluster_annotation_fine
proj$cluster_annotation_level1[match(dr_df$cell_id, proj$cell_id)] = dr_df$cluster_annotation_level1

proj$cluster_annotation_level2 = proj$cluster_annotation_fine
proj$cluster_annotation_level2[match(dr_df$cell_id, proj$cell_id)] = dr_df$cluster_annotation_level2

proj$cluster_annotation_level3 = proj$cluster_annotation_fine
proj$cluster_annotation_level3[match(dr_df$cell_id, proj$cell_id)] = dr_df$cluster_annotation_level3

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features = c('cluster_annotation_level1', 'cluster_annotation_level2', 'cluster_annotation_level3'))

proj$cluster_annotation_level1[proj$cluster_annotation_level1 == 'undef'] = 'undefined'
proj$cluster_annotation_level2[proj$cluster_annotation_level2 == 'undef'] = 'undefined'
proj$cluster_annotation_level3[proj$cluster_annotation_level3 == 'undef'] = 'undefined'

#saveArchRProject(proj)

vec = seq_along(unique(dr_df$cluster_annotation_level3))
myColors = c("darkblue", "lightblue")
myRangeFunction = colorRampPalette(myColors)
myColorRange = myRangeFunction(max(vec))
myColors = myColorRange[vec]
names(myColors) = unique(dr_df$cluster_annotation_level3)

col_map = jj_get_colours(dr_df$cluster_annotation_level3, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

pdf(paste0(storeFigPath, 'mouse_atlas_umap.pdf'), width = 12, height=10)
  jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation_level3', custom_colors = col_map)
  jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
  jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
dev.off()



# take b cell subset: cluster 2 could also be ILC, cluster 4 is also suspici --------

proj_b = proj[proj$Clusters_0.5 %in% paste0('C',1:6), ]
proj_b = archr_dim_reduction(proj_b, cluster_res = NULL)
#saveArchRProject(proj_b, outputDirectory='ArchRProject_B_subset', load=F)
dr_df = jj_get_reduction_coords(proj_b, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)

gg = plot_sigs(dr_df)

pdf(paste0(storeFigPath, 'mouse_atlas_b_tissue.pdf'), width = 6, height= 5)
print(add_umap_arrows(gg[[1]] + labs(color=''), theme_use = theme_minimal2))
dev.off()
pdf(paste0(storeFigPath, 'mouse_atlas_b_tisTreg_signatures.pdf'), width = 8, height= 5)
print(plot_grobs(gg[2:7], grobs=6, cols=3))
dev.off()



# take myeloid subset -----------------------------------------------------

proj_i = proj[proj$Clusters_0.5 %in% paste0('C',16:28), ]
proj_i = archr_dim_reduction(proj_i, cluster_res = NULL)
dr_df = jj_get_reduction_coords(proj_i, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
#saveArchRProject(proj_i, outputDirectory='ArchRProject_macrophage_dc_subset', load=F)

gg = plot_sigs(dr_df)

pdf(paste0(storeFigPath, 'mouse_atlas_i_tissue.pdf'), width = 6, height= 5)
print(add_umap_arrows(gg[[1]] + labs(color=''), theme_use = theme_minimal2))
dev.off()
pdf(paste0(storeFigPath, 'mouse_atlas_i_tisTreg_signatures.pdf'), width = 8, height= 5)
print(plot_grobs(gg[2:7], grobs=6, cols=3))
dev.off()



# take T/NK/ILC subset ----------------------------------------------------

proj_t = proj[proj$Clusters_0.5 %in% paste0('C',7:15), ]
proj_t = archr_dim_reduction(proj_t, cluster_res = NULL)
proj_t$Clusters_0.5_atlas = proj_t$Clusters_0.5
proj_t = archr_clustering(proj_t)
proj_t = addImputeWeights(proj_t)
dr_df = jj_get_reduction_coords(proj_t, redname = 'UMAP')

# sig_se = getMatrixFromProject(proj_t, 'tisTregMatrix')
# z_score_mat = t(assays(sig_se)[['z']])
# z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
# stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))
# proj_t = add_df_to_cellcoldata(proj_t, z_score_mat)

jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
jj_plot_features(reduction = dr_df, meta_features = 'core_tisTreg_sig', pt.size = 0.5, colorScale = 'viridis', cap_top = 'q95', cap_bottom = 'q05')
jj_plot_features(reduction = dr_df, meta_features = 'late_precursor_sig', pt.size = 0.5, colorScale = 'viridis', cap_top = 'q95', cap_bottom = 'q05')
jj_plot_features(reduction = dr_df, meta_features = c('Clusters_0.5_atlas'), pt.size = 0.5, colorScale = 'viridis', cap_top = 'q95', cap_bottom = 'q05')

gg = plot_sigs(dr_df)

pdf(paste0(storeFigPath, 'mouse_atlas_t_tissue.pdf'), width = 6, height= 5)
print(add_umap_arrows(gg[[1]] + labs(color=''), theme_use = theme_minimal2))
dev.off()
pdf(paste0(storeFigPath, 'mouse_atlas_t_tisTreg_signatures.pdf'), width = 8, height= 5)
print(plot_grobs(gg[2:7], grobs=6, cols=3))
dev.off()

dr_df = get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features='cluster_annotation_fine', label = T, box_col = 'white', pt.size = 0.5)

#saveArchRProject(proj_t, outputDirectory='ArchRProject_t_nk_ilc_subset_corrected', load=T, overwrite = T, dropCells = T)

proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features='cluster_annotation_fine', label = F, pt.size = 0.5)
jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F, pt.size = 0.5, custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))

marker_list = list(
  NK = c('Ncr1', 'Klrb1c', 'Eomes'),
  ILC1 = c('Kit','Ncr1', 'Tbx21', 'Irf8'),
  ILC2=c('Il4','Il5','Il13','Kit','Rora'),
  ILC3=c('Kit','Rorc','Il17a'),
  CD8_Tnaive=c('Ccr7','Sell','Cd8a','Cd8b1'),
  CD8_Teff = c('Cd8a','Cd8b1','Pdcd1','Tox'),
  CD8_Tmem= c('Cd8a','Cd8b1','Ccr7','Ifng'),
  CD4_Tnaive= c('Ccr7','Sell','Cd4'),
  Treg = c('Foxp3','Ikzf2','Il10'),
  tisTregST2 = c('Foxp3','Ikzf2','Areg','Klrg1','Il1rl1','Tigit'),
  pTreg = c('Foxp3','Rorc'),
  Th1 = c('Ifng','Tbx21','Cd4'),
  Th2 = c('Cd3e','Il4','Il1rl1'),
  Th17 = c('Cd3e','Rorc','Il17a'),
  Th_Il21 = c('Il21', 'Maf')
)
marker_genes = jj_list_to_df(marker_list)
marker_genes$name = mixed_relevel(marker_genes$name)
gene_mat = get_gene_mat(proj)

#cell type coarse annotation
proj$cluster_annotation_level1 = from_to(vec = proj$Clusters_1.2, 
                                             old_new_map_vec = c(
                                               C1='non_T', 
                                               C2='non_T', 
                                               C3='CD4_T', #mystery CD4 T cell
                                               C4='CD4_T', 
                                               C5='non_T',
                                               C6='non_T',  
                                               C7='T', 
                                               C8='CD4_T', 
                                               C9='CD4_T', 
                                               C10='non_T', 
                                               C11='non_T', 
                                               C12='non_T', 
                                               C13='T', 
                                               C14='CD4_T', 
                                               C15='CD8_T',
                                               C16='T', 
                                               C17='T', 
                                               C18='CD8_T',
                                               C19='CD4_T',
                                               C20='CD4_T', 
                                               C21='CD4_T',
                                               C22='CD4_T', # 
                                               C23='CD4_T',
                                               C24='CD4_T',
                                               C25='CD4_T', 
                                               C26='CD4_T',
                                               C27='CD8_T',
                                               C28='CD4_T',
                                               C29='CD4_T',
                                               C30='CD4_T'
                                             ))

#cell type fine annotation
proj$cluster_annotation_level2 = from_to(vec = proj$Clusters_1.2, 
                                       old_new_map_vec = c(
                                         C1='ILC2', 
                                         C2='ILC2', 
                                         C3='CD4_other', 
                                         C4='Th2', 
                                         C5='ILC2', 
                                         C6='ILC3',  
                                         C7='gd_T', 
                                         C8='Th17', 
                                         C9='Th17', 
                                         C10='NK', 
                                         C11='NK', 
                                         C12='ILC1', 
                                         C13='gd_T', 
                                         C14='NKT', 
                                         C15='CD8_Tmem',
                                         C16='NKT', 
                                         C17='undef', 
                                         C18='CD8_Tnaive',
                                         C19='CD4_Tnaive',
                                         C20='Treg', 
                                         C21='Treg',
                                         C22='Treg', 
                                         C23='Treg',
                                         C24='undef',
                                         C25='Treg', 
                                         C26='Th_Il21',
                                         C27='CD8_Teff',
                                         C28='Th1',
                                         C29='Th1',
                                         C30='Th1'
                                       ))

#tisTreg categorization
proj$cluster_annotation_level3 = from_to(vec = proj$Clusters_1.2, 
                                         old_new_map_vec = c(
                                           C1='ILC2', 
                                           C2='ILC2', 
                                           C3='CD4_other', 
                                           C4='Th2', 
                                           C5='ILC2',
                                           C6='ILC3',  
                                           C7='gd_T', 
                                           C8='Th17', 
                                           C9='Th17', 
                                           C10='NK', 
                                           C11='NK', 
                                           C12='ILC1', 
                                           C13='gd_T', 
                                           C14='NKT', 
                                           C15='CD8_Tmem',
                                           C16='NKT', 
                                           C17='undef',
                                           C18='CD8_Tnaive',
                                           C19='CD4_Tnaive',
                                           C20='tisTregST2', 
                                           C21='pTreg',
                                           C22='tisTregST2', 
                                           C23='tisTregST2',
                                           C24='undef',
                                           C25='Treg', 
                                           C26='Th_Il21',
                                           C27='CD8_Teff',
                                           C28='Th1',
                                           C29='Th1',
                                           C30='Th1'
                                         ))

dr_df = jj_get_reduction_coords(proj, 'UMAP')

dr_df$Clusters_1.2 = mixsort_factor(dr_df$Clusters_1.2)
dr_df$cluster_annotation_level1 = mixsort_factor(dr_df$cluster_annotation_level1)
dr_df$cluster_annotation_level2 = mixsort_factor(dr_df$cluster_annotation_level2)
dr_df$cluster_annotation_level3 = mixsort_factor(dr_df$cluster_annotation_level3)

feature_vec = structure(marker_genes$name, names=marker_genes$feature)
col_annot = HeatmapAnnotation(marker=feature_vec, col = list(marker=jj_get_jj_colours(levels(dr_df$cluster_annotation_level3))))
h1 = jj_feature_heatmap(gene_mat, features_use = marker_genes$feature, grouping_vector = dr_df$Clusters_1.2,
                   cluster_columns = F, cluster_rows = F, top_annot = col_annot,row_annot = row_annot) #scale is fine since it is column wise
h1_b = jj_feature_heatmap(gene_mat, features_use = unique(marker_genes$feature), grouping_vector = dr_df$Clusters_1.2,
                        cluster_columns = T, cluster_rows = T, row_annot = row_annot)

convert_color('darkblue')
c_annot_df = jj_summarize_dataframe(dplyr::select(dr_df, starts_with('cluster_annotation_level')), 
                                    dr_df$Clusters_1.2, method = 'mode')
colnames(c_annot_df) = gsub('cluster_annotation_', '', colnames(c_annot_df))
#row_annot = rowAnnotation(celltype_l1 = c_annot_df$clust, celltype_l2=
#                          col=list(celltype_l1=c('CD4_T'="blue",'CD8_T'="green",'gd_T'="purple",'non_T'='red','undef'='#CCCCCC')))
row_annot = HeatmapAnnotation(df = c_annot_df, which = 'row',
                              col = list(level1 = c('CD4_T'="blue",'CD8_T'="green",'T'='cyan', 'gd_T'="purple",'non_T'='red','undef'='#CCCCCC'),
                                         level2 = jj_get_jj_colours(levels(c_annot_df$level2)),
                                         level3 = jj_get_jj_colours(levels(c_annot_df$level3))))
h2 = jj_feature_heatmap(t(dr_df), features_use = modules_plot, grouping_vector = mixed_relevel(dr_df$Clusters_1.2),
                   cluster_columns = F, cluster_rows = T, row_annot = row_annot)

gglist = list()
gglist[1] = jj_plot_features(reduction = dr_df, meta_features = 'Clusters_1.2', return_gg_object = T)
gglist[2] = jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level1', return_gg_object = T, custom_colors = c('CD4_T'="blue",'CD8_T'="green",'T'='cyan','gd_T'="purple",'non_T'='red','undef'='#CCCCCC'))
gglist[3] = jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level2', return_gg_object = T, custom_colors = jj_get_jj_colours(levels(dr_df$cluster_annotation_level2)))
gglist[4] = jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level3', return_gg_object = T, custom_colors = jj_get_jj_colours(levels(dr_df$cluster_annotation_level3)))


proj = addModuleScore(proj, features = marker_list, useMatrix = 'GeneScoreMatrix')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
(modules_plot = grep('Module', colnames(dr_df), value = T))

pdf(paste0(storeFigPath, 'mouse_t_nk_ilc_marker_signatures.pdf'), width = 12, height=10)
jj_arrange_ggplots(gglist)
print(h1)
print(h1_b)
print(h2)
jj_arrange_ggplots(jj_plot_features(reduction = dr_df, meta_features = modules_plot, return_gg_object = T, cap_top = 'auto'),
                   nplots = 4, cols =2)
dev.off()

#differences
comparison_df = tibble::tribble(~name, ~g1, ~g2,
  'CD4_other_C3_vs_tisTreg_colon_C20', 'C3',  'C20',
  'CD4_other_C3_vs_tisTreg_all_C20_C22_C23', 'C3', 'C20,C22,C23',
  'ILC1_C10_C11_vs_NK_cells_C12', 'C10,C11', 'C12',
  'ILC2_C1_vs_ILC2_C2', 'C1', 'C2',
  'ILC2_C1_vs_tisTreg_colon_C20','C1','C20',
  'ILC2_C1_vs_tisTreg_all_C20_C22_C23', 'C1', 'C20,C22,C23',
  'ILC2_C2_vs_tisTreg_colon_C20','C2','C20',
  'ILC2_C2_vs_tisTreg_all_C20_C22_C23', 'C2', 'C20,C22,C23',
  'ILC_C1_vs_ILC_C12', 'C1', 'C12', 
  'ILC_C1_vs_ILC_C6', 'C1', 'C6',
  'ILC_C2_vs_ILC_C12', 'C2', 'C12',
  'ILC_C2_vs_ILC_C6', 'C2', 'C6', 
  'ILC_C6_vs_ILC_C12', 'C6', 'C12', 
  'ILC3_C6_vs_Th17_C8','C6', 'C8',
  'pTreg_colon_C20_vs_tisTreg_C21','C20','C21',
  'Il21_Th_C26_vs_C29', 'C26', 'C29',
  'Il21_Th_C26_vs_C8_C9', 'C26', 'C8,C9',
  'Mystery_C3_vs_C1', 'C3','C1',
  'Mystery_C3_vs_C8', 'C3','C8',
  'Mystery_C3_vs_C2', 'C3','C2',
  'Mystery_C3_vs_C4', 'C3','C4',
  'Mystery_C3_vs_C28_C29', 'C3','C28,C29')
#commonalities
comparison_df = tibble::tribble(~name, ~g1, ~g2,
                                'ILC2_C1_vs_ILC1_ILC3_NK_C12_C6_C10_C11', 'C1',  'C12,C6,C10,C11',
                                'CD4_other_C3_vs_Th17_Th1_C8_C9_C29_C30', 'C3', 'C8,C9,C28,C29,C30') #

comparison_df = tibble::tribble(~name, ~g1, ~g2,
                                'ILC2_C2_vs_ILC2_C1', 'C2',  'C1') #

comparison_df$name_simple = gsub(',','_',paste(comparison_df$g1, comparison_df$g2, sep='_vs_'))

### cluster comparisons, find markers
markers_list = list()
(use_matrix = c('PeakMatrix', 'GeneScoreMatrix')[1])
for(i in 1:nrow(comparison_df)){
  message(i)
  markers_list[[i]] <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = use_matrix, 
    useGroups = unlist(strsplit(pull(comparison_df[i,], g1), ',')), 
    bgdGroups = unlist(strsplit(pull(comparison_df[i,], g2), ',')),
    groupBy = "Clusters_1.2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
}

names(markers_list) = comparison_df$name

markers_list = list()
for(i in unique(proj$cluster_annotation_level0)){
  message(i)
  if(i == 'undefined') next
  markers_list[[i]] <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = use_matrix, 
    useGroups = i, 
    groupBy = "cluster_annotation_level0",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
}
marker_df_list = lapply(markers_list,  function(x) archr_get_markers_as_df(x, proj, annotate_closest_gene = F, cutOff = "FDR <= 0.01 & Log2FC >= 0.58"))
names(marker_df_list) = unique(proj$cluster_annotation_level0)
marker_df_list$undefined = NULL
sapply(marker_df_list, nrow)
#jj_save_excel(marker_df_list,file_name = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-08-08-update_0809/t_nk_ilc_markers.xlsx')

#fig1 marker peaks intersection
marker_id_list = lapply(marker_df_list, function(x) x$feature)
jj_plot_upsetr(marker_id_list)

marker_df_list$ILC %>% head


pgr = getPeakSet(proj)
peak_df = as.data.frame(unname(pgr))
peak_df$group = 'union peaks (254545)'
jj_plot_categorical_by_group(peak_df,
                             feature_column =  'peakType',
                             group_column = 'group', 
                             absolute_numbers = F,
                             flip_coordinates = T, add_text = T, text_size = 3) + 
  labs(y='n peaks', x = '')

rownames(peak_df) = with(peak_df, paste(seqnames, start, end, sep = '-'))
peak_df = peak_df[, c('nearestGene', 'distToGeneStart', 'nearestTSS', 'distToTSS', 'peakType')]
for(i in 1:length(marker_df_list)){
  markers_annotate = marker_df_list[[i]]
  rownames(markers_annotate) =  with(markers_annotate, paste(seqnames, start, end, sep = '-'))
  peak_annot_use = peak_df[match(rownames(markers_annotate), rownames(peak_df)), ]
  stopifnot(identical(rownames(markers_annotate), rownames(peak_annot_use)))
  marker_df_list[[i]] = cbind(markers_annotate, peak_annot_use)
  
}
for(i in 1:length(marker_df_list)){
  marker_df_list[[i]]$feature = rownames(marker_df_list[[i]])
  rownames(marker_df_list[[i]]) = NULL
}
#jj_save_excel(marker_df_list,file_name = "/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-11-16-percentage_signature/fig1_all_markers.xlsx")
#jj_save_excel(marker_df_list,file_name = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-08-08-update_0809/c2_vs_c1_marker_peaks.xlsx')

marker_df_list = jj_load_excel("/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-11-16-percentage_signature/fig1_all_markers.xlsx")
marker_gr_list = lapply(marker_df_list, function(x) convert_granges(x$feature))
marker_df_all = do.call(rbind, marker_df_list)
subset_df = marker_df_all %>% dplyr::select(comparison, peakType)
nmarkers = as.data.frame(table(subset_df$comparison)) %>% arrange(Freq) %>% pull(Var1) %>% as.character
subset_df$comparison = factor(subset_df$comparison, levels = nmarkers)
jj_plot_categorical_by_group(subset_df, 'peakType', 'comparison', flip_coordinates = T, absolute_numbers = T) + 
  labs(y='n peaks', x = 'Cell type')


heatfeatures = unique(marker_df_all$feature[marker_df_all$Log2FC > 2]) # 83319 marker peaks with 0.58, 60000 with 2
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df = dr_df[!dr_df$cluster_annotation_level0 == 'undefined', ]
pmat_use = pmat_use[, !proj$cluster_annotation_level0 == 'undefined']
#should not plot such big heatmaps (massive ram usage)
ht = jj_plot_heatmap(pmat_use,
                     features_use = heatfeatures,
                     group_vec = dr_df$cluster_annotation_level0, 
                     clustering_method_rows = "ward.D2",  clustering_method_columns = "ward.D2",
                     cluster_rows = T, cluster_columns = T,
                     show_row_names = T, show_column_names = F, show_row_dend=F, show_column_dend = F,
                     use_raster=T, raster_quality = 5, raster_device = 'png')
#png(paste0(storeFigPath, 'heatmap_all_markers.png'), res = 300, width = 720, height = 480 )
pdf(paste0(storeFigPath, 'heatmap_all_markers.pdf'), width = 12, height = 5)
draw(ht)
dev.off()

#marker_df_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-08-08-update_0809/t_nk_ilc_markers.xlsx')
marker_top100_list = lapply(marker_df_list, function(x) head(x$name, 100))

proj = ArchR::addModuleScore(proj,
                             useMatrix = 'GeneScoreMatrix', 
                             features = marker_top100_list)

dr_df = jj_get_reduction_coords(proj, 'UMAP')
modules_plot = grep('Module', colnames(dr_df), value = T)
pdf(paste0(storeFigPath, 'marker_gene_top100_modules.pdf'), width = 12, height=10)
jj_arrange_ggplots(jj_plot_features(reduction = dr_df, meta_features = modules_plot,
                                    return_gg_object = T, cap_top = 'auto'),
                   nplots = 6, cols =3)
dev.off()


##compare periphery for each celltype against spleen
proj$peripheral_spleen = ifelse(proj$Tissue == 'Spleen', 'Spleen', 'Peripheral')
for(i in unique(proj$cluster_annotation_level0)){
  message(i)
  proj_sub = proj[proj$cluster_annotation_level0 == i, ]
  markers_list[[i]] <- getMarkerFeatures(
    ArchRProj = proj_sub, 
    useMatrix = use_matrix, 
    useGroups = 'Peripheral', 
    groupBy = "peripheral_spleen",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
}


# plot peak signatures ----------------------------------------------------

marker_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-08-08-update_0809/t_nk_ilc_marker_peaks.xlsx')
marker_list2 = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-08-08-update_0809/t_nk_ilc2_signature_ilc3_signature_marker_peaks.xlsx')
markers_all_list = c(marker_list, marker_list2)
mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr_list = lapply(c(0.25, 0.5, 1, 1.5), function(x) StringToGRanges(mouse_cd4_tisTreg_df$feature[mouse_cd4_tisTreg_df$avg_logFC > x]))
names(mouse_cd4_tisTreg_gr_list) = paste0('scATAC_tisTreg_logFC', c(0.25, 0.5, 1, 1.5))
mouse_c11_tisTreg = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-03-16-scATAC_tisTreg_signature_subset/scATAC_tisTreg_signature_C11_subset.xlsx')
c11_tisTreg = lapply(mouse_c11_tisTreg, function(x) StringToGRanges(x$feature))
names(c11_tisTreg) = c('scATAC_tisTreg_C11', 'scATAC_tisTreg_generic')
for(i in seq_along(markers_all_list)){
  mdf = markers_all_list[[i]]
  #mdf = mdf[mdf$Log2FC>1,]
}
markers_all_gr = lapply(markers_all_list, function(x) StringToGRanges(with(x,paste(seqnames, start,end,sep='-'))))
markers_all_gr = c(mouse_cd4_tisTreg_gr_list,c11_tisTreg, markers_all_gr)

proj = archr_add_peak_signatures(proj, markers_all_gr, signature_name = 'sig')
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df[is.na(dr_df)] = 0
gg = jj_plot_features(reduction = dr_df,
                      pt.size = 0.5,
                      meta_features = grep('z_', colnames(dr_df), value=T),
                      cont_or_disc = 'c',
                      colorScale = 'viridis',
                      cap_top = 'q99', cap_bottom = 'q01',
                      custom_theme = theme_minimal(), return_gg_object = T)
pdf(paste0(storeFigPath, 'mouse_t_nk_ilc_signature_scores.pdf'), width = 12, height=10)
jj_arrange_ggplots(gg, 4,2)
dev.off()


# get Treg subset for a higher resolution ---------------------------------

#see 'manuscript_figures.R'

# -------------------------------------------------------------------------



 nbmarkers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C1_C2'), 
  bgdGroups = c('C3'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

table(proj$Clusters_1) %>% sort(., decreasing = T)

marker_df_list = list()

#NK cells (C1, C2) vs. ILC1 (C3)
proj$comparison = proj$Clusters_1
proj$comparison[proj$comparison %in% c('C1', 'C2')] = 'C1_C2'
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C1_C2'), 
  bgdGroups = c('C3'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C1_C2_vs_C3']] = get_markers_as_df(markers, proj)

#ILC3 (C5) vs. TH17 (C8, C7) 
proj$comparison = proj$Clusters_1
proj$comparison[proj$comparison %in% c('C7', 'C8')] = 'C7_C8'
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C5'), 
  bgdGroups = c('C7_C8'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C5_vs_C7_C8']] = get_markers_as_df(markers, proj)

#ILC2 (C13) vs ILC2 (C12)
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C13'), 
  bgdGroups = c('C12'),
  groupBy = "Clusters_1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C13_vs_C12']] = get_markers_as_df(markers, proj)

#ILC2 (C12) vs. CD4 (C14)
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C12'), 
  bgdGroups = c('C14'),
  groupBy = "Clusters_1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C12_vs_C14']] = get_markers_as_df(markers, proj)
jj_save_excel(marker_df_list, paste0(storeFigPath, 'ILC_vs_colon_tisTreg.xlsx'))

marker_peaks = lapply(marker_df_list, function(x) x$feature)
plot_upsetr(marker_peaks)

marker_peaks = lapply(marker_df_list, function(x) x$feature[x$Log2FC>0])
plot_upsetr(marker_peaks)
#ILC2 (C12) vs tisTregs (C11, C9, C10) ###### C10 sind eigentlich zwei Cluster (oben und unten)  - wird deutlich bei Areg und Klrg1, Areg+ und Klrg1+ sind tisTreg aus Colon, Areg- und Klrg1- sind converted Treg aus colon; bitte den tisTreg Anteil von C10 verwenden

#ab res 1.2 splitted C10 in die besagten Cluster auf
library(clustree)
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df = dr_df[, colnames(dr_df) %in% paste0('Clusters_', seq(0,2, 0.1))]
clustree(dr_df[,1:20], prefix = "Clusters_", show_axis=T)
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df$Clusters_1.2 = sorted_factor(dr_df$Clusters_1.2)
jj_plot_features(reduction=dr_df, meta_features = 'Clusters_1.2', label=T, box_col = 'white')


proj$comparison = proj$Clusters_1
proj$comparison[proj$comparison %in% c('C11', 'C9')] = 'C9_C10tTreg_C11'
proj$comparison[proj$Clusters_1.2 %in% c('C20')] = 'C9_C10tTreg_C11'

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C12'), 
  bgdGroups = c('C9_C10tTreg_C11'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C12_vs_C9_C10tTreg_C11']] = get_markers_as_df(markers, proj)

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C20'), 
  bgdGroups = c('C2'),
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C20_vs_C2']] = get_markers_as_df(markers, proj, cutOff = "FDR <= 0.05")
proj$comparison = proj$Clusters_1.2
proj$comparison[proj$comparison %in% c('C1', 'C3')] = 'C1_C3'

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C20'), 
  bgdGroups = c('C1_C3'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C20_vs_C1_C3']] = get_markers_as_df(markers, proj, cutOff = "FDR <= 0.05")

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C2'), 
  bgdGroups = c('C1_C3'),
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C2_vs_C1_C3']] = get_markers_as_df(markers, proj, cutOff = "FDR <= 0.05")

#Spannend wären auch die Gemeinsamkeiten der einzelnen Vergleiche, 
#z.B. was haben ILC2 (C12) und tisTregs (C11, C9, C10*) gemeinsam im Vergleich zum Rest.
proj$comparison = proj$Clusters_1
proj$comparison[proj$comparison %in% c('C11', 'C9', 'C12')] = 'C9_C10tTreg_C11_C12'
proj$comparison[proj$Clusters_1.2 %in% c('C20')] = 'C9_C10tTreg_C11_C12'

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C9_C10tTreg_C11_C12'), 
  bgdGroups = NULL,
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df_list[['C9_C10tTreg_C11_C12_vs_all']] = get_markers_as_df(markers, proj)

#jj_save_excel(list_of_df = marker_df_list, file_name = paste0(storeFigPath, 'mouse_t_nk_ilc_cluster_comparison.xlsx'))


signature_list=list(C1_C2_vs_C3=StringToGRanges(marker_df_list$C1_C2_vs_C3$feature), 
                    C5_vs_C7_C8 = StringToGRanges(marker_df_list$C5_vs_C7_C8$feature))
feat_avail = ArchR::getPeakSet(proj)
assay_sp_mat = assays(getMatrixFromProject(proj, 'PeakMatrix', binarize = T))[[1]]
rownames(assay_sp_mat) = paste(seqnames(feat_avail), ranges(feat_avail), sep='-')
library(BSgenome.Mmusculus.UCSC.mm10)
z_score_df = AddChromatinModuleOwn(sparse_peak_mat = assay_sp_mat, features_gr_list = signature_list, genome = BSgenome.Mmusculus.UCSC.mm10)
dr_df = dr_df %>% dplyr::select(!ends_with('sig')) %>% dplyr::left_join(rownames_to_column(z_score_df, 'cell_id'), by = 'cell_id')
plot_grobs(jj_plot_features(reduction=dr_df, 
                            meta_features = c('C1_C2_vs_C3','C5_vs_C7_C8'), 
                            colorScale = 'viridis', pt.size = 1,
                            cap_top = 'q99', cap_bottom = 'q01', 
                            return_gg_object = T))

### make bigwigs für Clusters_1.2
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
dr_df = get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features='Clusters_1.2', label = F, pt.size = 0.5)
res = getGroupBW(
  ArchRProj = proj,
  groupBy = "Clusters_1.2",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 1000,
  ceiling = 4
)

# derive new tisTreg signature --------------------------------------------

jj_plot_features(reduction=dr_df, meta_features='Clusters_1.2', label = T, pt.size = 0.5)
#use clusters 20,22,23 as tisTregST2 cells
#naive CD4: C19 
proj$comparison = proj$Clusters_1.2
proj$comparison[proj$comparison %in% c('C20', 'C22', 'C23')] = 'C20_22_23_tisTregST2'
proj$comparison[proj$comparison %in% c('C19')] = 'C19_naive_CD4'
proj$comparison[proj$comparison %in% c('C25')] = 'C25_naive_Treg'

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C20_22_23_tisTregST2'), 
  bgdGroups = 'C19_naive_CD4',
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = c('C20_22_23_tisTregST2'), 
  bgdGroups = 'C25_naive_Treg',
  groupBy = "comparison",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
marker_df = get_markers_as_df(markers, proj)
#write_csv(marker_df, paste0(storeFigPath, 'mouse_atlas_C20_22_23_tisTregST2_vs_C19_naive_CD4.csv'))
#write_csv(marker_df, paste0(storeFigPath, 'mouse_atlas_C20_22_23_tisTregST2_vs_C25_naive_Treg.csv'))

#marker_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-02-15-tisTreg_signature_new/mouse_atlas_C20_22_23_tisTregST2_vs_C19_naive_CD4.csv')
marker_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-03-16-scATAC_tisTreg_signature_subset/mouse_atlas_C20_22_23_tisTregST2_vs_C25_naive_Treg.csv')

volcano_plot_simple(marker_df,logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'symbol', use_text = F, labs_range = c(0, 10), pt.size = 0.5,
                    markers_highlight = c('Batf', 'Il1rl1', 'Areg', 'Ccr8','Pparg', 'Bach2'), col_by_highlight = T) + 
  labs(x='Log2FC', y='FDR') + guides(shape='none')


marker_gr = StringToGRanges(marker_df$feature)
marker_gr2 = StringToGRanges(marker_df2$feature)
olap1 = a_regions_that_have_minOverlap_with_b_regions(marker_gr, marker_gr2, minOverlap = 1e-10)
anyDuplicated(olap1$selectedRegions)
anyDuplicated(olap1$selectedCoordinates)
marker_df_keep = marker_df[match(olap1$selectedRegions, marker_df$feature), ]
identical(marker_df_keep$feature, olap1$selectedRegions)
marker_df2_keep = marker_df2[match(olap1$selectedCoordinates, marker_df2$feature), ]
identical(marker_df2_keep$feature, olap1$selectedCoordinates)
marker_combined_df = bind_cols(marker_df_keep, marker_df2_keep)
logfc_logfc_plot_simple(as.data.frame(marker_combined_df), 'Log2FC...2', 'Log2FC...10', 'symbol...5', colour_aes = NULL, labs_range = NULL)
#olap1 = olap1[!duplicated(olap1$selectedRegions), ]


marker_gr = StringToGRanges(marker_df$feature)

### compare with previous tisTreg signatures
mouse_cd4_tisTreg_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_tisTreg_df[mouse_cd4_tisTreg_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = StringToGRanges(mouse_cd4_tisTreg_df$feature)

olap = a_regions_that_have_minOverlap_with_b_regions(mouse_cd4_tisTreg_gr, marker_gr, minOverlap = 1e-10)
length(unique(olap$selectedRegions)) / length(mouse_cd4_tisTreg_gr) #4947/7655 0.64

olap = a_regions_that_have_minOverlap_with_b_regions(marker_gr, mouse_cd4_tisTreg_gr, minOverlap = 1e-10)
length(unique(olap$selectedRegions)) / length(marker_gr) # 6387/11504 0.55

core_sig_gr = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
core_sig_gr = GenomicRanges::makeGRangesFromDataFrame(core_sig_gr$core_tisTreg_sig)

olap = a_regions_that_have_minOverlap_with_b_regions(marker_gr, core_sig_gr, minOverlap = 1e-10)
length(unique(olap$selectedRegions)) / length(marker_gr) # 781/11504 0.067

olap = a_regions_that_have_minOverlap_with_b_regions(core_sig_gr, marker_gr, minOverlap = 1e-10)
length(unique(olap$selectedRegions)) / length(core_sig_gr) # 717/2267 0.316

#get midpoints of granges
gr.mid = function(x){
  start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
  return(x)
}
core_mid_gr = gr.mid(core_sig_gr)
mouse_cd4_tisTreg_mid_gr = gr.mid(mouse_cd4_tisTreg_gr)
marker_mid_gr = gr.mid(marker_gr)

m_core = m_sc = vector()
for(i in 1:10){
  m_core[i] = mean(as.data.frame(distanceToNearest(marker_mid_gr, core_mid_gr[sample(1:length(core_mid_gr), 1000, replace = F)]))$distance)
  m_sc[i] = mean(as.data.frame(distanceToNearest(marker_mid_gr, mouse_cd4_tisTreg_mid_gr[sample(1:length(mouse_cd4_tisTreg_mid_gr), 1000, replace = F)]))$distance)
}
mean(m_core) #1116809
mean(m_sc) #916846
916846/1116809

#dist_mat[dist_mat>1e6] = 1e6
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(0,1, 1e5, 5e5, 1e6), c("green","blue", "white","red", "grey"))

cmap_list = list()
for(chrom in seqlevels(marker_gr)){
  message(chrom)
  #distance heatmap
  core_use = core_mid_gr[seqnames(core_mid_gr) == chrom]
  strand(core_use) = '*'
  core_use = sort(core_use)
  
  mouse_use =  mouse_cd4_tisTreg_gr[seqnames(mouse_cd4_tisTreg_gr) == chrom]
  strand(mouse_use) = "*"
  mouse_use = sort(mouse_use)
  
  marker_use = marker_mid_gr[seqnames(marker_mid_gr) == chrom]
  strand(marker_use) = '*'
  marker_use = sort(marker_use)
  
  dist_core_mat = matrix(nrow=length(marker_use), ncol = length(core_use))
  for(i in seq_along(core_use)){
    dist_core_mat[, i] = GenomicRanges::distance(core_use[i], marker_use, select = 'all')
  }
  
  dist_sc_mat = matrix(nrow=length(marker_use), ncol = length(mouse_use))
  for(i in seq_along(mouse_use)){
    dist_sc_mat[, i] = GenomicRanges::distance(mouse_use[i], marker_use, select = 'all')
  }
  
  cmap_list[[chrom]] = ComplexHeatmap::Heatmap(dist_core_mat, name='distance', row_title = 'immune atlas tisTreg', column_title = 'core tisTreg', 
                          col = col_fun, cluster_rows = F, cluster_columns = F,  heatmap_legend_param = list(
                            at = c(0, 1, 1e5,5e5, 1e6), labels = c('0 (green)', '1 (blue)', '1e5 (white)', '5e5 (red)', '1e6+ (grey)'))) + 
    ComplexHeatmap::Heatmap(dist_sc_mat, name='distance', row_title = 'immune atlas tisTreg', column_title = 'mouse scATAC Cd4 tisTreg', 
                            col = col_fun, cluster_rows = F, cluster_columns = F, heatmap_legend_param = list(
                              at = c(0, 1, 1e5,5e5, 1e6), labels = c('0 (green)', '1 (blue)', '1e5 (white)', '5e5 (red)', '1e6+ (grey)')))
}

pdf(paste0(storeFigPath, 'tisTreg_overlap_heatmaps.pdf'), width = 10)
for(i in seq_along(cmap_list)){
  draw(cmap_list[[i]], column_title=seqlevels(marker_gr)[i], column_title_gp = gpar(fontsize = 16))
}
dev.off()

### genomic ranges heatmap
library(circlize)
library(GenomicRanges)
chr_df = read.chromInfo(species='mm10')$df
chr_df = chr_df[!chr_df$chr %in% paste0("chrY"), ]
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_gr
library(EnrichedHeatmap)
chr_window = makeWindows(chr_gr, w = 1e4)
chr_window

bed1 = generateRandomBed(nr = 1000, nc = 10) # generateRandomBed() is from circlize package
# convert to a GRanes object
gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))

num_mat = average_in_window(window = chr_window, gr = gr1, v = bed1[, -(1:3)])

calc_wighted_sum = function(x, w, gw) sum(x*w/gw)
calc_sum = function(x, w, gw) sum(x)

nmat1 = average_in_window(window = chr_window, gr = core_sig_gr, method=calc_sum)
nmat2 = average_in_window(window = chr_window, gr = mouse_cd4_tisTreg_gr, method=calc_sum)
nmat3 = average_in_window(window = chr_window, gr = marker_gr, method=calc_sum)
num_mat = cbind(nmat1, nmat2, nmat3)
colnames(num_mat) = c('core tisTreg', 'scATAC Cd4 tisTreg', 'Immune atlas tisTreg')

chr = as.vector(seqnames(chr_window))
chr_level = gtools::mixedsort(unique(chr))
chr = factor(chr, levels = chr_level)

Heatmap(t(num_mat), name = "mat", col = colorRamp2(c(0, 0.1, 0.2), c("darkblue", "red", "yellow")),
        column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE,na_col = 'white',
        #row_split = subgroup, cluster_row_slices = FALSE,
        row_title = "tisTregST2 signatures",
        #left_annotation = rowAnnotation(subgroup = subgroup, show_annotation_name = FALSE,
        #                                annotation_legend_param = list(
        #                                  subgroup = list(direction = "horizontal", title_position = "lefttop", nrow = 1))),
        column_title_gp = gpar(fontsize = 10), border = TRUE,
        column_gap = unit(4, "points"),
        column_title = ifelse(1:20 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
        #heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop")
)

chr_show = 'chrX'
mat_use = t(num_mat[as.character(seqnames(chr_window))==chr_show, ])
Heatmap(mat_use, 
        name = "peaks overlapping",
        col = colorRamp2(c(0, 0.5*max(mat_use, na.rm = T), max(mat_use, na.rm = T)), c("darkblue", "red", "yellow")),
        cluster_columns = FALSE, show_row_dend = FALSE, na_col = 'white',
        row_title = "tisTregST2 signatures",
                  column_title_gp = gpar(fontsize = 10), border = TRUE,
                  column_title = chr_show
)


mouse_cd4_tisTreg_df

marker_df$symbol


# treg subset analysis ----------------------------------------------------
#use T/NK/ILC subset (ArchRProject_t_nk_ilc_subset_corrected)
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
proj = proj[proj$Clusters_1 %in% c('C9', 'C10','C11','C18', 'C19'), ]  #C18 is naive CD4 + naive Treg
proj = archr_dim_reduction(proj)

dr_df = get_reduction_coords(proj, 'UMAP')

jj_plot_features(reduction=dr_df, meta_features='Clusters_0.5', label = T, pt.size = 0.5, return_gg_object = T)
jj_plot_features(reduction=dr_df, meta_features='Tissue', pt.size = 0.5, custom_colors = jj_get_colours(dr_df$Tissue, colour_csv = '/abi/data2/simonma/projects/scATAC/scripts/colour_map.csv'))

gg = plot_sigs(dr_df, features_plot = c("fat_treg_sig","skin_treg_sig", "core_tisTreg_sig","LNnaiveTreg_sig","early_precursor_sig", "late_precursor_sig"))
plot_grobs(gg[-1], grobs = 6, cols = 3)

proj = addImputeWeights(proj) #needs to be recalculated when cells are dropped
gene_mat = get_gene_mat(proj)
features_plot = c('Cd4','Foxp3','Sell','Klrg1', 'Nfil3','Ccr8')
gg = archr_plot_markers(proj, features_plot, use_magic = F)


# goi_mat = t(as.matrix(msGetGOIRNAMat(gene_mat, goi = features_plot, normFold = F)))
# dr_df = cbind(dr_df, goi_mat)
# colnames(dr_df)[duplicated(colnames(dr_df))]
# gg = jj_plot_features(reduction=dr_df, meta_features=features_plot,
#                       cap_top = 'auto', pt.size = 1, return_gg_object = T, colorScale = 'viridis')
plot_grobs(gg, grobs = 6, cols = 3)

annot_df = dr_df %>% dplyr::select(Clusters_1, Tissue) %>% 
  dplyr::group_by_all() %>% dplyr::summarise(n=n()) %>% 
  dplyr::ungroup() %>% dplyr::group_by(Clusters_1) %>%  dplyr::filter(n == max(n)) %>%
  dplyr::rename(Group=Tissue, cluster=Clusters_1) %>% dplyr::select(cluster, Group)
annot_df = data.frame(Group=c('tisTregST2', 'tisTregST2', 'pTreg', 'tisTregST2', 'tisTregST2 progenitor', 'naive'), cluster=gtools::mixedsort(unique(dr_df$Clusters_1)))
ggplot_dotplot(gene_mat, features_plot, dr_df$Clusters_1, group_annot_df = annot_df)

ggplot_dotplot(gene_mat, features_plot, 
               dr_df$Clusters_1,
               group_annot_df = annot_df, 
               column_colors = jj_get_colours(dr_df$Tissue, colour_csv = '/abi/data2/simonma/projects/scATAC/scripts/colour_map.csv'))

#saveArchRProject(proj, outputDirectory='ArchRProject_treg_corrected', load=T, overwrite = T, dropCells = T)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR,'ArchRProject_treg_corrected'))

#26.01. more genes and shiny app cluster dotplot
gene_mat = get_gene_mat(proj)

marker_genes = c('Eomes','ICOS','IL5','IL13','IL9','IL23a','IL23R','IL11','IL32','CXCR5',
                 'CXCL13','CCR5','CCR3','CCD4','CXCR4','CCR6','CSF2','KIT','AHR','IL7R',
                 'LTA','LTB','Thy1','TGFB1')
marker_genes <- switch(which_case(rownames(gene_mat)),
                       "lowercase"=tolower(marker_genes),
                       "titlecase"=tools::toTitleCase(tolower(marker_genes)),
                       "uppercase"=toupper(marker_genes))
marker_genes[!marker_genes %in% rownames(gene_mat)]
marker_genes = marker_genes[marker_genes %in% rownames(gene_mat)]

dr_df <- get_reduction_coords(proj, 'UMAP')
dr_df = bind_features_with_dr_df(gene_mat, features = marker_genes, dr_df = dr_df)
gg = plotFeatures(seurat_obj = NULL, 
                  reduction = dr_df,
                  pt.size = 0.5,
                  meta_features = marker_genes, #my_title = cell_type, #cell_type,
                  cont_or_disc = 'c', colorScale = 'viridis', topCap = 'auto', order=F,
                  custom_theme = theme_minimal(), return_gg_object = T)

#put the gene label to the title and remove gene name next to color scale
gg = lapply(gg, function(x) x + labs(color='', title=x$labels$colour) + theme_minimal(base_size = 8))

pdf(paste0(storeFigPath, dataset_use, 't_nk_ilc_subset_markers.pdf'), width=12, height = 10)
  jj_plot_features(reduction=dr_df,
                   meta_features=c('Clusters_1'), 
                   cap_top = 'auto', 
                   custom_colors = jj_get_jj_colours(gtools::mixedsort(unique(dr_df$Clusters_1))),
                   pt.size = 0.5, label = T, box_col = 'white')
  ggplot_dotplot(gene_mat, marker_genes, proj$Clusters_1, scale_data = T)
  plot_grobs(gg, grobs = 9, cols=3)
dev.off()
lowerRes(paste0(storeFigPath, dataset_use, 't_nk_ilc_subset_markers.pdf'), remove_original=T)

fr_df = get_fraction_expressing(gene_mat, genes_plot = rownames(gene_mat), proj$Clusters_1)
mean_df = get_mean_expressing(gene_mat, genes_plot = rownames(gene_mat), proj$Clusters_1, scale_data = F)
mean_df2 = get_mean_expressing(gene_mat, genes_plot = rownames(gene_mat), proj$Clusters_1, scale_data = T)
mean_df = mean_df %>% 
  dplyr::left_join(dplyr::rename(mean_df2, scaled_count=count), by = c('Gene', 'cluster')) %>% 
  dplyr::left_join(fr_df, by = c('Gene', 'cluster'))
#write_csv(mean_df, '/abi/data2/simonma/projects/imm_cell_atlas/analysis/mouse_t_nk_ilc_subset_dotplot_counts.csv')



# chromvar ----------------------------------------------------------------
proj = proj_t
if(pconfig$GENOME=='mm10'){
  signature_list = read_rds('/abi/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
  signature_list = signature_list[-c(4)]
  signature_list = lapply(signature_list, makeGRangesFromDataFrame)
  proj = archr_add_peak_signatures(proj, signature_list, signature_name = 'tisTreg_mouse1')
  dr_df = get_reduction_coords(proj, 'UMAP')
  gg = jj_plot_features(reduction = dr_df,
                        pt.size = 0.5,
                        meta_features = grep('z_', colnames(dr_df), value=T),
                        cont_or_disc = 'c',
                        colorScale = 'viridis',
                        cap_top = 'q95', cap_bottom = 'q05',
                        custom_theme = theme_minimal(), return_gg_object = T)
  
  pdf(paste0(storeFigPath, 'tisTreg_signatures.pdf'), width=12, height=8)
  plot_grobs(gg, grobs = 4, cols = 2)
  dev.off()
  #early and late precursor sig cause error 
  # Error in Matrix::sparseMatrix(i = queryHits(overlapRegions), j = match(names(allPositions),  : 
  #                                                                          NA's in (i,j) are not allowed
  #getPeakAnnotation(proj, 'tisTreg')
  #stores list of 3 anntoations: 
  #$Name: name of the signatures
  #$Positions: GRangesList with Granges for each signature
  #$Matches: SummarizedExperiment with logical sparse assay 'matches' that contains cells in rows and signatures in columns (True = Match of peak and signature)
}else{
  signature_list_hg19 = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-04-30-tconv_signature_substraction/tisTreg_vs_Tconv_filtering_plot_signatures.csv')
  signature_list_hg19 = apply(signature_list_hg19, 2, function(x) StringToGRanges(na.omit(x)))
  cons_sig = readLines('/abi/data2/simonma/projects/scATAC/analysis/2021-01-24-human_hcc/tis_treg_conserved_peaks643_human.txt')
  signature_list_hg19 = list(nr1 = signature_list_hg19$nr1, conserved_tisTreg_sig = StringToGRanges(cons_sig))
  #liftover from hg19 to hg38
  #http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
  chain = rtracklayer::import.chain('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/hg19ToHg38.over.chain')
  signature_list = list()
  for(i in seq_along(signature_list_hg19)){
    lift_res = rtracklayer::liftOver(signature_list_hg19[[i]], chain = chain)  
    lift_res = lift_res[sapply(lift_res, length) == 1] #exclude regions with more than 1 match per original region
    signature_list[[i]] = do.call(c, lift_res)
  }
  names(signature_list) = names(signature_list_hg19)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = signature_list, name = "tisTreg", force=T)
  
  proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "tisTreg",
    force = TRUE
  )
  
  plotVarDev <- getVarDeviations(proj, plot = TRUE, name = "tisTregMatrix")
  plotVarDev
  
  signatures <- getFeatures(proj, useMatrix = "tisTregMatrix")
  signatures <- sort(grep("z:", signatures, value = TRUE))
  signatures
  
  sig_se = getMatrixFromProject(proj, 'tisTregMatrix')
  z_score_mat = t(assays(sig_se)[['z']])
  dr_df = get_reduction_coords(proj, redname='UMAP')
  dr_df = dr_df %>% dplyr::select(-ends_with('sig')) %>% 
    dplyr::left_join(rownames_to_column(as.data.frame(z_score_mat), 'cellNames'), by = 'cellNames')
  rownames(dr_df) = dr_df$cellNames
  
  #plot(dr_df$colon_treg_sig.x, dr_df$colon_treg_sig.y)
  features_plot = colnames(z_score_mat)
  gg = plotFeatures(seurat_obj = NULL, 
                    reduction = dr_df,
                    pt.size = 0.5,
                    meta_features = features_plot,
                    cont_or_disc = 'c', 
                    colorScale = 'viridis', 
                    topCap = 'auto',
                    custom_theme = theme_minimal(), return_gg_object = T)
  
  pdf(paste0('/abi/data2/simonma/projects/imm_cell_atlas/analysis/2021-10-25-human_normal_donor11/tisTreg_signatures.pdf'), width=12, height=8)
  gg
  dev.off()
}


# plot markers ------------------------------------------------------------

marker_genes = sort(unique(c(
  unlist(list(
    Tcells= c('Cd3e','Cd247','Cd3d','Cd4','Cd8A','Cd8B1'), #CD3z=CD247
    Tregs = c('FoxP3','Il2ra','Ikzf2'),
    TisTregs=c('Il1rl1','Klrg1','Areg','Tox'),
    Bcells=c('Cd19','Ms4a1'),
    Macrophages=c('Itgam','Mrc1','Spi1','Csf1r'),
    DCs=c('Itgax','Siglech','Bst2','Xcr1','Clec9A','Sirpa','Irf8'),
    NK=c('Klrb1c','Ncr1'),
    ILCs=c('Il7r','Il2ra'), #(negative for CD3e, cd19…)
    Mastcells_Basophil=c('Itgam','Fcer1a','Hnmt','Areg') #fcer1 replaced with fcer1a
  )), 
    'Tbx21', 'Rorc',  'Rora',  'Gata3',  'Irf4',  'Maf',  'Bcl6',  'Tob2',  'Il4',  'Il2',  'Ifng',  'Ifna1',  'Ifnb1',  'Il10',
    'Tnf',  'Gzmb',  'Il21',  'Il22',  'Il17a',  'Batf',  'Il27',  'IL4I1',  'IDO1',  'Pdcd1',  'Tigit',  'CD274',  'Cd40',
    'CD40lg',  'Lilrb4a',  'Lgals1',  'Il6',  'Il18bp',  'Lgals9',  'Il7',  'Il33', 'Nfil3',
    'Arg1','Marco','Nos2','Pdcd1lg2','Ccr7','Ccr2','Itgae','Cd207','Epcam',
    'Cd209a','Ly75','Ly6c1','Ly6g','Ptprc','Fcgr1'#,'H2-Ab1'
  )))


res = archr_plot_markers(proj = proj, marker_genes = marker_genes)
pdf(paste0(storeFigPath, dataset_use, '_markers_magic.pdf'), width=12, height = 10)
  plot_grobs(res, grobs = 9, cols = 3)
dev.off()
lowerRes(paste0(storeFigPath, dataset_use, '_markers_magic.pdf'), remove_original = T)

gene_mat = ArchR::getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
rnames = rowData(gene_mat)
gene_mat = assays(gene_mat)[[1]]
rownames(gene_mat) = rnames$name
dr_df = get_reduction_coords(proj, redname='UMAP')
gene_mat = gene_mat[, match(rownames(dr_df), colnames(gene_mat))]
stopifnot(identical(rownames(dr_df), colnames(gene_mat)))
names_match(dr_df, gene_mat)


marker_genes = c('Foxp3', 'Il1rl1', 'Tigit', 'Pdcd1', 'Il2ra')
marker_genes = c('Nt5e', 'Cd40') #c('Cd19','Sdc1','Prdm1','Cd1d1','Havcr2','Pdcd1')

#marker_df = read_csv(pconfig$CELLTYPE_MARKERS)
#marker_genes = marker_df$`official gene symbol`
marker_genes <- switch(which_case(rownames(gene_mat)),
                       "lowercase"=tolower(marker_genes),
                       "titlecase"=tools::toTitleCase(tolower(marker_genes)),
                       "uppercase"=toupper(marker_genes))
marker_genes_avail = marker_genes %in% rownames(gene_mat)
marker_genes[!marker_genes_avail]
grep('Tim', rownames(gene_mat), value=T, ignore.case = T)

#without magic:
dr_df = jj_bind_features_with_dr_df(gene_mat, features = marker_genes, dr_df = dr_df)

deduplicate_df = function(df){
  df = df[, !duplicated(colnames(df))]
  return(df)
}
dr_df = deduplicate_df(dr_df)
marker_genes = gsub('-', '_', marker_genes)
jj_plot_features(seurat_obj = NULL, 
                  reduction = dr_df,
                  pt.size = 0.5,
                  meta_features = marker_genes, #'Foxp3',
                  cont_or_disc = 'c', colorScale = 'wbr', cap_top = 'auto', order=F,
                  custom_theme = theme_minimal())
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", useGroups = c('C23'), #if bgdGroups is NULL, all other clusters are taken as background
  groupBy = "Clusters_1.1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markers, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList$C23

jj_plot_features(seurat_obj = NULL, 
                 reduction = dr_df,
                 pt.size = 0.5,
                 meta_features = 'Clusters_1.1', my_title = '', #cell_type,
                 cont_or_disc = 'd', colorScale = 'wbr', cap_top = 'auto', order=F,
                 custom_theme = theme_minimal())
jj_plot_features(seurat_obj = NULL, 
                 reduction = dr_df,
                 pt.size = 0.5,
                 meta_features = 'singler_label', my_title = '', #cell_type,
                 cont_or_disc = 'd', colorScale = 'wbr', cap_top = 'auto', order=F,
                 custom_theme = theme_minimal())


# tissue signatures -------------------------------------------------------

'/abi/data2/simonma/projects/scATAC/analysis/2021-12-02-thesis_plots/mouse_normal_cd4_skin_vs_all_peaks.csv'


# find markers ------------------------------------------------------------

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", useGroups = c('C11'), #if bgdGroups is NULL, all other clusters are taken as background
  groupBy = "Clusters_0.5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


# cd8 subset analysis -----------------------------------------------------
proj = proj[proj$Clusters_1 %in% c('C17', 'C21','C27'), ] #c21 could be mix of cd4/cd8 effectors, c25?



# other plots -------------------------------------------------------------

pmat_se = getMatrixFromProject(proj, 'PeakMatrix', binarize = T)
pset = getPeakSet(proj)
sort(table(gsub('C[0-9]+\\._\\.','', pset$GroupReplicate)))
pmat = assays(pmat_se)[[1]]
pmat = pmat[, match(proj$cellNames, colnames(pmat))]
stopifnot(identical(colnames(pmat), proj$cellNames))

dim(pmat)
max(pmat)
dr_df = get_reduction_coords(proj)
dr_df$n_peaks = Matrix::colSums(pmat)

ggplot(dr_df, aes(x=nFrags, y=n_peaks, color=Sample)) + geom_point() +
  scale_color_manual(values=msPickSampleColors(dr_df$Sample)) #+ facet_wrap(.~Sample)


# mouse normal CD4 annotation ---------------------------------------------

#mouse normal cd4
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal_cd4.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_mouse_normal_CD4'))

#upper part of '6' positive for Rorc
proj$annotation = from_to(proj$snn_harmony_res.1.7, old_new_map_vec = c(
  '0' = 'Treg_naive',
  '1' = 'CD4_Tnaive', #lung
  '2' = 'Tconv', #potentially more naive Th1 (or memory Th1), 4th highest th1 score, naive label
  '3' = 'Treg_naive',
  '4' = 'CD4_Tnaive',
  '5' = 'Treg_naive', #not as naive as 0,3,14
  '6' = 'Th1', #colon Th1?, corresponds to C20 in imm atlas
  '7' = 'CD4_Tnaive',
  '8' = 'CD4_Tnaive',
  '9' = 'CD4_Tnaive', #colon naive T cells
  '10' = 'tisTregST2', #colon, but only upper left part
  '11' = 'tisTregST2_prog',
  '12' = 'Th1',#skin 
  '13' = 'CD4_Tnaive', #fat, probably two fat naive clusters because of two different pools of mice (1-5 and 11-14)
  '14' = 'Treg_naive', #potentially homotypic doublets
  '15' = 'CD4_Tnaive',
  '16' = 'tisTregST2', #skin
  '17' = 'Th1', #fat
  '18' = 'tisTregST2_prog', #Nfil3,Klrg1,Foxp3 positive as Cluster 11
  '19' = 'pTreg', #Rorc+Foxp3
  '20' = 'Th17', #Rorc, maps correctly in mouse atlas
  '21' = 'Th2', #highest Th2 score and maps to ILC2 cells in mouse atlas, highest Gata3
  '22' = 'CD4_Tnaive', #present in lung, naive Treg and naive Tconv, equivalent to cluster 17 in imm cell atlas
  '23' = 'tisTregST2', #VAT
  '24' = 'Tconv', #between fat 13 (naive T) and fat 2 (Tconv)
  '25' = 'undef' #potential myeloid interacting cells
))

#lower part: doublet? has high count
proj$annotation[proj$Clusters_0.5 == 'C15'] = 'undef'

dr_df = jj_get_reduction_coords(proj, 'UMAP')
pheno_df = read_rds(pick_content('mouse_normal_CD4', 'seurat_pheno'))
dr_df = dr_df %>% dplyr::left_join(pheno_df[, c(1,2,3)], by = 'cells') %>% dplyr::select(UMAP_1, UMAP_2, everything())

jj_plot_features(reduction=dr_df, meta_features = 'annotation')

proj_use = proj[proj$annotation %in% c('Treg_naive', 'tisTregST2_prog', 'tisTregST2'), ]
proj_use = archr_dim_reduction(proj_use)
proj_use = archr_clustering(proj_use, 0.5)
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features = 'Clusters_0.5')
# proj_use = addMonocleTrajectory(proj_use, groupBy = 'annotation', useGroups = c('Treg_naive', 'tisTregST2_prog', 'tisTregST2'))
# so = CreateSeuratObject(counts = feat_mat, meta.data = meta_df)
# proj_use <- addTrajectory(
#   ArchRProj = proj_use, 
#   name = "tisTregST2", 
#   groupBy = "annotation",
#   trajectory =  c('Treg_naive', 'tisTregST2_prog', 'tisTregST2'), 
#   embedding = "UMAP", 
#   force = TRUE
# )
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
gg = jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'), return_gg_object = T)[[1]]
cont_df = jj_get_contour_lines(dr_df, 'annotation', .95)
pdf(paste0(storeFigPath, 'mouse_cd4_treg_subset_annotation_umap.pdf'), width =8, height=6)
gg + geom_path(data = cont_df, aes(x=x, y=y, group=cont_group, linetype = annotation), #, colour=snn_harmony_res.1), 
               size=1)
dev.off()

gmat = get_gene_mat(proj_use)
library(monocle3)
so = CreateSeuratObject(counts = gmat, meta.data = dr_df)
so = FindVariableFeatures(so) 
so = ScaleData(so, features = VariableFeatures(so))
gmat_use = GetAssayData(so, slot = 'scale.data')
#gmat_use <- gmat[rownames(gmat) %in% VariableFeatures(so), ]
pd_df <- droplevels(dr_df)
fd_df <- data.frame(gene_short_name=rownames(gmat_use),
                    row.names = rownames(gmat_use))
rownames(gmat_use) <- fd_df$gene_short_name

#MAKE SURE monocle3 is loaded and monocle is unloaded! Otherwise error later in plot_cells
cds <- new_cell_data_set(gmat_use,
                         cell_metadata = pd_df,
                         gene_metadata = fd_df)

#performs log-normalization + PCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap with cosine metric
cds <- reduce_dimension(cds)
#cds <- partitionCells(cds)

reducedDims(cds)$UMAP <- dr_df[, 1:2]
#colnames(reducedDims(cds)$UMAP) <- c('V1', 'V2') #otherwise error in orderCells plotly picker
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = T)
cds <- order_cells(cds)#, root_cells = 'AAACGAAAGACTAATG-7') #root pickec at -5, -3
cds@colData$sample_name = NULL
#write_rds(cds, paste0(bigFilesDir, 'mouse_cd4_treg_subset_monocle3_cds.RDS'))
cds = read_rds(paste0(bigFilesDir,  'mouse_cd4_treg_subset_monocle3_cds.RDS'))
pdf(paste0(storeFigPath, 'mouse_cd4_treg_subset_trajectory_umap.pdf'), width =10, height=6)
  monocle3::plot_cells(cds, color_cells_by = 'pseudotime', label_leaves = F, label_branch_points = F, label_roots = F, cell_size = 0.5) + theme_minimal()
dev.off()

identical(colnames(cds), proj_use$cellNames)
proj_use$pseudotime = pseudotime(cds, reduction_method = 'UMAP')
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
#jj_plot_features(reduction = dr_df, meta_features = 'pseudotime')
#ggplot(dr_df, aes(x = pseudotime, y = core_tisTreg_sig, colour=annotation)) + geom_point(size = 0.5)

pmat = get_peak_mat(proj_use)
mouse_cd4_tisTreg_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-10-31-tisTreg_sig_subsets/mouse_scATAC_tisTreg_sig_subsets.csv')
mouse_cd4_tisTreg_df = mouse_cd4_tisTreg_df[mouse_cd4_tisTreg_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature[])
emb_df = getReducedDims(proj_use)
polap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, #dr_df[, 1:2],
                                  signature_gr = mouse_cd4_tisTreg_gr, nFrags_vec = proj_use$nFrags, count_thres = 2e5, verbose = T)
proj_use$signature_pct_overlap = polap_df$signature_pct_overlap
proj_use$dataset_max_pct_overlap = polap_df$dataset_max_pct_overlap
proj_use$signature_n_overlap = polap_df$signature_n_overlap
proj_use$n_cells_aggregated = polap_df$n_cells_aggregated
proj_use$nFrags_aggregated = polap_df$nFrags_aggregated
#saveArchRProject(proj_use, 'mouse_cd4_treg_subset_monocle3_proj.RDS')
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
pdf(paste0(storeFigPath, 'mouse_cd4_treg_signature_overlap_umap.pdf'), width =10, height=6)
  jj_plot_features(reduction = dr_df, meta_features = 'signature_pct_overlap')
dev.off()
pdf(paste0(storeFigPath, 'mouse_cd4_treg_correlation.pdf'), width =8, height=6)
ggplot(dr_df, aes(x = pseudotime, y = signature_pct_overlap)) + geom_point(aes(colour=annotation), size = 1) + theme_minimal() +
  geom_smooth(method='lm', formula= y~x, colour='black')
  #scale_color_manual(values = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
dev.off()


doubScores <- addDoubletScores(
  input = getArrowFiles(proj),
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#Doublet Enrichments - These represent the enrichment of simulated doublets nearby each single cell 
#compared to the expected if we assume a uniform distribution.

#s1 = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/QualityControl/MD-scATAC_71_Spleen/MD-scATAC_71_Spleen-Doublet-Summary.rds')
#h(doubScores[[1]]$doubletEnrich)
#table(proj$Sample)
doubScores = doubScores[match(names(proj@sampleMetadata), names(doubScores) )]

doublet_df = data.frame(DoubletEnrichment = unlist(lapply(doubScores, function(x) x$doubletEnrich)),
                        DoubletScore =  unlist(lapply(doubScores, function(x) x$doubletScore)),
                        cellNames = gsub('^.*\\.(.*)$', '\\1', names(unlist(lapply(doubScores, function(x) x$doubletEnrich)))))
#doublet_df = doublet_df[match(proj$cellNames, doublet_df$cellNames), ]                        
identical(rownames(proj), doublet_df$cellNames) #only identical when doubletEnrichment is returned for all samples
proj$DoubletEnrichment = doublet_df$DoubletEnrichment
proj$DoubletScore = doublet_df$DoubletScore

#when doublet estimation fails
doubScoresKeep = doubScores[!sapply(doubScores, function(x) all(x[[1]] == -1))]
doublet_df = data.frame(DoubletEnrichment = unlist(lapply(doubScoresKeep, function(x) x$doubletEnrich)),
                        DoubletScore =  unlist(lapply(doubScoresKeep, function(x) x$doubletScore)),
                        cellNames = gsub('^.*\\.(.*)$', '\\1', names(unlist(lapply(doubScoresKeep, function(x) x$doubletEnrich)))))
dr_df = as.data.frame(proj@cellColData)
dr_df = dr_df %>% dplyr::left_join(doublet_df, by = 'cellNames')
rownames(dr_df) = dr_df$cellNames
proj@cellColData = as(dr_df, 'DataFrame')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_summarize_dataframe(dr_df[, c('snn_harmony_res.1.7', 'DoubletEnrichment')], summarize_by_vec = dr_df$snn_harmony_res.1.7) %>%
  arrange(desc(DoubletEnrichment))
dr_df$DoubletEnrichment[is.na(dr_df$DoubletEnrichment)] = -1
jj_plot_features(reduction=dr_df, meta_features = 'DoubletEnrichment', cap_top = 'auto')

#add the single sample umap coordinates to the project
ss_umap_files = Sys.glob(file.path(pconfig$ARCHR_DIR, 'QualityControl/*/*Doublet*.rds'))
ss_umap_list = lapply(ss_umap_files, function(x){
  ss = read_rds(x)
  return(ss$originalDataUMAP[, c('X1', 'X2')])
})
ss_umap_df = do.call('rbind', ss_umap_list)
ss_umap_df = ss_umap_df[match(rownames(proj), rownames(ss_umap_df)), ]
stopifnot(identical(rownames(proj), rownames(ss_umap_df)))
proj$UMAP_1_ss = ss_umap_df$X1
proj$UMAP_2_ss = ss_umap_df$X2



# mouse normal cd8 annotation ---------------------------------------------

#mouse normal cd8
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal_cd8.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_mouse_normal_CD8'))


proj$annotation = from_to(proj$Clusters_0.5, old_new_map_vec = c(
  'C1' = 'CD8_Tnaive',
  'C2' = 'CD8_Tnaive', 
  'C3' = 'CD8_Tnaive',
  'C4' = 'CD8_Tnaive',
  'C5' = 'CD8_Tnaive',
  'C6' = 'CD8_Tcm', 
  'C7' = 'CD8_Tcm', 
  'C8' = 'CD8_Tcm', 
  'C9' = 'CD8_Tmem/MAIT', 
  'C10' = 'NK/MAIT', 
  'C11' = 'CD8_Trm', 
  'C12' = 'CD8_Trm', 
  'C13' = 'CD8_Trm_Tem', 
  'C14' = 'CD8_Trm_Tem'  
))
#Il7r (Cd127 is low in all 'effector/memory' clusters)
singler_pred = read_csv('/omics/groups/OE0436/internal/msimon/scATAC/mouse_normal_cd8_gene_activity_matrix_archr_singler_ref_monacoimmunedata_only_T_label.fine_predictions.csv')
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df$label = singler_pred$labels
jj_plot_features(reduction=dr_df, meta_features = 'label', label = T)
ref_df2 = jj_summarize_sparse_mat(gmat, proj$Clusters_0.5)

genes_plot = c('Cd3d','Cd3e','Cd8a','Cd8b1','Cd4','Ptprc','Sell','Ccr7','Lef1','Il7r','Itgae',
               'Cd69','Ikzf2','Klrg1','Tbx21','Cx3cr1','Cxcr3','Gzma','Prf1','Itga1','Id2','Id3',
               'Klra1','Il2rb','Cd27','Klrb1','Mki67')
res = archr_plot_markers(gmat, genes_plot, iW, reduction = dr_df)
jj_plot_dotplot(gmat, features_use = genes_plot, group_vec=proj$Clusters_0.5, scale_data = T)
jj_plot_heatmap(gmat, features_use = genes_plot, group_vec = proj$Clusters_0.5, scale_data = T)

ref_df2 = jj_summarize_sparse_mat(gmat, proj$annotation)
ref_all = cbind(ref_df, ref_df2)
rownames(ref_all) = as.character(rownames(ref_all))
cdf =  data.frame(label.fine=colnames(ref_all), label.main=colnames(ref_all))
rownames(cdf) = cdf$label.fine
se = SummarizedExperiment(assays = SimpleList(logcounts = ref_all),
                          colData = cdf)
#write_rds(se, '/omics/groups/OE0436/internal/msimon/scATAC/singler/singler_ref_mouse_normal_cd4_cd8_se.RDS')
