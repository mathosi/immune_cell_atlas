library(ArchR)
library(tidyverse)
library(jj)
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
dr_df = jj_get_reduction_coords(proj, 'UMAP')

### 1A Tissue umap
gg = jj_plot_features(dr_df, features='Tissue', pt_size = 0.5, 
                      custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                      return_gg_object = T) 

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tissue.pdf'),  width = 10, height = 8)
png(paste0(storeFigPath, 'mouse_imm_cell_atlas_tissue.png'), width = 7, height = 6, res = 400, units = 'in')
gg
dev.off()

### S1C Clusters umap
dr_df$Clusters_0.5 = mixsort_factor(dr_df$Clusters_0.5)
gg = jj_plot_features(dr_df, features='Clusters_0.5', pt_size = 0.5, 
                                            return_gg_object = T, label_type = 'geom_label', fill_colors = 'white') 

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_clusters.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()

### 1B Annotation umap
# gg = jj_plot_features(reduction=dr_df, meta_features='cluster_annotation_level3', pt.size = 0.5, 
#                       custom_colors = jj_get_colours(dr_df$cluster_annotation_level3, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
#                                             return_gg_object = T) 

gg = jj_plot_features(dr_df, features='cluster_annotation_level0', pt_size = 0.5,
                      custom_colors = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_annotation.pdf'),  width = 10, height = 8)
png(paste0(storeFigPath, 'mouse_imm_cell_atlas_annotation.png'), width = 7, height = 6, res = 400, units = 'in')
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()


# S1C nFrags umap
gg = jj_plot_features(dr_df, features='nFrags', pt_size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_nFrags.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()

# S1C sample umap 
gg = jj_plot_features(dr_df, features='Sample', pt_size = 0.5,
                      return_gg_object = T)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_sample.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()

# ### Not included: singler prediction
# dr_df$singler_prediction = from_to(vec= dr_df$singler_label, old_new_map_vec=c(
#   'B cells'= 'B cell',
#   'Basophils'= 'Basophil',
#   'Eosinophils'= 'Eosinophil',
#   'Macrophages'= 'Macrophage',
#   'Mast cells'= 'Mast cell',
#   'Monocytes'= 'Monocyte',
#   'Neutrophils'= 'Neutrophil',
#   'NK cells'= 'NK cell',
#   'T cells'= 'T cell',
#   'Tgd'= 'gd_T'
# ))
# gg = jj_plot_features(dr_df, features='singler_prediction', pt_size = 0.5, return_gg_object = T, 
#                       custom_colors = jj_get_colours(dr_df$singler_prediction, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) 
# 
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_singler.pdf'),  width = 10, height = 8)
# gg[[1]] + labs(colour='SingleR prediction') 
# dev.off()


# 2E Mouse atlas marker genes
genes_plot = c( 'Entpd1', 'Maf', 'Irf4', 'Pparg','Col15a1', 'Vps8')
gg = archr_plot_markers(proj, genes_plot)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_magic_marker_umap.pdf'),  width = 7.5, height = 6)
gg
dev.off()



# #get the treg subset and subcluster to refine the Treg annotation (just calculate once)
# proj_treg = proj[proj$cluster_annotation_level3 %in% c('Treg', 'tisTregST2', 'pTreg'), ]
# proj_treg = archr_dim_reduction(proj_treg)
# 
# proj_treg = loadArchRProject('ArchR_Project_Treg_subset')
# proj_treg <- addIterativeLSI(
#     ArchRProj = proj_treg,
#     useMatrix = "TileMatrix",
#     name = "IterativeLSI",
#     force = T,
#     iterations = 2,
#     clusterParams = list( #See Seurat::FindClusters
#       resolution = c(0.2),
#       sampleCells = 10000,
#       n.start = 10
#     ),
#     varFeatures = 25000,
#     dimsToUse = 1:20 #lowered from 30
#   )
# proj_treg <- addUMAP(
#     ArchRProj = proj_treg,
#     reducedDims = "IterativeLSI",
#     name = "UMAP",
#     force = T,
#     nNeighbors = 30,
#     minDist = 0.5,
#     metric = "cosine"
# )
# proj_treg = archr_clustering(proj_treg)
# signature_list = read_rds('//omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
# signature_list = signature_list[-c(4)]
# signature_list = lapply(signature_list, makeGRangesFromDataFrame)
# proj_treg = archr_add_peak_signatures(proj_treg, signature_list, signature_name = 'tisTreg_mouse1')
# proj_treg$treg_annotation = from_to(vec= proj_treg$Clusters_1.2, old_new_map_vec=c(
#   'C1'= 'tisTregST2_prog',
#   'C2'= 'Treg_naive',
#   'C3'= 'tisTregST2', #skin
#   'C4'= 'tisTregST2', #skin
#   'C5'= 'tisTregST2', #skin
#   'C6'= 'pTreg',
#   'C7'= 'tisTregST2', #VAT
#   'C8' = 'tisTregST2', #colon
#   'C9' = 'tisTregST2_prog'
# ))

# #saveArchRProject(proj_treg, outputDirectory = 'ArchR_Project_Treg_subset')
proj_treg = loadArchRProject('ArchR_Project_Treg_subset')
dr_df_treg = jj_get_reduction_coords(proj_treg, 'UMAP')
gmat = get_gene_mat(proj_treg)
##dr_df$tnkilc_annot = dr_df_tnkilc$cluster_annotation_level3[match(dr_df$cellNames, dr_df_tnkilc$cellNames)]

markers_use =  c('Foxp3','Cd4','Sell','Lef1','Klrg1','Nfil3','Ccr8','Il1rl1','Batf','Rorc','Ikzf2')
# jj_arrange_ggplots(jj_plot_features(reduction = dr_df_treg, meta_features = grep('sig', colnames(dr_df_treg), value = T),
#                                     return_gg_object = T, pt.size = 1), 9, 3)
# res = archr_plot_markers(proj_treg, markers_use)
# jj_arrange_ggplots(res, nplots = 9, cols = 3)

# #not included: mouse atlas treg subset umaps
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_treg_subset_annotation.pdf'),  width = 8, height = 6)
# jj_plot_features(dr_df_treg, features = 'Tissue', pt_size = 1, custom_colors = jj_get_colours(dr_df_treg$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
# jj_plot_features(reduction = dr_df_treg, meta_features = 'Clusters_1.2', pt.size = 1)
# #jj_plot_features(reduction = dr_df_treg, meta_features = 'cluster_annotation_level3', pt.size = 1)
# jj_plot_features(reduction = dr_df_treg, meta_features = 'treg_annotation', pt.size = 1, custom_colors = jj_get_colours(dr_df_treg$treg_annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'))
# jj_plot_heatmap(gmat, markers_use, proj_treg$Clusters_1.2)
# dev.off()

##scATAC tisTreg signature overlap (just calculate once)
# library(genomic_region_tools)
# source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')
# 
# mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
# mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
# mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
# emb_df = getReducedDims(proj_treg, reducedDims = 'IterativeLSI')
# dr_df_treg = jj_get_reduction_coords(proj_treg, 'UMAP')
# pmat = get_peak_mat(proj_treg)
# olap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df,
#                                 nFrags_vec = proj_treg$nFrags, signature_gr = mouse_cd4_tisTreg_gr,
#                                 verbose = T, k = 100, count_thres = 2e5)
# #write_rds(olap_df, paste0(storeFigPath, 'mouse_atlas_treg_subset_scATAC_tisTreg_sig_olap_df.RDS'))

# olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-29-te_analysis/mouse_atlas_treg_subset_scATAC_tisTreg_sig_olap_df.RDS')
# dr_df_treg$pct_overlap = olap_df$signature_pct_overlap
# #not included: tisTreg signature percentage overlap
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_treg_subset_tisTreg_signature_overlap_umap.pdf'),  width = 8, height = 6)
# jj_plot_features(dr_df_treg, features = 'pct_overlap', pt_size = 1, return_gg_object = T)[[1]] + 
#   labs(colour = '% overlap') 
# dev.off()
# 
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_treg_subset_tisTreg_signature_overlap_boxplot.pdf'),  width = 6, height = 5)
# jj_plot_numeric_by_group(dr_df_treg, 'pct_overlap', group_column = 'treg_annotation', 
#                          order = T, flip_coordinates = T, type = 'boxplot',
#                          custom_colors = jj_get_colours(dr_df$treg_annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
#   theme(legend.position = 'none') + labs(y='% overlap', x='Cell type')
# dev.off()

## add monocle trajectory

# #execute once
# gmat = get_gene_mat(proj_treg)
# library(monocle3)
# so = CreateSeuratObject(counts = gmat[, dr_df_treg$treg_annotation !='pTreg'], meta.data = dr_df_treg[dr_df_treg$treg_annotation !='pTreg', ])
# so = FindVariableFeatures(so) 
# so = ScaleData(so, features = VariableFeatures(so))
# gmat_use = GetAssayData(so, slot = 'scale.data')
# gmat_use <- gmat[rownames(gmat) %in% VariableFeatures(so), dr_df_treg$treg_annotation !='pTreg']
# pd_df <- droplevels(dr_df_treg[dr_df_treg$treg_annotation !='pTreg', ])
# fd_df <- data.frame(gene_short_name=rownames(gmat_use),
#                     row.names = rownames(gmat_use))
# rownames(gmat_use) <- fd_df$gene_short_name
# #MAKE SURE monocle3 is loaded and monocle is unloaded! Otherwise error later in plot_cells
# cds <- new_cell_data_set(gmat_use,
#                         cell_metadata = pd_df,
#                         gene_metadata = fd_df)
# # 
# # #performs log-normalization + PCA
# cds <- preprocess_cds(cds, num_dim = 50)
# # #umap with cosine metric
# cds <- reduce_dimension(cds)
# # 
# reducedDims(cds)$UMAP <- dr_df_treg[dr_df_treg$treg_annotation !='pTreg', 1:2]
# colnames(reducedDims(cds)$UMAP) <- c('V1', 'V2') #otherwise error in orderCells plotly picker
# cds <- cluster_cells(cds)
# # 
# # ## Step 5: Learn a graph
# cds <- learn_graph(cds, use_partition = T)
# cds <- order_cells(cds)#, root_cells = 'AAACGAAAGACTAATG-7') #root pickec at -5, -3
# cds@colData$sample_name = NULL
# #write_rds(cds, paste0(bigFilesDir, 'mouse_imm_cell_atlas_treg_subset_monocle3_cds.RDS'))
# #
# identical(colnames(cds), proj_treg$cellNames)
# proj_treg$pseudotime[match(cds@colData$cellNames, proj_treg$cellNames)] = pseudotime(cds, reduction_method = 'UMAP')
# dr_df_treg = jj_get_reduction_coords(proj_treg, 'UMAP')

# #not included: monocle trajectory of mouse atlas treg subset
# gg = jj_plot_features(reduction=dr_df_treg, meta_features = 'Tissue', pt.size = 1, 
#                       custom_colors = jj_get_colours(dr_df_treg$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'), return_gg_object = T)[[1]]
# cont_df = jj_get_contour_lines(dr_df_treg, 'treg_annotation', .95)
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_treg_subset_annotation_umap.pdf'), width =8, height=6)
# gg + geom_path(data = cont_df, aes(x=x, y=y, group=cont_group, linetype = treg_annotation), #, colour=snn_harmony_res.1), 
#                size=1) + labs(linetype = 'Cell type') + scale_linetype_manual(values = c(tisTregST2='solid', tisTregST2_prog='dotted', Treg_naive='dashed', pTreg='twodash'))
# dev.off()
# 
# gg = jj_plot_features(reduction=dr_df_treg, meta_features = 'pseudotime', pt.size = 1,colorScale = 'bry', return_gg_object = T)[[1]]
# gg2 =   monocle3::plot_cells(cds, color_cells_by = 'pseudotime', label_leaves = F, label_branch_points = F,trajectory_graph_color = 'black',  label_roots = F, cell_size = 0.5) + 
#   theme_minimal() + coord_fixed()
# gg$layers[[2]] = gg2$layers[[2]]
# 
# pdf(paste0(storeFigPath, 'mouse_atlas_treg_subset_trajectory_umap.pdf'), width =8, height=6)
# gg + labs(colour='Pseudotime')
# dev.off()

#saveArchRProject(proj_treg)


#### #full datset (add annotation to full dataset once)
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# #add in annotation from t nk ilc subset
# dr_df_tnkilc = jj_get_reduction_coords(loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected'))
# dr_df$cluster_annotation_fine = dr_df$cluster_annotation_level0
# dr_df$cluster_annotation_fine[match(dr_df_tnkilc$cellNames, dr_df$cellNames)] = dr_df_tnkilc$cluster_annotation_level3
# dr_df$cluster_annotation_fine[dr_df$cluster_annotation_fine == 'undef'] = 'undefined'
# dr_df$cluster_annotation_fine[match(dr_df_treg$cellNames, dr_df$cellNames)] = dr_df_treg$treg_annotation
# #write_csv(dr_df[, c('cellNames', 'cluster_annotation_fine')], '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
# cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
# proj$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
# proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'
# proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'CD4_other'] = 'Th17_Areg'
# saveArchRProject(proj)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(dr_df, features='cluster_annotation_fine', pt_size = 0.5, 
                      custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T, label_type =  'geom_label') 

# # #include labels later on
# gg = add_label(gg_obj = gg[[1]], text =  'tisTregST2', text_coord =  c(-13, -11), pointer_coord =  c(-8.5, -11))
# gg = add_label(gg, 'ILC2', c(-13, -4), c(-9, -4))
# gg = add_label(gg, 'ILC1', c(-7, -1), c(-7, -4))
# gg = add_label(gg, 'ILC3', c(-13, -4.5), c(-10, -4.5))
# gg = add_label(gg, 'ILC2', c(-13, -4.5), c(-10, -4.5))

# S2A mouse atlas with fine resolution in the t/nk/ilc subset
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_annotation_tnkilc_fine.pdf'),  width = 10, height = 8)
gg
dev.off()

# ### not included: cell types per tissue
# 
# dr_df$Tissue2 = factor(dr_df$Tissue, levels = c('Skin', 'VAT','Colon','Spleen'))
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_annotation_barplot.pdf'),  width = 6, height = 3.5)
# jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
#                              flip_coordinates = T, 
#                              custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) + 
#   labs(fill = 'Cell type', x='Tissue', y = 'Fraction')
# dev.off()
# 
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_annotation_barplot_absolute.pdf'),  width = 5, height = 6)
# jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
#                              flip_coordinates = F,absolute_numbers = T, add_text = T, text_size = 3, 
#                              custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
#   labs(fill = 'Cell type', x='Tissue', y = 'n cells')
# 
# dev.off()


# S2D cell type tissue distribution
dr_df$cluster_annotation_levelx = factor(dr_df$cluster_annotation_level0,
                                         levels = rev(gtools::mixedsort(unique(dr_df$cluster_annotation_level0))))
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tissue_by_cell_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'Tissue2', 'cluster_annotation_levelx', flip_coordinates = T,
                             custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  labs(x = 'Cell type', y = 'Fraction', fill = 'Tissue')
dev.off()

### marker dotplot
proj_use = proj[!proj$cluster_annotation_level0 %in% c('undefined'), ]
gmat = get_gene_mat(proj_use)
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')

# features_plot = list(
#   Tcells= c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b1'), #Cd247, CD3z=CD247
#   Tregs = c('Foxp3','Il2ra','Ikzf2',   'Kit','Klrg1','Thy1','ID2'),
#   T_naive = c('Sell', 'Ccr7', 'Lef1','Il7r'),
#   TisTregs=c('Il1rl1','Klrg1','Areg','Tox', 'Nfil3', 'Ccr8'),
#   Bcells=c('Cd19','Ms4a1', 'Cd27','Jchain'),
#   Macrophages=c('Itgam','Mrc1','Spi1','Csf1r','S100a9'),
#   DCs=c('Itgax','Siglech','Bst2','Xcr1','Clec9a','Sirpa','Irf8'),
#   NK=c('Klrb1c','Ncr1', 'Gzma', 'Gzmb','Prf1'),
#   ILCs=c('Il7r','Il2ra'),
#   Mastcells_Basophil=c('Itgam','Hnmt','Areg')
# ) 

features_plot = list(
  Tcells= c('Cd3e','Cd3d'),#,'Cd4','Cd8a','Cd8b1'), #CD3z=CD247
  Bcells=c('Cd19','Ms4a1'),
  Macrophages=c('Mrc1','Csf1r'),
  Monocyte = c('Csf1r','Cx3cr1'),
  DCs=c('Siglech','Clec9a'),
  NK=c('Gzma','Prf1', "Klrb1"),
  ILCs=c('Il7r','Il2ra','Kit'),
  plasma = c('Sdc1','Bst2'),
  neutrophil = c('Ly6g', 'Itgam','S100a9'),
  Mastcells_Eo_Basophil=c('Itgam',"Syne1","Il4"),
  t_nk_ilc = c('Klrg1','Thy1','Id2')
) 
#gene_avail('Id2', rownames(gmat))
# gg = jj_plot_heatmap(obj = gmat, group_vec = proj_use$cluster_annotation_level0, 
#                      features_use = unique(unname(unlist(features_plot))), scale_data = T, 
#                      row_annot = rowAnnotation(annot = sort(unique(proj_use$cluster_annotation_level0)), col=list(annot = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) ))
# gg
# 
# gg = jj_plot_dotplot(obj = gmat, group_vec = proj_use$cluster_annotation_level0, scale_data = T,
#                       features_use = unique(unname(unlist(features_plot))),
#                        )
# gg

library(Seurat)
library(cowplot)
# 1C mouse atlas marker gene dotplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_marker_heatmap.pdf'),  width = 10, height = 5)
seurat_dotplot(gmat=gmat, metadf = dr_df,
               features = unique(unname(unlist(features_plot))),
               group_column = 'cluster_annotation_level0')
dev.off()


### bZIP TF gene activity
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
proj_use = proj[dr_df$cluster_annotation_fine != 'undefined', ]
dr_df = dr_df[dr_df$cluster_annotation_fine != 'undefined', ]
gmat = get_gene_mat(proj_use)

# ##get bZIP genes from homer library
# motifPositions <- getPositions(proj, name = 'Motif')
# bzip_factors = grep('bzip', names(motifPositions), value = T, ignore.case = T)
# bzip_factors = gsub('(.*)\\.bZIP.*','\\1', bzip_factors)
genes_plot = c('Atf1','Atf2','Atf3','Atf4','Atf7', 'Bach1', 'Bach2', 'Batf', 'Batf2','Batf3','Cebpa', 'Cebpb',
               'Cebpd','Cebpg', 'Creb1','Crebzf', 'Fosl1','Fosl2', "Jun","Junb", 'Jund','Mafa', 'Maff', 'Mafk', 'Nfe2l2', 'Nfe2',
               "Nfatc4", "Nfatc1", "Nfatc2", "Nfatc3", "Ddit3", "Irf4", "Irf8")

gmat_use = gmat[rownames(gmat) %in% genes_plot, ]
gmat_use = gmat_use[rowVars(gmat_use) > 0, ]

proj_use$cluster_annotation_fine[proj_use$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'
dr_df$cluster_annotation_fine[dr_df$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'

# S2I mouse atlas bZIP gene activity dotplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_bZIP_genes_dotplot.pdf'),  width = 10, height = 7)
seurat_dotplot(gmat=gmat_use, metadf = dr_df,
               features = rownames(gmat_use),
               group_column = 'cluster_annotation_fine')
# jj_plot_heatmap(obj = gmat_use, group_vec = dr_df$cluster_annotation_fine, 
#                 features_use = rownames(gmat_use), scale_data = T)
dev.off()

# immune cell atlas qc plots ----------------------------------------------

# distribution of fragments -----------------------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
dr_df = as.data.frame(proj@cellColData)
### fragment size distribution
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% Fragments', color='Sample') + 
  scale_color_manual(values = jj_get_jj_colours(fsize_df$group))

# S1B mouse atlas fragment size distribution
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_fragment_size_distribution.pdf'), width = 8, height = 3)
p1
dev.off()

# S1A mouse atlas TSS-nFrags scatterplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_nFrags_TSS_scatterplot.pdf'), width = 10, height = 8)
  nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6)
dev.off()

# #not included: doublet enrichment umap
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# gg = jj_plot_features(reduction=dr_df, meta_features='DoubletEnrichment', pt.size = 0.5, 
#                       return_gg_object = T, cap_top = 'q95')
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_doublet_enrichment.pdf'),  width = 10, height = 8)
# gg
# dev.off()


# correlation heatmap -----------------------------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
pmat = get_peak_mat(proj)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#sum_mat = jj_summarize_sparse_mat(pmat, dr_df$tissue_celltype)
dr_df$tissue_celltype = paste(dr_df$Tissue, dr_df$cluster_annotation_level0, sep ='_')
dr_df$tissue_celltype = replace_if(dr_df$tissue_celltype, count_below = 100, 'undefined')
pmat_use = pmat[, !grepl('undefined', dr_df$tissue_celltype)]
dr_df = dr_df[!grepl('undefined', dr_df$tissue_celltype), ]
htmat = jj_cluster_correlation_heatmap(pmat_use, group_a = dr_df$tissue_celltype, plot_complex_heatmap = F)
col_file = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'
column_ha = HeatmapAnnotation(Tissue = sapply(strsplit(colnames(htmat), '_'), '[[', 1),
                              Celltype = sapply(strsplit(colnames(htmat), '_'), '[[', 2),
                              col = list(Tissue = jj_get_colours( sapply(strsplit(colnames(htmat), '_'), '[[', 1), '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                                         Celltype = jj_get_colours(sapply(strsplit(colnames(htmat), '_'), '[[', 2), col_file, comment_char = '$')))

row_ha = rowAnnotation(Tissue = sapply(strsplit(colnames(htmat), '_'), '[[', 1),
                       Celltype = sapply(strsplit(colnames(htmat), '_'), '[[', 2),
                       col = list(Tissue = jj_get_colours( sapply(strsplit(colnames(htmat), '_'), '[[', 1), '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                                  Celltype = jj_get_colours(sapply(strsplit(colnames(htmat), '_'), '[[', 2), col_file, comment_char = '$')),
                       show_legend = FALSE)

#orderuse = hclust(dist(htmat, method = 'euclidean'), method = 'complete')$order
# cnames = colnames(htmat)[order_use]
#better version than triangular heatmap
ht = Heatmap(htmat, rect_gp = gpar(type = "none"), column_dend_side = 'top', show_row_dend = F, row_dend_side = 'right',show_column_names = F,
             top_annotation = column_ha, right_annotation = row_ha, name = 'Cor',
             cell_fun = function(j, i, x, y, w, h, fill) {
               #if(i <= j){
               if(as.numeric(x) + 1e-6 >= 1 - as.numeric(y)) {
                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
               }
               if(i == j){ #i = row index, j = column index, x = column position in units, y = row position in units, h = rect size (subtract the value to move the lavel to bottom in y direction)
                 grid.text(label = colnames(htmat)[i], x = x, y = y-0.7*h, rot = 90, just = 'right')
               }
               # }else if(i == j+1){
               #   grid.rect(x, y, w, h, gp = gpar(fill = 'green', lwd = 0))
               #   grid.text(sprintf("%s", colnames(htmat)[j]), x, y, rot=90, just = 'right', gp = gpar(fontsize = 6))
               # }
             })

# S1J Mouse atlas tissue-cell type correlation heatmap
pdf(paste0(storeFigPath, 'celltype_tissue_correlation_heatmap.pdf'), width = 12, height = 10)
draw(ht, padding =  unit(c(45, 2, 2, 2), 'mm'))
dev.off()


# # markers per cell type -------------------------------------------------
# 
# markerGenes <- getMarkerFeatures(
#   ArchRProj = proj,
#   groupBy = "cluster_annotation_level3",
#   useMatrix = "GeneScoreMatrix",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# ### marker gene activity heatmap with the possibility to highlight specifig genes
# markerList <- getMarkers(markerGenes)#, cutOff = "FDR <= 0.001 & Log2FC >=1")
# marker_df = archr_get_markers_as_df(markerGenes, proj, cutOff = "FDR <= 0.001 & Log2FC >=1", annotate_closest_gene = F)
# table(marker_df$comparison)
# heatmapGenes <- plotMarkerHeatmap(
#   seMarker = markerGenes, 
#   cutOff = "FDR <= 0.001 & Log2FC >= 1", returnMatrix =T,
#   transpose = T, labelMarkers = unique(unname(unlist(features_plot))), #c("Cd19", "Cd79a" )
# ) #problem: need to pick the exact label .1, .2 etc for the peak which is the relevant one
# library(ComplexHeatmap)
# draw(heatmapGenes, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# 
# highlight_cols = c('Cd79a', 'Cd19', 'S100a9', 'Csf1r','Cd8a', 'Foxp3')
# highlight_cols = unique(unname(unlist(features_plot)))
# highlight_pos = match(highlight_cols, colnames(heatmapGenes))
# #colnames(pmat_use)[highlight_pos] = highlight_cols
# colAnno = columnAnnotation(mark = anno_mark(at = highlight_pos, labels = highlight_cols, 
#                                             labels_gp = gpar(fontsize = 8), padding = 0))
# 
# ht = Heatmap(heatmapGenes, show_column_names = F, cluster_columns = F, name = 'scaled\ngene activity',
#              cluster_rows = F, use_raster = T, show_column_dend = F, top_annotation = colAnno)
# draw(ht)


# number of DA peaks in the major subsets ---------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
# markers_se <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = 'PeakMatrix',  
#   groupBy = "cluster_annotation_level0",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )

#write_rds(markers_se, paste0(bigFilesDir, 'mouse_atlas_cluster_annotation_level0_all_markers.RDS'))
markers_se = read_rds(paste0(bigFilesDir, 'mouse_atlas_cluster_annotation_level0_all_markers.RDS'))
markers_use_se = markers_se[, !colnames(markers_se) %in% 'undefined']
#transpose = T results in wrong clustering of cell types, nLabel = 0 is not possible...
# 1D: mouse atlas marker peaks heatmap
pdf(paste0(storeFigPath, 'mouse_marker_peak_heatmap.pdf'), width =6, height=6)
plotMarkerHeatmap(markers_use_se, 
                  cutOff = "FDR <= 0.01 & Log2FC >=0.5",  
                  transpose = F, clusterCols = T, binaryClusterRows = T, #clutering not working with transpose = T
                  nLabel = 1) #83575
dev.off()


marker_df = archr_get_markers_as_df(markers_use_se, proj, cutOff = "FDR <= 0.01 & Log2FC >=0.5")
#write_csv(marker_df, paste0(storeFigPath, 'table_1d_mouse_immune_cell_atlas_marker_peaks.csv'))
marker_list = split(marker_df, marker_df$comparison)
marker_list = lapply(marker_list, '[[', 'feature')
#S1E Mouse atlas marker peaks upset plot
pdf(paste0(storeFigPath, 'mouse_marker_peak_upset_plot.pdf'), width =8, height=6)
  jj_plot_upsetr(marker_list)
dev.off()

pgr = getPeakSet(proj)
peak_df = as.data.frame(unname(pgr))
peak_df$comparison = sprintf('Union peaks')

# jj_plot_categorical_by_group(peak_df,
#                              feature_column =  'peakType',
#                              group_column = 'comparison', 
#                              absolute_numbers = F,
#                              flip_coordinates = T, add_text = T, text_size = 3) + 
#   labs(y='n peaks', x = '')

peak_df$feature = with(peak_df, paste(seqnames, start, end, sep = '-'))
marker_df = marker_df %>% dplyr::left_join(peak_df[, c('feature','nearestGene', 'distToGeneStart', 'nearestTSS', 'distToTSS', 'peakType')], by = 'feature')
subset_df = marker_df %>% dplyr::select(comparison, peakType)
subset_df = rbind(subset_df, peak_df[, c('comparison', 'peakType')])
subset_df = subset_df %>% dplyr::group_by(comparison) %>% add_count() %>% ungroup() %>% 
  dplyr::mutate(comparison = sprintf('%s (%i)', comparison, n))

nmarkers = as.data.frame(table(subset_df$comparison)) %>% arrange(Freq) %>% pull(Var1) %>% as.character
subset_df$comparison = factor(subset_df$comparison, levels = nmarkers)
# S1F mouse atlas marker peaks genomic region annotation barplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_peak_type_barplot.pdf'),  width = 6, height = 3.5)
  jj_plot_categorical_by_group(subset_df, 'peakType', 'comparison', flip_coordinates = T, absolute_numbers = F, add_text = F, text_size = 3) + 
  labs(y='n peaks', x = 'Cell type')
dev.off()

#enriched TF motifs in the marker peaks
motifsUp <- peakAnnoEnrichment(
  seMarker = markers_use_se,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(motifsUp, n = 5, transpose = T)
# S1G Mouse atlas enriched TF heatmap
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_cell_type_markers_enr_motifs_heatmap.pdf'),  width = 7, height =5)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


# homer transcription factor z-scores -------------------------------------

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
binding_domains = binding_domains[!binding_domains == '?']
binding_domains[binding_domains=='T-box'] = 'T.box'
plotVarDev$data$bindingDomain = 'Other'
for(i in binding_domains){
  plotVarDev$data$bindingDomain[grepl(i, plotVarDev$data$name)] = i  
}
table(plotVarDev$data$bindingDomain)
plotVarDev$data$name = gsub('(.*)_[0-9]+','\\1', plotVarDev$data$name)
#cat(paste(names(jj_get_jj_colours(plotVarDev$data$bindingDomain)), jj_get_jj_colours(plotVarDev$data$bindingDomain), sep = ','),sep='\n')

# # not included: mouse atlas TF deviation rank plot
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tf_deviation_rankplot.pdf'),  width = 8, height =5)
# ggplot() + geom_point(data = plotVarDev$data, aes(x=rank, y=combinedVars, colour=bindingDomain), size = 1.5) + 
#   geom_text_repel(data = plotVarDev$data[1:20, ], aes(x=rank, y=combinedVars, label=name)) + 
#   theme_minimal() + labs(x = 'Rank') + scale_colour_manual(values = jj_get_colours(plotVarDev$data$bindingDomain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
#   labs(y = 'Combined variation', colour = 'Binding domain')
# dev.off()

proj$peripheral_tissue = ifelse(proj$Tissue == 'Spleen', 'spleen',  'peripheral')
proj$cluster_annotation_level00 = from_to(proj$cluster_annotation_level0, old_new_map_vec = c(Macrophage = "Macrophage/Monocyte", Monocyte = "Macrophage/Monocyte"))
proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level00, sep='__')

#sig_se = get_archr_mat(proj, 'MotifMatrix')
sig_se = getMatrixFromProject(proj, 'MotifMatrix')
z_score_mat = t(assays(sig_se)[['z']])
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))

# chromvar_summarized = t(jj_summarize_sparse_mat(Matrix::t(z_score_mat), proj$peripheral_cell_type))
# Heatmap(chromvar_summarized, split = gsub('(.*)__.*', '\\1', rownames(chromvar_summarized)))
# 
# own_se = SummarizedExperiment(assays = assays(markers_use_se), 
#                               metadata = list(Params = list(useMatrix ='PeakMatrix')), #metadata(markers_use_se),
#                               rowData = rowData(markers_use_se))#list(Params = list(useMatrix ='PeakMatrix')))

# rowData
# DataFrame with 254545 rows and 4 columns
# seqnames     idx     start       end
# <Rle> <array>   <array>   <array>
#   1          chr1       1   4470163   4470663

# motifsUp <- peakAnnoEnrichment(
#   seMarker = own_se,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.58"
# )

comparisons_keep = dr_df %>% dplyr::group_by(cluster_annotation_level00) %>% dplyr::count(peripheral_tissue) %>% 
  tidyr::spread(peripheral_tissue, n) %>% 
  dplyr::filter(peripheral > 200 & spleen > 200 & !cluster_annotation_level00 %in%  c('undefined')) %>% 
  dplyr::pull(cluster_annotation_level00)
# "B cell",  "DC", "Macrophage/Monocyte", "NK cell", "T cell" 

# markers_list = list()
# for(i in unique(proj$cluster_annotation_level00)){
#   message(i)
#   if(!i %in%  comparisons_keep) next
#   markers_list[[i]] <- getMarkerFeatures(
#     ArchRProj = proj[proj$cluster_annotation_level00 == i, ], 
#     useMatrix = use_matrix, 
#     useGroups = 'peripheral', 
#     groupBy = "peripheral_tissue",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon"
#   )
# }
# head(assays(markers_list$`T cell`)[[1]])
# 
# #combine the markers from the individual comparisons in one assay, which is put in a summarized experiment
# new_ass_list = list()
# for(i in seq_along(assays(markers_list[[1]]))){
#   new_name = names(assays(markers_list[[1]]))[i]
#   ass_list = list()
#   for(j in seq_along(markers_list)){
#     ass_list[[j]] = assays(markers_list[[j]])[[i]]
#   }
#   ass_df = Reduce(cbind, ass_list)
#   colnames(ass_df) = names(markers_list)
#   new_ass_list[[new_name]] = ass_df[, colnames(ass_df) %in% comparisons_keep]
# }
# own_se = SummarizedExperiment(assays = new_ass_list, 
#                               metadata = list(Params = list(useMatrix ='PeakMatrix')), #metadata(markers_use_se),
#                               rowData = rowData(markers_list[[1]]))
# 
# 
# marker_df = archr_get_markers_as_df(own_se, proj,  cutOff = "FDR <= 0.01 & Log2FC >= 0.58")
# table(marker_df$comparison)
# # comparisons_keep = comparisons_keep[comparisons_keep %in% names(which(table(marker_df$comparison) > 100))]
# 
# #write_rds(own_se, paste0(bigFilesDir, 'mouse_atlas_peripheral_vs_spleen_marker_peaks_se.RDS'))
# plotMarkerHeatmap(own_se, cutOff = "FDR <= 0.01 & Log2FC >=0.58", transpose = T, nLabel = 1) #83575
# 
# motifsUp <- peakAnnoEnrichment(
#   seMarker = own_se,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.58"
# )
# heatmapEM <- plotEnrichHeatmap(motifsUp, n = 100, transpose = T)
# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# motif_list = motifs_select = list()
# for(i in comparisons_keep){
#   message(i)
#   if(!i %in%  comparisons_keep) next
#   dr_df_use = dr_df[dr_df$cluster_annotation_level00 == i, ]
#   z_score_mat_use = t(z_score_mat[dr_df$cluster_annotation_level00 == i, ])
#   z_score_mean_mat = as.data.frame(jj_summarize_sparse_mat(z_score_mat_use, dr_df_use$peripheral_tissue))
#   z_score_mean_mat$diff = NA
#   for(j in 1:nrow(z_score_mean_mat)){
#     z_score_mean_mat[j, 'diff'] = diff(c(z_score_mean_mat[j, 'spleen'], z_score_mean_mat[j, 'peripheral']))
#   }
#   z_score_mean_mat = z_score_mean_mat#[order(z_score_mean_mat$diff, decreasing = T), ]
#   motif_list[[i]] <- z_score_mean_mat#[z_score_mean_mat$diff > 1, ]
#   motifs_select[[i]] = rownames(z_score_mean_mat)[z_score_mean_mat$diff > 1]
# }
# jj_plot_upsetr(motifs_select)
# jj_make_overlap_df(motifs_select)
# length(unique(unlist(motifs_select)))
# motifs_plot = unique(unlist(motifs_select))

da_list = list()
for(i in comparisons_keep){
  da_list[[i]] = getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = 'MotifMatrix',
    groupBy = 'peripheral_cell_type',
    useGroups = paste0('peripheral__', i),
    bgdGroups = paste0('spleen__', i),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
}
# #res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0.05 & FDR <= 0.01'))
# res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0 & FDR <= 0.01'))
# 
# da_tf_list = lapply(res_se, function(x) x[[1]][, 'name'])
# jj_plot_upsetr(da_tf_list)
# jj_make_overlap_df(da_tf_list)
# motifs_plot = unique(unlist(motifs_select))
# length(motifs_plot)
# 
# motif_df = jj_initialize_df(ncol = length(motif_list), nrow = nrow(motif_list[[1]]), col.names = names(motif_list), row.names= rownames(motif_list[[1]]))
# for(i in seq_along(motif_list)){
#   motif_df[, i] = motif_list[[i]]$diff
# }
# col_fun = colorRamp2(c(-5, 0, 5), c("darkblue", "white", "red"))
# h1 <- Heatmap(t(motif_df[rownames(motif_df) %in% motifs_plot, ]), col = col_fun,
#               name = 'z-score difference\nperipheral - spleen',
#               column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 8))
# h1


#####
res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0.05 & FDR <= 0.01'))

da_tf_list = lapply(res_se, function(x) x[[1]][, 'name'])
jj_plot_upsetr(da_tf_list)
motif_olap_df =jj_make_overlap_df(da_tf_list)
motifs_plot = unique(unlist(da_tf_list))
motifs_plot = na.omit(unique(c(unlist(motif_olap_df$overlaps[,!(colnames(motif_olap_df$overlaps) %in% c('a','b','c','d','e'))]))))
length(motifs_plot)

res_se2 = lapply(da_list, function(x) getMarkers(x, cutOff = 'FDR <= 1'))
diff_list = lapply(res_se2, function(x) as.data.frame(x[[1]][, c('name','MeanDiff')]))
res_df = reduce(diff_list, full_join, by = "name") %>% column_to_rownames('name')
colnames(res_df) = names(diff_list)

plot_mat = t(res_df[rownames(res_df) %in% motifs_plot, ])
colnames(plot_mat) = gsub('(.*)_[0-9]+','\\1',colnames(plot_mat))

heatmap_binding_domain = rep('Other', ncol(plot_mat))
for(i in binding_domains){
  heatmap_binding_domain[grepl(i, colnames(plot_mat))] = i
}
ha = columnAnnotation('Binding Domain'= heatmap_binding_domain, col = list('Binding Domain' = jj_get_colours(heatmap_binding_domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))


library(circlize)
Heatmap(plot_mat) #define range based on automatic colorscale
col_fun = colorRamp2(c(-0.1, 0, 0.1, 0.5), c("darkblue", "white", "red", "darkred"))
h1 <- Heatmap(plot_mat, col = col_fun,
              name = 'MeanDiff\nperipheral - spleen',
              column_names_gp = gpar(fontsize = 8), 
              row_names_gp = gpar(fontsize = 8),
              top_annotation = ha)
# 1E: mouse atlas peripheral vs spleen TF deviation heatmap
pdf(paste0(storeFigPath, 'mouse_peripheral_tf_heatmap.pdf'), width =10, height=3)
png(paste0(storeFigPath, 'mouse_peripheral_tf_heatmap.png'), width = 8, height = 3, res = 400, units = 'in')
h1
dev.off()


# da_list = list()
# for(i in comparisons_keep){
#   da_list[[i]] = find_markers_single_cell(
#     feat_mat = t(z_score_mat),
#     meta_df = dr_df,
#     ident = 'peripheral_cell_type',
#     ident.1 = paste0('peripheral__', i),
#     ident.2 = paste0('spleen__', i),
#     only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
# }
# da_tf_list = lapply(da_list, function(x) rownames(x)[x$FDR < 0.01])


# # using archrs hypergeometric test ------------------------------------------------
# 
# dr_df = as.data.frame(cbind(dr_df, z_score_mat))
# jj_plot_features(reduction = dr_df, meta_features = 'RORgt.NR_246', cap_top = 'q99', cap_bottom = 'q01', 
#                  facet_by = 'peripheral_tissue', background_cells = T)
# jj_plot_numeric_by_group(dr_df[dr_df$cluster_annotation_level00 %in% comparisons_keep, ], 'Tcf3.HMG_288', 'peripheral_cell_type', order = T, flip_coordinates = T)
# 
# chromvar_summarized = t(jj_summarize_sparse_mat(Matrix::t(z_score_mat), proj$peripheral_cell_type))
# library(ComplexHeatmap)
# h1 <- Heatmap(chromvar_summarized,
#               name = 'chromvar deviation', 
#               column_names_gp = gpar(fontsize = 3))#,
# #top_annotation = HeatmapAnnotation(cluster = anno_text(x = colnames(chromvar_summarized))))
# pdf(paste0(storeFigPath, dataset_use, 'chromvar_heatmap.pdf'), width = 8, height=25)
# h1
# #+ rowAnnotation(link = anno_mark(at = 1:nrow(chromvar_mat),
# #labels = rownames(chromvar_summarized), 
# #labels_gp = gpar(fontsize = 10))) 
# dev.off()
# 
# #proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
# motifsUp <- peakAnnoEnrichment(
#   seMarker = markers_use_se,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.58"
# )
# heatmapEM <- plotEnrichHeatmap(motifsUp, n = 6, transpose = T)
# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")


# motif comparison --------------------------------------------------------
#https://bioconductor.org/packages/release/bioc/html/universalmotif.html

library(universalmotif)
#library(MotifDb)
library(chromVARmotifs)

#same data as used by 'addMotifAnnotations' from ArchR
data("homer_pwms")
motifs <- homer_pwms
obj <- ArchR:::.summarizeChromVARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary
motifs[['AP.1.bZIP_1']]
cor_mat = compare_motifs(motifs[names(motifs) %in% motifs_plot], method = 'PCC')
#Heatmap(cor_mat, name = 'Pcor',column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6))
clust_res = hclust(dist(cor_mat))
clust_res$order
# S1I TF motif sequence comparison heatmap
pdf(paste0(storeFigPath, 'mouse_motif_sequence_triangular_heatmap.pdf'), width =8, height=7)
triangular_heatmap(cor_mat[clust_res$order, clust_res$order], text_format = NULL,
                   col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), name = 'Pcor', label_size = 10)
dev.off()
# library(cowplot)
# library(MotifDb)
# motifs2 <- convert_motifs(MotifDb[1:10])
# motifs2 <- convert_motifs(motifs[grep('bZIP', names(motifs), value= T)])[1:15]
# ## Get the tree: make sure it's a horizontal type layout
# tree <- motif_tree(motifs2, layout = "rectangular", linecol = "none")
# ## Now, make sure we order our list of motifs to match the order of tips:
# mot.names <- sapply(motifs2, function(x) x["name"])
# names(motifs2) <- mot.names
# new.order <- tree$data$label[tree$data$isTip]
# new.order <- rev(new.order[order(tree$data$y[tree$data$isTip])])
# motifs2 <- motifs2[new.order]
# ## Plot the two together (finessing of margins and positions may be required):
# plot_grid(nrow = 1, rel_widths = c(1, -0.15, 1),
#           tree + xlab(""), NULL,
#           view_motifs(motifs2, names.pos = "right") +
#             ylab(element_blank()) +
#             theme(
#               axis.line.y = element_blank(),
#               axis.ticks.y = element_blank(),
#               axis.text.y = element_blank(),
#               axis.text = element_text(colour = "white")
#             )
# )

# TF footprints -----------------------------------------------------------

library(BSgenome.Mmusculus.UCSC.mm10)
motifPositions <- getPositions(proj, name = 'Motif')
grep('NFkB.p65', names(motifPositions), value = T)
motifs <- c("Fra1.bZIP_90", "Atf3.bZIP_12", "BATF.bZIP_20", 'NFkB.p65.RHD_208')[3] #only show BATF for now
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
proj_use = proj[proj$peripheral_cell_type %in% c('spleen__T cell', 'peripheral__T cell'), ] #, 'spleen__Monocyte', 'peripheral__Monocyte'), ]
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)
#plot(gg$BATF.bZIP_20)

proj_use = proj[proj$peripheral_cell_type %in% c('spleen__Macrophage/Monocyte', 'peripheral__Macrophage/Monocyte'), ] #, 'spleen__Monocyte', 'peripheral__Monocyte'), ]
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg2 = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)

plist = list(`T cell BATF.bZIP` = gg$BATF.bZIP_20, `Macrophage/Monocyte BATF.bZIP` = gg2$BATF.bZIP_20)
# 1F mouse atlas BATF TF footprints
pdf(paste0(storeFigPath, 'tf_footprints.pdf'), width =6, height=5)
do.call(cowplot::plot_grid, c(list(ncol = 2, labels = c('T cell', 'Macrophage/Monocyte'), label_size=10),plist))
dev.off()


# scATAC tisTreg signature percentage overlap -----------------------------

#devtools::install_github('mathosi/genomic_region_tools')
library(genomic_region_tools)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)

config_list = get_config()
dataset_use = 'mouse_normal'
pconfig = yaml::read_yaml(config_list[[dataset_use]])
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)

# emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
# pmat = get_peak_mat(proj)
# olap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df,
#                                  nFrags_vec = proj$nFrags, signature_gr = mouse_cd4_tisTreg_gr,
#                                  verbose = T, k = 100, count_thres = 2e5)
#write_rds(olap_df, paste0(storeFigPath, 'mouse_atlas_scATAC_tisTreg_sig_olap_df.RDS'))
#olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-11-16-percentage_signature/mouse_atlas_scATAC_tisTreg_sig_olap_df.RDS')
olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-homer2/mouse_atlas_scATAC_tisTreg_sig_olap_df.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$pct_overlap = olap_df$signature_pct_overlap

# 2A: Mouse atlas tisTreg signature percentage overlap
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_overlap_umap.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = 'pct_overlap', return_gg_object = T)[[1]] + labs(colour = '% overlap') 
dev.off()
# celltypes_keep = dr_df$cluster_annotation_fine %>% table %>% .[. > 200] %>% names
# dr_df = dr_df[dr_df$cluster_annotation_fine %in% celltypes_keep, ]

# 2B: Mouse atlas tisTreg signature percentage overlap boxplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_overlap_boxplot.pdf'),  width = 6, height = 5)
jj_plot_numeric_by_group(dr_df[!dr_df$cluster_annotation_fine == 'undefined', ], 'pct_overlap', group_column = 'cluster_annotation_fine', 
                         order = T, flip_coordinates = T, type = 'boxplot',
                         custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  theme(legend.position = 'none') + labs(y='% overlap', x='Cell type')
dev.off()

# scATAC tisTreg signature subsets ----------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
#contains the scATAC tisTreg regions as peakset
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets', pconfig$DONOR)) 
# proj = addPeakSet(proj, peakSet = mouse_cd4_tisTreg_gr, force=T)
# proj <- addPeakMatrix(proj,binarize = T)
# #saveArchRProject(proj)
cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/treg_fine_annotation.csv')
cluster_annotation_fine_df$cluster_annotation_fine[cluster_annotation_fine_df$cluster_annotation_fine == 'Th_Il21'] = 'Tfh_like'
stopifnot(identical(proj$cellNames, cluster_annotation_fine_df$cellNames))
proj$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine

pmat = get_peak_mat(proj)
summarize_vec = proj$cluster_annotation_fine
mean_mat = jj_summarize_sparse_mat(pmat, summarize_vec)
mean_mat = mean_mat[, !colnames(mean_mat) %in% c('undefined')]
scaled_mat = t(scale(t(mean_mat)))
htmat = t(scaled_mat)
library(ComplexHeatmap)

ht = Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, gap = unit(3,'mm'),
             show_column_names = F, column_split = 5, name = 'scaled\naccessibility')
ht = draw(ht)

# 2C Mouse atlas tisTreg signature subset heatmap
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_subset_heatmap.pdf'), width =8, height=5)
ht
dev.off()


# htmat = htmat[, match(mouse_cd4_tisTreg_df$feature, colnames(htmat)) ]
# markers = c('Nfil3','Entpd1','Havcr2','Irf4','Maf','Ctla4','Ikzf2','Batf','Irf8',
#             'Areg','Rorc','Il1rl1','Gata3','Klrg1')
# 
# stopifnot(identical(mouse_cd4_tisTreg_df$feature, colnames(htmat)))
# ids_highlight = which(mouse_cd4_tisTreg_df$symbol %in% markers & mouse_cd4_tisTreg_df$distance < 2000)
# mouse_cd4_tisTreg_df$symbol_use = mouse_cd4_tisTreg_df$symbol
# mouse_cd4_tisTreg_df$symbol_use[!mouse_cd4_tisTreg_df$symbol_use %in% markers] = ' '
# ha = columnAnnotation(`Closest gene` = anno_mark(at = ids_highlight, labels = mouse_cd4_tisTreg_df$symbol[ids_highlight]))
# #ha = columnAnnotation(foo1 = mouse_cd4_tisTreg_df$symbol_use, col = list(foo1 = jj_get_jj_colours( mouse_cd4_tisTreg_df$symbol_use)))
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_subset_heatmap_annotated.pdf'), width = 12, height = 8)
# Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, gap = unit(3,'mm'),
#         show_column_names = F, column_split = 5, name = 'scaled\naccessibility', top_annotation = ha)
# dev.off()

corder = column_order(ht)
subset_list = list(
  c1_myeloid =  colnames(htmat)[corder[[1]]],
  c2_cd8_nk =  colnames(htmat)[corder[[2]]],
  c3_treg =  colnames(htmat)[corder[[3]]],
  c4_plasma =  colnames(htmat)[corder[[4]]],
  c5_ilc_th17 = colnames(htmat)[corder[[5]]]
)
#jj_plot_upsetr(subset_list)
# Heatmap(htmat[, colnames(htmat) %in% subset_list$c2_cd8_nk], show_row_names=T, show_column_dend = F, show_row_dend = F, 
#         show_column_names = F)

mouse_cd4_tisTreg_df = mouse_cd4_tisTreg_df %>% dplyr::select(p_val:distance)
for(i in names(subset_list)){
  mouse_cd4_tisTreg_df[, i] = mouse_cd4_tisTreg_df$feature %in% subset_list[[i]]
}
#write_csv(mouse_cd4_tisTreg_df, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/mouse_scATAC_tisTreg_sig_subsets.csv')

mouse_cd4_tisTreg_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/mouse_scATAC_tisTreg_sig_subsets.csv')

mouse_cd4_tisTreg_table_2c = mouse_cd4_tisTreg_df %>% dplyr::select(-comparison) %>% 
  dplyr::relocate(feature, p_val, p_val_adj, avg_logFC, pct.1, pct.2)
#write_csv(mouse_cd4_tisTreg_table_2c, paste0(storeFigPath, 'table_2c_tisTregST2_signature_subsets.csv'))

# c3peaks = mouse_cd4_tisTreg_df$feature[mouse_cd4_tisTreg_df$c3_treg]
# c1peaks = mouse_cd4_tisTreg_df$feature[mouse_cd4_tisTreg_df$c1_myeloid]
# proj = archr_add_peak_signatures(proj, signature_list = list(c1 = convert_granges(c1peaks), c3 = convert_granges(c3peaks)), signature_name = 'tisTreg_subsets')
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_plot_features(reduction = dr_df, meta_features = c('z_c1', 'z_c3'), cap_top = 'q99', cap_bottom = 'q01')
# 
# mean_expr = apply(mean_mat[subset_list$c1_myeloid, colnames(mean_mat) %in% c('Macrophage','DC','Monocyte')],1,mean)
# mean_expr_other = apply(mean_mat[subset_list$c1_myeloid, !colnames(mean_mat) %in% c('Macrophage','DC','Monocyte')],1,mean)
# sd_expr = apply(mean_mat[subset_list$c1_myeloid, c('Macrophage','DC','Monocyte')],1,sd)
# table(mean_expr > 3*sd_expr + mean_expr_other)
# res_narrow = names(mean_expr)[mean_expr > 3*sd_expr + mean_expr_other]
# diff_res = mean_expr - mean_expr_other
# diff_res = sort(diff_res, decreasing = T)
# diff_res = head(diff_res[names(diff_res) %in% res_narrow], 50)
# gene_vec = table(mouse_cd4_tisTreg_df[mouse_cd4_tisTreg_df$feature %in% names(diff_res) & mouse_cd4_tisTreg_df$distance < 2000 , 'symbol']) %>% sort
# 
# 
# proj$peripheral_tissue = ifelse(proj$Tissue == 'Spleen', 'spleen',  'peripheral')
# proj$cluster_annotation_level00 = from_to(proj$cluster_annotation_level0, old_new_map_vec = c(Macrophage = "Macrophage/Monocyte", Monocyte = "Macrophage/Monocyte"))
# proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level00, sep='__')
# dr_df = jj_get_reduction_coords(proj, redname='UMAP')
# comparisons_keep = dr_df %>% dplyr::group_by(cluster_annotation_level00) %>% dplyr::count(peripheral_tissue) %>% 
#   tidyr::spread(peripheral_tissue, n) %>% 
#   dplyr::filter(peripheral > 200 & spleen > 200 & !cluster_annotation_level00 %in%  c('undefined')) %>% 
#   dplyr::pull(cluster_annotation_level00)
# groups_select = paste0(rep(c('spleen__', 'peripheral__'), each = length(comparisons_keep)), comparisons_keep)
# p = plotBrowserTrack(proj, geneSymbol = 'Rptor', groupBy = 'peripheral_cell_type', useGroups = groups_select)
# plot(p[[1]])

## homer results dotplot

mouse_cd4_tisTreg_df = as.data.frame(read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/mouse_scATAC_tisTreg_sig_subsets.csv'))

groups = c('c1_myeloid','c2_cd8_nk', 'c3_treg', 'c4_plasma', 'c5_ilc_th17')
# for(i in groups){
#   jj_save_bed(convert_granges(mouse_cd4_tisTreg_df[mouse_cd4_tisTreg_df[, i], 'feature']), file_name = paste0(storeFigPath, i, '.bed'))
# }
#specify the paths in /scATAC/scripts/config_files/config_downstream_analysis.yaml for findMotifsGenome
#option -mask was used (in contrast to TE regions only)

#for the subsets use:
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-homer2/'
#for only tisTreg and naive Treg parts, use:
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-11-finalizing_te_analysis/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
   #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_list[[i]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
sapply(homer_list, nrow)
homer_motif_list = lapply(homer_list, '[[', 1)
jj_plot_upsetr(homer_motif_list)

# #TFs that are significant only in one subset
# diff_list = list()
# for(i in seq_along(homer_motif_list)){
#   diff_list[[i]] = setdiff(homer_motif_list[[i]], unlist(homer_motif_list[-i]))
# }
# names(diff_list) = names(homer_motif_list)

homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]

# res = as.matrix(pivot_wider(gsea_res_comb, id_cols = motif, names_from =  name, values_from = enrichment, values_fill = NA))
# rnames = res[, 1]
# res = res[, -1]
# res = apply(res, 2, as.numeric)
# rownames(res) = rnames
# Heatmap(res)

gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
#for subsets use
gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'c1_myeloid'= '1',
  'c2_cd8_nk'= '2',
  'c3_treg'= '3',
  'c4_plasma'= '4',
  'c5_ilc_th17'= '5'
))
#for only tisTreg/naive Treg parts use:
gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'mouse_cd4_naiveTreg_peaks'= 'naive Treg',
  'mouse_cd4_tisTreg_peaks'= 'tisTregST2'
))

# 2D Mouse atlas tisTreg signature subset homer known motif heatmap (and S7H)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_subset_homer.pdf'), width =5, height=8)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_complete_homer.pdf'), width =5, height=8)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'Signature subset', y = '', fill = 'Enrichment') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# #exclude TFs that are unique for one subset
# gsea_res_comb = gsea_res_comb[!gsea_res_comb$motif %in% unlist(diff_list), ]
# #only unique
# gsea_res_comb = gsea_res_comb[gsea_res_comb$motif %in% unlist(diff_list), ]

#gsea_res_comb$id <- factor(gsea_res_comb$id, levels = cell_type_names)
#http://homer.ucsd.edu/homer/motif/motifDatabase.html

#gsea_res_comb$motif <- factor(gsea_res_comb$motif, levels = rev(sort(unique(gsea_res_comb$motif))))
# gsea_res_comb$rank[gsea_res_comb$rank > 10] = 10
# levels_use = gsea_res_comb %>% dplyr::group_by(motif) %>% dplyr::summarise(rank_sum = (6-n())*10 + sum(rank), 
#                                                                            enr_sum = sum(enrichment)) %>% 
#   dplyr::arrange( enr_sum) %>% pull(motif) #desc(rank_sum),
# gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
# gsea_res_comb$rank = as.character(gsea_res_comb$rank)
# gsea_res_comb$rank[gsea_res_comb$rank == '10'] = '>=10'
# gsea_res_comb$rank = factor(gsea_res_comb$rank, levels = c(as.character(1:9), '>=10'))
# 
# library(circlize)
# myColors <- colorRamp2(c(1,10), colors = c('darkblue', 'grey60'))
# col_scale = myColors(1:10)
# names(col_scale) = levels(gsea_res_comb$rank) 
# pdf(paste0(storeFigPath, 'homer_motifs_on_tisTreg_sig_subsets_dotplot.pdf'), width =6, height=8)
# ggplot(gsea_res_comb) +
#   aes(x = name, y = motif, colour = rank, size = enrichment) +
#   scale_colour_manual(values=col_scale) +
#   #scale_size_continuous(range = c(0.5, 3)) +
#   scale_size_continuous(range = c(.5, 3), 
#                         breaks = seq(2, 4.5, by = 0.5)) + 
#   geom_point() +
#   guides(colour = guide_legend(title = 'rank', order = 1),
#          size = guide_legend(title = "enrichment", order = 2, override.aes = list(size =  seq(.5, 3, length.out = 6))),
#          shape = "none") +
#   xlab("Signature subset") +
#   ylab("Homer motif") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 1, vjust = 0.5),
#         axis.text.y = element_text(colour = "black"),
#         axis.ticks = element_line(colour = "black"))
# dev.off()


# parse homer de novo motif results for tisTreg signature subsets ---------

#use R 4.2.0
library(ggplot2)
library(ggseqlogo)
library(flextable)
library(patchwork)
set_flextable_defaults(font.family = '') #standard font not found when saving as pdf

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-homer2/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)

for(i in homer_subsets){
  pwm_list = pwm_gg_list = match_df_list = match_df_gg = combined_gg = list()
  for(j in seq(10)){
    similar_motifs = readLines(sprintf('%s%s/homerResults/motif%s.info.html',homer_res, i, j))
    #get the statistics for the de novo motif
    motif_pval = gsub('.*(1e-[0-9]+).*', '\\1', grep('<TD>p-value:</TD>', similar_motifs, value = T))
    motif_target = gsub('.*?([0-9]+\\.[0-9]+\\%).*', '\\1', grep('Percentage of Target Sequences with motif', similar_motifs, value = T))
    motif_bg = gsub('.*?([0-9]+\\.[0-9]+\\%).*', '\\1', grep('Percentage of Background Sequences with motif', similar_motifs, value = T))
    title_use = sprintf('Rank: %i, P-value: %s, Target: %s, Background: %s', j, motif_pval, motif_target, motif_bg)
    #rows are the positions specific probabilities for each nucleotide (A/C/G/T)
    #read the pwm
    pwm_list[[j]] = t(read.delim(sprintf('%s%s/homerResults/motif%s.motif',homer_res, i, j),
                                 header = F, skip = 1, sep = '\t', col.names = c('A', 'C', 'G', 'T')))
    #plot sequence logo
    pwm_gg_list[[j]] = ggplot() + geom_logo(pwm_list[[j]], method = 'probability' ) + theme_logo() + ggtitle(title_use) + theme(plot.title = element_text(size = 10), text = element_text(size=10)) #axis.title=element_text(size=10), axis.text=element_text(size=12))
    #get top known motif matches and their match scores
    match_scores = as.numeric(gsub('^.*([0,1]\\.[0-9]+).*$', '\\1', grep('<TR><TD>Score:</TD><TD>', similar_motifs, value = T)))
    match_names = gsub('<H4>(.*)</H4>','\\1', grep('<H4>', similar_motifs, value = T))
    #match_names = stringr::str_pad(match_names, side = 'right', width = 70)
    match_df_list[[j]] = data.frame(`Known motif matches` = match_names, `Match score` = match_scores)
    match_df_gg[[j]] = gen_grob(flextable(match_df_list[[j]], cwidth = 1)) #ggtexttable(match_df_list[[j]]) #ggplot() + theme_void() + annotate(geom = "table", x= 1, y = 1,  table.hjust=0, label = list(match_df_list[[j]]))
    combined_gg[[j]] = pwm_gg_list[[j]] +  match_df_gg[[j]] # + plot_layout(widths=c(1,2))
  }
  all_gg = combined_gg[[1]] / combined_gg[[2]] / combined_gg[[3]] / combined_gg[[4]] / combined_gg[[5]]
  # S2G S2H S2I Mouse atlas tisTreg signature subset de novo homer results
  pdf(paste0(storeFigPath, i, '_top5_de_novo.pdf'), width = 9, height = 9)
  print(all_gg)
  dev.off()
} 

# t/nk/ilc subset umap ----------------------------------------------------

proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

# cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/treg_fine_annotation.csv')
# cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cluster_annotation_fine %in% c('tisTregST2_prog', 'Treg_naive', 'tisTregST2','pTreg'), ]
# #jj_compare_vectors(cluster_annotation_fine_df$cellNames, proj$cellNames[proj$cluster_annotation_level3 == 'Treg'])
# cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cellNames %in% proj$cellNames[proj$cluster_annotation_level3 == 'Treg'], ]
# proj$cluster_annotation_fine = proj$cluster_annotation_level3
# proj$cluster_annotation_fine[match(cluster_annotation_fine_df$cellNames, proj$cellNames)] = cluster_annotation_fine_df$cluster_annotation_fine
# 
# #jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level3')
# proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'undef'] = 'undefined'
# proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'
# proj$cluster_annotation_fine[proj$cluster_annotation_fine == 'CD4_other'] = 'Th17_Areg'
# saveArchRProject(proj)

dr_df = jj_get_reduction_coords(proj, 'UMAP')

dr_df$Clusters_1.2 = mixsort_factor(dr_df$Clusters_1.2)
gg = jj_plot_features(reduction=dr_df, meta_features='Clusters_1.2',  
                      return_gg_object = T, label_type = 'geom_label', label_col = 'white')
# S1K Mouse tnkilc subset clusters
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tnkilc_clusters.pdf'),  width = 8, height = 6)
gg
dev.off()

### Tissue umap
gg = jj_plot_features(dr_df, features='Tissue', 
                      custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                      return_gg_object = T) 
# 1G Mouse tnkilc subset tissue umap
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tnkilc_tissue.pdf'),  width = 8, height = 6)
png(paste0(storeFigPath, 'mouse_imm_cell_atlas_tnkilc_tissue.png'), width = 7, height = 6, res = 400, units = 'in')
gg
dev.off()

#jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_level3')
dr_df$cluster_annotation_fine[dr_df$cluster_annotation_fine == 'undef'] = 'undefined'
dr_df$cluster_annotation_fine[dr_df$cluster_annotation_fine == 'Th_Il21'] = 'Tfh-like'
gg = jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_fine', label_type = 'geom_label',
                 custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                 return_gg_object = T)[[1]] + labs(colour = 'Cell type')
# 1H Mouse tnkilc subset cell type umap  
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tnkilc_annotation.pdf'),  width = 8, height = 6)
gg  
dev.off()

# ### not included: singler prediction
# dr_df$singler_prediction = replace_if(dr_df$singler_label_fine, count_below = 20, 'other')
# cols_use = jj_get_jj_colours(dr_df$singler_prediction)
# gg = jj_plot_features(reduction=dr_df, meta_features='singler_prediction', 
#                       pt.size = 0.5, return_gg_object = T, label_type = 'geom_label_repel', 
#                       custom_colors = cols_use, label_subset = names(cols_use)[!names(cols_use) %in% c('other','B cells')], 
#                       label_col = cols_use)
# 
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_t_nk_ilc_singler.pdf'),  width = 12, height = 8)
# gg[[1]] + labs(colour='SingleR prediction') 
# dev.off()


### marker genes plot
genes_plot = c( 'Entpd1', 'Maf', 'Irf4', 'Pparg','Col15a1', 'Vps8')
genes_plot = c('Areg','Batf', 'Bcl6', 'Cd3e', 'Cd8a', 'CTLA4','Cxcr5',
               'Foxp3', 'Gata3','Gzmb', 'Ifng', 'Ikzf2','Il10', 'Il17a',
               'Il1rl1','Il21', 'Il2ra','Il4','Klrb1c', 'Klrg1','Maf','Rora',
               'Rorc','Tbx21','Tox','Il5','Il13')
genes_plot = c( 'Emp1', 'Krt28', 'Ctsh')
gg = archr_plot_markers(proj, genes_plot)
# S2B 3L Mouse tnkilc subset marker gene umaps
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tnkilc_subset_magic_marker_umap.pdf'),  width = 7.5, height = 6)
gg
dev.off()

# ### not included: marker heatmap and dotplot
# 
# gmat = get_gene_mat(proj)
# 
# genes_plot = tools::toTitleCase(tolower(c('Areg','Batf','CCR7','Il7r','CD3e','CD4','CD44','Ctla4','Itga1','CD8a','Cxcr5',
#                'Eomes','Foxp3','Gata3','Gzmb','Icos','Ifng','Ikzf2','Il10','Il13','Il17a','Il1rl1',
#                'Il21','Il2ra','Il4','Il5','Klrb1c','Klrg1','Maf','Pdcd1','Ncr1','Rorc','Sell',
#                'Tbx21','Rora','Tox', 'Tra','Trb','Trd')))
# genes_plot[!genes_plot %in% rownames(gmat)]
# 
# # custom_regions = c(Tcra = 'chr14-52427967-54224198', Tcrb='chr6-40891296-41558371', Tcrd = 'chr14-53946073-54159198')
# # custom_regions = sort(convert_granges(custom_regions))
# # strand(custom_regions) = '+'
# # fmat = getFeatureMatrix(proj, custom_regions)
# # rownames(fmat) = names(custom_regions)
# # identical(colnames(gmat), colnames(fmat))
# gmat_use = gmat[rownames(gmat) %in% genes_plot, ]
# #gmat_use = rbind(gmat_use, fmat)
# gmat_use = gmat_use[, dr_df$cluster_annotation_fine != 'undefined']
# gg = jj_plot_heatmap(obj = gmat_use, group_vec = dr_df$cluster_annotation_fine[dr_df$cluster_annotation_fine !='undefined'], 
#                      features_use = rownames(gmat_use), scale_data = T)
# 
# 
# pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_t_nk_ilc_marker_heatmap.pdf'),  width = 10, height = 5)
# gg
# seurat_dotplot(gmat=gmat_use, metadf = dr_df[dr_df$cluster_annotation_fine != 'undefined', ],
#                features = rownames(gmat_use),
#                group_column = 'cluster_annotation_fine')
# dev.off()


# mouse normal cd4 tisTreg signature introduction -------------------------

#execute once
##mouse normal cd4
# pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal_cd4.yaml')
# addArchRGenome(pconfig$GENOME)
# bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
# setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
# proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_mouse_normal_CD4'))
# proj_use = proj[proj$annotation %in% c('Treg_naive', 'tisTregST2_prog', 'tisTregST2'), ]
# proj_use = archr_dim_reduction(proj_use)
# proj_use = archr_clustering(proj_use, 0.5)
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
# dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
# gmat = get_gene_mat(proj_use)
# library(monocle3)
# so = CreateSeuratObject(counts = gmat, meta.data = dr_df)
# so = FindVariableFeatures(so) 
# so = ScaleData(so, features = VariableFeatures(so))
# gmat_use = GetAssayData(so, slot = 'scale.data')
# #gmat_use <- gmat[rownames(gmat) %in% VariableFeatures(so), ]
# pd_df <- droplevels(dr_df)
# fd_df <- data.frame(gene_short_name=rownames(gmat_use),
#                     row.names = rownames(gmat_use))
# rownames(gmat_use) <- fd_df$gene_short_name
# 
# #MAKE SURE monocle3 is loaded and monocle is unloaded! Otherwise error later in plot_cells
# cds <- new_cell_data_set(gmat_use,
#                          cell_metadata = pd_df,
#                          gene_metadata = fd_df)
# 
# #performs log-normalization + PCA
# cds <- preprocess_cds(cds, num_dim = 50)
# #umap with cosine metric
# cds <- reduce_dimension(cds)
# #cds <- partitionCells(cds)
# 
# reducedDims(cds)$UMAP <- dr_df[, 1:2]
# #colnames(reducedDims(cds)$UMAP) <- c('V1', 'V2') #otherwise error in orderCells plotly picker
# cds <- cluster_cells(cds)
# 
# ## Step 5: Learn a graph
# cds <- learn_graph(cds, use_partition = T)
# cds <- order_cells(cds)#, root_cells = 'AAACGAAAGACTAATG-7') #root pickec at -5, -3
# cds@colData$sample_name = NULL
# #write_rds(cds, paste0(bigFilesDir, 'mouse_cd4_treg_subset_monocle3_cds.RDS'))
#
# identical(colnames(cds), proj_use$cellNames)
# proj_use$pseudotime = pseudotime(cds, reduction_method = 'UMAP')
# dr_df = jj_get_reduction_coords(proj_use, 'UMAP')
# pmat = get_peak_mat(proj_use)
# mouse_cd4_tisTreg_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-10-31-tisTreg_sig_subsets/mouse_scATAC_tisTreg_sig_subsets.csv')
# mouse_cd4_tisTreg_df = mouse_cd4_tisTreg_df[mouse_cd4_tisTreg_df$comparison == '16_23_versus_0_3_14', ]
# mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature[])
# emb_df = getReducedDims(proj_use)
# polap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, #dr_df[, 1:2],
#                                   signature_gr = mouse_cd4_tisTreg_gr, nFrags_vec = proj_use$nFrags, count_thres = 2e5, verbose = T)
# proj_use$signature_pct_overlap = polap_df$signature_pct_overlap
# proj_use$signature_pct_overlap = polap_df$signature_pct_overlap
# proj_use$signature_n_overlap = polap_df$signature_n_overlap
# proj_use$n_cells_aggregated = polap_df$n_cells_aggregated
# proj_use$nFrags_aggregated = polap_df$nFrags_aggregated
# #saveArchRProject(proj_use, 'mouse_cd4_treg_subset_monocle3_proj.RDS')

setwd('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal_cd4/')
proj_use = loadArchRProject('mouse_cd4_treg_subset_monocle3_proj.RDS')
cds = read_rds(paste0(bigFilesDir,  'mouse_cd4_treg_subset_monocle3_cds.RDS'))
dr_df = jj_get_reduction_coords(proj_use, 'UMAP')

gg = jj_plot_features(dr_df, features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'), return_gg_object = T)[[1]]
cont_df = jj_get_contour_lines(dr_df, 'annotation', .95)
# S2E Mouse CD4 treg subset annotation umap
pdf(paste0(storeFigPath, 'mouse_cd4_treg_subset_annotation_umap.pdf'), width =8, height=6)
gg + geom_path(data = cont_df, aes(x=x, y=y, group=cont_group, linetype = annotation), #, colour=snn_harmony_res.1), 
               size=1) + labs(linetype = 'Cell type')
dev.off()

gg = jj_plot_features(reduction=dr_df, meta_features = 'pseudotime', colorScale = 'bry', return_gg_object = T)[[1]]
gg2 =   monocle3::plot_cells(cds, color_cells_by = 'pseudotime', label_leaves = F, label_branch_points = F,trajectory_graph_color = 'black',  label_roots = F, cell_size = 0.5) + 
  theme_minimal() + coord_fixed()
gg$layers[[2]] = gg2$layers[[2]]

# S2D Mouse CD4 treg subset trajectory umap
pdf(paste0(storeFigPath, 'mouse_cd4_treg_subset_trajectory_umap.pdf'), width =8, height=6)
gg + labs(colour='Pseudotime')
dev.off()

# S2C Mouse CD4 treg subset tisTreg signature percentage overlap umap
pdf(paste0(storeFigPath, 'mouse_cd4_treg_signature_overlap_umap.pdf'), width =8, height=6)
  jj_plot_features(dr_df, features = 'signature_pct_overlap', return_gg_object = T)[[1]] + labs(colour='% Overlap')
dev.off()
# S2F Mouse CD4 treg subset pseudotime-percentage overlap correlation scatterplot
pdf(paste0(storeFigPath, 'mouse_cd4_treg_correlation.pdf'), width =8, height=6)
ggplot(dr_df, aes(x = pseudotime, y = signature_pct_overlap)) + geom_point(aes(colour=annotation), size = 1) + theme_minimal() +
  geom_smooth(method='lm', formula= y~x, colour='black') + labs(x='Pseudotime', y = '% Overlap', colour='Cell type') + 
  annotate('text', x = 2.5, y = 70, label = paste0('Pcor: ', round(cor(dr_df$pseudotime, dr_df$signature_pct_overlap),2)))
#scale_color_manual(values = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
dev.off()


# cell type markers -------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
#contains the scATAC tisTreg regions as peakset
###
#full atlas
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected', pconfig$DONOR))
# #or full atlas, but with scATAC signature as peak matrix
# proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets', pconfig$DONOR))
# cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/treg_fine_annotation.csv')
# stopifnot(identical(proj$cellNames, cluster_annotation_fine_df$cellNames))
# proj$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
# #or t/nk/ilc subset only
# proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
# cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-23-treg_subsets/treg_fine_annotation.csv')
# cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cluster_annotation_fine %in% c('tisTregST2_prog', 'Treg_naive', 'tisTregST2','pTreg'), ]
# cluster_annotation_fine_df = cluster_annotation_fine_df[cluster_annotation_fine_df$cellNames %in% proj$cellNames[proj$cluster_annotation_level3 == 'Treg'], ]
# proj$cluster_annotation_fine = proj$cluster_annotation_level3
# proj$cluster_annotation_fine[match(cluster_annotation_fine_df$cellNames, proj$cellNames)] = cluster_annotation_fine_df$cluster_annotation_fine

###

rm(markers_peaks_se, markers_se)
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = c('tisTregST2','Th_Il21', 'Macrophage', 'CD4_other', 'B cell', 'Plasma cell', 'CD4_Tnaive'), #without macrophage for t/nk/ilc subset
  bgdGroups = NULL,
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix',#'PeakMatrix', 
  useGroups = c('tisTregST2','Th_Il21', 'Macrophage', 'CD4_other','B cell', 'Plasma cell', 'CD4_Tnaive'),
  bgdGroups = NULL,
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_peaks_list = markers_list = list()
markers_peaks_list[['many']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
markers_list[['many']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)

dr_df_tnkilc = jj_get_reduction_coords(loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected'))
c1_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C1']
c2_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C2']
#cd4_other_colon = dr_df_tnkilc$cellNames[dr_df_tnkilc$Tissue == 'Colon' & dr_df_tnkilc$Clusters_1.2 == 'C3'] #except for 1 cell colon anyway
th17_colon = dr_df_tnkilc$cellNames[dr_df_tnkilc$Tissue == 'Colon' & dr_df_tnkilc$Clusters_1.2 == 'C8']
# proj$cluster_annotation_fine[proj$cellNames %in% c1_bc] = 'ILC2_C1'
# proj$cluster_annotation_fine[proj$cellNames %in% c2_bc] = 'ILC2_C2'
# proj$cluster_annotation_fine[proj$cellNames %in% cd4_other_colon] = 'CD4_other_colon'
# proj$cluster_annotation_fine[proj$cellNames %in% th17_colon] = 'Th17_colon'
#'ILC2_C1'
rm(markers_peaks_se, markers_se)
proj$cluster_annotation_use = proj$cluster_annotation_fine
proj$cluster_annotation_use[proj$cellNames %in% c1_bc ] = 'ILC2_C1'
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'ILC2_C1', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix', 
  useGroups = 'ILC2_C1', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_peaks_list[['ILC2_C1']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
markers_list[['ILC2_C1']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
#'ILC2_C2',
rm(markers_peaks_se, markers_se)
proj$cluster_annotation_use = proj$cluster_annotation_fine
proj$cluster_annotation_use[proj$cellNames %in% c2_bc ] = 'ILC2_C2'

markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'ILC2_C2', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix', 
  useGroups = 'ILC2_C2', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_peaks_list[['ILC2_C2']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
markers_list[['ILC2_C2']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)

##Th17_colon
rm(markers_peaks_se, markers_se)
proj$cluster_annotation_use = proj$cluster_annotation_fine
proj$cluster_annotation_use[proj$cellNames %in% th17_colon ] = 'Th17_colon'
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'Th17_colon', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix', 
  useGroups = 'Th17_colon', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_use",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_peaks_list[['Th17_colon']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
markers_list[['Th17_colon']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)

## for tisTreg vs Th2
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'Th2', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix', 
  useGroups = 'Th2', 
  bgdGroups = NULL,
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_peaks_list[['Th2']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
markers_list[['Th2']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)


# ## for tisTreg vs pTreg
# proj$cluster_annotation_use = proj$cluster_annotation_fine
# proj$cluster_annotation_use[proj$cluster_annotation_use == 'tisTregST2' & proj$Tissue == 'Colon'] = 'Colon_tisTregST2'
# proj$cluster_annotation_use[proj$cluster_annotation_use == 'pTreg' & proj$Tissue == 'Colon'] = 'Colon_pTreg'
# markers_peaks_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'PeakMatrix', 
#   useGroups = 'Colon_pTreg', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'GeneScoreMatrix', 
#   useGroups = 'Colon_pTreg', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_peaks_list[['Colon_pTreg']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
# markers_list[['Colon_pTreg']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# rm(markers_peaks_se, markers_se)
# markers_peaks_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'PeakMatrix', 
#   useGroups = 'Colon_tisTregST2', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'GeneScoreMatrix', 
#   useGroups = 'Colon_tisTregST2', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_peaks_list[['Colon_tisTregST2']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
# markers_list[['Colon_tisTregST2']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# 


marker_df = do.call(rbind, markers_list)
marker_peaks_df = do.call(rbind, markers_peaks_list)
#write_rds(list(ga=marker_df, peaks=marker_peaks_df), '/omics/groups/OE0436/internal/msimon/scATAC/markers_tisTreg_ilc2_tfh_th17areg.RDS')

marker_df_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_tisTreg_ilc2_tfh_th17areg.RDS')
marker_df_list$ga = rbind(marker_df_list$ga, marker_df) #add in th2 results
marker_df_list$peaks = rbind(marker_df_list$peaks, marker_peaks_df)
marker_df = marker_df_list$ga

marker_df_sig = marker_df[marker_df$FDR <= 0.01 & marker_df$Log2FC > 0.5, ]
#jj_save_excel(split(marker_df_sig, marker_df_sig$comparison), paste0(storeFigPath, 'atlas_marker_genes_fdr001_log2fc05.xlsx'))

# markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_tisTreg_ilc2_tfh_th17areg.RDS')
# markers_se = markers_list$ga
# markers_peaks_se = markers_list$peaks
# marker_df = archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 0.01 & Log2FC >=0.5", annotate_closest_gene = F)
# table(marker_df$comparison)
# #jj_save_excel(split(marker_df, marker_df$comparison), paste0(storeFigPath, 'atlas_marker_genes_fdr001_log2fc05.xlsx'))
# marker_list = split(marker_df, marker_df$comparison)
# jj_plot_upsetr(lapply(marker_list, '[[', 'name'))
# #volcano_df = dplyr::inner_join(marker_list$tisTregST2, marker_list$ILC2, by = 'name', suffix = c('.tisTregST2', '.ILC2'))
# #jj_fc_fc_plot(volcano_df, 'Log2FC.tisTregST2', 'Log2FC.ILC2', 'name', labs_range = c(0, 4), marker_thres = 2, use_text = T)
# 
# marker_df = archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# marker_peaks_df = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)

markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_tisTreg_ilc2_tfh_th17areg.RDS')
marker_df = markers_list$ga

cdf = tibble(
  group1 = 'tisTregST2',
  group2 = c('ILC2_C1','Th_Il21','CD4_other','Macrophage', 'B_cell', 'CD4_Tnaive', 'ILC2_C2', 'Th2'),
  markers_highlight = list( c('Batf', 'Areg','Il1rl1','Gata3','Klrg1'), 
                            c('Ctla4', 'Ikzf2','Maf','Batf','Irf8'),
                            c('Batf', 'Areg','Rorc'),
                            c('Nfil3','Entpd1','Havcr2','Irf4','Maf'),
                            c('Cd19','Ms4a1','Foxp3','Klrg1','Il1rl1'),
                            c('Cd28','Ccr4','Ikzf4'),
                            c('Batf', 'Areg','Il1rl1','Gata3','Klrg1'),
                            c('Batf', 'Areg','Il1rl1','Gata3','Klrg1', 'Rora','Ccr8','Tigit'))
)
cdf = cdf %>% add_row(group1='ILC2_C1', group2='ILC2_C2', markers_highlight = list(c('Gata3','Il12r','Klrg1','Il2ra','Il7r','Rora','Areg','Il4','Il13')))
cdf = cdf %>% add_row(group1='CD4_other', group2='Th17_colon', markers_highlight = list(c('Rora', 'Il17a','Rorc','Il17f','Il23r','Il22','Il21','Ccl6','Batf','Cd3e','Cd4')))
cdf = cdf %>% add_row(group1='tisTregST2', group2='Plasma_cell', markers_highlight = list(c('Ctla4','Tigit', 'Il1b','Csf1r','Ccl4','Clec5a')))
cdf = cdf %>% add_row(group1='B_cell', group2='Plasma_cell', markers_highlight = list(c('Cd19','Cd79a','Cd79b','Sdc1','Sell','Ly6d', 'Il5ra', 'Mzb1','Nkg7','Clec4d','Ccl2','Ikzf2','Gzme')))


marker_df$comparison[marker_df$comparison == 'B cell'] = 'B_cell'
marker_peaks_df$comparison[marker_peaks_df$comparison == 'B cell'] = 'B_cell'
marker_df$comparison[marker_df$comparison == 'Plasma cell'] = 'Plasma_cell'
marker_peaks_df$comparison[marker_peaks_df$comparison == 'Plasma cell'] = 'Plasma_cell'
#highlight_thres = 3
gg_list = list()
for(i in 1:nrow(cdf)){
  gg_list[[i]] = fc_fc_plot_broad(marker_df, #marker_peaks_df, 
                                  group1 = unlist(cdf[i, 'group1']), 
                   group2 = unlist(cdf[i, 'group2']), markers_highlight = unlist(cdf[i, 'markers_highlight']),
                   symbol_column = 'name', log2FC_thres = 0.5, labs_range = c(-5,5,-5,5))
}
# 3C, 3D, 3E, 3J, 4B, 4C, 4F fc-fc scatterplots tisTreg vs others
pdf(paste0(storeFigPath, 'tisTregST2_other_cell_types_fc_fc_plots.pdf'),  width = 10, height = 8)
gg_list
dev.off()


table_3_4_cell_type_marker_genes_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-21-te_browsertracks/atlas_marker_genes_fdr001_log2fc05.xlsx')
names(table_3_4_cell_type_marker_genes_list) = from_to(vec= names(table_3_4_cell_type_marker_genes_list), old_new_map_vec=c(
  'CD4_other'= 'Th17_Areg',
  'Th_Il21'= 'Tfh-like'
))
table_3_4_cell_type_marker_genes_list = table_3_4_cell_type_marker_genes_list[c('tisTregST2','ILC2_C1', 'ILC2_C2','Th17_colon','Th17_Areg', 'Tfh-like', 'B cell')]
#jj_save_excel(table_3_4_cell_type_marker_genes_list, paste0(storeFigPath, 'table_3_4_cell_type_marker_genes.xlsx'))

# marker_df_list = split(marker_df, marker_df$comparison)
# marker_df_list = marker_df_list[names(marker_df_list) %in% c('B_cell', 'Plasma_cell')]
# marker_df_use = marker_df_list[[1]] %>% dplyr::inner_join(marker_df_list[[2]], by = 'name', suffix=c("B_cell",'Plasma_cell'))
# jj_fc_fc_plot(marker_df_use, 'Log2FCB_cell', 'Log2FCPlasma_cell', 'name', labs_range = c(0,5), marker_thres = 1.5, use_text = T)

#points on the axis: only significant and above minimum log2FC threshold in one of the comparisons
# ggplot(volcano_df, aes_string(x =  lfc_g1, y = lfc_g2))  + theme_minimal() + 
#   annotate("rect", xmin = labs_range[1], xmax = labs_range[2]-0.5, ymin = -0.5, ymax = 0.5,
#            alpha = 1,fill = "grey90") + 
#   annotate("rect", xmin = -0.5, xmax = 0.5, ymin = labs_range[3], ymax = labs_range[4]-0.5,
#            alpha = 1,fill = "grey90") + 
#   annotate(geom = "polygon", x = c(labs_range[2]-0.5, labs_range[2]-0.5, labs_range[2]), y = c(-0.5, 0.5, 0), fill = "grey90", alpha = 1 ) + 
#   annotate(geom = "polygon", x = c(-0.5, 0.5, 0), y = c(labs_range[4]-0.5, labs_range[4]-0.5, labs_range[4]), fill = "grey90", alpha = 1 ) + 
#   geom_point(size = 1, aes(shape=has_marker_peak), show.legend = T) + 
#   scale_shape_manual(values=c(1, 19)) + 
#   coord_fixed(xlim = labs_range[1:2], ylim = labs_range[3:4], 
#                 expand = F, clip = "off") + 
#   scale_x_continuous(breaks=c(-rev(seq(0.5,5,1)), seq(0.5,5,1))) + 
#   scale_y_continuous(breaks=c(-rev(seq(0.5,5,1)), seq(0.5,5,1)))  + 
#   geom_point(data = volcano_df[volcano_df$label_use != '', ], aes(shape = has_marker_peak), size = 1.5, colour = 'red', show.legend = F) +
#   ggrepel::geom_text_repel(aes_string(label = "label_use"), max.overlaps = 500) 


# volcano_df = volcano_df %>% dplyr::left_join(marker_peaks_count_df, by = c('name' = 'symbol'))
# volcano_df$norm_marker_peaks[is.na(volcano_df$norm_marker_peaks)] = 0  
# volcano_df$norm_marker_peaks[volcano_df$norm_marker_peaks > 1] = 1  
# #points on the axis: only significant and above minimum log2FC threshold in one of the comparisons


#write_csv(marker_df_all, paste0(storeFigPath, 'ILC2_C1_vs_tisTreg_marker_peaks.csv'))
jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR',labs_range = c(-8,8), 
                symbol_column = 'symbol', marker_thres = c(Inf, 30), pt.size = 0.8)



# # compare tissue celltype against naive spleen cell type ------------------
# 
# #results in only few marker genes
# # #tissue macrophage vs spleen monocyte
# # #tisTregST2 vs naive treg spleen
# # #ILC2 vs CD4_Tnaive
# # #CD4 other vs CD4_Tnaive
# # #Th_il21 vs CD4_Tnaive
# # 
# # proj$peripheral_tissue = ifelse(proj$Tissue == 'Spleen', 'spleen',  'peripheral')
# # #proj$cluster_annotation_level00 = from_to(proj$cluster_annotation_fine, old_new_map_vec = c(Macrophage = "Macrophage/Monocyte", Monocyte = "Macrophage/Monocyte"))
# # proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_fine, sep='__')
# # 
# comparisons_df = data.frame(
#   g1 = paste0('peripheral__', c('tisTregST2','ILC2_C1','Th_Il21','CD4_other','Macrophage')),
#   g2 = paste0('spleen__', c('Treg_naive','CD4_Tnaive','CD4_Tnaive','CD4_Tnaive','Monocyte'))
# )
# mdf_list = list()
# for(i in 1:nrow(comparisons_df)){
#   rm(markers_se)
#   markers_se <- getMarkerFeatures(
#     ArchRProj = proj,
#     useMatrix = 'GeneScoreMatrix',#'PeakMatrix',
#     useGroups = comparisons_df$g1[i],
#     bgdGroups = comparisons_df$g2[i],
#     groupBy = "peripheral_cell_type",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon"
#   )
#   mdf_list[[comparisons_df$g1[i]]] = archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# }
# 
# mdf = do.call(rbind, mdf_list)
# # 
# # cdf = tibble(
# #   group1 = paste0('peripheral__', 'tisTregST2'),
# #   group2 = paste0('peripheral__', c('ILC2_C1','Th_Il21','CD4_other','Macrophage')),
# #   markers_highlight = list( c('Batf', 'Areg','Il1rl1','Gata3','Klrg1'), 
# #                             c('Ctla4', 'Ikzf2','Maf','Batf','Irf8'),
# #                             c('Batf', 'Areg','Rorc'),
# #                             c('Nfil3','Entpd1','Havcr2','Irf4','Maf'))
# # )
# # gg_list = list()
# # for(i in 1:nrow(cdf)){
# #   gg_list[[i]] = fc_fc_plot_broad(mdf, group1 = unlist(cdf[i, 'group1']), 
# #                                   group2 = unlist(cdf[i, 'group2']), markers_highlight = unlist(cdf[i, 'markers_highlight']),
# #                                   symbol_column = 'name', log2FC_thres = 0.5, labs_range = c(-5,5,-5,5))
# # }
# # #jj_arrange_ggplots(gg_list)
# # pdf(paste0(storeFigPath, 'tisTregST2_other_cell_types_peripheral_vs_spleen_fc_fc_plots.pdf'),  width = 10, height = 8)
# # gg_list
# # dev.off()

# ILC2 C1 vs ILC2 C2 ------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

proj_use = proj[proj$cluster_annotation_fine == 'ILC2' & proj$Tissue == 'Colon', ]
#proj_use = proj_use[proj_use$Tissue == 'Colon'] #no big difference -> not a big problem that some skin/fat cells are included
rm(markers_peaks_se, markers_se)
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'PeakMatrix', 
  useGroups = 'C1',
  bgdGroups = 'C2',
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'GeneScoreMatrix',#'PeakMatrix', 
  useGroups = 'C1',
  bgdGroups = 'C2',
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markers_TF_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'MotifMatrix', 
  useGroups = 'C1',
  bgdGroups = 'C2',
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#write_rds(list(ga=markers_se, peaks=markers_peaks_se), '/omics/groups/OE0436/internal/msimon/scATAC/markers_ILC2_C1_vs_C2.RDS')
marker_df = archr_get_markers_as_df(proj = proj_use, markers = markers_se, cutOff = "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj_use, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
#jj_save_excel(list(marker_genes=marker_df, marker_peaks = marker_peaks_df), paste0(storeFigPath, 'ILC2_C1_vs_C2_markers.xlsx'))
cat(marker_df$name[marker_df$Log2FC >= 1], sep='\n')

markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_ILC2_C1_vs_C2.RDS')
markers_se = markers_list$ga
markers_peaks_se = markers_list$peaks

marker_df = archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)

#markers_highlight = c('Itga4','Tbx21','Slamf7','Il10ra','Pdcd1','Runx3', 'Itga1', 'Tox','Ncr1','Batf','Ifng','Pparg', 'Klrb1c','Il1rl1','Areg','Ereg','Ccr8')
markers_highlight = c('Areg','Il4','Il13','Il6')

# 4G ILC2 subset volcano plot
pdf(paste0(storeFigPath, 'ILC2_volcano_plot.pdf'),  width = 8, height = 6)
#jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', 
#                labs_range = c(-4,4), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight)
#jj_volcano_plot(marker_peaks_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'symbol', 
#                labs_range = c(-6,6), pt.size = 1, marker_thres = c(Inf,10), markers_highlight = NULL)
jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', 
                labs_range = c(-4,4), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight,
                highlight.pt.size = 2, markers_highlight_col = structure(rep('red',length(markers_highlight)), names = markers_highlight),
                col_by_highlight = T, group_names = c('ILC2_C1','ILC2_C2'), arrow_pos = list(x=0.25, xend = 1, y = 45)) + guides(colour="none")

dev.off()

#compare chromvar TF activities (differential chromvar, differential TF)
marker_tf_df = archr_get_markers_as_df(proj = proj, markers = markers_TF_se, cutOff =  "FDR <= 1", annotate_closest_gene = F)
marker_tf_df$name = sapply(strsplit(marker_tf_df$name, '_'),'[[',1)
markers_highlight = marker_tf_df$name[abs(marker_tf_df$MeanDiff) >=0.075]
# S4M ILC2 TF chromvar score volcano plot
pdf(paste0(storeFigPath, 'ILC2_chromvar_volcano_plot.pdf'),  width = 8, height = 6)
jj_volcano_plot(marker_tf_df, logfc_column = 'MeanDiff', pval_column = 'FDR', symbol_column = 'name', 
                labs_range = c(-0.25, 0.25), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight,
                highlight.pt.size = 2, markers_highlight_col = structure(rep('red',length(markers_highlight)), names = markers_highlight),
                col_by_highlight = T, group_names = c('ILC2_C1','ILC2_C2'), arrow_pos = list(x=0.05, xend=0.15, y=40))+ guides(colour="none", shape = 'none')
dev.off()
# 
# p <- plotEmbedding(
#   ArchRProj = proj, 
#   colorBy = "MotifMatrix", 
#   name = c('z:Tbet.T.box_282', 'z:GATA.Zf..IR3_105'), 
#   embedding = "UMAP",
#   imputeWeights = getImputeWeights(proj)
# )
# plot(p$`z:Tbet.T.box_282`)
# plot(p$`z:GATA.Zf..IR3_105`)


# # TF comparison ILC2_C1 vs ILC2_C2 ----------------------------------------
# 
# pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
# addArchRGenome(pconfig$GENOME)
# bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
# setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
# proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
# 
# proj_use = proj[proj$cluster_annotation_fine == 'ILC2' & proj$Tissue == 'Colon', ]
# #proj_use = proj_use[proj_use$Tissue == 'Colon'] #no big difference -> not a big problem that some skin/fat cells are included
# 
# 
# res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0.05 & FDR <= 0.01'))
# 
# da_tf_list = lapply(res_se, function(x) x[[1]][, 'name'])
# jj_plot_upsetr(da_tf_list)
# motif_olap_df =jj_make_overlap_df(da_tf_list)
# motifs_plot = unique(unlist(da_tf_list))
# motifs_plot = na.omit(unique(c(unlist(motif_olap_df$overlaps[,!(colnames(motif_olap_df$overlaps) %in% c('a','b','c','d','e'))]))))
# length(motifs_plot)
# 
# res_se2 = lapply(da_list, function(x) getMarkers(x, cutOff = 'FDR <= 1'))
# diff_list = lapply(res_se2, function(x) as.data.frame(x[[1]][, c('name','MeanDiff')]))
# res_df = reduce(diff_list, full_join, by = "name") %>% column_to_rownames('name')
# colnames(res_df) = names(diff_list)
# 
# plot_mat = t(res_df[rownames(res_df) %in% motifs_plot, ])
# colnames(plot_mat) = gsub('(.*)_[0-9]+','\\1',colnames(plot_mat))
# 
# heatmap_binding_domain = rep('Other', ncol(plot_mat))
# for(i in binding_domains){
#   heatmap_binding_domain[grepl(i, colnames(plot_mat))] = i
# }
# ha = columnAnnotation('Binding Domain'= heatmap_binding_domain, col = list('Binding Domain' = jj_get_colours(heatmap_binding_domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))
# 
# 
# 
# Heatmap(plot_mat) #define range based on automatic colorscale
# col_fun = colorRamp2(c(-0.1, 0, 0.1, 0.5), c("darkblue", "white", "red", "darkred"))
# h1 <- Heatmap(plot_mat, col = col_fun,
#               name = 'MeanDiff\nperipheral - spleen',
#               column_names_gp = gpar(fontsize = 8), 
#               row_names_gp = gpar(fontsize = 8),
#               top_annotation = ha)
# pdf(paste0(storeFigPath, 'mouse_peripheral_tf_heatmap.pdf'), width =10, height=3)
# h1
# dev.off()


# cd4_other vs colon_th17 -------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

proj_use = proj[proj$cluster_annotation_fine %in% c('CD4_other','Th17') & proj$Tissue == 'Colon', ]
rm(markers_peaks_se, markers_se)
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'PeakMatrix', 
  useGroups = 'C3',
  bgdGroups = 'C8',
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'GeneScoreMatrix',#'PeakMatrix', 
  useGroups = 'C3',
  bgdGroups = 'C8',
  groupBy = "Clusters_1.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#write_rds(list(ga=markers_se, peaks=markers_peaks_se), '/omics/groups/OE0436/internal/msimon/scATAC/markers_CD4_other_C3_vs_colon_Th17_C8.RDS')
markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_CD4_other_C3_vs_colon_Th17_C8.RDS')
markers_se = markers_list$ga
markers_peaks_se = markers_list$peaks

marker_df = archr_get_markers_as_df(proj = proj_use, markers = markers_se, cutOff = "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj_use, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
#jj_save_excel(list(marker_genes=marker_df, marker_peaks = marker_peaks_df), paste0(storeFigPath, 'CD4_other_C3_vs_colon_Th17_C8_markers.xlsx'))
cat(marker_df$name[marker_df$Log2FC >= 1], sep='\n')

# markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_CD4_other_C3_vs_colon_Th17_C8.RDS')
# markers_se = markers_list$ga
# markers_peaks_se = markers_list$peaks

marker_df = archr_get_markers_as_df(proj = proj_use, markers = markers_se, cutOff = "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)

#markers_highlight = c('Rora','Il17a')
markers_highlight = c('Areg','Il4','Il13','Il6')

# 3K Th17_Areg vs Th17 volcano plot
pdf(paste0(storeFigPath, 'CD4_other_Th17_volcano_plot.pdf'),  width = 8, height = 6)
#jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', 
#                labs_range = c(-4,4), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight, highlight.pt.size = 5, col_vec = c("black", "black", "black","red"))
#jj_volcano_plot(marker_peaks_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'symbol', 
#                labs_range = c(-6,6), pt.size = 1, marker_thres = c(Inf,10), markers_highlight = NULL)
jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', 
                labs_range = c(-4,4), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight,
                highlight.pt.size = 2, markers_highlight_col = structure(rep('red',length(markers_highlight)), names = markers_highlight),
                col_by_highlight = T, group_names = c('CD4_other','Th17'), arrow_pos = c(0.25, 1, 45))+ guides(colour="none")
dev.off()



# tisTregST2 colon vs pTreg -----------------------------------------------

# ## for tisTreg vs pTreg
# proj$cluster_annotation_use = proj$cluster_annotation_fine
# proj$cluster_annotation_use[proj$cluster_annotation_use == 'tisTregST2' & proj$Tissue == 'Colon'] = 'Colon_tisTregST2'
# proj$cluster_annotation_use[proj$cluster_annotation_use == 'pTreg' & proj$Tissue == 'Colon'] = 'Colon_pTreg'
# markers_peaks_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'PeakMatrix', 
#   useGroups = 'Colon_pTreg', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'GeneScoreMatrix', 
#   useGroups = 'Colon_pTreg', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_peaks_list[['Colon_pTreg']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
# markers_list[['Colon_pTreg']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# rm(markers_peaks_se, markers_se)
# markers_peaks_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'PeakMatrix', 
#   useGroups = 'Colon_tisTregST2', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_se <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = 'GeneScoreMatrix', 
#   useGroups = 'Colon_tisTregST2', 
#   bgdGroups = NULL,
#   groupBy = "cluster_annotation_fine",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markers_peaks_list[['Colon_tisTregST2']] = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
# markers_list[['Colon_tisTregST2']] =  archr_get_markers_as_df(proj = proj, markers = markers_se, cutOff = "FDR <= 1", annotate_closest_gene = F)
# 
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

proj_use = proj[proj$cluster_annotation_fine %in% c('tisTregST2','pTreg') & proj$Tissue == 'Colon', ]
rm(markers_peaks_se, markers_se)
markers_peaks_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'PeakMatrix', 
  useGroups = 'tisTregST2',
  bgdGroups = 'pTreg',
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_se <- getMarkerFeatures(
  ArchRProj = proj_use,
  useMatrix = 'GeneScoreMatrix',#'PeakMatrix', 
  useGroups = 'tisTregST2',
  bgdGroups = 'pTreg',
  groupBy = "cluster_annotation_fine",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#write_rds(list(ga=markers_se, peaks=markers_peaks_se), '/omics/groups/OE0436/internal/msimon/scATAC/markers_colon_tisTregST2_vs_pTreg.RDS')
markers_list = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/markers_colon_tisTregST2_vs_pTreg.RDS')
markers_se = markers_list$ga
markers_peaks_se = markers_list$peaks

marker_df = archr_get_markers_as_df(proj = proj_use, markers = markers_se, cutOff = "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj_use, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)
#jj_save_excel(list(marker_genes=marker_df, marker_peaks = marker_peaks_df), paste0(storeFigPath, 'markers_colon_tisTregST2_vs_pTreg_markers.xlsx'))
cat(marker_df$name[marker_df$Log2FC >= 1], sep='\n')

marker_df = archr_get_markers_as_df(proj = proj_use, markers = markers_se, cutOff = "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = F)
marker_peaks_df = archr_get_markers_as_df(proj = proj, markers = markers_peaks_se, cutOff =  "FDR <= 0.01 & abs(Log2FC) >=0.5", annotate_closest_gene = T)

#markers_highlight = c('Rora','Il17a')
markers_highlight = c('Areg','Il4','Il13','Il6')

# 3B: tisTreg vs pTreg volcano plot
pdf(paste0(storeFigPath, 'colon_tisTreg_vs_ptreg_volcano_plot.pdf'),  width = 8, height = 6)
#jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', 
#                labs_range = c(-4,4), pt.size = 1, marker_thres = Inf, markers_highlight = markers_highlight, highlight.pt.size = 5, col_vec = c("black", "black", "black","red"))
#jj_volcano_plot(marker_peaks_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'symbol', 
#                labs_range = c(-6,6), pt.size = 1, marker_thres = c(Inf,10), markers_highlight = NULL)
jj_volcano_plot(marker_df, logfc_column = 'Log2FC', pval_column = 'FDR', symbol_column = 'name', add_thres_line = F, 
                labs_range = c(-4,4), pt.size = 1, marker_thres = c(2, 5), #markers_highlight = markers_highlight,
                highlight.pt.size = 2, #markers_highlight_col = structure(rep('red',length(markers_highlight)), names = markers_highlight),
                col_by_highlight = F, group_names = c('tisTregST2', 'pTreg'), arrow_pos = list(x=0.25, xend=1, y=15)) + labs(title = 'Colon') + theme(plot.title = element_text(hjust = 0.5))#+ guides(colour="none") 
dev.off()

table_3_4_vs_markers_list = list()
table_3_4_vs_markers_list[['tisTregST2_vs_pTreg']] = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-08-02-table/markers_colon_tisTregST2_vs_pTreg_markers.xlsx')[['marker_genes']]
table_3_4_vs_markers_list[['ILC2_C1_vs_ILC2_C2']] = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-01-snapshot_for_titan_fig2/ILC2_C1_vs_C2_markers.xlsx')[['marker_genes']]
table_3_4_vs_markers_list[['Th17_vs_Th17_Areg']] = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-07-fc_fc_plots/CD4_other_C3_vs_colon_Th17_C8_markers.xlsx')[['marker_genes']]
jj_save_excel(table_3_4_vs_markers_list, paste0(storeFigPath, 'table_3_4_one_vs_one_cell_type_comparison_markers.xlsx'))


# browser tracks ----------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
#contains the scATAC tisTreg regions as peakset
###
#full atlas
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected', pconfig$DONOR))

# browsertracks:
# th21tistreg: Ikzf2, Ctla4, Maf, tisTreg haut fett colon, th21 other main organ. vgl naive cd4, spleen b cell
# cd4other: areg, il17a, gleiche vgl
# ilc: alle 4 top marker, vgl same
# ilcc1 and c2: alle die rot sind (eigene). Vgl alle ilcs, keine tisTreg or naive cd4
# th17tistreg: Rorc, Il1a, Il17f, Cd3, Il22
# volcano plot: Il71a, Rorc (nur die beiden Vergleiche)
# volcano plot ilc c1 vs c2, tbx21,ifng, ncr1,areg, il4, il13, keine vgl

dr_df_tnkilc = jj_get_reduction_coords(loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected'))
c1_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C1']
c2_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C2']

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_categorical_by_group(dr_df, 'Tissue', 'cluster_annotation_fine', add_text = T, absolute_numbers = T, text_size = 2, flip_coordinates = T)
proj$tissue_cell_type = proj$cluster_annotation_fine
#tfh-like are almost exclusive in VAT
proj$tissue_cell_type[proj$cluster_annotation_fine=='tisTregST2'] = paste(proj$cluster_annotation_fine[proj$cluster_annotation_fine=='tisTregST2'], proj$Tissue[proj$cluster_annotation_fine=='tisTregST2'], sep = '__')
proj$tissue_cell_type[proj$cluster_annotation_fine=='B cell'] = paste(proj$cluster_annotation_fine[proj$cluster_annotation_fine=='B cell'], proj$Tissue[proj$cluster_annotation_fine=='B cell'], sep = '__')
proj$tissue_cell_type[proj$cluster_annotation_fine=='CD4_Tnaive'] = paste( proj$cluster_annotation_fine[proj$cluster_annotation_fine=='CD4_Tnaive'], proj$Tissue[proj$cluster_annotation_fine=='CD4_Tnaive'], sep = '__')
proj$tissue_cell_type[proj$cluster_annotation_fine=='Th17'] = paste( proj$cluster_annotation_fine[proj$cluster_annotation_fine=='Th17'], proj$Tissue[proj$cluster_annotation_fine=='Th17'], sep = '__')
proj$tissue_cell_type[proj$cluster_annotation_fine=='ILC2' & proj$cellNames %in% c1_bc] = 'ILC2__C1'
proj$tissue_cell_type[proj$cluster_annotation_fine=='ILC2' & proj$cellNames %in% c2_bc] = 'ILC2__C2'

my_cols = jj_get_jj_colours(proj$tissue_cell_type)
my_cols_names = names(my_cols)
my_cols_ref = jj_get_colours(proj$cluster_annotation_fine, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')
my_cols = unname(my_cols_ref[match(sapply(strsplit(names(my_cols), '__'), '[[', 1), names(my_cols_ref))])
names(my_cols) = my_cols_names

gannot = getGeneAnnotation(proj)$genes

# library(colorspace)
# for(i in seq_along(my_cols)){
#   if(grepl('Colon', names(my_cols)[i])){
#     my_cols[i] = darken(my_cols[i], 0.1)
#   }else if(grepl('Skin',  names(my_cols)[i])){
#     my_cols[i] = lighten(my_cols[i], 0.1)
#   }else if(grepl('VAT',  names(my_cols)[i])){
#     my_cols[i] = lighten(my_cols[i], 0.01)
#   }else if(grepl('Spleen',  names(my_cols)[i])){
#     my_cols[i] = darken(my_cols[i], 0.01)
#   }
# }
#anyDuplicated(unname(my_cols))

adjust_ranges = function(gr, new_ranges){
  for(i in new_ranges$symbol){
    if(i %in% gr$symbol){
      ranges(gr)[gr$symbol == i] = ranges(new_ranges)[new_ranges$symbol == i]
    }
  }
  gr
}
region_strings = c(Areg = 'chr5-91121604-91150517',
                   Il17a = 'chr1-20724552-20735413',
                   Ctla4 = 'chr1-60897703-60915868',
                   Maf = 'chr8-115596703-115848045')
new_ranges = genomic_region_tools::convert_granges(region_strings)
new_ranges$symbol = names(region_strings)


proj$custom = proj$tissue_cell_type
proj$custom[proj$Clusters_0.5 == 'C5'] = 'C5'
proj$custom[proj$Clusters_0.5 == 'C20'] = 'C20'
proj$custom[proj$Clusters_0.5 == 'C24'] = 'C24'
groups_select = c('tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon','C5','C20','C24', 'CD4_Tnaive__Spleen')
stopifnot(all(groups_select %in% proj$custom))
genes_plot = c('Batf', 'Entpd1', 'Maf', 'Irf4', 'Vps8', 'Col15a1', 'Pparg')
my_cols = c(my_cols, c(C20='red',C24='darkred',C5='green'))

gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     groupBy = 'custom',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
# 2H 3M S3E 4E 4J Mouse atlas browser tracks
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_additional.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()

### Tfh-like
table(proj$tissue_cell_type)
groups_select = c('Tfh-like','tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon','B cell__Spleen', 'CD4_Tnaive__Spleen')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Ikzf2','Ctla4','Maf')

gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_tfh_like_tisTregST2.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()

###CD4 other
groups_select = c('CD4_other','tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon','B cell__Spleen', 'CD4_Tnaive__Spleen')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Areg','Il17a')
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)


p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_CD4_other_tisTregST2.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()


###Th17
groups_select = c('Th17__Skin','Th17__VAT','Th17__Colon','tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon','B cell__Spleen', 'CD4_Tnaive__Spleen')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Rorc','Il1a','Il17f','Cd3e','Il22','Il17a')
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_Th17_tisTregST2.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()


### ILC2_C1
groups_select = c('ILC2__C1','tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon','B cell__Spleen', 'CD4_Tnaive__Spleen')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Klrg1','Il1rl1','Areg','Gata3')
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)


pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_ilc2_c1_tisTregST2.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()


### ILC2_C1 vs ILC2_C2
groups_select = c('ILC2__C1','ILC2__C2','ILC1','ILC3','tisTregST2__Skin','tisTregST2__VAT','tisTregST2__Colon')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Klrg1','Il13','Il4','Rora','Gata3','Il2ra','Il7r','Areg','Tbx21','Ifng','Il1rl1')  
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)

pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_ilc2_c1_ilc2_c2.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()

### ILC2_C1 vs ILC2_C2 differences
groups_select = c('ILC2__C1','ILC2__C2')
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Tbx21','Ifng','Ncr1','Il13','Il4','Areg')
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)

pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_ilc2_c1_ilc2_c2_differences.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()

### th17 vs cd4_other differences
groups_select = c('tisTregST2__Colon','Th17__Colon','CD4_other') #cd4 other is colon already
stopifnot(all(groups_select %in% proj$tissue_cell_type))
genes_plot = c('Rorc', 'Il17a','Areg')
gannot_use = gannot[0]
for(i in genes_plot){
  ruse = gannot[gannot$symbol %in% i]
  stopifnot(length(ruse)==1)
  ruse = ruse + 0.1*width(ruse)
  gannot_use =  c(gannot_use, ruse)
}
gannot_use = adjust_ranges(gannot_use, new_ranges)

p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_Th17_cd4_other_differences.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()


### ptreg vs tistreg 
region_strings = c(
  'Ccl20-Daw1' = 'chr1:83065353-83222701',
  Itm2b = 'chr14:73353416-73394058',
  Inpp4 = 'chr8:81123194-82168801',
  Emp1 = 'chr6:135347162-135394789',
  Krt28 = 'chr11:99358504-99382092',
  Foxp3 = 'chrX:7577676-7597243',
  Ceacam12 = 'chr7:18063818-18079392',
  Swap70 = 'chr7:110217010-110291678',
  Ctsh = 'chr9:90052109-90078877'
)

groups_select = c('pTreg', 'tisTregST2__Colon', 'tisTregST2__Skin','tisTregST2__VAT', 'CD4_Tnaive__Spleen') 
stopifnot(all(groups_select %in% proj$tissue_cell_type))

#proj$tissue_cell_type = factor( proj$tissue_cell_type, levels = )
library(genomic_region_tools)
gannot_use = convert_granges(region_strings)
p = plotBrowserTrack(proj,
                     region = gannot_use,
                     #geneSymbol = genes_plot[1],
                     groupBy = 'tissue_cell_type',
                     useGroups = groups_select, 
                     tileSize = 50,
                     pal = my_cols)
pdf(paste0(storeFigPath, 'mouse_atlas_browser_tracks_ptreg_tistreg.pdf'), width = 8, height = 6)
sapply(p, plot)
dev.off()

# library("ggplotify")
# gg = as.ggplot(p$Ikzf2)
# gg + scale_colour_manual(values = )

# Figure 3 T/NK/ILC subset comparisons ------------------------------------

proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
dr_df = jj_get_reduction_coords(proj, 'UMAP')

jj_plot_features(reduction = dr_df, meta_features = 'Clusters_1.2')
jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_fine', label = T,
                 custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'))

proj$tisTregST2  = proj$Clusters_1.2
proj$tisTregST2[proj$tisTregST2 %in%  c('C20','C22','C23')] = 'tisTregST2'
cell_groups = unique(proj$tisTregST2)

markers_list = list()

#get differences
markers_list[['C1']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'C1', 
  bgdGroups = c('tisTregST2'),
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_list[['tisTregST2']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'tisTregST2', 
  bgdGroups = c('C1'),
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#get commonalities
markers_list[['C1_vs_all']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'C1', 
  bgdGroups = cell_groups[!cell_groups %in% c('tisTregST2', 'C1')],
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markers_list[['tisTregST2_vs_all']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'PeakMatrix', 
  useGroups = 'tisTregST2', 
  bgdGroups = cell_groups[!cell_groups %in% c('tisTregST2', 'C1')],
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


marker_df = archr_get_markers_as_df(proj = proj, markers = markers_list[['C1']], cutOff = "FDR <= 0.01 & Log2FC >=0.58")
marker_df2 = archr_get_markers_as_df(proj = proj, markers = markers_list[['tisTregST2']], cutOff = "FDR <= 0.01 & Log2FC >=0.58")
marker_df2$Log2FC = -marker_df2$Log2FC
marker_df_all = rbind(marker_df, marker_df2)
table(marker_df_all$comparison)
#write_csv(marker_df_all, paste0(storeFigPath, 'ILC2_C1_vs_tisTreg_marker_peaks.csv'))
jj_volcano_plot(marker_df_all, logfc_column = 'Log2FC', pval_column = 'FDR',labs_range = c(-8,8), 
                                 symbol_column = 'symbol', marker_thres = c(Inf, 30), pt.size = 0.8)
#highlight = marker_df_all %>% dplyr::filter( (-log10(FDR) > 30 | abs(Log2FC) > 8) & distance < 2000) %>% pull(symbol)

go_files = c(
  hallmark="/omics/groups/OE0436/internal/msimon/scATAC/mh.all.v2022.1.Mm.symbols.gmt",
  curated_sets='/omics/groups/OE0436/internal/msimon/scATAC/m2.all.v2022.1.Mm.symbols.gmt',
  go_sets='/omics/groups/OE0436/internal/msimon/scATAC/m5.all.v2022.1.Mm.symbols.gmt'
  #celltype_signatures='/omics/groups/OE0436/internal/msimon/scATAC/m8.all.v2022.1.Mm.symbols.gmt'
)
go_list <- lapply(go_files, clusterProfiler::read.gmt)

marker_list = list(c1 = marker_df, tisTregST2 = marker_df2)
go_res_list = list()
for(i in seq_along(marker_list)){
  message(i, '/', length(marker_list))
  go_res_list[[i]] = lapply(go_list, function(x) {
    enricher(
      gene = marker_list[[i]]$symbol[marker_list[[i]]$distance < 2000],
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = x)
  }
  )
}
names(go_res_list) = names(marker_list)
dp_list = list()

for(j in names(go_list)){
  hallmark_list = lapply(go_res_list, function(x) x[[j]]@result)
  for(i in seq_along(hallmark_list)){
    halltemp = hallmark_list[[i]]
    halltemp$padj_BH = halltemp$qvalue
    halltemp$pathway = halltemp$ID
    halltemp$NES = sapply(strsplit(halltemp$GeneRatio, '/'), function(x) as.integer(x[1]) / as.integer(x[2]))
    halltemp = halltemp[halltemp$padj_BH < 0.05, ]
    halltemp = halltemp[halltemp$NES > 0.05, ]
    hallmark_list[[i]] = halltemp
  }
  hallmark_list = hallmark_list[sapply(hallmark_list, nrow) > 0]
  dp_list[[j]] = generate_gsea_dot_plot(hallmark_list, cell_type_names = names(hallmark_list))
  
}

tf_markers_list = list()
tf_markers_list[['C1']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'MotifMatrix', 
  useGroups = 'C1', 
  bgdGroups = c('tisTregST2'),
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
tf_markers_list[['tisTregST2']] <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'MotifMatrix', 
  useGroups = 'tisTregST2', 
  bgdGroups = c('C1'),
  groupBy = "tisTregST2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

tf_marker_df = archr_get_markers_as_df(proj = proj, 
                                       markers = tf_markers_list[['C1']], 
                                       cutOff = "FDR <= 0.01 & MeanDiff >=0", annotate_closest_gene = F)
tf_marker_df2 = archr_get_markers_as_df(proj = proj, markers = tf_markers_list[['tisTregST2']],
                                        cutOff = "FDR <= 0.01 & MeanDiff >=0", annotate_closest_gene = F)
tf_marker_df2$MeanDiff = -tf_marker_df2$MeanDiff
tf_marker_df_all = rbind(tf_marker_df, tf_marker_df2)
table(tf_marker_df_all$comparison)
#write_csv(tf_marker_df_all, paste0(storeFigPath, 'ILC2_C1_vs_tisTreg_marker_peaks.csv'))
jj_volcano_plot(tf_marker_df_all, logfc_column = 'MeanDiff', pval_column = 'FDR',labs_range = c(-0.1,0.4), 
                symbol_column = 'name', marker_thres = c(Inf, 50), pt.size = 0.8)
plotEmbedding(proj, 'UMAP', colorBy = 'MotifMatrix', name = c('z:bZIP.IRF.bZIP.IRF_29'))


## commonalities
marker_df3 = archr_get_markers_as_df(proj = proj, markers = markers_list[['C1_vs_all']], cutOff = "FDR <= 0.01 & Log2FC >=0.58")
marker_df4 = archr_get_markers_as_df(proj = proj, markers = markers_list[['tisTregST2_vs_all']], cutOff = "FDR <= 0.01 & Log2FC >=0.58")
jj_plot_upsetr(list(ILC2_C1=marker_df3$feature, tisTregST2=marker_df4$feature))
jj_plot_upsetr(list(ILC2_C1=marker_df3$symbol, tisTregST2=marker_df4$symbol))
marker_df_all = rbind(marker_df3, marker_df4)
table(marker_df_all$comparison)
#write_csv(marker_df_all, paste0(storeFigPath, 'ILC2_C1_and_tisTreg_vs_all_marker_peaks.csv'))
markers_common = marker_df3 %>% dplyr::inner_join(marker_df4, by = 'feature',suffix=c('.C1','.tisTregST2'))
jj_fc_fc_plot(markers_common, logfc_column1 = 'Log2FC.C1', logfc_column2 = 'Log2FC.tisTregST2',
              labs_range = c(0,8,0,8), marker_thres = Inf, symbol_column = 'symbol.C1', use_text = T, 
              markers_highlight = c('Areg','Rora','Ctsh','Hhip','Tox','Klrg1','Il2ra','Arid4b','Irak4','Ccr1','Irf6'))
jj_fc_fc_plot(markers_common, logfc_column1 = 'Log2FC.C1', logfc_column2 = 'Log2FC.tisTregST2',
              labs_range = c(0,8,0,8), marker_thres = 4, symbol_column = 'symbol.C1', use_text = T)


# mouse areg gfp ----------------------------------------------------------


pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_areg_gfp.yaml')

#print config variables
pconfig

set.seed(1)
#addArchRThreads(threads = 1) 
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, pconfig$ARCHR_PROJECT))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gmat = get_gene_mat(proj)

### S6D: Mouse Areg-GFP Tissue umap
gg = jj_plot_features(dr_df, features='Tissue', pt_size = 0.5,
                      custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                      return_gg_object = T)

pdf(paste0(storeFigPath, 'mouse_areg_gfp_tissue.pdf'),  width = 10, height = 8)
gg
dev.off()

# ###not included: Clusters umap
# dr_df$Clusters_1.5 = mixsort_factor(dr_df$Clusters_1.5)
# gg = jj_plot_features(reduction=dr_df, meta_features='Clusters_1.5', pt.size = 0.5, 
#                       return_gg_object = T, label_type = 'geom_label', label_col = 'white') 
# 
# pdf(paste0(storeFigPath, 'mouse_areg_gfp_clusters.pdf'),  width = 10, height = 8)
# gg
# dev.off()

### Annotation umap

# ref = loadArchRProject(path = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_filtered_no_doublets_corrected')
# cluster_annotation_fine_df = readr::read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
# ref$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
# gmat_ref = get_gene_mat(ref)
# gmat_sum = jj_summarize_sparse_mat(gmat_ref, ref$cluster_annotation_fine)
# preds = SingleR::SingleR(test = gmat, ref = gmat_sum, labels = colnames(gmat_sum))
# #write_rds(preds, paste0(bigFilesDir, 'mouse_areg_gfp_singler_pred_mouse_atlas_ref.RDS'))

preds = read_rds(paste0(bigFilesDir, 'mouse_areg_gfp_singler_pred_mouse_atlas_ref.RDS'))
dr_df$predictions = from_to(preds$labels, c(Th_Il21 = 'Tfh_like'))
dr_df$predictions = replace_if(dr_df$predictions, count_below = 20, 'undefined')

label_subset = unique(dr_df$predictions)[!unique(dr_df$predictions) %in% c('undefined','Plasma cell','B cell', 'DC', 'Monocyte','Macrophage')]
gg = jj_plot_features(reduction=dr_df, meta_features='predictions', pt.size = 0.5, label_type = 'geom_label_repel', label_subset = label_subset,
                      custom_colors = jj_get_colours(dr_df$predictions, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

# not included
pdf(paste0(storeFigPath, 'mouse_areg_gfp_mouse_atlas_singler.pdf'),  width = 10, height = 8)
gg
dev.off()

dr_df$annotation = from_to(vec= dr_df$Clusters_1.5, old_new_map_vec=c(
  'C1'= 'DC',
  'C2'= 'Macrophage/Monocyte',
  'C3'= 'Macrophage/Monocyte',
  'C4'= 'Macrophage/Monocyte',
  'C5'= 'ILC2',
  'C6'= 'ILC2',
  'C7'= 'CD4_other',
  'C8'= 'Tfh_like',
  'C9'= 'Tfh_like',
  'C10'= 'tisTregST2_prog',
  'C11'= 'tisTregST2_prog',
  'C12'= 'tisTregST2_prog',
  'C13'= 'Tnaive',
  'C14'= 'tisTregST2_prog',
  'C15'= 'tisTregST2_prog',
  'C16'= 'ILC1/NK/gd_T',
  'C17'= 'ILC2',
  'C18'= 'tisTregST2',
  'C19'= 'tisTregST2',
  'C20'= 'tisTregST2',
  'C21'= 'pTreg',
  'C22'= 'Th17',
  'C23'= 'B cell',
  'C24'= 'Plasma cell'
))
label_subset = unique(dr_df$predictions)[!unique(dr_df$predictions) %in% c('undefined','Plasma cell','B cell', 'DC', 'Monocyte','Macrophage')]
gg = jj_plot_features(reduction=dr_df, meta_features='annotation', pt.size = 0.5, label_type = 'geom_label_repel', #label_subset = label_subset,
                      custom_colors = jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 
# not included
pdf(paste0(storeFigPath, 'mouse_areg_gfp_annotation.pdf'),  width = 10, height = 8)
gg
dev.off()


#nFrags umap
gg = jj_plot_features(reduction=dr_df, meta_features='nFrags', pt.size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
pdf(paste0(storeFigPath, 'mouse_areg_gfp_nFrags.pdf'),  width = 10, height = 8)
gg
dev.off()

#sample umap 
dr_df = dr_df[sample(1:nrow(dr_df)), ] #shuffle to reduce overplotting problem
gg = jj_plot_features(reduction=dr_df, meta_features='Sample', pt.size = 0.5, 
                      return_gg_object = T)
pdf(paste0(storeFigPath, 'mouse_areg_gfp_sample.pdf'),  width = 10, height = 8)
gg
dev.off()

### singler prediction
dr_df$singler_prediction = from_to(vec= dr_df$singler_label, old_new_map_vec=c(
  'B cells'= 'B cell',
  'Basophils'= 'Basophil',
  'Eosinophils'= 'Eosinophil',
  'Macrophages'= 'Macrophage',
  'Mast cells'= 'Mast cell',
  'Monocytes'= 'Monocyte',
  'Neutrophils'= 'Neutrophil',
  'NK cells'= 'NK cell',
  'T cells'= 'T cell',
  'Tgd'= 'gd_T'
))
dr_df$singler_prediction = replace_if(dr_df$singler_prediction, count_below = 20, 'undefined')
gg = jj_plot_features(reduction=dr_df, meta_features='singler_prediction', pt.size = 0.5, return_gg_object = T, 
                      custom_colors = jj_get_colours(dr_df$singler_prediction, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) 

pdf(paste0(storeFigPath, 'mouse_areg_gfp_monaco_singler.pdf'),  width = 10, height = 8)
gg[[1]] + labs(colour='SingleR prediction') 
dev.off()



genes_plot = c( 'Il17a','Il17f', 'Areg', 'Rorc', 'Batf', 'Cd4')
gg = archr_plot_markers(proj, genes_plot)
cont_df = jj_get_contour_lines(dr_df[dr_df$predictions == 'CD4_other', ], 'predictions', .7)

pdf(paste0(storeFigPath, 'mouse_areg_gfp_magic_marker_umap.pdf'),  width = 8, height = 6)
gg
lapply(gg, function(x) x + geom_path(data = cont_df, aes(x=x, y=y, group=cont_group, linetype=predictions),
               size=1) + labs(linetype = 'Cell type'))
dev.off()

### cell types per tissue

dr_df$Tissue2 = factor(dr_df$Tissue, levels = c('Skin', 'VAT','Colon','Spleen'))
pdf(paste0(storeFigPath, 'mouse_areg_gfp_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'annotation', 'Tissue2', 
                             flip_coordinates = T, 
                             custom_colors =jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) + 
  labs(fill = 'Cell type', x='Tissue', y = 'Fraction')
dev.off()

pdf(paste0(storeFigPath, 'mouse_areg_gfp_annotation_barplot_absolute.pdf'),  width = 5, height = 6)
jj_plot_categorical_by_group(dr_df, 'annotation', 'Tissue2', 
                             flip_coordinates = F,absolute_numbers = T, add_text = T, text_size = 3, 
                             custom_colors =jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
  labs(fill = 'Cell type', x='Tissue', y = 'n cells')

dev.off()


#cell type tissue distribution
dr_df$cluster_annotation_levelx = factor(dr_df$annotation,
                                         levels = rev(gtools::mixedsort(unique(dr_df$annotation))))
pdf(paste0(storeFigPath, 'mouse_areg_gfp_tissue_by_cell_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'Tissue2', 'cluster_annotation_levelx', flip_coordinates = T,
                             custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  labs(x = 'Cell type', y = 'Fraction', fill = 'Tissue')
dev.off()


features_plot = list(
  Tcells= c('Cd3e','Cd3d','Cd4','Cd8a'),
  Bcells=c('Cd19','Ms4a1'),
  Macrophages=c('Mrc1','Csf1r'),
  Monocyte = c('Csf1r'),
  DCs=c('Siglech','Clec9a'),
  NK=c('Gzma','Prf1', "Klrb1"),
  ILCs=c('Il7r','Il2ra','Kit'),
  plasma = c('Sdc1','Bst2'),
  naive = c('Sell'),
  tisTreg = c('Ccr8','Ikzf2','Il1rl1','Klrg1','Nfil3','Foxp3'),
  pTreg = c('Rorc'),
  tfh_like = c('Il21'),
  cd4_other = c('Areg', 'Il17a')
) 
#gene_avail('Id2', rownames(gmat))

dr_df_use = dr_df
# dr_df_use$predictions = from_to(vec= dr_df_use$predictions, old_new_map_vec=c(
#   'CD8_Teff'= 'T cell',
#   'CD8_Tmem'= 'T cell',
#   'gd_T'= 'T cell',
#   'ILC3'= 'undefined',
#   'Mast cell/Eosinophil/Basophil'= 'undefined',
#   'NKT'= 'T cell'
# ))
#dr_df_use = dr_df_use[!dr_df_use$predictions == 'undefined', ]

pdf(paste0(storeFigPath, 'mouse_areg_gfp_marker_heatmap.pdf'),  width = 10, height = 5)
seurat_dotplot(gmat=gmat[, colnames(gmat) %in% dr_df_use$cellNames], metadf = dr_df_use,
               features = unique(unname(unlist(features_plot))),
               group_column = 'annotation')
dev.off()


proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
dr_df = as.data.frame(proj@cellColData)
### fragment size distribution
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% Fragments', color='Sample') + 
  scale_color_manual(values = jj_get_jj_colours(fsize_df$group))

pdf(paste0(storeFigPath, 'mouse_areg_gfp_fragment_size_distribution.pdf'), width = 8, height = 3)
p1
dev.off()

pdf(paste0(storeFigPath, 'mouse_areg_gfp_nFrags_TSS_scatterplot.pdf'), width = 10, height = 8)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6)
dev.off()


# colon subset ------------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_areg_gfp.yaml')

#print config variables
pconfig

set.seed(1)
#addArchRThreads(threads = 1) 
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

# proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, pconfig$ARCHR_PROJECT))
# proj = proj[proj$Tissue == 'Colon', ]
# proj = archr_dim_reduction(proj)
#saveArchRProject(proj, outputDirectory = 'ArchRProject_filtered_no_doublets_peaks_colon_subset')
proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_areg_gfp/ArchRProject_filtered_no_doublets_peaks_colon_subset')

dr_df = jj_get_reduction_coords(proj, 'UMAP')
gmat = get_gene_mat(proj)

### Clusters umap
dr_df$Clusters_0.4 = mixsort_factor(dr_df$Clusters_0.4)
gg = jj_plot_features(dr_df, features='Clusters_0.4', pt_size = 0.5, 
                      return_gg_object = T) 
#not included
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_clusters.pdf'),  width = 8, height = 6)
gg
dev.off()

### Annotation umap

# ref = loadArchRProject(path = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_filtered_no_doublets_corrected')
# cluster_annotation_fine_df = readr::read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
# ref$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
# gmat_ref = get_gene_mat(ref)
# gmat_sum = jj_summarize_sparse_mat(gmat_ref, ref$cluster_annotation_fine)
# preds = SingleR::SingleR(test = gmat, ref = gmat_sum, labels = colnames(gmat_sum))
# #write_rds(preds, paste0(bigFilesDir, 'mouse_areg_gfp_singler_pred_mouse_atlas_ref.RDS'))

preds = read_rds(paste0(bigFilesDir, 'mouse_areg_gfp_singler_pred_mouse_atlas_ref.RDS'))
preds = as.data.frame(preds)
preds = preds[match(dr_df$cellNames, rownames(preds)), ]
stopifnot(identical(rownames(preds), dr_df$cellNames))
dr_df$predictions = from_to(preds$labels, c(Th_Il21 = 'Tfh_like'))
dr_df$predictions = replace_if(dr_df$predictions, count_below = 20, 'undefined')

label_subset = unique(dr_df$predictions)[!unique(dr_df$predictions) %in% c('undefined','Plasma cell','B cell', 'DC', 'Monocyte','Macrophage')]
gg = jj_plot_features(reduction=dr_df, meta_features='predictions', pt.size = 0.5, label_type = 'geom_label_repel', label_subset = label_subset,
                      custom_colors = jj_get_colours(dr_df$predictions, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_mouse_atlas_singler.pdf'),  width = 6, height = 6)
gg
dev.off()

dr_df$annotation = from_to(vec= dr_df$Clusters_0.4, old_new_map_vec=c(
  'C1'= 'Macrophage/Monocyte',
  'C2'= 'DC',
  'C3'= 'Plasma cell',
  'C4'= 'B cell',
  'C5'= 'Th17_Areg',
  'C6'= 'ILC2',
  'C7'= 'ILC2',
  'C8'= 'pTreg',
  'C9'= 'tisTregST2',
  'C10'= 'Tfh_like',
  'C11'= 'ILC1',
  'C12'= 'Th17'
))
#label_subset = unique(dr_df$predictions)[!unique(dr_df$predictions) %in% c('undefined','Plasma cell','B cell', 'DC', 'Monocyte','Macrophage')]
gg = jj_plot_features(dr_df, features='annotation', pt_size = 0.5, label_type = 'geom_label_repel', #label_subset = label_subset,
                      custom_colors = jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

# S4E Mouse Areg-GFP colon ssubset cell type umap
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_annotation.pdf'),  width = 6, height = 6)
gg
dev.off()


#nFrags umap
gg = jj_plot_features(reduction=dr_df, meta_features='nFrags', pt.size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_nFrags.pdf'),  width = 6, height = 6)
gg
dev.off()

#sample umap 
dr_df = dr_df[sample(1:nrow(dr_df)), ] #shuffle to reduce overplotting problem
gg = jj_plot_features(reduction=dr_df, meta_features='Sample', pt.size = 0.5, 
                      return_gg_object = T)
pdf(paste0(storeFigPath, 'mouse_areg_gfp_sample.pdf'),  width = 6, height = 6)
gg
dev.off()

### singler prediction
dr_df$singler_prediction = from_to(vec= dr_df$singler_label, old_new_map_vec=c(
  'B cells'= 'B cell',
  'Basophils'= 'Basophil',
  'Eosinophils'= 'Eosinophil',
  'Macrophages'= 'Macrophage',
  'Mast cells'= 'Mast cell',
  'Monocytes'= 'Monocyte',
  'Neutrophils'= 'Neutrophil',
  'NK cells'= 'NK cell',
  'T cells'= 'T cell',
  'Tgd'= 'gd_T'
))
dr_df$singler_prediction = replace_if(dr_df$singler_prediction, count_below = 20, 'undefined')
gg = jj_plot_features(reduction=dr_df, meta_features='singler_prediction', pt.size = 0.5, return_gg_object = T, 
                      custom_colors = jj_get_colours(dr_df$singler_prediction, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) 

pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_monaco_singler.pdf'),  width = 6, height = 6)
gg[[1]] + labs(colour='SingleR prediction') 
dev.off()



genes_plot = c( 'Il17a','Il17f', 'Areg', 'Rorc', 'Batf', 'Cd4', 'Cd247', 'Cd3e', 'Cd3d')
gg = archr_plot_markers(proj, genes_plot)
# get_contour_lines = function(reduction_df, grouping_var, cont_thres){
#   reduction_df <- reduction_df[!is.na(reduction_df[, grouping_var]), 
#   ]
#   ls <- jj:::getContourLines(reduction_df, grouping_var, cont_thres)
#   ls = ls[sapply(ls, function(x) length(x) > 0)]
#   cont_df <- jj:::getContourCoordinates(ls)
#   names(cont_df)[5] <- grouping_var
#   return(cont_df)
# }
# cont_df = get_contour_lines(dr_df, 'predictions', cont_thres = 0.5)

# S4F Mouse areg gfp colon subset marker genes
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_magic_marker_umap.pdf'),  width = 6, height = 6)
gg
# lapply(gg, function(x) x + geom_path(data = cont_df[cont_df$predictions == 'CD4_other', ], aes(x=x, y=y, group=cont_group, linetype=predictions),
#                                      size=1) + labs(linetype = 'Cell type'))
dev.off()

### cell types per tissue

# not included:
dr_df$Tissue2 = factor(dr_df$Tissue, levels = c('Skin', 'VAT','Colon','Spleen'))
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'annotation', 'Tissue2',
                             flip_coordinates = T,
                             custom_colors =jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
  labs(fill = 'Cell type', x='Tissue', y = 'Fraction')
dev.off()

pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_annotation_barplot_absolute.pdf'),  width = 5, height = 6)
jj_plot_categorical_by_group(dr_df, 'annotation', 'Tissue2',
                             flip_coordinates = F,absolute_numbers = T, add_text = T, text_size = 3,
                             custom_colors =jj_get_colours(dr_df$annotation, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
  labs(fill = 'Cell type', x='Tissue', y = 'n cells')

dev.off()


features_plot = list(
  Tcells= c('Cd3e','Cd3d','Cd4','Cd8a'),
  Bcells=c('Cd19','Ms4a1'),
  Macrophages=c('Mrc1','Csf1r'),
  Monocyte = c('Csf1r'),
  DCs=c('Siglech','Clec9a'),
  ILCs=c('Il7r','Il2ra','Kit'),
  plasma = c('Sdc1','Bst2'),
  tisTreg = c('Ccr8','Ikzf2','Il1rl1','Klrg1','Nfil3','Foxp3'),
  pTreg = c('Rorc'),
  tfh_like = c('Il21'),
  cd4_other = c('Areg', 'Il17a')
) 
#gene_avail('Id2', rownames(gmat))

dr_df_use = dr_df
# dr_df_use$predictions = from_to(vec= dr_df_use$predictions, old_new_map_vec=c(
#   'CD8_Teff'= 'T cell',
#   'CD8_Tmem'= 'T cell',
#   'gd_T'= 'T cell',
#   'ILC3'= 'undefined',
#   'Mast cell/Eosinophil/Basophil'= 'undefined',
#   'NKT'= 'T cell'
# ))
#dr_df_use = dr_df_use[!dr_df_use$predictions == 'undefined', ]

#not included
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_marker_heatmap.pdf'),  width = 8, height = 5)
seurat_dotplot(gmat=gmat[, colnames(gmat) %in% dr_df_use$cellNames], metadf = dr_df_use,
               features = unique(unname(unlist(features_plot))),
               group_column = 'annotation')
dev.off()


proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
proj = proj[proj$Sample == 'MD_scATAC_62_Colon_Areg_GFP', ]
dr_df = as.data.frame(proj@cellColData)
### fragment size distribution
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% Fragments', color='Sample') + 
  scale_color_manual(values = jj_get_jj_colours(fsize_df$group))

#not included
pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_fragment_size_distribution.pdf'), width = 6, height = 3)
p1
dev.off()

pdf(paste0(storeFigPath, 'mouse_areg_gfp_colon_nFrags_TSS_scatterplot.pdf'), width = 6, height = 5)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6)
dev.off()

# gmat = get_gene_mat(proj)
# 
# p2 = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_filtered_no_doublets_corrected')
# refmat = get_gene_mat(p2)
# refmat = jj_summarize_sparse_mat(refmat, summarize_by_vec = p2$cluster_annotation_level3)
# refmat = refmat[, !colnames(refmat) %in% 'undefined']
# #write_rds(refmat, paste0(bigFilesDir, 'mouse_atlas_gene_ref_mat.RDS'))
# library(SingleR)
# pred_df = SingleR::SingleR(test = gmat, ref = refmat, labels = colnames(refmat))
# #write_rds(pred_df, paste0(bigFilesDir, 'mouse_areg_gfp_singler_pred.RDS'))
# plotScoreHeatmap(pred_df)
# proj$singler_pred = pred_df$labels
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# table(dr_df$singler_pred)
# dr_df$singler_pred = replace_if(dr_df$singler_pred, count_below = 75, 'Other')
# jj_plot_features(reduction = dr_df, meta_features = 'singler_pred')
# jj_plot_features(reduction = dr_df, meta_features = 'singler_pred', use_pointdensity = T, facet_by = 'singler_pred', background_cells = T)
# 
# #https://www.nature.com/articles/s41467-018-07492-4/figures/2
# tfh_sig = c('Bach2','Ccr7','Dapl1','Ets1','Ifnar1','Il17ra',
#             'Klf13','Klf2','Lef1','Myc','S1pr1','Sell')
# favail = getFeatures(proj)
# tfh_genes = 
# proj = ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = list(tfh=tfh_sig), name = 'Tfh_signature')
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_plot_features(reduction = dr_df, meta_features = 'Tfh_signature.tfh', cap_top = 'q99')
# 
# tfh_df = jj_load_excel('~/1-s2.0-S1074761319304509-mmc2.xlsx')[[1]]
# tfh_genes = tfh_df$GeneNames[tfh_df$`log2(Fold_change).normalized` > 2]
# 
# tfh_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/./analysis/2020-05-12-monocle_ccr8_and_treg/seurat_cd4_satpathy_Tfh_vs_Th17_Naive_CD4_T_Activated_CD4_T_Th1_scATAC_raw_T_cell_annot_simple_LR_0.1minpct_0.25logfc_n_peaks,patient_regressed.csv')
# tfh_df = tfh_df[tfh_df$comparison == 'Tfh_versus_Th17_Naive_CD4_T_Activated_CD4_T_Th1', ]
# tfh_df = tfh_df[tfh_df$avg_logFC > 0.5, ]
# chain = rtracklayer::import.chain('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/hg19ToMm10.over.chain')
# lift_res = rtracklayer::liftOver(human_peaks_gr, chain = chain)  
# lift_res = lift_res[sapply(lift_res, length) != 0]
# lift_res2 = lapply(lift_res, function(x) x[1])
# lift_gr = do.call(c, lift_res2)
# 
# stubbington_list = jj_load_excel('~/13062_2015_45_MOESM1_ESM.xls')
# th1 = readxl::read_xls('~/13062_2015_45_MOESM1_ESM.xls', sheet = 2, skip = 3)$`Gene name`
# th2 = readxl::read_xls('~/13062_2015_45_MOESM1_ESM.xls', sheet = 3, skip = 3)$`Gene name`
# th17 = readxl::read_xls('~/13062_2015_45_MOESM1_ESM.xls', sheet = 4, skip = 3)$`Gene name`
# th_list = list(th1=th1, th2=th2, th17=th17)
# th_list = lapply(th_list, function(x) x[x %in% favail])
# proj = ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = th_list, name = 'Th', nBin = 50, nBgd = 200)
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_plot_features(reduction = dr_df, meta_features = c('Th.th1','Th.th2','Th.th17'), cap_top = 'q99')
# 
# 
# #Table S2
# monaco = jj_load_excel('~/1-s2.0-S2211124719300592-mmc3.xlsx', startRow = 2)
# fc_df = monaco$`FoldChange TPM`
# fdr_df = monaco$`FDR TPM`
# cols_keep = c('Gene.name','T.CD4.Naive','Tregs','Tfh','Th1','Th17','Th2')#,'T.CD4.Naive+Tregs')
# fdr_df=fdr_df[, colnames(fdr_df) %in% cols_keep]
# fc_df=fc_df[, colnames(fc_df) %in% cols_keep]
# sig_list = list()
# for(i in cols_keep[-1]){
#   print(i)
#   test_vec = fdr_df[, i]
#   rest_df = fdr_df[, !colnames(fdr_df) %in% c(i,'Gene.name')]
#   thres = 0.2
#   keep1 = 0
#   # while(sum(keep1) < 100){
#     keep1 = test_vec < thres
#     # thres = thres + 0.1
#     
#   # }
#   keep2= apply(rest_df, 1, min) > thres
#   print(table(keep1))
#   print(table(keep1 & keep2))
#   #fdr_df$Gene.name[keep1 & keep2]
#   sig_list[[i]] = fdr_df$Gene.name[keep1]
# }
# for(i in cols_keep[-1]){
#   print(i)
#   test_vec = fc_df[, i]
#   rest_df = fc_df[, !colnames(fc_df) %in% c(i,'Gene.name')]
#   thres = 2
#   keep1 = test_vec > thres
#   keep2= apply(rest_df, 1, max) < thres
#   print(table(keep1))
#   print(table(keep1 & keep2))
#   sig_list[[i]] = fc_df$Gene.name[keep1 & keep2]
# }
# 
# for(i in cols_keep[-1]){
#   print(i)
#   test_vec = fc_df[, i]
#   rest_df = fc_df[, !colnames(fc_df) %in% c(i,'Gene.name')]
#   thres = 2
#   keep_df = data.frame(Gene.name = fc_df$Gene.name)
#   keep_df$keep= test_vec -  apply(rest_df, 1, max)
#   keep_df = keep_df[keep_df$keep > 0.5, ]
#   keep_df = keep_df[order(keep_df$keep, decreasing = T), ]
#   sig_list[[i]] = keep_df$Gene.name #head(keep_df$Gene.name, 200)
# }
# 
# 
# sapply(sig_list, length)
# #sig_list = lapply(sig_list, function(x) tools::toTitleCase(tolower(x)))
# #http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
# gene_conversion <- read_tsv('~/HOM_MouseHumanSequence.rpt.1')
# gene_conversion = gene_conversion %>% dplyr::select(`DB Class Key`, `Common Organism Name`, Symbol) %>% 
#   dplyr::rename(Key=`DB Class Key`, Species=`Common Organism Name`)
# gene_conversion$Species[gene_conversion$Species == 'mouse, laboratory'] = 'mouse'
# # head(gene_conversion)
# # # A tibble: 6 x 3
# # Key Species Symbol
# # <dbl> <chr>   <chr> 
# # 1 42532293 mouse   Coa4  
# # 2 42532293 human   COA4  
# # 3 42532294 mouse   Mvb12a
# # 4 42532294 human   MVB12A
# #write_rds(gene_conversion, paste0(bigFilesDir, 'HOM_MouseHumanSequence_tidy.RDS'))
# convert_genes = function(gene_vec, conf_df, map_to_mouse = TRUE, mode='many_to_many', verbose=F){
#   return_list = list()
#   is_mouse = conf_df$Species == 'mouse'
#   for(i in gene_vec){
#     if(i %in% conf_df$Symbol){
#       key_use = conf_df$Key[conf_df$Symbol == i]
#       #if(length(key_use)>1) warning('Multiple keys for one input gene found')
#       if(map_to_mouse){
#         new_symbols = conf_df$Symbol[conf_df$Key %in% key_use & is_mouse]
#       }else{
#         new_symbols = conf_df$Symbol[conf_df$Key %in% key_use & !is_mouse]
#       }
#       if(verbose) message(i,": ", paste(new_symbols, collapse = ', '))
#       return_list[[i]] = new_symbols
#     }else{
#       if(verbose) message(i, ' not found')
#       return_list[[i]] = NA
#     }
#   }
#   return_list
# }
# convert_genes(c('GZMB','GZMA', 'dfs'), conf_df = gene_conversion)
# sig_list_converted = lapply(sig_list, function(x) na.omit(unlist(convert_genes(x, conf_df = gene_conversion))))
# sapply(sig_list, length)
# sapply(sig_list_converted, length)
# 
# sig_list_converted = lapply(sig_list_converted, function(x) x[x %in% favail])
# sapply(sig_list_converted, length)
# proj = ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = sig_list_converted, name = 'Th')
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_arrange_ggplots(nplots=6, cols =3, jj_plot_features(reduction = dr_df, meta_features = paste0('Th.', names(sig_list)), cap_top = 'q99',cap_bottom = 'q01', return_gg_object = T))
# 
# # occurences = table(gene_conversion$`DB Class Key`)
# # occurences_keep = names(occurences)[occurences == 2]
# # gene_conversion_df = gene_conversion[gene_conversion$`DB Class Key` %in% occurences_keep, ]
# # anyDuplicated(gene_conversion_df$Symbol)
# # gene_conversion_df = gene_conversion_df[!duplicated(gene_conversion_df$Symbol), ]
# # res = gene_conversion_df %>% dplyr::select(`DB Class Key`,`Common Organism Name`, Symbol) %>%
# #   dplyr::group_by(`Common Organism Name`) %>% 
# #   mutate(row = row_number()) %>% 
# #   pivot_wider(
# #     id_cols = 'DB Class Key', 
# #     names_from = 'Common Organism Name', 
# #     values_from = 'Symbol'
# #   ) %>% 
# #   dplyr::select(-`DB Class Key`) %>% 
# #   dplyr::rename(mouse=`mouse, laboratory`)
# 
# #tibbit 2019
# #https://www.sciencedirect.com/science/article/pii/S1074761319302341?via%3Dihub
# #wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131937/suppl/GSE131937_Naive_CD4_vs_Th2.bed.gz
# # tibbit = read_tsv('~/GSE131937_Naive_CD4_vs_Th2.bed', col_names = c('seqnames', 'start','end','strand','diff'))
# # tibbit = tibbit[order(tibbit$diff, decreasing = T),]
# # tibbit = tibbit[tibbit$diff >5, ]
# # tibbit_gr = convert_granges(tibbit)
# # chain = rtracklayer::import.chain('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm9ToMm10.over.chain')
# # lift_res = rtracklayer::liftOver(tibbit_gr, chain = chain)  
# # table(sapply(lift_res, length))
# # lift_res = lift_res[sapply(lift_res, length) == 1]
# # lift_gr = do.call(c, lift_res)
# # proj = archr_add_peak_signatures(proj, signature_list = list(th2=lift_gr, th2_small=lift_gr[1:500]), signature_name = 'tibbit')
# # dr_df = jj_get_reduction_coords(proj, 'UMAP')
# # jj_plot_features(reduction = dr_df, meta_features = c('z_th2','z_th2_small'), cap_top = 'q99',cap_bottom = 'q01')
# 
# 
# # th17_sig_df = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mmc5.xlsx')[['Gene_correlations_with_voronoi']]
# # th17_sig_df$Voronoi.cells %>% table
# # th17_sig_df[th17_sig_df$Voronoi.cells == 'Th17/Th1-like effector' & th17_sig_df$`Up/down` == 'plus', 'Gene']
# 
# #A Validated Regulatory Network for Th17 Cell Specification
# #Table S1. Literature-Curated Validation List for Genes with Critical Influence for Th17 Development or Function, 
# th17_sig_df = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/1-s2.0-S0092867412011233-mmc2.xlsx')
# th17_sig = tools::toTitleCase(tolower(th17_sig_df$Sheet1$Gene[th17_sig_df$Sheet1$Effect == 'positive']))
# favail = getFeatures(proj)
# th17_sig = th17_sig[th17_sig %in% favail]
# proj = ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = list(th17 = th17_sig), name = 'th17')
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_plot_features(reduction = dr_df, meta_features = 'th17.th17', cap_top = 'q99')
# 
# 
# #Table S4
# tibbit = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/1-s2.0-S1074761319302341-mmc5.xlsx')[[3]]
# tibbit_sig = tibbit$gene[tibbit$p_val_adj < 0.001 & tibbit$avg_logFC > 1]
# tibbit_sig = tibbit_sig[tibbit_sig %in% favail]
# proj = ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = list(th2 = tibbit_sig), name = 'tibbit')
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# jj_plot_features(reduction = dr_df, meta_features = 'tibbit.th2', cap_top = 'q99')

# human donor 11 ----------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor11.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
dr_df = jj_get_reduction_coords(proj, 'UMAP')

### Tissue umap
gg = jj_plot_features(dr_df, features='Tissue', pt_size = 0.5, 
                      custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                      return_gg_object = T) 
# S5N Human donor 11 tissue umap
pdf(paste0(storeFigPath, 'human_d11_tissue.pdf'),  width = 10, height = 8)
gg
dev.off()


### Annotation umap

#undefined at -8,-4 are most likely fibroblasts
dr_df$cluster_annotation_level0 = from_to(vec= dr_df$cluster_annotation_level1, old_new_map_vec=c(
  'B'= 'B cell',
  'CD4_T'= 'T cell',
  'CD8_T'= 'T cell',
  'DC'= 'DC',
  'Monocyte'= 'Macrophage/Monocyte',
  'Neutrophil/Basophil'= 'Granulocyte',
  'NK'= 'NK cell',
  'other'= 'undefined',
  'T'= 'T cell'
))
proj$cluster_annotation_level0 = dr_df$cluster_annotation_level0
#saveArchRProject(proj)

gg = jj_plot_features(dr_df, features='cluster_annotation_level0', pt_size = 0.5, 
                      custom_colors = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

# S5O Human donor 11 cell type annotation umap
pdf(paste0(storeFigPath, 'human_d11_annotation.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()


### Clusters umap
#clusters used for subsetting
dr_df$Clusters_0.6 = mixsort_factor(dr_df$Clusters_0.6)
gg = jj_plot_features(dr_df, features='Clusters_0.6', pt_size = 0.5, 
                      return_gg_object = T, label_type = 'geom_text', fill_colors = 'black') 

# S5P Human donor 11 clusters umap
pdf(paste0(storeFigPath, 'human_d11_clusters.pdf'),  width = 10, height = 8)
gg
dev.off()


gg = jj_plot_features(dr_df, features='nFrags', pt_size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
# S5K nFrags umap
pdf(paste0(storeFigPath, 'human_d11_nFrags.pdf'),  width = 10, height = 8)
gg
dev.off()

# S5M sample umap 
gg = jj_plot_features(dr_df, features='Sample', pt_size = 0.5, 
                      return_gg_object = T)
pdf(paste0(storeFigPath, 'human_d11_sample.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()


### not included: singler prediction 
dr_df$singler_prediction = from_to(vec= dr_df$singler_label, old_new_map_vec=c(
  'B cells'= 'B cell',
  'Basophils'= 'Basophil',
  'CD4+ T cells'= 'T cell',
  'CD8+ T cells'= 'T cell',
  'Dendritic cells'= 'DC',
  'Monocytes'= 'Monocyte',
  'Neutrophils'= 'Neutrophil',
  'NK cells'= 'NK cell',
  'Progenitors'= 'Progenitor',
  'T cells'= 'T cell'
))
gg = jj_plot_features(reduction=dr_df, meta_features='singler_prediction', pt.size = 0.5, return_gg_object = T, 
                      custom_colors = jj_get_colours(dr_df$singler_prediction, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) 

pdf(paste0(storeFigPath, 'human_d11_singler.pdf'),  width = 10, height = 8)
gg[[1]] + labs(colour='SingleR prediction') 
dev.off()



### marker dotplot
proj_use = proj[!proj$cluster_annotation_level0 %in% c('undefined'), ]
gmat = get_gene_mat(proj_use)
dr_df = jj_get_reduction_coords(proj_use)

features_plot = list(
  Tcells= c('Cd3e','Cd3d'),#,'Cd4','Cd8a','Cd8b1'), #CD3z=CD247
  t_nk_ilc = c('Klrg1','Id2','Il7r'),
  NK=c('Gzma','Prf1', "Klrb1"),
  Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
  Monocyte = c('Csf1r','Cx3cr1','Cd14'),
  DCs=c('LAMP3','HLA-DRA','CXCR3','FLT3','CD74','SIGLEC1'),#siglech not availaable
  neutrophil = c('Itgam','S100a9'), #'Ly6g' not available
  Mastcells_Eo_Basophil=c('Itgam',"Il4")
) 
features_plot = lapply(features_plot, toupper)

# S5Q Human donor 11 marker gene dotplot
pdf(paste0(storeFigPath, 'human_d11_marker_heatmap.pdf'),  width = 10, height = 5)
seurat_dotplot(gmat=gmat, metadf = dr_df,
               features = unique(unname(unlist(features_plot))),
               group_column = 'cluster_annotation_level0')
dev.off()

### not included: cell types per tissue

dr_df = dr_df[dr_df$cluster_annotation_level0 != 'undefined', ]
dr_df$Tissue2 = factor(dr_df$Tissue, levels = c('Skin', 'Fat','Blood'))
pdf(paste0(storeFigPath, 'human_d11_annotation_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
                             flip_coordinates = T, 
                             custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) + 
  labs(fill = 'Cell type', x='Tissue', y = 'Fraction')
dev.off()

pdf(paste0(storeFigPath, 'human_d11_annotation_barplot_absolute.pdf'),  width = 5, height = 6)
jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
                             flip_coordinates = F,absolute_numbers = T, add_text = T, text_size = 3, 
                             custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
  labs(fill = 'Cell type', x='Tissue', y = 'n cells')

dev.off()


#cell type tissue distribution
dr_df$cluster_annotation_levelx = factor(dr_df$cluster_annotation_level0,
                                         levels = rev(gtools::mixedsort(unique(dr_df$cluster_annotation_level0))))
pdf(paste0(storeFigPath, 'human_d11_tissue_by_cell_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'Tissue2', 'cluster_annotation_levelx', flip_coordinates = T,
                             custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  labs(x = 'Cell type', y = 'Fraction', fill = 'Tissue')
dev.off()

# number of DA peaks in the major subsets ---------------------------------


markers_se <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = 'PeakMatrix',  
  groupBy = "cluster_annotation_level0",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#write_rds(markers_se, paste0(bigFilesDir, 'human_atlas_donor11_cluster_annotation_level0_all_markers.RDS'))
rm(markers_se)
markers_se = read_rds(paste0(bigFilesDir, 'human_atlas_donor11_cluster_annotation_level0_all_markers.RDS'))
markers_use_se = markers_se[, !colnames(markers_se) %in% 'undefined']
#transpose = T results in wrong clustering of cell types, nLabel = 0 is not possible...

# S5R Human donor 11 marker peak heatmap
pdf(paste0(storeFigPath, 'human_d11_marker_peak_heatmap.pdf'), width =6, height=6)
plotMarkerHeatmap(markers_use_se, 
                  cutOff = "FDR <= 0.01 & Log2FC >=0.5",  
                  transpose = F, clusterCols = T, binaryClusterRows = T, #clutering not working with transpose = T
                  nLabel = 1) #33653
dev.off()


marker_df = archr_get_markers_as_df(markers_use_se, proj, cutOff = "FDR <= 0.01 & Log2FC >=0.5")
write_csv(marker_df, paste0(storeFigPath, 'table_S5r_human_donor11_marker_peaks.csv'))
marker_list = split(marker_df, marker_df$comparison)
marker_list = lapply(marker_list, '[[', 'feature')
pdf(paste0(storeFigPath, 'human_d11_peak_upset_plot.pdf'), width =8, height=6)
jj_plot_upsetr(marker_list)
dev.off()

pgr = getPeakSet(proj)
peak_df = as.data.frame(unname(pgr))
peak_df$comparison = sprintf('Union peaks')

peak_df$feature = with(peak_df, paste(seqnames, start, end, sep = '-'))
marker_df = marker_df %>% dplyr::left_join(peak_df[, c('feature','nearestGene', 'distToGeneStart', 'nearestTSS', 'distToTSS', 'peakType')], by = 'feature')
subset_df = marker_df %>% dplyr::select(comparison, peakType)
subset_df = rbind(subset_df, peak_df[, c('comparison', 'peakType')])
subset_df = subset_df %>% dplyr::group_by(comparison) %>% add_count() %>% ungroup() %>% 
  dplyr::mutate(comparison = sprintf('%s (%i)', comparison, n))

nmarkers = as.data.frame(table(subset_df$comparison)) %>% arrange(Freq) %>% pull(Var1) %>% as.character
subset_df$comparison = factor(subset_df$comparison, levels = nmarkers)
#not included
pdf(paste0(storeFigPath, 'human_d11_peak_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(subset_df, 'peakType', 'comparison', flip_coordinates = T, absolute_numbers = F, add_text = F, text_size = 3) + 
  labs(y='n peaks', x = 'Cell type')
dev.off()

#enriched TF motifs in the marker peaks
motifsUp <- peakAnnoEnrichment(
  seMarker = markers_use_se,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(motifsUp, n = 5, transpose = T)
# not included
pdf(paste0(storeFigPath, 'human_d11_cell_type_markers_enr_motifs_heatmap.pdf'),  width = 7, height =5)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# human donor 11 qc -------------------------------------------------------

### distribution of fragments
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
dr_df = as.data.frame(proj@cellColData)
### fragment size distribution
#takes very long!
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% Fragments', color='Sample') + 
  scale_color_manual(values = jj_get_jj_colours(fsize_df$group))

### frags-tss scatter

pdf(paste0(storeFigPath, 'human_d11_fragment_size_distribution.pdf'), width = 8, height = 3)
p1
dev.off()

pdf(paste0(storeFigPath, 'human_d11_nFrags_TSS_scatterplot.pdf'), width = 10, height = 8)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6)
dev.off()

#doublet enrichment umap
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(reduction=dr_df, meta_features='DoubletEnrichment', pt.size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
pdf(paste0(storeFigPath, 'human_donor11_doublet_enrichment.pdf'), width = 10, height = 8)
gg
dev.off()

# combined donor 11+13 tss-nfrags scatterplot

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

### distribution of fragments
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
dr_df = as.data.frame(proj@cellColData)
dr_df$Sample = mixsort_factor(dr_df$Sample)

### frags-tss scatter (to be included)
pdf(paste0(storeFigPath, 'human_d11_d13_nFrags_TSS_scatterplot.pdf'), width = 10, height = 10)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6, 
                   facet_ncol = 3, text_size = 5)
dev.off()



# homer transcription factor z-scores -------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor11.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
binding_domains = binding_domains[!binding_domains == '?']
binding_domains[binding_domains=='T-box'] = 'T.box'
plotVarDev$data$bindingDomain = 'Other'
for(i in binding_domains){
  plotVarDev$data$bindingDomain[grepl(i, plotVarDev$data$name)] = i  
}
table(plotVarDev$data$bindingDomain)
plotVarDev$data$name = gsub('(.*)_[0-9]+','\\1', plotVarDev$data$name)

#not included
pdf(paste0(storeFigPath, 'human_d11_tf_deviation_rankplot.pdf'),  width = 8, height =5)
ggplot() + geom_point(data = plotVarDev$data, aes(x=rank, y=combinedVars, colour=bindingDomain), size = 1.5) + 
  geom_text_repel(data = plotVarDev$data[1:20, ], aes(x=rank, y=combinedVars, label=name)) + 
  theme_minimal() + labs(x = 'Rank') + scale_colour_manual(values = jj_get_colours(plotVarDev$data$bindingDomain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  labs(y = 'Combined variation', colour = 'Binding domain')
dev.off()

proj$peripheral_tissue = ifelse(proj$Tissue == 'Blood', 'blood',  'peripheral')
proj$cluster_annotation_level00 = proj$cluster_annotation_level0 #from_to(proj$cluster_annotation_level0, old_new_map_vec = c(Macrophage = "Macrophage/Monocyte", Monocyte = "Macrophage/Monocyte"))
proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level00, sep='__')

#sig_se = get_archr_mat(proj, 'MotifMatrix')
sig_se = getMatrixFromProject(proj, 'MotifMatrix')
z_score_mat = t(assays(sig_se)[['z']])
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))

comparisons_keep = dr_df %>% dplyr::group_by(cluster_annotation_level00) %>% dplyr::count(peripheral_tissue) %>% 
  tidyr::spread(peripheral_tissue, n) %>% 
  dplyr::filter(peripheral > 200 & blood > 200 & !cluster_annotation_level00 %in%  c('undefined')) %>% 
  dplyr::pull(cluster_annotation_level00)
# "Monocyte" "NK cell"  "T cell"  

da_list = list()
for(i in comparisons_keep){
  da_list[[i]] = getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = 'MotifMatrix',
    groupBy = 'peripheral_cell_type',
    useGroups = paste0('peripheral__', i),
    bgdGroups = paste0('blood__', i),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
}

res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0.05 & FDR <= 0.01'))

da_tf_list = lapply(res_se, function(x) x[[1]][, 'name'])
jj_plot_upsetr(da_tf_list)
motif_olap_df =jj_make_overlap_df(da_tf_list)
motifs_plot = unique(unlist(da_tf_list))
motifs_plot = na.omit(unique(c(unlist(motif_olap_df$overlaps[,!(colnames(motif_olap_df$overlaps) %in% c('a','b','c'))]))))
length(motifs_plot)

res_se2 = lapply(da_list, function(x) getMarkers(x, cutOff = 'FDR <= 1'))
diff_list = lapply(res_se2, function(x) as.data.frame(x[[1]][, c('name','MeanDiff')]))
res_df = reduce(diff_list, full_join, by = "name") %>% column_to_rownames('name')
colnames(res_df) = names(diff_list)

plot_mat = t(res_df[rownames(res_df) %in% motifs_plot, ])
colnames(plot_mat) = gsub('(.*)_[0-9]+','\\1',colnames(plot_mat))

heatmap_binding_domain = rep('Other', ncol(plot_mat))
for(i in binding_domains){
  heatmap_binding_domain[grepl(i, colnames(plot_mat))] = i
}
ha = columnAnnotation('Binding Domain'= heatmap_binding_domain, col = list('Binding Domain' = jj_get_colours(heatmap_binding_domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))

Heatmap(plot_mat) #define range based on automatic colorscale
col_fun = colorRamp2(c(-0.5, 0, 0.2, 0.6), c("darkblue", "white", "red", "darkred"))
h1 <- Heatmap(plot_mat, col = col_fun,
              name = 'MeanDiff\nperipheral - spleen',
              column_names_gp = gpar(fontsize = 8), 
              row_names_gp = gpar(fontsize = 8),
              top_annotation = ha)
# S5S Human donor 11 chromvar z-score heatmap peripheral vs blood
pdf(paste0(storeFigPath, 'human_d11_peripheral_tf_heatmap.pdf'), width =10, height=3)
h1
dev.off()

# TF footprints -----------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
motifPositions <- getPositions(proj, name = 'Motif')
grep('BATF.bZIP_20', names(motifPositions), value = T)
motifs <- c("Fra1.bZIP_90", "Atf3.bZIP_12", "BATF.bZIP_20", 'NFkB.p65.RHD_208')[3] #only show BATF for now
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
proj_use = proj[proj$peripheral_cell_type %in% c('blood__T cell', 'peripheral__T cell'), ] #, 'spleen__Monocyte', 'peripheral__Monocyte'), ]
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)
#plot(gg$BATF.bZIP_20)

proj_use = proj[proj$peripheral_cell_type %in% c('blood__Macrophage/Monocyte', 'peripheral__Macrophage/Monocyte'), ] 
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg2 = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)

plist = list(`T cell BATF.bZIP` = gg$BATF.bZIP_20, `Macrophage/Monocyte BATF.bZIP` = gg2$BATF.bZIP_20)
# S5T Human donor 11 TF footprints
pdf(paste0(storeFigPath, 'human_d11_tf_footprints.pdf'), width =6, height=5)
do.call(cowplot::plot_grid, c(list(ncol = 2, labels = c('T cell', 'Macrophage/Monocyte'), label_size=10),plist))
dev.off()


# human donor 11 tisTreg signature subsets --------------------------------

## not included:
# liftover_sig_gr = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_tisTreg_sig_liftover_hg38_gr.RDS')
# 
# config_list = get_config()
# dataset_use = 'human_normal_donor11'
# pconfig = yaml::read_yaml(config_list[[dataset_use]])
# setwd(pconfig$ARCHR_DIR)
# # #contains the scATAC tisTreg regions as peakset
# #proj = loadArchRProject(pconfig$ARCHR_PROJECT)
# #proj = addPeakSet(proj, peakSet = liftover_sig_gr, force=T)
# #saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_filtered_no_doublets_peaks_donor11_human_tisTreg_sig')
# #proj <- addPeakMatrix(proj,binarize = T)
# #saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_filtered_no_doublets_peaks_donor11_human_tisTreg_sig')
# 
# proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_filtered_no_doublets_peaks_donor11_human_tisTreg_sig')
# 
# proj = proj[!proj$cluster_annotation_level3 == 'other', ]
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# dr_df$cluster_annotation_fine = from_to(vec= dr_df$cluster_annotation_level3, old_new_map_vec=c(
#   'CD4 memory T cell'= 'CD4_Tmem',
#   'CD4 naive T cell'= 'CD4_Tnaive',
#   'CD8 memory T cell'= 'CD8_Tmem',
#   'CD8 naive T cell'= 'CD8_Tnaive',
#   'Memory B cell'= 'B cell',
#   'Monocyte'= 'Macrophage/Monocyte',
#   'Naive B cell'= 'B cell',
#   'Neutrophil/Basophil'= 'Granulocyte',
#   'pDC'= 'DC',
#   'Plasmablast'= 'Plasma cell',
#   'Tfh/Th1'= 'CD4_Tmem',
#   'Tfh/Th1/Th2'= 'CD4_Tmem',
#   'Th1/Th17'= 'CD4_Tmem',
#   'Th17 cell'= 'CD4_Tmem'
# ))
# 
# pmat = get_peak_mat(proj)
# stopifnot(identical(colnames(pmat), proj$cellNames))
# proj_peaks_gr = ArchR::getPeakSet(proj)
# names(proj_peaks_gr) = paste(seqnames(proj_peaks_gr), ranges(proj_peaks_gr), sep='-')
# 
# 
# mean_mat = jj_summarize_sparse_mat(pmat,  dr_df$cluster_annotation_fine)
# scaled_mat = t(scale(t(mean_mat)))
# htmat = t(scaled_mat)
# library(ComplexHeatmap)
# 
# ht = Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, gap = unit(3,'mm'),
#              show_column_names = F, column_split = 9, name = 'scaled\naccessibility')
# ht = draw(ht)
# pdf(paste0(storeFigPath, 'human_d11_tisTreg_liftover_signature_subset_heatmap.pdf'), width =8, height=5)
# ht
# dev.off()
# 
# corder = column_order(ht)
# subset_list = list(
#   c1_mono_dc =  colnames(htmat)[corder[[1]]],
#   c2_granulocyte =  colnames(htmat)[corder[[2]]],
#   c3_tisTreg =  colnames(htmat)[corder[[3]]],
#   c4_plasma =  colnames(htmat)[corder[[4]]],
#   c5_t = colnames(htmat)[corder[[5]]],
#   c689_maitcd8nk =  colnames(htmat)[c(corder[[6]], corder[[8]], corder[[9]])],
#   c7_b = colnames(htmat)[corder[[7]]]
# )
# 
# 
# #jj_plot_upsetr(subset_list)
# Heatmap(htmat[, colnames(htmat) %in% subset_list$c3_tisTreg], show_row_names=T, show_column_dend = F, show_row_dend = F,
#         show_column_names = F)
# 
# homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-15-human_d11_homer/'
# #save subset regions as bed
# #sapply(seq_along(subset_list), function(x) jj_save_bed(convert_granges(subset_list[[x]]), paste0(homer_res, names(subset_list)[x], '.bed')))
# #run homer using findMotfsGenome from Snakefile_runHomer.py
# 
# homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
# homer_list = list()
# for(i in homer_subsets){
#   homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
#   homer_df$motif = homer_df$`Motif Name`
#   homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
#   homer_df$rank = 1:nrow(homer_df)
#   if(nrow(homer_df)>20){
#     homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
#   }else{
#     homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
#   }
#   #1= size 2, 20+ = size 0.2
#   homer_df$name = i
#   homer_df = as.data.frame(homer_df)
#   homer_list[[i]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
# }
# sapply(homer_list, nrow)
# homer_motif_list = lapply(homer_list, '[[', 1)
# jj_plot_upsetr(homer_motif_list)
# 
# # #TFs that are significant only in one subset
# # diff_list = list()
# # for(i in seq_along(homer_motif_list)){
# #   diff_list[[i]] = setdiff(homer_motif_list[[i]], unlist(homer_motif_list[-i]))
# # }
# # names(diff_list) = names(homer_motif_list)
# 
# homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# # Combine the cell-type-specific data sets.
# gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
# rownames(gsea_res_comb) = NULL
# gsea_res_comb$rank = NULL
# gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]
# 
# # res = as.matrix(pivot_wider(gsea_res_comb, id_cols = motif, names_from =  name, values_from = enrichment, values_fill = NA))
# # rnames = res[, 1]
# # res = res[, -1]
# # res = apply(res, 2, as.numeric)
# # rownames(res) = rnames
# # Heatmap(res)
# 
# gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))
# 
# levels_use = gsea_res_comb %>% 
#   dplyr::group_by(motif) %>% 
#   dplyr::summarise(enr_sum = sum(enrichment)) %>% 
#   dplyr::arrange( enr_sum) %>%
#   pull(motif)
# gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
# gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
#   'c1_mono_dc'= '1',
#   'c2_granulocyte'= '2',
#   'c3_tisTreg'= '3',
#   'c4_plasma'= '4',
#   'c5_t'= '5',
#   'c7_b'= '7',
#   'c689_maitcd8nk'= '689'
# ))
# 
# pdf(paste0(storeFigPath, 'human_d11_tisTreg_signature_subset_homer.pdf'), width =5, height=8)
# ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "red") +
#   coord_fixed() + theme_minimal() + 
#   scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
#   labs(x = 'Signature subset', y = '', fill = 'Enrichment')
# dev.off()


# human donor 13 ----------------------------------------------------------

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor13.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
# proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
# cells_exclude = fractionsByCoordinates(reduction = dr_df,x.coord = c(7,14), y.coord = c(-3,-9))
# proj = proj[!proj$cellNames %in% cells_exclude$cellNames, ]
# proj = archr_dim_reduction(proj, cluster_res = NULL)
# proj = proj[!proj$cluster_annotation_level3 == 'other', ] #just 1 cell
# dr_df = jj_get_reduction_coords(proj, 'UMAP')
# proj = archr_clustering(proj)
# proj$cluster_annotation_level3 = from_to(vec= dr_df$Clusters_1.2, old_new_map_vec=c(
#   'C1'= 'Granulocyte',
#   'C2'= 'Granulocyte',
#   'C3'= 'Macrophage/Monocyte',
#   'C4'= 'Macrophage/Monocyte',
#   'C5'= 'Macrophage/Monocyte',
#   'C6'= 'DC',
#   'C7'= 'Macrophage/Monocyte',
#   'C8'= 'DC',
#   'C9'= 'DC',
#   'C10'= 'DC',
#   'C11'= 'CD4_Tmem',
#   'C12'= 'CD4_Tnaive',
#   'C13'= 'CD8_Tmem',
#   'C14'= 'NK cell',
#   'C15'= 'NK cell',
#   'C16'= 'CD8_Tmem',
#   'C17'= 'tisTreg',
#   'C18'= 'CD4_Tmem',
#   'C19'= 'CD4_Tmem',
#   'C20'= 'Plasma cell',
#   'C21'= 'B cell',
#   'C22'= 'pDC'
# ))
# #saveArchRProject(proj)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
dr_df = jj_get_reduction_coords(proj, 'UMAP')

### Tissue umap
gg = jj_plot_features(dr_df, features='Tissue', pt_size = 0.5,
                      custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                      return_gg_object = T) 
# 5A Human donor 13 Tissue umap
pdf(paste0(storeFigPath, 'human_d13_tissue.pdf'),  width = 10, height = 8)
gg
dev.off()

### Annotation umap
# dr_df$cluster_annotation_level0 = from_to(vec= dr_df$cluster_annotation_level3, old_new_map_vec=c(
#   'B cell'= 'B cell',
#   'CD4 T cell'= 'T cell',
#   'DC'= 'DC',
#   'effector CD8 T cell'= 'T cell',
#   'Monocytes'= 'Macrophage/Monocyte',
#   'Myeloid progenitor/Neutrophil (netosis?)'= 'Neutrophil',
#   'naive B cell'= 'B cell',
#   'naive CD4 T cell'= 'T cell',
#   'Neutrophil'= 'Neutrophil',
#   'NK cell'= 'NK cell',
#   'other'= 'undefined', #progenitors, hematopoietic stem cells...
#   'Treg'= 'T cell'
# ))
proj$cluster_annotation_level0 = from_to(vec= proj$cluster_annotation_level3, old_new_map_vec=c(
  'CD4_Tmem'= 'T cell',
  'CD4_Tnaive'= 'T cell',
  'CD8_Tmem'= 'T cell',
  'pDC'= 'DC',
  'Plasma cell'= 'B cell',
  'tisTreg'= 'T cell'
))

dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(dr_df, features='cluster_annotation_level0', pt_size = 0.5, 
                      custom_colors = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$'),
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

pdf(paste0(storeFigPath, 'human_d13_annotation.pdf'),  width = 10, height = 8)
gg
dev.off()

cols_use = jj_get_colours(dr_df$cluster_annotation_level3, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')
cols_use['pDC'] = 'pink'
gg = jj_plot_features(dr_df, features='cluster_annotation_level3', pt_size = 0.5, 
                      custom_colors = cols_use,
                      return_gg_object = T)[[1]] + labs(colour='Cell type') 

# 5B Human donor 13 Tissue umap
pdf(paste0(storeFigPath, 'human_d13_annotation_fine.pdf'),  width = 10, height = 8)
gg
dev.off()


### Clusters umap
#clusters used for subsetting
# dr_df$Clusters_0.5 = mixsort_factor(dr_df$Clusters_0.5)
# gg = jj_plot_features(reduction=dr_df, meta_features='Clusters_0.5', pt.size = 0.5, 
#                       return_gg_object = T, label_type = 'geom_label', label_col = 'white') 

#clusters used for annotation
dr_df$Clusters_1.2 = mixsort_factor(dr_df$Clusters_1.2)
gg = jj_plot_features(dr_df, features='Clusters_1.2', pt_size = 0.5, 
                       return_gg_object = T, label_type = 'geom_text', fill_colors = 'black') 

# S5C Human donor 13 clusters umap
pdf(paste0(storeFigPath, 'human_d13_clusters.pdf'),  width = 10, height = 8)
gg
dev.off()


#nFrags umap
gg = jj_plot_features(dr_df, features='nFrags', pt_size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
# S5A Human donor 13 nFrags umap
pdf(paste0(storeFigPath, 'human_d13_nFrags.pdf'),  width = 10, height = 8)
gg
dev.off()



#sample umap 
gg = jj_plot_features(dr_df, features='Sample', pt_size = 0.5,
                      return_gg_object = T)
# S5B Human donor 13 Sample umap
pdf(paste0(storeFigPath, 'human_d13_sample.pdf'),  width = 10, height = 8)
#add_umap_arrows(gg[[1]], theme_use = theme_minimal2)
gg
dev.off()

### singler prediction
dr_df$singler_prediction = from_to(vec= dr_df$singler_monaco, old_new_map_vec=c(
  'B cells'= 'B cell',
  'Basophils'= 'Basophil',
  'CD4+ T cells'= 'T cell',
  'CD8+ T cells'= 'T cell',
  'Dendritic cells'= 'DC',
  'Monocytes'= 'Monocyte',
  'Neutrophils'= 'Neutrophil',
  'NK cells'= 'NK cell',
  'Progenitors'= 'Progenitor',
  'T cells'= 'T cell'
))
gg = jj_plot_features(dr_df,features='singler_prediction', pt_size = 0.5, return_gg_object = T, 
                      custom_colors = jj_get_colours(dr_df$singler_prediction, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) 
#S5D Human donor 13 singler prediction
pdf(paste0(storeFigPath, 'human_d13_singler.pdf'),  width = 10, height = 8)
gg[[1]] + labs(colour='SingleR prediction') 
dev.off()

### marker dotplot
# features_plot = list(
#   Tcells= c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b'), #CD3z=CD247
#   Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
#   Monocyte = c('Csf1r','Cx3cr1','Cd14'),
#   DCs=c('LAMP3','HLA-DRA','CXCR3','FLT3','CD74','SIGLEC1'),#siglech not availaable
#   NK=c('Gzma','Prf1', "Klrb1",'Gzmb'),
#   neutrophil = c('Itgam','S100a9')
# ) 
#human d11:
# features_plot = list(
#   Tcells= c('Cd3e','Cd3d'),#,'Cd4','Cd8a','Cd8b1'), #CD3z=CD247
#   Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
#   Monocyte = c('Csf1r','Cx3cr1','Cd14'),
#   DCs=c('LAMP3','HLA-DRA','CXCR3','FLT3','CD74','SIGLEC1'),#siglech not availaable
#   NK=c('Gzma','Prf1', "Klrb1"),
#   neutrophil = c('Itgam','S100a9'), #'Ly6g' not available
#   Mastcells_Eo_Basophil=c('Itgam',"Il4"),
#   t_nk_ilc = c('Klrg1','Id2','Il7r')
# )
features_plot = list(
  Tcells= c('Cd3e','Cd3d'),
  Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
  Monocyte = c('Csf1r','Cx3cr1','Cd14'),
  DCs=c('LAMP3','HLA-DRA','FLT3','CD74','SIGLEC1'),#siglech not availaable
  NK=c('Gzma','Prf1', "Klrb1"),
  neutrophil = c('Itgam','S100a9'),
  t_nk_ilc = c('Klrg1','Id2','Il7r')
)

#donor 11 genes:
# features_plot = list(
#   Tcells= c('Cd3e','Cd3d'),
#   t_nk_ilc = c('Klrg1','Id2','Il7r'),
#   NK=c('Gzma','Prf1', "Klrb1"),
#   Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
#   Monocyte = c('Csf1r','Cx3cr1','Cd14'),
#   DCs=c('LAMP3','HLA-DRA','CXCR3','FLT3','CD74','SIGLEC1'),#siglech not availaable
#   neutrophil = c('Itgam','S100a9'), #'Ly6g' not available
#   Mastcells_Eo_Basophil=c('Itgam',"Il4")
# ) 

features_plot = lapply(features_plot, toupper)

dr_df = jj_get_reduction_coords(proj)
gmat = get_gene_mat(proj)

# 5D: Human donor 13 marker gene dotplot
pdf(paste0(storeFigPath, 'human_d13_marker_heatmap.pdf'),  width = 10, height = 5)
seurat_dotplot(gmat=gmat, metadf = dr_df,
               features = unique(unname(unlist(features_plot))),
               group_column = 'cluster_annotation_level0')
dev.off()

# ### not included: cell types per tissue
# 
# dr_df$Tissue2 = factor(dr_df$Tissue, levels = c('Skin', 'Fat','Blood'))
# pdf(paste0(storeFigPath, 'human_d13_annotation_barplot.pdf'),  width = 6, height = 3.5)
# jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
#                              flip_coordinates = T, 
#                              custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) + 
#   labs(fill = 'Cell type', x='Tissue', y = 'Fraction')
# dev.off()
# 
# pdf(paste0(storeFigPath, 'human_d13_annotation_barplot_absolute.pdf'),  width = 5, height = 6)
# jj_plot_categorical_by_group(dr_df, 'cluster_annotation_level0', 'Tissue2', 
#                              flip_coordinates = F,absolute_numbers = T, add_text = T, text_size = 3, 
#                              custom_colors =jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$') ) +
#   labs(fill = 'Cell type', x='Tissue', y = 'n cells')
# 
# dev.off()


#cell type tissue distribution
dr_df$cluster_annotation_levelx = factor(dr_df$cluster_annotation_level0,
                                         levels = rev(gtools::mixedsort(unique(dr_df$cluster_annotation_level0))))
# S5J Human donor 13 cell type distribution in tissues
pdf(paste0(storeFigPath, 'human_d13_tissue_by_cell_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(dr_df, 'Tissue2', 'cluster_annotation_levelx', flip_coordinates = T,
                             custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  labs(x = 'Cell type', y = 'Fraction', fill = 'Tissue')
dev.off()

# 5C: marker umaps
genes_plot = c( 'FOXP3', 'GATA3', 'TBX21', 'RORC')
#genes_plot = c('CTLA4', 'IKZF2', 'IL2RA', 'GITR', 'BCL6', 'MAF')
gg = archr_plot_markers(proj, genes_plot)
pdf(paste0(storeFigPath, 'human_donor13_magic_marker_umap.pdf'),  width = 7.5, height = 6)
gg
dev.off()

# number of DA peaks in the major subsets ---------------------------------


markers_se <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = 'PeakMatrix',  
  groupBy = "cluster_annotation_level0",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#write_rds(markers_se, paste0(bigFilesDir, 'human_atlas_donor13_cluster_annotation_level0_all_markers.RDS'))
rm(markers_se)
markers_se = read_rds(paste0(bigFilesDir, 'human_atlas_donor13_cluster_annotation_level0_all_markers.RDS'))
markers_use_se = markers_se[, !colnames(markers_se) %in% 'undefined']
#transpose = T results in wrong clustering of cell types, nLabel = 0 is not possible...

#5E Human donor 13 marker peak heatmap
pdf(paste0(storeFigPath, 'human_d13_marker_peak_heatmap.pdf'), width =6, height=6)
plotMarkerHeatmap(markers_use_se, 
                  cutOff = "FDR <= 0.01 & Log2FC >=0.5",  
                  transpose = F, clusterCols = T, binaryClusterRows = T, #clutering not working with transpose = T
                  nLabel = 1) #33653
dev.off()


marker_df = archr_get_markers_as_df(markers_use_se, proj, cutOff = "FDR <= 0.01 & Log2FC >=0.5")
#write_csv(marker_df, paste0(storeFigPath, 'table_5e_human_donor13_marker_peaks.csv'))
marker_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-08-02-table/table_5e_human_donor13_marker_peaks.csv')

marker_list = split(marker_df, marker_df$comparison)
marker_list = lapply(marker_list, '[[', 'feature')
# S5H Human donor 13 upset plot
pdf(paste0(storeFigPath, 'human_d13_peak_upset_plot.pdf'), width =8, height=6)
jj_plot_upsetr(marker_list)
dev.off()

pgr = getPeakSet(proj)
peak_df = as.data.frame(unname(pgr))
peak_df$comparison = sprintf('Union peaks')

peak_df$feature = with(peak_df, paste(seqnames, start, end, sep = '-'))
marker_df = marker_df %>% dplyr::left_join(peak_df[, c('feature','nearestGene', 'distToGeneStart', 'nearestTSS', 'distToTSS', 'peakType')], by = 'feature')
subset_df = marker_df %>% dplyr::select(comparison, peakType)
subset_df = rbind(subset_df, peak_df[, c('comparison', 'peakType')])
subset_df = subset_df %>% dplyr::group_by(comparison) %>% add_count() %>% ungroup() %>% 
  dplyr::mutate(comparison = sprintf('%s (%i)', comparison, n))

nmarkers = as.data.frame(table(subset_df$comparison)) %>% arrange(Freq) %>% pull(Var1) %>% as.character
subset_df$comparison = factor(subset_df$comparison, levels = nmarkers)

pdf(paste0(storeFigPath, 'human_d13_peak_type_barplot.pdf'),  width = 6, height = 3.5)
jj_plot_categorical_by_group(subset_df, 'peakType', 'comparison', flip_coordinates = T, absolute_numbers = F, add_text = F, text_size = 3) + 
  labs(y='n peaks', x = 'Cell type')
dev.off()

#enriched TF motifs in the marker peaks
motifsUp <- peakAnnoEnrichment(
  seMarker = markers_use_se,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(motifsUp, n = 5, transpose = T)
# S5G Enriched TFs in the marker peaks of human donor 13
pdf(paste0(storeFigPath, 'human_d13_cell_type_markers_enr_motifs_heatmap.pdf'),  width = 7, height =5)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# human tisTreg signature -------------------------------------------------

library(genomic_region_tools)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

# human_tisTreg_sig_df = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-04-30-tconv_signature_substraction/tisTreg_vs_Tconv_filtering.csv')
# human_tisTreg_sig_df = human_tisTreg_sig_df$tisTreg_vs_Treg_not_in_Tconv_vs
human_tisTreg_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-02-08-cluster_comparison_3_vs_7_human_cd4/diff_results_7_versus_3_seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.csv')
human_tisTreg_sig_df = human_tisTreg_sig_df[human_tisTreg_sig_df$comparison == '3_versus_7', ]
human_tisTreg_sig_gr = convert_granges(human_tisTreg_sig_df$feature)

chainfile = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/hg19ToHg38.over.chain'
liftover_sig_gr = liftover_granges(human_tisTreg_sig_gr, chainfile) #5 peaks could not be transferred
#write_rds(liftover_sig_gr, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_tisTreg_sig_liftover_hg38_gr.RDS')
liftover_sig_gr = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_tisTreg_sig_liftover_hg38_gr.RDS')

config_list = get_config()
dataset_use = 'human_normal_donor13'
pconfig = yaml::read_yaml(config_list[[dataset_use]])
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)

proj = archr_add_peak_signatures(proj, list(c3_vs_c7 = liftover_sig_gr), signature_name = 'lifover_tisTreg_sig')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df[is.na(dr_df)] = 0
gg = jj_plot_features(reduction = dr_df,
                      pt.size = 0.5,
                      meta_features = grep('z_', colnames(dr_df), value=T),
                      cont_or_disc = 'c',
                      colorScale = 'viridis',
                      cap_top = 'q99', cap_bottom = 'q01',
                      custom_theme = theme_minimal(), return_gg_object = T)

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, 
                                 nFrags_vec = proj$nFrags, signature_gr = liftover_sig_gr,
                                 verbose = T, k = 100, count_thres = 2e5)
#write_rds(olap_df, paste0(storeFigPath, 'human_d13_c3_vs7_scATAC_tisTreg_sig_olap_df.RDS'))
#write_rds(olap_df, paste0(storeFigPath, 'human_d11_nr1_scATAC_tisTreg_sig_olap_df.RDS'))
#write_rds(olap_df, paste0(storeFigPath, 'human_d11_c3_vs7_scATAC_tisTreg_sig_olap_df.RDS'))
#write_rds(olap_df, paste0(storeFigPath, 'human_d13_c3_vs_c7_scATAC_tisTreg_sig_olap_df.RDS'))
#olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_d11_nr1_scATAC_tisTreg_sig_olap_df.RDS')
#olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_d11_nr1_scATAC_tisTreg_sig_olap_df.RDS')
#olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis//2023-02-27-d13/human_d13_c3_vs7_scATAC_tisTreg_sig_olap_df.RDS')
olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-01-peripheral_markers/human_d13_c3_vs_c7_scATAC_tisTreg_sig_olap_df.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$pct_overlap = olap_df$signature_pct_overlap

# S5F Human donor 13 tisTreg signature liftover percentage overlap
pdf(paste0(storeFigPath, 'human_d13_tisTreg_signature_overlap_umap.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = 'pct_overlap', return_gg_object = T)[[1]] + labs(colour = '% overlap') 
dev.off()
# celltypes_keep = dr_df$cluster_annotation_fine %>% table %>% .[. > 200] %>% names
# dr_df = dr_df[dr_df$cluster_annotation_fine %in% celltypes_keep, ]


dr_df$cluster_annotation_fine = from_to(vec= dr_df$cluster_annotation_level3, old_new_map_vec=c(
  'pDC' = 'DC'
))

# 5H Human donor 13 tisTreg signature percentage overlap boxplot
pdf(paste0(storeFigPath, 'human_d13_tisTreg_signature_overlap_boxplot.pdf'),  width = 6, height = 5)
jj_plot_numeric_by_group(dr_df[!dr_df$cluster_annotation_fine == 'undefined', ], 'pct_overlap', group_column = 'cluster_annotation_fine', 
                         order = T, flip_coordinates = T, type = 'boxplot',
                         custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  theme(legend.position = 'none') + labs(y='% overlap', x='Cell type')
dev.off()


# human tisTreg signature deconvolution -----------------------------------

liftover_sig_gr = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-02-09-human_analysis/human_tisTreg_sig_liftover_hg38_gr.RDS')

config_list = get_config()
dataset_use = 'human_normal_donor13'
pconfig = yaml::read_yaml(config_list[[dataset_use]])
setwd('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_filtered_no_doublets_peaks_donor13_human_tisTreg_sig')
# #contains the scATAC tisTreg regions as peakset
# proj = addPeakSet(proj, peakSet = liftover_sig_gr, force=T)
# proj <- addPeakMatrix(proj,binarize = T)
# saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_filtered_no_doublets_peaks_donor13_human_tisTreg_sig')

proj = loadArchRProject()
pmat = get_peak_mat(proj)
stopifnot(identical(colnames(pmat), proj$cellNames))
proj_peaks_gr = ArchR::getPeakSet(proj)
names(proj_peaks_gr) = paste(seqnames(proj_peaks_gr), ranges(proj_peaks_gr), sep='-')
#stopifnot(identical(convert_granges(liftover_sig_gr), rownames(pmat)))

dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$cluster_annotation_fine = from_to(vec= dr_df$cluster_annotation_level3, old_new_map_vec=c(
  'pDC' = 'DC'
))

mean_mat = jj_summarize_sparse_mat(pmat,  dr_df$cluster_annotation_fine)
scaled_mat = t(scale(t(mean_mat)))
htmat = t(scaled_mat)
library(ComplexHeatmap)

ht = Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, gap = unit(3,'mm'),
             show_column_names = F, column_split = 6, name = 'scaled\naccessibility')
ht = draw(ht)
# 5I Human donor 13 human tisTreg signature subset heatmap
pdf(paste0(storeFigPath, 'human_d13_tisTreg_liftover_signature_subset_heatmap.pdf'), width =8, height=5)
ht
dev.off()

corder = column_order(ht)
subset_list = list(
  c1_mono_dc =  colnames(htmat)[corder[[1]]],
  c2_tmem_nk =  colnames(htmat)[corder[[2]]],
  c3_plasma =  colnames(htmat)[corder[[3]]],
  c45_tisTreg =  colnames(htmat)[c(corder[[4]], corder[[5]])],
  c6_granulocyte = colnames(htmat)[corder[[6]]]
)


#jj_plot_upsetr(subset_list)
Heatmap(htmat[, colnames(htmat) %in% subset_list$c45_tisTreg], show_row_names=T, show_column_dend = F, show_row_dend = F,
       show_column_names = F)

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-15-human_d13_homer/'
#save subset regions as bed
#sapply(seq_along(subset_list), function(x) jj_save_bed(convert_granges(subset_list[[x]]), paste0(homer_res, names(subset_list)[x], '.bed')))
#run homer using findMotfsGenome from Snakefile_runHomer.py

subset_files = list.files('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-15-human_d13_homer/', pattern = '.bed', full.names = T)
subset_list = lapply(subset_files, jj_load_bed)
names(subset_list) = gsub('\\.bed', '', basename(subset_files))
liftover_sig_df = annotate_surrounding_genes(liftover_sig_gr,gene_annotation = getGeneAnnotation(proj)$genes)
for(i in seq_along(subset_list)){
  liftover_sig_df[, names(subset_list)[i]] = liftover_sig_df$query_region %in% subset_list[[i]]$name
}
liftover_sig_df = liftover_sig_df %>% dplyr::rename(feature = query_region) %>% 
  dplyr::select(feature, closest_gene_id:closest_distance, c1_mono_dc:c6_granulocyte)
write_csv(liftover_sig_df, paste0(storeFigPath, 'table_5i_human_tisTreg_signature_subsets.csv'))


homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_list[[i]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
sapply(homer_list, nrow)
homer_motif_list = lapply(homer_list, '[[', 1)
jj_plot_upsetr(homer_motif_list)

# #TFs that are significant only in one subset
# diff_list = list()
# for(i in seq_along(homer_motif_list)){
#   diff_list[[i]] = setdiff(homer_motif_list[[i]], unlist(homer_motif_list[-i]))
# }
# names(diff_list) = names(homer_motif_list)

homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]

# res = as.matrix(pivot_wider(gsea_res_comb, id_cols = motif, names_from =  name, values_from = enrichment, values_fill = NA))
# rnames = res[, 1]
# res = res[, -1]
# res = apply(res, 2, as.numeric)
# rownames(res) = rnames
# Heatmap(res)

gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'c1_mono_dc'= '1',
  'c2_tmem_nk'= '2',
  'c3_plasma'= '3',
  'c45_tisTreg'= '4/5',
  'c6_granulocyte'= '6'
))

# 5J Human donor 13 tisTreg signature subsets homer known motif heatmap
pdf(paste0(storeFigPath, 'human_d13_tisTreg_signature_subset_homer.pdf'), width =5, height=8)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'Signature subset', y = '', fill = 'Enrichment')
dev.off()


# # not included: d11 d13 tisTreg signature subset overlaps -------------------------------
# 
# d13 = list.files('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-15-human_d13_homer/', pattern = '.bed', full.names = T)
# subset_list_d13 = lapply(d13, function(x) pull(jj_load_bed(x), 'name'))
# names(subset_list_d13) = gsub('.bed', '_d13', basename(d13))
# 
# d11 = list.files('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-15-human_d11_homer/', pattern = '.bed', full.names = T)
# subset_list_d11 = lapply(d11, function(x) pull(jj_load_bed(x), 'name'))
# names(subset_list_d11) = gsub('.bed', '_d11', basename(d11))
# 
# 
# upset_df = UpSetR::fromList(c(subset_list_d13, subset_list_d11))
# upset_df_d13= upset_df[, grep('_d13', colnames(upset_df))]
# upset_df_d11= upset_df[, grep('_d11', colnames(upset_df))]
# olap_df = jj_initialize_df(ncol = ncol(upset_df_d13), nrow = ncol(upset_df_d11), init =  0, col.names = colnames(upset_df_d13), row.names = colnames(upset_df_d11))
# for(i in seq_along(upset_df_d13)){
#   for(j in seq_along(upset_df_d11)){
#     olap_df[j, i] = sum(upset_df_d13[, i] == 1 & upset_df_d11[, j] == 1)
#   }
# }
# 
# olap_gg_df = olap_df %>% rownames_to_column('Group_D11') %>%  reshape2::melt(id = 'Group_D11', value.name = 'Peak_overlap') %>% dplyr::rename(Group_D13 = variable)
# 
# min_max_val = range(olap_gg_df$Peak_overlap)
# ggplot(olap_gg_df, aes(x=Group_D11, y = Group_D13, fill=Peak_overlap)) + 
#   geom_bin2d(stat='identity') + geom_text(aes(label=Peak_overlap, colour=Peak_overlap), size=3) + 
#   viridis::scale_fill_viridis() + theme_minimal() + #scale_color_stepsn(colours= c('white', rep('black',5)), guide = 'none')+
#   #scale_x_continuous(breaks = seq(min(res_df$fragment_cutoff), max(res_df$fragment_cutoff), by = 1000)) + 
#   #scale_y_continuous(breaks = seq(min(res_df$tss_cutoff), max(res_df$tss_cutoff), by = 1))  +
#   # scale_color_stepsn(
#   #   colours = c("white", "white", "black", "black", "black"), guide = 'none'
#   # )
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   binned_scale(aesthetics = "color",
#                scale_name = "stepsn", 
#                palette = function(x) c("white", "black"),
#                breaks = (min_max_val[2]-min_max_val[1])/2 + min_max_val[1],
#                limits = c(min_max_val[1], min_max_val[2]),
#                show.limits = TRUE, 
#                guide = "none"
#   )


# human donor 13 qc -------------------------------------------------------

### distribution of fragments
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
dr_df = as.data.frame(proj@cellColData)
### fragment size distribution
#takes very long!
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% Fragments', color='Sample') + 
  scale_color_manual(values = jj_get_jj_colours(fsize_df$group))

### frags-tss scatter

pdf(paste0(storeFigPath, 'human_fragment_size_distribution.pdf'), width = 8, height = 4)
p1
dev.off()

pdf(paste0(storeFigPath, 'human_nFrags_TSS_scatterplot.pdf'), width = 18, height = 12)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 3000, tss_cutoff = 6)
dev.off()

#doublet enrichment umap
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(reduction=dr_df, meta_features='DoubletEnrichment', pt.size = 0.5, 
                      return_gg_object = T, cap_top = 'q95')
pdf(paste0(storeFigPath, 'human_donor13_doublet_enrichment.pdf'), width = 10, height = 8)
gg
dev.off()


# homer transcription factor z-scores -------------------------------------


#sig_se = get_archr_mat(proj, 'MotifMatrix')
sig_se = getMatrixFromProject(proj, 'MotifMatrix')
z_score_mat = t(assays(sig_se)[['z']])
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))

proj$peripheral_tissue = ifelse(proj$Tissue == 'Blood', 'blood',  'peripheral')
proj$cluster_annotation_level00 = dr_df$cluster_annotation_level0 #from_to(dr_df$cluster_annotation_level0, old_new_map_vec = c(Monocyte = "Monocyte/DC", DC = "Monocyte/DC"))
proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level00, sep='__')
dr_df = jj_get_reduction_coords(proj, redname='UMAP')

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
binding_domains = binding_domains[!binding_domains == '?']
binding_domains[binding_domains=='T-box'] = 'T.box'
plotVarDev$data$bindingDomain = 'Other'
for(i in binding_domains){
  plotVarDev$data$bindingDomain[grepl(i, plotVarDev$data$name)] = i  
}
table(plotVarDev$data$bindingDomain)
plotVarDev$data$name = gsub('(.*)_[0-9]+','\\1', plotVarDev$data$name)

pdf(paste0(storeFigPath, 'human_d13_tf_deviation_rankplot.pdf'),  width = 8, height =5)
ggplot() + geom_point(data = plotVarDev$data, aes(x=rank, y=combinedVars, colour=bindingDomain), size = 1.5) + 
  geom_text_repel(data = plotVarDev$data[1:20, ], aes(x=rank, y=combinedVars, label=name)) + 
  theme_minimal() + labs(x = 'Rank') + scale_colour_manual(values = jj_get_colours(plotVarDev$data$bindingDomain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  labs(y = 'Combined variation', colour = 'Binding domain')
dev.off()

#sig_se = get_archr_mat(proj, 'MotifMatrix')
sig_se = getMatrixFromProject(proj, 'MotifMatrix')
z_score_mat = t(assays(sig_se)[['z']])
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))

comparisons_keep = dr_df %>% dplyr::group_by(cluster_annotation_level00) %>% dplyr::count(peripheral_tissue) %>% 
  tidyr::spread(peripheral_tissue, n) %>% 
  dplyr::filter(peripheral > 200 & blood > 200 & !cluster_annotation_level00 %in%  c('undefined')) %>% 
  dplyr::pull(cluster_annotation_level00)

da_list = list()
for(i in comparisons_keep){
  da_list[[i]] = getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = 'MotifMatrix',
    groupBy = 'peripheral_cell_type',
    useGroups = paste0('peripheral__', i),
    bgdGroups = paste0('blood__', i),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
}

res_se = lapply(da_list, function(x) getMarkers(x, cutOff = 'MeanDiff >= 0.05 & FDR <= 0.01'))

da_tf_list = lapply(res_se, function(x) x[[1]][, 'name'])
jj_plot_upsetr(da_tf_list)
motif_olap_df =jj_make_overlap_df(da_tf_list)
motifs_plot = unique(unlist(da_tf_list))
motifs_plot = na.omit(unique(c(unlist(motif_olap_df$overlaps[,!(colnames(motif_olap_df$overlaps) %in% c('a','b','c'))]))))
length(motifs_plot)

res_se2 = lapply(da_list, function(x) getMarkers(x, cutOff = 'FDR <= 1'))
diff_list = lapply(res_se2, function(x) as.data.frame(x[[1]][, c('name','MeanDiff')]))
res_df = purrr::reduce(diff_list, full_join, by = "name") %>% column_to_rownames('name')
colnames(res_df) = names(diff_list)

plot_mat = t(res_df[rownames(res_df) %in% motifs_plot, ])
colnames(plot_mat) = gsub('(.*)_[0-9]+','\\1',colnames(plot_mat))

heatmap_binding_domain = rep('Other', ncol(plot_mat))
for(i in binding_domains){
  heatmap_binding_domain[grepl(i, colnames(plot_mat))] = i
}
ha = columnAnnotation('Binding Domain'= heatmap_binding_domain, col = list('Binding Domain' = jj_get_colours(heatmap_binding_domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))

Heatmap(plot_mat) #define range based on automatic colorscale
col_fun = colorRamp2(c(-0.05, 0, 0.2, 0.4), c("darkblue", "white", "red", "darkred"))
h1 <- Heatmap(plot_mat, col = col_fun,
              name = 'MeanDiff\nperipheral - spleen',
              column_names_gp = gpar(fontsize = 8), 
              row_names_gp = gpar(fontsize = 8),
              top_annotation = ha)
# 5F Human donor 13 cromvar z-score heatmap peripheral vs blood
pdf(paste0(storeFigPath, 'human_d13_peripheral_tf_heatmap.pdf'), width =10, height=4)
h1
dev.off()

# TF footprints -----------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
motifPositions <- getPositions(proj, name = 'Motif')
grep('BATF.bZIP_20', names(motifPositions), value = T)
motifs <- c("Fra1.bZIP_90", "Atf3.bZIP_12", "BATF.bZIP_20", 'NFkB.p65.RHD_208')[3] #only show BATF for now
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
proj_use = proj[proj$peripheral_cell_type %in% c('blood__T cell', 'peripheral__T cell'), ] #, 'spleen__Monocyte', 'peripheral__Monocyte'), ]
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)
#plot(gg$BATF.bZIP_20)

proj_use = proj[proj$peripheral_cell_type %in% c('blood__Macrophage/Monocyte', 'peripheral__Macrophage/Monocyte'), ] 
proj_use <- addGroupCoverages(ArchRProj = proj_use, groupBy = "peripheral_cell_type")
seFoot <- getFootprints(
  ArchRProj = proj_use, 
  positions = motifPositions[markerMotifs], 
  groupBy = "peripheral_cell_type"
)
gg2 = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_use, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, plot = F
)

plist = list(`T cell BATF.bZIP` = gg$BATF.bZIP_20, `Macrophage/Monocyte BATF.bZIP` = gg2$BATF.bZIP_20)
# 5G Human donor 13 TF foottprints
pdf(paste0(storeFigPath, 'human_d13_tf_footprints.pdf'), width =6, height=5)
do.call(cowplot::plot_grid, c(list(ncol = 2, labels = c('T cell', 'Monocyte'), label_size=10),plist))
dev.off()



# te quantification validation --------------------------------------------
#use pbmc scATAC dataset also used in fig S11h of the paper
proj = loadArchRProject('/omics/groups/OE0436/internal/Regensburg/scATAC_public/10x_10k_pbmc_human/ArchRProject_filtered/')
comb_norm_df = prepare_comb_norm_df(proj = proj, 
                                    config = '/omics/groups/OE0436/internal/msimon/scATAC/transposable_elements/human_10x_pbmc_comb_df_binary2_20230425.RDS', 
                                    normalize = T)
features_plot = toupper(c('LTR4','LTR26','LTR8B','LTR71A','LTR10C','HERV17_int','MER89','MER48','MER57B2','HERV9N_int',
                          'LTR16D1','LTR33C','LTR10E','MLT1D_int','MER72B','MER54B','HERVP71A_int','LTR24','LTR27B',
                          'LTR18A','LTR26D','LTR1C','LTR10B','LTR47B','MER57E2','HERV1_LTRc','LTR12','LTR12_',
                          'HERVS71_int','LTR21B','LTR10A','LTR43B','LTR2B','LTR77','MER49','LTR2','LTR13_','LTR43',
                          'LOR1a','LTR9','LTR1A1','MER57C2','MER65B','LTR53B','LTR06','HERV1_LTRe'))

names(features_plot) = c(rep('Monocytes', 26), rep('T_NK_cells', 6), rep('B_cells', 14))
colnames(comb_norm_df) = gsub('-', '_', toupper(colnames(comb_norm_df)))
features_plot_avail = features_plot[features_plot %in% colnames(comb_norm_df)]
features_plot[!features_plot %in% colnames(comb_norm_df)]

proj$cell_type = proj$singler_label
cells_keep = proj$cell_type %in% c('Monocytes', 'B cells', 'NK cells', 'CD4+ T cells', 'CD8+ T cells')
identical(proj$cellNames, rownames(comb_norm_df))
#ggplot_dotplot(t(comb_norm_df), genes_plot = features_plot_avail, group_vec = proj$singler_label, scale_data = T, cluster_rows = F)
te_heatmap = jj_plot_heatmap(obj = t(comb_norm_df[cells_keep, ]), features_use = features_plot_avail, 
                group_vec = factor(proj$cell_type[cells_keep], levels = c('CD8+ T cells', 'NK cells', 'CD4+ T cells', 'B cells', 'Monocytes')),
                cluster_columns = F, cluster_rows = F)
te_heatmap@top_annotation@anno_list$Feature@color_mapping@name = 'PMID 33674594 annotation'
te_heatmap@matrix_legend_param$title = 'Scaled mean accessibility'

#S7B TE quantification validation
pdf(paste0(storeFigPath, 'human_pbmc_scatac_figS11_validation.pdf'), width = 10, height=4)
draw(te_heatmap,  annotation_legend_side = 'bottom', heatmap_legend_side ='bottom', merge_legend = TRUE)
dev.off()

# te analysis mouse atlas -------------------------------------------------

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
proj = proj[proj$cluster_annotation_level0 != 'undefined', ]
#proj = proj[proj$Tissue == 'Colon', ]
proj$cluster_annotation_level0 = replace_if(proj$cluster_annotation_level0, 100, 'undefined')
proj = proj[!proj$cluster_annotation_level0 %in% 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')
comb_norm_df = prepare_comb_norm_df(pconfig, proj)
stopifnot(identical(rownames(dr_df), rownames(comb_norm_df)))
# proj$peripheral_tissue = ifelse(proj$Tissue == 'Blood', 'blood',  'peripheral')
# proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level1, sep='__')


# ###not included: Heatmap per tissue
# sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$Tissue)
# sum_mat_scaled = scale(t(sum_mat))
# 
# pdf(paste0(storeFigPath, 'mouse_te_tissue_heatmap.pdf'), width =6, height=3)
# Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
#         cluster_rows = T)
# dev.off()

### Heatmap per cell type
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$cluster_annotation_level0)
sum_mat_scaled = scale(t(sum_mat))
# S7F Mouse atlas global TE heatmap per cell type
pdf(paste0(storeFigPath, 'mouse_te_celltype_heatmap.pdf'), width =8, height=4)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = T)
dev.off()

proj_colon = proj[proj$Tissue == 'Colon',]
proj_colon$cluster_annotation_level0 = replace_if(proj_colon$cluster_annotation_level0, 100, 'undefined')
proj_colon = proj_colon[!proj_colon$cluster_annotation_level0 %in% 'undefined', ]
dr_df = jj_get_reduction_coords(proj_colon, 'UMAP')
comb_norm_df = prepare_comb_norm_df(pconfig, proj_colon)
stopifnot(identical(rownames(dr_df), rownames(comb_norm_df)))
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj_colon$cluster_annotation_level0)
sum_mat_scaled = scale(t(sum_mat))
# 6E Mouse atlas colon subset te heatmap per cell type
pdf(paste0(storeFigPath, 'mouse_colon_te_celltype_heatmap.pdf'), width =8, height=4)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = T)
dev.off()

### TE UMAP
#all te insertions
comb_df = prepare_comb_norm_df(pconfig, proj, normalize = F)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
sobj = CreateSeuratObject(assay = 'RNA', counts = Matrix::t(comb_df), meta.data = dr_df)
sobj = NormalizeData(sobj)
#sobj = FindVariableFeatures(sobj, selection.method = 'mean.var.plot', dispersion.cutoff = c(0.5, Inf))
#VariableFeaturePlot(sobj)
sobj = ScaleData(sobj, features = rownames(sobj), vars.to.regress = 'nFrags') #cor nFrags with nCount_RNA = 0.97
sobj = RunPCA(sobj, features = rownames(sobj))
#sobj <- JackStraw(sobj, num.replicate = 100)
#sobj <- ScoreJackStraw(sobj, dims = 1:20)
sobj = RunUMAP(sobj, dims = 1:20)
#write_rds(sobj, paste0(bigFilesDir, 'mouse_atlas_te_seurat.RDS')) #23.02.2023
dr_df = jj_get_reduction_coords(sobj, 'umap')
#pca_res = prcomp(comb_norm_df[, !colVars(comb_norm_df)==0], center = T, scale. = T)
#umap_res = uwot::umap(comb_norm_df, pca=5, pca_center = T, scale = T)
#dr_df = cbind(umap_res, as.data.frame(proj@cellColData))

col_map = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

# not included: umap based on all TE insertions
# pdf(paste0(storeFigPath, 'mouse_atlas_te_umap.pdf'), width = 10, height=8)
# jj_plot_features(dr_df, features = 'cluster_annotation_level0', custom_colors = col_map)
# jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
#                  do_shuffle = T)
# jj_plot_features(reduction=dr_df, meta_features = 'nFrags', cap_top = 'q95')
# #jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
# dev.off()


#only distal tes
proj = proj[proj$Tissue == 'Colon', ] #only for figure S7E: subset to colon
comb_df = prepare_comb_norm_df(config = '/omics/groups/OE0436/internal/msimon/scATAC/transposable_elements/corrected/mouse_imm_cell_atlas_te_df_binary_distal_only.RDS',
                               proj, normalize = F)
dr_df = as.data.frame(proj@cellColData)
stopifnot(identical(rownames(dr_df), rownames(comb_df)))
sobj = CreateSeuratObject(assay = 'RNA', counts = Matrix::t(comb_df), meta.data = dr_df)
sobj = NormalizeData(sobj)
sobj = ScaleData(sobj, features = rownames(sobj), vars.to.regress = 'nFrags') #cor nFrags with nCount_RNA = 0.97
sobj = RunPCA(sobj, features = rownames(sobj))
#sobj <- JackStraw(sobj, num.replicate = 100)
#sobj <- ScoreJackStraw(sobj, dims = 1:20)
sobj = RunUMAP(sobj, dims = 1:20)
#write_rds(sobj, paste0(bigFilesDir, 'mouse_atlas_te_seurat.RDS')) #23.02.2023
dr_df = jj_get_reduction_coords(sobj, 'umap')
col_map = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

#6C 6D S7E Mouse atlas (colon subset) TE umap
pdf(paste0(storeFigPath, 'mouse_atlas_distal_te_umap.pdf'), width = 10, height=8)
jj_plot_features(dr_df, features = 'cluster_annotation_level0', custom_colors = col_map)
jj_plot_features(dr_df, features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
jj_plot_features(dr_df, features = 'nFrags', cap_top = 'q95')
#jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
dev.off()


### T/NK/ILC subset
config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$cluster_annotation_fine = replace_if(dr_df$cluster_annotation_fine, count_below = 100, 'undefined')

proj$cluster_annotation_fine = dr_df$cluster_annotation_fine
jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation_fine')
proj = proj[proj$cluster_annotation_fine != 'undefined', ]
comb_norm_df = prepare_comb_norm_df(pconfig, proj)

### Heatmap per tissue
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), summarize_by_vec = proj$Tissue)
sum_mat_scaled = scale(t(sum_mat))
colnames(sum_mat_scaled)[!complete.cases(t(sum_mat_scaled))]
sum_mat_scaled = t(t(sum_mat_scaled)[complete.cases(t(sum_mat_scaled)), ])

#not included
pdf(paste0(storeFigPath, 'mouse_tnkilc_te_tissue_heatmap.pdf'), width =6, height=3)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = T)
dev.off()

### Heatmap per cell type
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), summarize_by_vec = proj$cluster_annotation_fine)
sum_mat_scaled = scale(t(sum_mat))
colnames(sum_mat_scaled)[!complete.cases(t(sum_mat_scaled))]
sum_mat_scaled = t(t(sum_mat_scaled)[complete.cases(t(sum_mat_scaled)), ])

#6F tnkilc subset TE heatmap per celltype
pdf(paste0(storeFigPath, 'mouse_tnkilc_te_celltype_heatmap.pdf'), width =8, height=4)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = T)
dev.off()

### tisTreg signature significant TE heatmap T/NK/ILC subset
config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#dr_df$cluster_annotation_fine = replace_if(dr_df$cluster_annotation_fine, count_below = 100, 'undefined')

proj = proj[proj$cluster_annotation_fine != 'undefined', ]
dr_df = jj_get_reduction_coords(proj, 'UMAP')

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 

#only tisTreg-sig overlapping sites
comb_df = read_rds('/omics/groups/OE0436/internal/msimon/scATAC/transposable_elements/mouse_cd4_tisTreg_sig_olap_tes_04072023.RDS')
dr_df = as.data.frame(proj@cellColData)
comb_df = comb_df[match(rownames(dr_df), rownames(comb_df)), ]
stopifnot(identical(rownames(dr_df), rownames(comb_df)))
comb_norm_df = apply(comb_df, 2, function(x) x / proj$nFrags * 10000)

#all TE integration sites (S8A)
comb_norm_df = prepare_comb_norm_df(pconfig, proj)


comb_norm_df = comb_norm_df[, colnames(comb_norm_df) %in% olaps_sig_df$name]
stopifnot(identical(rownames(dr_df), rownames(comb_norm_df)))
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$cluster_annotation_fine)
sum_mat_scaled = scale(t(sum_mat))
family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
colcsv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'
col_list_use = list(Class = jj_get_colours(family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)],colour_csv = colcsv, comment_char = '$'),
                    Family = jj_get_colours(family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)], colour_csv = colcsv, comment_char = '$'))
column_ha = HeatmapAnnotation(Class = family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)], 
                              Family = family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)],
                              col = col_list_use)

mean_scaled_scores = apply(sum_mat_scaled, 1, mean)
row_ha = rowAnnotation(Mean =anno_barplot(mean_scaled_scores))

# S8A tnkilc TE integrations in the tisTreg signature heatmap per cell type
pdf(paste0(storeFigPath, 'mouse_tnkilc_tisTreg_sig_fisher_signif_te_celltype_heatmap.pdf'), width =14, height=5)
pdf(paste0(storeFigPath, 'mouse_tnkilc_tisTreg_sig_fisher_signif_te_only_tisTreg_sig_olap_celltype_heatmap.pdf'), width =14, height=5)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = T,
        cluster_rows = T, top_annotation = column_ha, right_annotation = row_ha)

dev.off()

#ggplot(data = data.frame(celltype = reorder(names(mean_scaled_scores), mean_scaled_scores), mean_scaled_mean = mean_scaled_scores)) + 
#  geom_bar(aes(x = celltype,  y = mean_scaled_mean), stat = 'identity') + coord_flip()
comb_norm_scaled = scale(comb_norm_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$tisTreg_te_scaled_mean = apply(comb_norm_scaled,1,mean)

pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_tisTreg_sig_fisher_signif_te_umap.pdf'),  width = 8, height = 6)
pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_tisTreg_sig_fisher_signif_te_only_tisTreg_sig_olap_celltype_umap.pdf'),width = 8, height = 6)
jj_plot_features(dr_df, 'tisTreg_te_scaled_mean', cap_top = 'q99', cap_bottom = 'q01')
dev.off()

# 6H S8B only tisTreg-sig overlapping sites tisTreg TE scaled mean boxplot
# pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_tisTreg_sig_fisher_signif_te_boxplot.pdf'),  width = 4, height = 4)
# pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_tisTreg_sig_fisher_signif_te_only_tisTreg_sig_olap_boxplot.pdf'),  width = 4, height = 4)
jj_plot_numeric_by_group(dr_df, 'tisTreg_te_scaled_mean',
                         group_column = 'cluster_annotation_fine', order = T, flip_coordinates = T, type = 'boxplot',
                         custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) +
  labs(y = 'mean(scaled TE accessibility)',x='Celltype') + ylim(-1, 1) + theme(legend.position = 'none')
# dev.off()


### TE UMAP
comb_df = prepare_comb_norm_df(pconfig, proj, normalize = F)
sobj = CreateSeuratObject(assay = 'RNA', counts = Matrix::t(comb_df), meta.data = as.data.frame(proj@cellColData))
sobj = NormalizeData(sobj)
#sobj = FindVariableFeatures(sobj, selection.method = 'mean.var.plot', dispersion.cutoff = c(0.5, Inf))
#VariableFeaturePlot(sobj)
sobj = ScaleData(sobj, features = rownames(sobj), vars.to.regress = 'nFrags') #cor nFrags with nCount_RNA = 0.97
sobj = RunPCA(sobj, features = rownames(sobj))
#sobj <- JackStraw(sobj, num.replicate = 100)
#sobj <- ScoreJackStraw(sobj, dims = 1:20)
sobj = RunUMAP(sobj, dims = 1:20)
#write_rds(sobj, paste0(bigFilesDir, 'mouse_atlas_te_seurat.RDS')) #23.02.2023
dr_df = jj_get_reduction_coords(sobj, 'umap')
#pca_res = prcomp(comb_norm_df[, !colVars(comb_norm_df)==0], center = T, scale. = T)
#umap_res = uwot::umap(comb_norm_df, pca=5, pca_center = T, scale = T)
#dr_df = cbind(umap_res, as.data.frame(proj@cellColData))

col_map = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')
pdf(paste0(storeFigPath, 'mouse_atlas_tnkilc_te_umap.pdf'), width = 10, height=8)
jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation_fine', #custom_colors = col_map, 
                 do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                 do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'nFrags', cap_top = 'q95')
#jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
dev.off()



# homer known motifs for TE overlap with tisTreg signature subsets --------

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-03-peripheral_markers/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = T)
homer_res2 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/'
homer_subsets2 = list.dirs(homer_res2, recursive = F, full.names = T)
homer_subsets = c(homer_subsets, homer_subsets2)


#map homer names to binding domains
homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains_keep = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
mapping_df= data.frame(homer_name = homer_headers_names, domain = binding_domains)
mapping_df$domain[!mapping_df$domain %in% binding_domains_keep] = 'Other'
mapping_df$domain[mapping_df$domain == '?'] = 'Other'
mapping_df$domain[mapping_df$domain =='T-box'] = 'T.box'
mapping_df$domain[mapping_df$domain =="NR,IR3"] = 'NR'
mapping_df$domain[mapping_df$domain =="NR,DR1"] = 'NR'
mapping_df$domain[mapping_df$domain =="POU,Homeobox"] = 'Homeobox'
mapping_df$domain[mapping_df$domain =="Paired,Homeobox"] = 'Homeobox'
#write_csv(mapping_df, paste0(storeFigPath, 'homer_tf_domain_mapping.csv'))

homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  homer_df$domain = mapping_df$domain[match(homer_df$`Motif Name`, mapping_df$homer_name)]
  homer_df$domain_col = domain_cols[match(homer_df$domain, names(domain_cols))]
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = basename(i)
  homer_df = as.data.frame(homer_df)
  homer_list[[basename(i)]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
homer_motif_list = lapply(homer_list, '[[', 1)
jj_plot_upsetr(homer_motif_list)
homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif','domain','domain_col', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]
gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))


n_tistreg_no_olap = nrow(jj_load_bed('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/./2023-05-04-peripheral_markers/mouse_tisTreg_no_overlapping_tes.bed'))
n_naivetreg_no_olap = nrow(jj_load_bed('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/./2023-05-04-peripheral_markers/mouse_naiveTreg_no_overlapping_tes.bed'))
n_tistreg_olap = nrow(jj_load_bed('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/./2023-05-03-peripheral_markers/mouse_tisTreg_overlapping_tes.bed'))
n_naivetreg_olap = nrow(jj_load_bed('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/./2023-05-03-peripheral_markers/mouse_naive_Treg_overlapping_tes.bed'))


gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'mouse_naive_Treg_overlapping_tes'= sprintf('naive Treg sig. TE overlap (%i)', n_naivetreg_olap),
  'mouse_naiveTreg_no_overlapping_tes'= sprintf('naive Treg sig. no TE overlap (%i)', n_naivetreg_no_olap),
  'mouse_tisTreg_no_overlapping_tes'= sprintf('tisTreg sig. no TE overlap (%i)', n_tistreg_no_olap),
  'mouse_tisTreg_overlapping_tes'= sprintf('tisTreg sig. TE overlap (%i)', n_tistreg_olap)
))
levels_use = gsea_res_comb %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(enr_sum = sum(enrichment)) %>%
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb_wide = pivot_wider(gsea_res_comb, id_cols = c('name', 'motif'), values_from = 'enrichment',  values_fill = 1, values_fn = mean) %>% 
  as.data.frame %>% column_to_rownames('motif') %>% as.matrix
gsea_res_comb_wide = jj_plot_heatmap(gsea_res_comb_wide, features_use = rownames(gsea_res_comb_wide), group_vec = colnames(gsea_res_comb_wide), return_matrix = T)
gsea_res_comb$name = factor(gsea_res_comb$name, levels = rownames(gsea_res_comb_wide))
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = colnames(gsea_res_comb_wide))


gg = ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() +
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) +
  labs(x = 'Signature subset', y = '', fill = 'Enrichment') +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  annotate(
    xmin = 5,
    xmax = 5.5,
    ymin = seq(0.5, length(unique(gsea_res_comb$motif))-0.5, 1),
    ymax = seq(1.5, length(unique(gsea_res_comb$motif))+0.5, 1),
    geom = "rect",
    fill = gsea_res_comb$domain_col[!duplicated(gsea_res_comb$motif)]
  ) 

library(ggpubr)

ggfake = ggplot(gsea_res_comb, aes(x = domain, fill = domain)) + geom_bar() + 
  scale_fill_manual(values = domain_cols) + labs(fill = 'Domain')
domain_legend  = as_ggplot(get_legend(ggfake))

library(patchwork)
# 6I tisTreg signature peak overlap vs nonoverlapping part TE homer enrichment heatmap
pdf(paste0(storeFigPath, 'mouse_tisTreg_sig_te_olap_nonolap_subsets_homer_heatmap.pdf'),  width = 5, height = 10)
gg + domain_legend #+ plot_layout(widths = c(3, 1))
dev.off()



# te analysis human donor 11/13 -----------------------------------------------

#2 human donor 11
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor11.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
proj$cluster_annotation_level1 = replace_if(proj$cluster_annotation_level1, count_below = 100, 'other')
proj = proj[proj$cluster_annotation_level1 != 'other', ]

dr_df = jj_get_reduction_coords(proj, redname='UMAP')
comb_norm_df = prepare_comb_norm_df(pconfig, proj)
stopifnot(identical(rownames(dr_df), rownames(comb_norm_df)))
# proj$peripheral_tissue = ifelse(proj$Tissue == 'Blood', 'blood',  'peripheral')
# proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level1, sep='__')
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$Tissue)
sum_mat_scaled = scale(t(sum_mat))

sum_mat_cell_type = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$cluster_annotation_level0)
sum_mat_scaled_cell_type = scale(t(sum_mat_cell_type))

pdf(paste0(storeFigPath, 'human_d11_te_tissue_heatmap.pdf'), width =6, height=3)
Heatmap(sum_mat_scaled_cell_type, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = F)
dev.off()
#same picture if only taking eg. CD4_T cells only (and also with the same number of cells being sampled for each group)


pdf(paste0(storeFigPath, 'human_d11_te_celltype_heatmap.pdf'), width =8, height=4)
Heatmap(sum_mat_scaled_cell_type, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = F)
dev.off()


### TE UMAP
comb_df = prepare_comb_norm_df(pconfig, proj, normalize = F)
sobj = CreateSeuratObject(assay = 'RNA', counts = Matrix::t(comb_df), meta.data = as.data.frame(proj@cellColData))
sobj = NormalizeData(sobj)
#sobj = FindVariableFeatures(sobj, selection.method = 'mean.var.plot', dispersion.cutoff = c(0.5, Inf))
#VariableFeaturePlot(sobj)
sobj = ScaleData(sobj, features = rownames(sobj), vars.to.regress = 'nFrags') #cor nFrags with nCount_RNA = 0.97
sobj = RunPCA(sobj, features = rownames(sobj))
#sobj <- JackStraw(sobj, num.replicate = 100)
#sobj <- ScoreJackStraw(sobj, dims = 1:20)
sobj = RunUMAP(sobj, dims = 1:20)
#write_rds(sobj, paste0(bigFilesDir, 'mouse_atlas_te_seurat.RDS')) #23.02.2023
dr_df = jj_get_reduction_coords(sobj, 'umap')
#pca_res = prcomp(comb_norm_df[, !colVars(comb_norm_df)==0], center = T, scale. = T)
#umap_res = uwot::umap(comb_norm_df, pca=5, pca_center = T, scale = T)
#dr_df = cbind(umap_res, as.data.frame(proj@cellColData))
dr_df$cluster_annotation_level0 = from_to(vec= dr_df$cluster_annotation_level1, old_new_map_vec=c(
  'B'= 'B cell',
  'CD4_T'= 'T cell',
  'CD8_T'= 'T cell',
  'DC'= 'DC',
  'Monocyte'= 'Macrophage/Monocyte',
  'Neutrophil/Basophil'= 'Granulocyte',
  'NK'= 'NK cell',
  'other'= 'undefined',
  'T'= 'T cell'
))

col_map = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')
pdf(paste0(storeFigPath, 'human_d11_te_umap.pdf'), width = 10, height=8)
jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation_level0', custom_colors = col_map, do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                 do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'nFrags', cap_top = 'q95')
#jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
dev.off()



#2 human donor 13
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor13.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
#nothing removed by these steps
#proj$cluster_annotation_level0 = replace_if(proj$cluster_annotation_level0, count_below = 100, 'undefined') 
#proj = proj[proj$cluster_annotation_level0 != 'undefined', ]
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
comb_norm_df = prepare_comb_norm_df(pconfig, proj)
stopifnot(identical(rownames(dr_df), rownames(comb_norm_df)))
# proj$peripheral_tissue = ifelse(proj$Tissue == 'Blood', 'blood',  'peripheral')
# proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level1, sep='__')
sum_mat = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$Tissue)
sum_mat_scaled = scale(t(sum_mat))

sum_mat_cell_type = jj_summarize_sparse_mat(Matrix::t(comb_norm_df), proj$cluster_annotation_level0)
sum_mat_scaled_cell_type = scale(t(sum_mat_cell_type))

pdf(paste0(storeFigPath, 'human_d13_te_tissue_heatmap.pdf'), width =6, height=3)
Heatmap(sum_mat_scaled, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = F)
dev.off()
pdf(paste0(storeFigPath, 'human_d13_te_celltype_heatmap.pdf'), width =8, height=4)
Heatmap(sum_mat_scaled_cell_type, name = 'scaled\nmean accessibility', show_column_names = F,
        cluster_rows = F)
dev.off()


### TE UMAP
comb_df = prepare_comb_norm_df(pconfig, proj, normalize = F)
sobj = CreateSeuratObject(assay = 'RNA', counts = Matrix::t(comb_df), meta.data = as.data.frame(proj@cellColData))
sobj = NormalizeData(sobj)
#sobj = FindVariableFeatures(sobj, selection.method = 'mean.var.plot', dispersion.cutoff = c(0.5, Inf))
#VariableFeaturePlot(sobj)
sobj = ScaleData(sobj, features = rownames(sobj), vars.to.regress = 'nFrags') #cor nFrags with nCount_RNA = 0.97
sobj = RunPCA(sobj, features = rownames(sobj))
#sobj <- JackStraw(sobj, num.replicate = 100)
#sobj <- ScoreJackStraw(sobj, dims = 1:20)
sobj = RunUMAP(sobj, dims = 1:20)
#write_rds(sobj, paste0(bigFilesDir, 'mouse_atlas_te_seurat.RDS')) #23.02.2023
dr_df = jj_get_reduction_coords(sobj, 'umap')
#pca_res = prcomp(comb_norm_df[, !colVars(comb_norm_df)==0], center = T, scale. = T)
#umap_res = uwot::umap(comb_norm_df, pca=5, pca_center = T, scale = T)
#dr_df = cbind(umap_res, as.data.frame(proj@cellColData))

col_map = jj_get_colours(dr_df$cluster_annotation_level0, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')
pdf(paste0(storeFigPath, 'human_d13_te_umap.pdf'), width = 10, height=8)
jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation_level0', custom_colors = col_map, do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'),
                 do_shuffle = T)
jj_plot_features(reduction=dr_df, meta_features = 'nFrags', cap_top = 'q95')
#jj_plot_features(reduction=dr_df, meta_features = 'core_tisTreg_sig', cap_top = 'q95', cap_bottom = 'q05')
dev.off()


# te analysis -------------------------------------------------------------

### overview plot
#this includes all TEs whereas in my analysis only TEs on standard chromosomes (and wo Y are used), ie. 880 instead of 960
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)
family_df = family_df %>% .[!duplicated(.), ] %>% arrange(Class, Family)
#write_csv(family_df, paste0(storeFigPath, 'TE_family_class_annotation.csv'))
family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-07-te_analysis_summary/TE_family_class_annotation.csv')
family_df = family_df %>%
  dplyr::group_by(Family) %>% dplyr::add_count() %>% 
  dplyr::select(Class, Family, n) %>% 
  .[!duplicated(.), ] %>%
  dplyr::arrange(desc(n))
family_df$Class = factor(family_df$Class, levels = c('LINE','SINE','LTR','Retroposon'))
family_df = family_df %>% dplyr::arrange(Class)
family_df$Family = factor(family_df$Family, levels = rev(unique(family_df$Family)))
cols_use = jj_get_colours(family_df$Class, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

#S7C TE class and family overview
pdf(paste0(storeFigPath, 'TE_family_overview.pdf'), width = 6, height = 4)
ggplot(family_df, aes(x=Family, y = n, fill = Class)) + geom_bar(stat = 'identity') + coord_flip() + 
  theme_minimal() + labs(y = 'n (TE members)')  + scale_fill_manual(values = cols_use)
dev.off()

pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))

te_bed = read_te('mm10', return_gr = T)
###TE overlap with peaks (average number of peaks per TE?, or fraction of TE sites that have overlap with at least one peak)
pgr = getPeakSet(proj)
peak_df = as.data.frame(unname(pgr))
peak_df$comparison = sprintf('Union peaks')
#peaks that overlap with a TE
pgr$polap = genomic_region_tools::granges_overlap(pgr, te_bed, minOverlap = 0, return_type = 'logical')
#TEs that overlap with a peak
te_bed$polap = genomic_region_tools::granges_overlap(te_bed, pgr, minOverlap = 0, return_type = 'logical')

#peaks that overlap with a specific TE
te_names = unique(te_bed$name)
olap_vec = vector()
for(i in seq_along(te_names)){
  if(i %% 100 == 0 ) message(i)
 te_use = te_bed[ te_bed$name == te_names[i] ]
 olap_vec[i] = sum(genomic_region_tools::granges_overlap(pgr, te_use, minOverlap = 0, return_type = 'logical'))
}
names(olap_vec) = te_names
#write_csv(data.frame(name = names(olap_vec), peaks_olapping_TE=olap_vec), paste0(storeFigPath, 'mouse_atlas_peakset_te_overlap.csv'))
olap_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-04-14-fig1to4_additional_plots/mouse_atlas_peakset_te_overlap.csv')
olap_vec = olap_df$peaks_olapping_TE
names(olap_vec) = olap_df$name
olap_vec_pml = log10( (olap_vec +1) / (length(pgr) + 1))
olap_vec_pct = log10(olap_vec  / length(pgr) * 100)


library(tidyverse)
#wget -c http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz -O /omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.txt.gz
#cd /omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/
#zcat mm10.te.txt.gz | grep -E 'LINE|SINE|LTR|Retroposon' | cut -f5,6-8,11-13 >mm10.te.bed
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)
#write_csv(family_df, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
#family_df$Class[grep('\\?', family_df$Class)] = 'undefined'
family_col_df = family_df %>% dplyr::select(Class, Family) %>% .[!duplicated(.), ] %>% dplyr::arrange(Class, Family)
cols_use_family = jj_get_colours(family_col_df$Family, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')


sum_df = as.data.frame(te_bed@elementMetadata) %>% 
  dplyr::group_by(name) %>%
  dplyr::summarise(total=n(), polap=sum(polap)) %>% 
  dplyr::mutate(fraction=polap/total)

sum_df = sum_df %>% dplyr::left_join(family_df, by = c('name' = 'TE'))

tab_6a_df = sum_df %>% dplyr::relocate(name, Family, Class) %>% 
  dplyr::rename(n_insertions = total, n_peak_overlap = polap, class = Class, family=Family)
#write_csv(tab_6a_df, paste0(storeFigPath, 'table_6a_te_overview.csv'))
#polap = TEs overlapping with any peak, peaks_olapping_TE = peak overlapping with any of the te insertions (by te)
sum_df = sum_df %>% dplyr::left_join(data.frame(percent_olap= olap_vec_pct, peaks_olapping_TE=olap_vec, TE = names(olap_vec_pml)), by = c('name' = 'TE'))
# hist(sum_df$fraction, breaks = 100)
cols_use = jj_get_colours(sum_df$Class, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

# 6A: TE count and peak overlap fraction overview
pdf(paste0(storeFigPath, 'te_abundance_and_overlap_overview.pdf'), width = 10, height = 10)
png(paste0(storeFigPath, 'te_abundance_and_overlap_overview.png'), width = 5, height = 5, res = 400, units = 'in')
# cols_use = jj_get_jj_colours(sum_df$Family)
# cols_use['other'] = 'grey90'
# cols_use['Repetitive sequence'] = 'red'
ggplot(sum_df[sum_df$total > 10, ], aes(y=fraction, x = log10(total))) + 
  #geom_histogram(bins = 300) +
  theme_minimal() + 
  geom_point(aes(colour=Class, shape = Class), size = 1.5, alpha = .8) + 
  scale_colour_manual(values = cols_use) +
  #geom_label(data = sum_df[sum_df$fraction > 0.08 & sum_df$total > 10, ], aes(label=name), size = 3) + 
  geom_text_repel(data = sum_df[sum_df$fraction > 0.1 & sum_df$total > 10, ], aes(label=name), max.overlaps = 500, size = 3) + 
  #geom_point(data = sum_df[sum_df$name == 'RLTR48A', ], aes(y=fraction, x = log10(total)), color = 'red') + 
  coord_cartesian(ylim=c(0, 0.25)) +
  labs(y='Fraction TE insertions with peak overlap' , x= 'log10( total TE insertions in mm10 )')

ggplot(sum_df[sum_df$total > 10, ], aes(y=fraction, x = log10(total))) + 
  #geom_histogram(bins = 300) +
  theme_minimal() + 
  geom_point(aes(colour=Family, shape = Class),  size = 3) + 
  scale_colour_manual(values =cols_use_family) +
  #geom_label(data = sum_df[sum_df$fraction > 0.08 & sum_df$total > 10, ], aes(label=name), size = 3) + 
  geom_text_repel(data = sum_df[sum_df$fraction > 0.1 & sum_df$total > 10, ], aes(label=name), max.overlaps = 500, size = 3) + 
  #geom_point(data = sum_df[sum_df$name == 'RLTR48A', ], aes(y=fraction, x = log10(total)), color = 'red') + 
  coord_cartesian(ylim=c(0, 0.25)) +
  labs(y='Fraction TE insertions with peak overlap' , x= 'log10( total TE insertions in mm10 )')

ggplot(sum_df[sum_df$total > 10, ], aes(y=percent_olap, x = log10(total))) + 
  #geom_histogram(bins = 300) +
  theme_minimal() + 
  geom_point(aes(colour=Family), size = 1.5, alpha = .8) + 
  scale_colour_manual(values = cols_use) +
  geom_text_repel(data = sum_df[sum_df$fraction > 0.1 & sum_df$total > 10, ], aes(label=name), max.overlaps = 500, size = 3) +
  labs(y='log10( % peaks overlapping with TE )' , x= 'log10( total TE insertions in mm10 )')

dev.off()


peak_anno_df = genomic_region_tools::granges_overlap(te_bed, pgr, minOverlap = 0, return_type = 'pairs')
peak_anno_df$name = te_bed$name[as.integer(peak_anno_df$a_index)]
peak_anno_df$peakType = pgr$peakType[as.integer(peak_anno_df$b_index)]
ggplot(peak_anno_df) + geom_bar(aes(x=forcats::fct_infreq(peakType))) + labs(x='Type of associated peaks', y = 'Count')
peak_anno_df$Group = 'All TE-overlapping peaks'
jj_plot_categorical_by_group(peak_anno_df, feature_column = 'peakType', group_column = 'Group', add_text = T)

sum_df2 = jj_plot_categorical_by_group(peak_anno_df, feature_column = 'peakType', group_column = 'name', return_df = T)
sum_df2[sum_df2$name=='RLTR48A', ]

tesplot = sum_df$name[sum_df$polap >=10]
jj_plot_categorical_by_group(peak_anno_df[peak_anno_df$name %in% tesplot, ], feature_column = 'peakType', group_column = 'name')


### peak-independent genomic feature overlaps of TEs

pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))

###TE overlap with peaks (average number of peaks per TE?, or fraction of TE sites that have overlap with at least one peak)
gannot = getGeneAnnotation(proj)
pgr = getPeakSet(proj)

te_bed = read_te('mm10', return_gr = T)
te_bed$polap = genomic_region_tools::granges_overlap(te_bed, pgr, minOverlap = 0, return_type = 'logical')

library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
peaks_list = list(all_te = te_bed, all_peaks = unname(pgr), te_in_peaks=te_bed[te_bed$polap])
peaks_anno_list = list()
for(i in names(peaks_list)){
  peaks_anno_list[[i]] <- annotatePeak(peaks_list[[i]],
                               tssRegion=c(-2000, 500),
                               TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene, 
                               annoDb="org.Mm.eg.db", 
                               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                             "Downstream", "Intergenic"))
}

peaks_anno_df = jj_initialize_df(ncol = length(peaks_anno_list), nrow = 6, init = 0, col.names = names(peaks_anno_list))
for(i in seq_along(peaks_anno_list)){
  peak_annot = peaks_anno_list[[i]]@detailGenomicAnnotation
  peak_annot_vec =  c(Promoter = sum(peak_annot$Promoter), #promoter sum
                      `5' UTR` = sum(peak_annot$fiveUTR & !peak_annot$Promoter), #5UTR
                      `3' UTR` = sum(peak_annot$threeUTR & !peak_annot$fiveUTR & !peak_annot$Promoter), #3UTR
                      Exon = sum(peak_annot$Exon & !peak_annot$threeUTR & !peak_annot$fiveUTR & !peak_annot$Promoter), #Exon
                      Intron = sum(peak_annot$Intron & !peak_annot$Exon & !peak_annot$threeUTR & !peak_annot$fiveUTR & !peak_annot$Promoter), #Intron
                      `Distal Intergenic` = sum(peak_annot$distal_intergenic & !peak_annot$Intron & !peak_annot$Exon & !peak_annot$threeUTR & !peak_annot$fiveUTR & !peak_annot$Promoter)) #distal intergenic
  peaks_anno_df[, i] = peak_annot_vec
  
}
rownames(peaks_anno_df) = names(peak_annot_vec)
peaks_anno_df['Total', ] = apply(peaks_anno_df, 2, sum)
#S7D Table of genomic annoations for TEs and peaks
#write.csv(peaks_anno_df, paste0(storeFigPath, 'mouse_atlas_te_genomic_annotations.csv'))

#library(RColorBrewer)
#6B TE genomic region annotation barplot
pdf(paste0(storeFigPath, 'mouse_atlas_te_genomic_annotations.pdf'), width = 6, height = 3)
cols_use = brewer.pal(6, 'RdYlBu')
names(cols_use) = c('0-1kb', '1-3kb' ,'3-5kb', '5-10kb','10-100kb' ,'>100kb')
plotDistToTSS(peaks_anno_list,
              title="Distribution of regions relative to TSS") +
  #viridis::scale_fill_viridis(discrete = TRUE) 
  scale_fill_manual(values = cols_use)

cols_use = c('#D73027', '#FC8D59', "lightgreen", "darkgreen", "purple", "violet", "burlywood2", "burlywood4", "#91BFDB", "#4575B4")
names(cols_use) = c("Promoter (<=1kb)", "Promoter (1-2kb)", "5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron",
                    "Other Intron", "Downstream (<=300)", "Distal Intergenic")
plotAnnoBar(peaks_anno_list, title = "Genomic region annotation") + 
  scale_fill_manual(values = cols_use)
dev.off()
#upsetplot(peakAnno.edb)


# tisTreg vs tisTreg precursors (not used) --------------------------------

# .libPaths(c("/home/simonma/R/x86_64-pc-linux-gnu-library/3.6",
#             "/omics/groups/OE0436/data2/simonma/R/x86_64-pc-linux-gnu-library/3.6"
# ))
# 
# library(Seurat)
# library(Signac)
# library(GenomeInfoDb)
# library(tidyverse)
# library(maltesFunctions)
# library(genomic_region_tools)
# library(jj)
# 
# choose_seurat_file('mouse_normal', 'CD4')
# contents_use <- pick_content(dataset_cell_type, c('seurat_file', 'fragment_path', 'ident_use'))
# 
# seurat_file <- contents_use[1]
# ident_use <- contents_use[3]
# 
# seurat_atac <- read_rds(seurat_file)
# dr_df = jj_get_reduction_coords(seurat_atac, 'umap')
# jj_plot_features(reduction = dr_df, meta_features = ident_use)
# Idents(seurat_atac) <- seurat_atac@meta.data[, ident_use]
# DefaultAssay(seurat_atac) = 'scATAC_raw'
# 
# mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
# mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
# mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
# 
# da_peaks <- FindMarkers(
#   object = seurat_atac,
#   ident.1 = c('16', '23'),
#   ident.2 = c('11', '5'),
#   min.pct = 0.1, 
#   test.use = 'LR',
#   latent.vars = 'n_peaks'
# )
# da_peaks$feature <- rownames(da_peaks)
# da_peaks <- da_peaks %>%
#   as_tibble %>% 
#   dplyr::filter(p_val_adj < 0.05) %>% 
#   dplyr::mutate(comparison = ifelse(avg_logFC > 0, comp1, comp2)) %>% 
#   dplyr::group_by(comparison) %>% 
#   dplyr::arrange(p_val_adj, .by_group = TRUE) #groups are ignored unless specified by .by_group = TRUE
# 
# gene.coords <- read_rds('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/1_res/mm10_genes_gr.RDS')
# clFeat <- as_tibble(Signac::ClosestFeature(regions=da_peaks$feature, annotation = gene.coords))
# da_peaks <- bind_cols(da_peaks, clFeat) 
# #write_csv(da_peaks, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-12-08-mouse_precursor_scatac_subset/diff_results_16_23_versus_5_11_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
# 
# jj_compare_vectors(mouse_cd4_scATAC_sig_df$feature, da_peaks$feature)
# 
# mouse_cd4_scATAC_sig_df$precursor_shared = !mouse_cd4_scATAC_sig_df$feature %in% da_peaks$feature
# #write_csv(mouse_cd4_scATAC_sig_df, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-12-08-mouse_precursor_scatac_subset/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_precursor_shared_annot.csv')
# 
# sig_list = list(precursor_shared = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14' & mouse_cd4_scATAC_sig_df$precursor_shared, ], 
#                 tisTreg_specific = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14' & !mouse_cd4_scATAC_sig_df$precursor_shared, ])
# sig_list = lapply(sig_list, function(x) convert_granges(x$feature))
# 
# pmat = GetAssayData(seurat_atac, assay = 'scATAC_raw', slot = 'data')
# pmat = pmat[rownames(pmat) %in% mouse_cd4_tisTreg_df$feature, ]
# pmat = pmat[!rownames(pmat) %in% mouse_cd4_scATAC_sig_df$feature[mouse_cd4_scATAC_sig_df$precursor_shared], ]
# summarize_vec = seurat_atac$snn_harmony_res.1.7
# mean_mat = jj_summarize_sparse_mat(pmat, summarize_vec)
# mean_mat = mean_mat[, colnames(mean_mat) %in% c('0','3','14','11','5','16','23','10')]
# htmat = scale(t(mean_mat))
# Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, 
#         show_column_names = F)


# nebula peaks plots ------------------------------------------------------

# plot signatures in T/NK/ILC subset --------------------------------------


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

dr_df = jj_get_reduction_coords(proj, 'UMAP')

sig_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-26-th1_signature/tnkilc_nebula_fragment_count_da_peaks_lfc0_fdr05_complete_nonbinary_peakset.xlsx')
#jj_save_excel(sig_list, paste0(storeFigPath, 'table_S6a_nebula_tnkilc_subset_da_peaks.xlsx'))

# #only exclusive peaks
# sig_list_plot = list()
# for(i in seq_along(sig_list)){
#   peaks_keep = sig_list[[i]] 
#   peaks_keep = peaks_keep[peaks_keep$exclusive & peaks_keep$logFC > 1, ]
#   sig_list_plot[[ names(sig_list)[i] ]] = peaks_keep
# }
# summary_df = sapply(sig_list_plot, nrow)
# summary_df = data.frame(variable = names(summary_df), n_peaks = summary_df)
# my_cols = rep('grey50', length(summary_df$variable))
# names(my_cols) = summary_df$variable
# my_cols2 = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')
# my_cols[names(my_cols2)] = my_cols2
# 
# pdf(paste0(storeFigPath, 'tnkilc_nebula_exclusive_peaks_summary.pdf'), width = 4, height = 4)
# ggplot(summary_df) + geom_bar(aes(x = reorder(variable, n_peaks), y = n_peaks, fill = variable), stat = 'identity') + coord_flip() +
#   scale_fill_manual(values = my_cols) + theme_minimal()  + theme(legend.position = 'none') + labs(x = '', y = 'n (exclusive peaks)')
# dev.off()

#significant peaks with logFC > 1
#difference with exclusive: if logFC is filtered after searching for exclusive, then the result is more stringent
sig_list_plot = list()
for(i in seq_along(sig_list)){
  peaks_keep = sig_list[[i]]
  peaks_keep = peaks_keep[peaks_keep$logFC > 1, ]
  sig_list_plot[[ names(sig_list)[i] ]] = peaks_keep
}

upset_list = lapply(sig_list_plot, '[[', 1)
pdf(paste0(storeFigPath, 'nebula_tnkilc_peaks_upset.pdf'),  width = 11, height = 6)
#
jj_plot_upsetr(upset_list, nintersects = 48)
dev.off()

summary_df = sapply(sig_list_plot, nrow)
summary_df = data.frame(variable = names(summary_df), n_peaks = summary_df)
my_cols = rep('grey50', length(summary_df$variable))
names(my_cols) = summary_df$variable
my_cols2 = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')
my_cols[names(my_cols2)] = my_cols2

#use exclusive peaks for the tissue signatures (skin 512, vat 355, colon 1129)
sig_list_keep = list()
for(i in names(sig_list_plot)){
  message(i,": ", nrow(sig_list_plot[[i]]))
  peaks_use = sig_list_plot[[i]]$feature
  others = unlist(lapply(sig_list_plot[names(sig_list_plot) != i], '[[', 'feature'))
  sig_list_keep[[i]] = sig_list_plot[[i]][!peaks_use %in% others, ]
  message(i, " unique: ", nrow(sig_list_keep[[i]]))
}
#save as bed for running homer analysis
tissue_sig_gr = lapply(sig_list_keep[c('Colon', 'VAT', 'Skin')], function(x) convert_granges(x$feature))
for(i in seq_along(tissue_sig_gr)){
  jj_save_bed(tissue_sig_gr[[i]], file_name = paste0(storeFigPath, names(tissue_sig_gr)[i], '.bed'))
}
proj = archr_add_peak_signatures(proj, signature_list = tissue_sig_gr, signature_name = 'nebula_atlas_tissue')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(dr_df, features = grep('z_', colnames(dr_df), value = T), cap_top = 'q99', cap_bottom = 'q01', return_gg_object = T)
pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_nebula_exclusive_peaks_chromvar_scores.pdf'),  width = 8, height = 6)
gg
dev.off()

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df_list = list()
for(i in names(tissue_sig_gr)){
  message(i)
  olap_df_list[[i]] = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, 
                                        nFrags_vec = proj$nFrags, signature_gr = tissue_sig_gr[[i]],
                                        verbose = T, k = 100, count_thres = 2e5)
}

#write_rds(olap_df_list, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_tnkilc_nebula_tissue_exlusive_pct_olap_list.RDS')
olap_df_list = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_tnkilc_nebula_tissue_exlusive_pct_olap_list.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#S6B Nebula tissue signature percentage overlap
pdf('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-12-finalizing_te_analysis/umap_mouse_tnkilc_nebula_exlcusive_peaks_pct_olap.pdf',  width = 8, height = 6)
for(i in names(olap_df_list)){
  dr_df[[i]] = olap_df_list[[i]]$signature_pct_overlap #signature_pct_overlap
  print(jj_plot_features(dr_df, features = i, return_gg_object = T, my_title = i)[[1]] + labs(colour = '% overlap'))
}
dev.off()


library(tidyverse)
library(clusterProfiler)
library(jj)

#M2_CP
go_collection = clusterProfiler::read.gmt("/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/msigdb/m2.cp.v2022.1.Mm.symbols.gmt")
go_res =  enricher(
    gene = unique(sig_list_plot$Colon$symbol),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE = go_collection)
#S6E dotplot is sorted by GeneRatio
pdf(paste0(storeFigPath, 'mouse_tnkilc_nebula_exclusive_colon_m2_cp_go_dotplot.pdf'),  width = 8, height = 6)
dotplot(go_res, showCategory = 15, font.size = 8 )
dev.off()
go_res_df = go_res@result

#over 400 genes in this set
pathway_genes_use = go_res_df$geneID[go_res_df$ID == 'REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM']
pathway_genes_use = unlist(strsplit(pathway_genes_use, '/'))
go_complete = go_collection$gene[go_collection$term=='REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM']
go_complete = go_complete[go_complete %in% rownames(gmat)]

proj= ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = list(cytokine_signaling=pathway_genes_use, cytokine_signaling_full=go_complete), name = 'go')
dr_df = jj_get_reduction_coords(proj, 'UMAP')

#S6F
pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_reactome_cytokine_signaling_in_immune_system.pdf'),  width = 8, height = 6)
jj_plot_features(dr_df, features = 'go.cytokine_signaling', cap_top = 'q99', cap_bottom = 'q01', return_gg_object = T)[[1]] +
  labs(colour = 'REACTOME\nCYTOKINE SIGNALING\nIN IMMUNE SYSTEM')
jj_plot_features(dr_df, features = 'go.cytokine_signaling_full', cap_top = 'q99', cap_bottom = 'q01', return_gg_object = T)[[1]] +
  labs(colour = 'Complete REACTOME\nCYTOKINE SIGNALING\nIN IMMUNE SYSTEM')
dev.off()


gmat = get_gene_mat(proj)


proj= ArchR::addModuleScore(proj, useMatrix = 'GeneScoreMatrix', features = list(cytokine_signaling=go_complete), name = 'go_complete')

pdf(paste0(storeFigPath, 'heatmap_mouse_tnkilc_reactome_cytokine_signaling_in_immune_system.pdf'),  width = 6, height = 7)
jj_plot_heatmap(gmat, pathway_genes_use, proj$Tissue, transpose = T)
dev.off()


### plot tissue signatures in mouse cd4 and areg gfp datasets

### mouse areg gfp
library(yaml)
pconfig = yaml::read_yaml(config_list$mouse_areg_gfp)
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)
proj = archr_add_peak_signatures(proj, signature_list = tissue_sig_gr, signature_name = 'nebula_atlas_tissue')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(dr_df, features = grep('z_', colnames(dr_df), value = T), cap_top = 'q99', cap_bottom = 'q01', return_gg_object = T)
gg2 = jj_plot_features(dr_df, features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv', comment_char = '#'), return_gg_object = T)
pdf(paste0(storeFigPath, 'umap_mouse_areg_gfp_nebula_exclusive_peaks_chromvar_scores.pdf'),  width = 8, height = 6)
gg
gg2
dev.off()

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df_list = list()
for(i in names(tissue_sig_gr)){
  message(i)
  olap_df_list[[i]] = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, 
                                             nFrags_vec = proj$nFrags, signature_gr = tissue_sig_gr[[i]],
                                             verbose = T, k = 100, count_thres = 2e5)
}

#write_rds(olap_df_list, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_areg_gfp_nebula_tissue_exlusive_pct_olap_list.RDS')
olap_df_list = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_areg_gfp_nebula_tissue_exlusive_pct_olap_list.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#S6D
pdf(paste0(storeFigPath, 'umap_mouse_areg_gfp_nebula_exlcusive_peaks_pct_olap.pdf'),  width = 8, height = 6)
for(i in names(olap_df_list)){
  dr_df[[i]] = olap_df_list[[i]]$signature_pct_overlap
  print(jj_plot_features(dr_df, features = i, return_gg_object = T, my_title = i)[[1]] + labs(colour = '% overlap'))
}
dev.off()

### mouse normal cd4
pconfig = yaml::read_yaml(config_list$mouse_normal_cd4)
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)
proj = archr_add_peak_signatures(proj, signature_list = tissue_sig_gr, signature_name = 'nebula_atlas_tissue')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
gg = jj_plot_features(dr_df, features = grep('z_', colnames(dr_df), value = T), cap_top = 'q99', cap_bottom = 'q01', return_gg_object = T)
gg2 = jj_plot_features(dr_df, features = 'Tissue', custom_colors = jj_get_colours(dr_df$Tissue, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'), return_gg_object = T)
pdf(paste0(storeFigPath, 'umap_mouse_normal_cd4_nebula_exclusive_peaks_chromvar_scores.pdf'),  width = 8, height = 6)
gg
gg2
dev.off()

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df_list = list()
for(i in names(tissue_sig_gr)){
  message(i)
  olap_df_list[[i]] = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df, 
                                             nFrags_vec = proj$nFrags, signature_gr = tissue_sig_gr[[i]],
                                             verbose = T, k = 100, count_thres = 2e5)
}

#write_rds(olap_df_list, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_normal_cd4_nebula_tissue_exlusive_pct_olap_list.RDS')
olap_df_list = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mouse_normal_cd4_nebula_tissue_exlusive_pct_olap_list.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#S6C
pdf(paste0(storeFigPath, 'umap_mouse_normal_cd4_nebula_exlcusive_peaks_pct_olap.pdf'),  width = 8, height = 6)
for(i in names(olap_df_list)){
  dr_df[[i]] = olap_df_list[[i]]$signature_pct_overlap
  print(jj_plot_features(dr_df, features = i, return_gg_object = T, my_title = i)[[1]] + labs(colour = '% overlap'))
}
dev.off()



# homer analysis on colon and skin tissue signatures ----------------------

mapping_df= read_csv("/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap/homer_tf_domain_mapping.csv")
domain_cols = jj_get_colours(mapping_df$domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = T)

homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  homer_df$domain = mapping_df$domain[match(homer_df$`Motif Name`, mapping_df$homer_name)]
  homer_df$domain_col = domain_cols[match(homer_df$domain, names(domain_cols))]
  #1= size 2, 20+ = size 0.2
  homer_df$name = basename(i)
  homer_df = as.data.frame(homer_df)
  homer_list[[basename(i)]] = homer_df # homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif','domain','domain_col', 'enrichment','q-value (Benjamini)',  'rank')])
homer_list_gsea_plot = homer_list_gsea_plot[c('Colon', 'Skin')]

plot_df = dplyr::full_join(homer_list_gsea_plot$Colon, homer_list_gsea_plot$Skin, by = 'motif', suffix = c('_Colon', '_Skin'))
plot_df$signif_Colon = plot_df$`q-value (Benjamini)_Colon` < 0.001
plot_df$signif_Skin = plot_df$`q-value (Benjamini)_Skin` < 0.001
plot_df$significant = 'n.s.'
plot_df$significant[plot_df$signif_Colon] = 'Colon'
plot_df$significant[plot_df$signif_Skin] = 'Skin'
plot_df$significant[plot_df$signif_Skin & plot_df$signif_Colon] = 'Both'
plot_df$significant = factor(plot_df$significant, levels = c('Both', 'Colon', 'Skin', 'n.s.'))
plot_df$motif = sapply(strsplit(plot_df$motif, split = '/'), '[[', 1) 
# plot_df$log_FDR = 0
# plot_df$log_FDR[plot_df$signif_Colon] = (-1) * log10(plot_df$`q-value (Benjamini)_Colon`[plot_df$signif_Colon])
# plot_df$log_FDR[plot_df$signif_Skin] = (-1) * log10(plot_df$`q-value (Benjamini)_Skin`[plot_df$signif_Skin])
# plot_df$log_FDR[plot_df$signif_Colon & plot_df$signif_Skin] = (-1) * log10(pmax(plot_df$`q-value (Benjamini)_Colon`[plot_df$signif_Colon & plot_df$signif_Skin], plot_df$`q-value (Benjamini)_Skin`[plot_df$signif_Colon & plot_df$signif_Skin]))
# plot_df$log_FDR[plot_df$log_FDR == Inf] = 5
# plot_df$log_FDR = min_max_normalize(plot_df$log_FDR) *3
cols_use = c('n.s.' = 'grey', 'Colon' = 'purple', 'Skin' = 'darkblue', 'Both' = 'red')

pdf(paste0(storeFigPath, 'colon_skin_nebula_signature_homer_tf_scatterplot.pdf'), width = 10, height = 9)
ggplot(plot_df, aes(x = enrichment_Colon, y = enrichment_Skin)) +
  geom_point(aes(colour = significant)) + #size = log_FDR
  geom_text_repel(data = plot_df[plot_df$significant != 'n.s.', ], aes(label = motif), size = 3, max.overlaps = 100) + 
  scale_colour_manual(values = cols_use) + theme_minimal() + 
  labs(x='Enrichment (Colon)', y = 'Enrichment (Skin)', colour = 'Significant', title = 'Homer known motif enrichment') + 
  theme(plot.title = element_text(hjust = 0.5)) #+ 
  #scale_size_area(na.value = 0.2)

dev.off()







# tisTreg signature TE overlap fisher test -------------------------------------

# calculate te enrichment in mouse cd4 tisTreg sig

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
seurat_atac = read_rds(pick_content('mouse_normal_CD4', 'seurat_file'))
mouse_cd4_universe_gr = convert_granges(rownames(GetAssayData(seurat_atac, assay = 'scATAC_raw')))

te_gr = read_te('mm10', return_gr = T)

mouse_cd4_universe_gr$tisTreg_sig_olap = granges_overlap(mouse_cd4_universe_gr, mouse_cd4_tisTreg_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)

### on the family basis instead of TE basis (replace te name with family name)
te_gr$family = family_df$Family[match(te_gr$name, family_df$TE)]
te_gr$name = te_gr$family
###

te_names = unique(te_gr$name)

olaps_df = jj_initialize_df(5, length(te_names), init = 0, col.names = c('both','te_only','tisTreg_sig_only','none', 'pval'), row.names=te_names)
for(i in 1:length(te_names)){
  
  te_test_gr = te_gr[te_gr$name == te_names[i]]
  mouse_cd4_universe_gr$te_olap = granges_overlap(mouse_cd4_universe_gr, te_test_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)
  if(is.null(mouse_cd4_universe_gr$te_olap)){
    mouse_cd4_universe_gr$te_olap = FALSE
  }
  present_both = sum(mouse_cd4_universe_gr$tisTreg_sig_olap & mouse_cd4_universe_gr$te_olap)
  present_te_only = sum(!mouse_cd4_universe_gr$tisTreg_sig_olap & mouse_cd4_universe_gr$te_olap)
  present_batf_only = sum(mouse_cd4_universe_gr$tisTreg_sig_olap & !mouse_cd4_universe_gr$te_olap)
  present_none = sum(!mouse_cd4_universe_gr$tisTreg_sig_olap & !mouse_cd4_universe_gr$te_olap)
  
  cont_df = data.frame('te_yes' = c(present_both, present_te_only),
                       'te_no' = c(present_batf_only, present_none),
                       row.names = c('tisTreg_sig_yes', 'tisTreg_sig_no'))
  
  test <- fisher.test(cont_df)
  olaps_df[i,] = c(present_both, present_te_only, present_batf_only, present_none, test$p.value)
}

olaps_df$fraction_tisTreg_sig_olap = olaps_df$both / length(mouse_cd4_tisTreg_gr)
olaps_df$fraction_other_peaks_olap = olaps_df$te_only / (length(mouse_cd4_universe_gr) - length(mouse_cd4_tisTreg_gr))
olaps_df$relative_enrichment = olaps_df$fraction_tisTreg_sig_olap / olaps_df$fraction_other_peaks_olap
olaps_df = olaps_df[order(olaps_df$relative_enrichment, decreasing = T), ]

te_n_df_mouse = as.data.frame(table(te_gr$name))
colnames(te_n_df_mouse) = c('name', 'n_te')
olaps_df = olaps_df %>% 
  rownames_to_column('name') %>% 
  dplyr::left_join(te_n_df_mouse, by = 'name') %>% 
  dplyr::mutate(fraction_olap_both = both / n_te)
olaps_df$fdr = p.adjust(olaps_df$pval, method = 'fdr')
olaps_df = olaps_df %>% 
  relocate(name, n_te, both, te_only, tisTreg_sig_only, none, pval, fdr, relative_enrichment) %>% 
  arrange(fdr)
olaps_df$odds_ratio = (olaps_df$both * olaps_df$none) / (olaps_df$te_only * olaps_df$tisTreg_sig_only)
View(olaps_df)

#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_naiveTreg_sig_vs_universe_te_overlap_df.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_naiveTreg_sig_vs_universe_te_overlap_df.csv'))
olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_df$odds_ratio = (olaps_df$both * olaps_df$none) / (olaps_df$te_only * olaps_df$tisTreg_sig_only)
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
olaps_sig_df$logfdr = -log10(olaps_sig_df$fdr)
family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
olaps_sig_df = olaps_sig_df %>% dplyr::left_join(family_df, by = c('name' = 'TE'))
family_col_df = family_df %>% dplyr::select(Class, Family) %>% .[!duplicated(.), ] %>% dplyr::arrange(Class, Family)
cols_use_family = jj_get_colours(family_col_df$Family, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

#6G Fisher test enrichment of TE insertions in the tisTreg signature
pdf(paste0(storeFigPath, 'tisTreg_sig_te_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
jj_plot_volcano(olaps_sig_df, logfc_column = 'odds_ratio',
                pval_column = 'fdr', symbol_column = 'name',  marker_thres = c(2, 10), labs_range = c(0,10, 0, 10)) +
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'Odds ratio') + 
  scale_colour_manual(
    values = c(">-2 & <2" = 'black', ">= 2" = 'blue'),
    breaks = c(">= 2")
  ) 
ggplot() + geom_point(data = olaps_sig_df, aes(x = odds_ratio, y = logfdr, colour = Family, shape = Class), size = 3) + 
  scale_colour_manual(values = cols_use_family) + theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_sig_df[olaps_sig_df$odds_ratio > 2 | olaps_sig_df$logfdr > 5, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + labs(x = 'Odds ratio', y = '-log10(FDR)')

dev.off()



###for family level (use the family result)
olaps_df = read_csv(paste0("/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/", 'mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv'))
olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
olaps_plot_df$Class = family_df$Class[match(olaps_plot_df$name, family_df$Family)]

#S7G tisTreg TE Family enrichment volcano plot
pdf(paste0(storeFigPath, 'tisTreg_sig_te_family_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

ggplot() + geom_point(data = olaps_plot_df[olaps_plot_df$both > 10, ], aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$both > 10, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

dev.off()

### comparison with nebula TEs

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
#writeLines(olaps_sig_df$name, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-19-footprints/tisTreg_sig_61_enriched_tes.txt'))

res_list = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_tnkilc_20230605.RDS')
res_list = get_nebula_df(nebula_res_df = res_list$summary)
nebula_list = convert_nebula_df_to_list(res_list, fdr_thres = 0.05, logFC_thres = 0)
da_te_list = lapply(nebula_list, '[[', 'feature')
tisTreg_da_te = da_te_list$cluster_annotation_finetisTregST2
gg = jj_plot_upsetr(list(nebula_fdr005_logFC0 = tisTreg_da_te, tisTreg_signature = olaps_sig_df$name))


#S8C
jj_plot_upsetr(lapply(convert_nebula_df_to_list(res_list, fdr_thres = 0.05, logFC_thres = 0.5),'[[', 'feature'))


#intersect(tisTreg_da_te, olaps_sig_df$name)
s1 = sample(1:880, length(tisTreg_da_te))
olap_found = vector()
for(i in 1:10000){
  s2 = sample(1:880, 61)
  olap_found[i] = sum(s1 %in% s2)
}
#summary(olap_found)
gg1 = jj_plot_hist(olap_found, thres = length(intersect(tisTreg_da_te, olaps_sig_df$name))) + 
  labs(x = 'n overlapping TEs', title = 'expected overlap distribution')

nebula_list = convert_nebula_df_to_list(res_list, fdr_thres = 0.05, logFC_thres = 0.5)
da_te_list = lapply(nebula_list, '[[', 'feature')
tisTreg_da_te2 = da_te_list$cluster_annotation_finetisTregST2
gg2 = jj_plot_upsetr(list(nebula_fdr005_logFC05 = tisTreg_da_te2, tisTreg_signature = olaps_sig_df$name))

#intersect(tisTreg_da_te, olaps_sig_df$name)
s1 = sample(1:880, length(tisTreg_da_te2))
olap_found = vector()
for(i in 1:10000){
  s2 = sample(1:880, 61)
  olap_found[i] = sum(s1 %in% s2)
}
#summary(olap_found)
gg3 = jj_plot_hist(olap_found, thres = length(intersect(tisTreg_da_te2, olaps_sig_df$name))) + 
  labs(x = 'n overlapping TEs', title = 'expected overlap distribution')
pdf(paste0(storeFigPath, 'nebula_tisTreg_TE_mouse_tisTregST2_signature_TE_overlaps.pdf'), width = 4, height = 3)
gg
gg1
gg2
gg3
dev.off()
quantile(olap_found, 0.95)

tes_plot = intersect(tisTreg_da_te, olaps_sig_df$name)

re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_tnkilc_20230605.RDS')
#combine dr_df with te matrix to plot tes in umap space
te_dr_df = cbind(dr_df, prepare_comb_norm_df(pconfig, proj, normalize = T))

te_heatmap_df = te_dr_df[, colnames(te_dr_df) %in% tes_plot]
### Heatmap per cell type
#use newest version of this function: there is a bug in the ordering in the old one!
sum_mat = jj_summarize_sparse_mat(Matrix::t(te_heatmap_df), summarize_by_vec = te_dr_df$cluster_annotation_fine)
sum_mat = sum_mat[, !colnames(sum_mat) %in% c('undefined')]
sum_mat_scaled = scale(t(sum_mat))
colnames(sum_mat_scaled)[!complete.cases(t(sum_mat_scaled))]
sum_mat_scaled = t(t(sum_mat_scaled)[complete.cases(t(sum_mat_scaled)), ])

family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
colcsv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'

col_list_use = list(Class = jj_get_colours(family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)],colour_csv = colcsv, comment_char = '$'),
                    Family = jj_get_colours(family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)], colour_csv = colcsv, comment_char = '$'))
column_ha = HeatmapAnnotation(Class = family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)], 
                              Family = family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)],
                              col = col_list_use)

ht = Heatmap(sum_mat_scaled,  name = 'scaled\nmean accessibility',
             show_column_names = T, top_annotation = column_ha, cluster_columns = T, cluster_rows = T)
pdf(paste0(storeFigPath, 'mouse_tnkilc_nebula_fdr005_tisTreg_sig_fisher_fdr005_30te_heatmap.pdf'), width = 10, height=7)
draw(ht)
dev.off()



# tnkilc nebula TE heatmap and dotplot ------------------------------------


config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/nebula_TE_tnkilc_20230605.RDS')
#combine dr_df with te matrix to plot tes in umap space
te_dr_df = cbind(jj_get_reduction_coords(proj, 'UMAP'), prepare_comb_norm_df(pconfig, proj, normalize = T))
res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_fine|Tissue', '', colnames(res_df))
jj_fc_fc_plot(res_df[res_df$Colon_FDR < 0.05 & res_df$Skin_FDR < 0.05, ], logfc_column1 = 'Colon_logFC', logfc_column2 = 'Skin_logFC', symbol_column = 'feature', marker_thres = 0.75, use_text = T)
# fc_fc_plot_broad(res_df, group1 = 'Colon_logFC', group2 = 'Skin_logFC', symbol_column = 'feature')

res_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0.5)
sapply(res_list, nrow)

da_te_list = lapply(res_list, '[[', 'feature')
da_te_list = da_te_list[!names(da_te_list) %in% c('VAT','Skin', 'Colon')]

pdf(paste0(storeFigPath, 'mouse_tnkilc_te_fdr005_logfc05_upset.pdf'), width =6, height=5)
#reduce nr of intersects until only intersections of size > 1 are shown
jj_plot_upsetr(da_te_list, nintersects = 21)
dev.off()

family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
colcsv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'

marker_te_groups = c('tisTregST2', 'ILC2', 'CD4_other', 'Tfh-like')
for(k in marker_te_groups){
  message(k)
  tes_plot =  da_te_list[[k]] # unique(unlist(lapply(da_te_list[c('Colon','Skin', 'VAT')], head, 10)))
  te_heatmap_df = te_dr_df[, colnames(te_dr_df) %in% tes_plot]
  ### Heatmap per cell type
  #use newest version of this function: there is a bug in the ordering in the old one!
  sum_mat = jj_summarize_sparse_mat(Matrix::t(te_heatmap_df), summarize_by_vec = te_dr_df$cluster_annotation_fine)
  sum_mat_scaled = scale(t(sum_mat))
  colnames(sum_mat_scaled)[!complete.cases(t(sum_mat_scaled))]
  sum_mat_scaled = t(t(sum_mat_scaled)[complete.cases(t(sum_mat_scaled)), ])
  
  col_list_use = list(Class = jj_get_colours(family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)],colour_csv = colcsv, comment_char = '$'),
                      Family = jj_get_colours(family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)], colour_csv = colcsv, comment_char = '$'))
  column_ha = HeatmapAnnotation(Class = family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)], 
                                Family = family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)],
                                col = col_list_use)
  
  ht = Heatmap(sum_mat_scaled,  name = 'scaled\nmean accessibility', show_column_names = T, top_annotation = column_ha)
  ht = draw(ht)
  celltype_order = rownames(sum_mat_scaled)[row_order(ht)]
  te_order = colnames(sum_mat_scaled)[column_order(ht)]
  # celltype_order = rownames(sum_mat_scaled)[hclust(dist(sum_mat_scaled, method = 'euclidean'), method = 'complete')$order]
  # te_order = colnames(sum_mat_scaled)[hclust(dist(t(sum_mat_scaled, method = 'euclidean')), method = 'complete')$order]
  # Heatmap(sum_mat_scaled[match(celltype_order, rownames(sum_mat_scaled)), match(te_order, colnames(sum_mat_scaled))], 
  #         cluster_rows = F, cluster_columns = F)
  
  dotplot_df = res_df %>% rownames_to_column('feature') %>% dplyr::filter(feature %in% tes_plot) 
  dotplot_df = dotplot_df %>% pivot_longer(cols = -1, names_to = 'name', values_to = 'value') %>% dplyr::mutate(info_type = ifelse(grepl('_FDR', name), 'FDR', 'logFC')) %>%
    dplyr::mutate(name = gsub('_FDR|_logFC', '', name)) %>% 
    spread(key = info_type, value = value)
  dotplot_df = dotplot_df %>% dplyr::filter(!name %in% c('VAT','Skin', 'Colon'))
  dotplot_df$logFDR <- (-1) * log10(dotplot_df$FDR)
  dotplot_df$signif = dotplot_df$FDR < 0.05
  dotplot_df$logFC_capped = jj_cap_vals(dotplot_df$logFC, cap_bottom = -1, cap_top = NULL)
  
  dotplot_df$feature = factor(dotplot_df$feature, levels = te_order)
  dotplot_df$name = factor(dotplot_df$name, levels = celltype_order)
  
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
  #S8D S8E S8F
  pdf(paste0(storeFigPath, sprintf('mouse_tnkilc_te_celltype_heatmap_%s_tes.pdf', k)), width = 1+0.22*length(tes_plot), height=6)
  print(ht)
  dev.off()
  pdf(paste0(storeFigPath, sprintf('mouse_tnkilc_te_celltype_dotplot_%s_tes.pdf', k)), width = 8 , height=1+0.2*length(tes_plot))
  print(nebula_dotplot)
  dev.off()
}


#63 TEs, 25 unique TEs
top_use = 10
top_te_list = lapply(res_list[!names(res_list) %in% c('Colon','Skin','VAT')], function(x) head(x$feature, top_use))
jj_plot_upsetr(top_te_list)
top_tes = unique(unlist(top_te_list))
length(top_tes)
te_heatmap_df = te_dr_df[, colnames(te_dr_df) %in% top_tes]
### Heatmap per cell type
#use newest version of this function: there is a bug in the ordering in the old one!
sum_mat = jj_summarize_sparse_mat(Matrix::t(te_heatmap_df), summarize_by_vec = te_dr_df$cluster_annotation_fine)
sum_mat_scaled = scale(t(sum_mat))
colnames(sum_mat_scaled)[!complete.cases(t(sum_mat_scaled))]
sum_mat_scaled = t(t(sum_mat_scaled)[complete.cases(t(sum_mat_scaled)), ])

col_list_use = list(Class = jj_get_colours(family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)],colour_csv = colcsv, comment_char = '$'),
                    Family = jj_get_colours(family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)], colour_csv = colcsv, comment_char = '$'))
column_ha = HeatmapAnnotation(Class = family_df$Class[match(colnames(sum_mat_scaled), family_df$TE)], 
                              Family = family_df$Family[match(colnames(sum_mat_scaled), family_df$TE)],
                              col = col_list_use)

ht = Heatmap(sum_mat_scaled,  name = 'scaled\nmean accessibility',
             show_column_names = T, top_annotation = column_ha, cluster_columns = T, cluster_rows = T)
ht = draw(ht)
celltype_order = rownames(sum_mat_scaled)[row_order(ht)]
te_order = colnames(sum_mat_scaled)[column_order(ht)]

dotplot_df = res_df %>% rownames_to_column('feature') %>% dplyr::filter(feature %in% top_tes) 
dotplot_df = dotplot_df %>% pivot_longer(cols = -1, names_to = 'name', values_to = 'value') %>% dplyr::mutate(info_type = ifelse(grepl('_FDR', name), 'FDR', 'logFC')) %>%
  dplyr::mutate(name = gsub('_FDR|_logFC', '', name)) %>% 
  spread(key = info_type, value = value)
dotplot_df = dotplot_df %>% dplyr::filter(!name %in% c('VAT','Skin', 'Colon'))
dotplot_df$logFDR <- (-1) * log10(dotplot_df$FDR)
dotplot_df$signif = dotplot_df$FDR < 0.05
dotplot_df$logFC_capped = jj_cap_vals(dotplot_df$logFC, cap_bottom = -1, cap_top = NULL)

dotplot_df$feature = factor(dotplot_df$feature, levels = te_order)
dotplot_df$name = factor(dotplot_df$name, levels = celltype_order)

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
pdf(paste0(storeFigPath, 'mouse_tnkilc_te_celltype_heatmap_top10_tes.pdf'), width = 10, height=7)
print(ht)
dev.off()
pdf(paste0(storeFigPath, 'mouse_tnkilc_te_celltype_dotplot_top10_tes.pdf'), width = 8 , height=8)
print(nebula_dotplot)
dev.off()

### tisTreg CD4 other overlap
res_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0)
#same as
res_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc0_fdr05.xlsx')
#jj_save_excel(res_list, paste0(storeFigPath, 'table_6_7_nebula_tnkilc_subset_da_transposable_elements.xlsx'))

da_te_list = lapply(res_list, '[[', 'feature')
shared = intersect(da_te_list$tisTregST2, da_te_list$CD4_other)
### tisTreg ILC2 overlap
shared2 = intersect(da_te_list$tisTregST2, da_te_list$ILC2)

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_df$odds_ratio = (olaps_df$both * olaps_df$none) / (olaps_df$te_only * olaps_df$tisTreg_sig_only)
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
olaps_sig_df$logfdr = -log10(olaps_sig_df$fdr)
olaps_sig_df$ilc2_overlap = olaps_sig_df$name %in% shared2
olaps_sig_df$th17_areg_overlap = olaps_sig_df$name %in% shared

upset_list = list(tisTreg_sig = olaps_sig_df$name, 
                  ILC2 =  shared2,
                  Th17_Areg = shared)

#S7J
pdf(paste0(storeFigPath, 'tisTreg_signature_TE_nebula_TE_overlap.pdf'), width = 10, height = 8)
jj_plot_upsetr(upset_list)
ggplot() + geom_point(data = olaps_sig_df, aes(x = odds_ratio, y = logfdr, colour = ilc2_overlap, shape = th17_areg_overlap), size = 3) + 
  theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_sig_df[olaps_sig_df$odds_ratio > 2 | olaps_sig_df$logfdr > 5, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)')
dev.off()
write_csv(olaps_sig_df, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-06-finalizing_te_analysis/mouse_cd4_tisTreg_sig_nebula_olaps.csv')

### plot selected TEs on umap
tes_quantify = c('RSINE1','RMER19C','RLTR48A', 'ERV3-16A3_LTR', 'MERV1_LTR', 'MLT2B2', 'MIR3')
gg = jj_plot_features(te_dr_df, features = tes_quantify, cap_top = 'auto', return_gg_object = T)
pdf(paste0(storeFigPath, 'umap_mouse_tnkilc_selected_tes.pdf'),  width = 8, height = 6)
gg
dev.off()

#16 TE overlap
writeLines(olaps_sig_df$name[olaps_sig_df$ilc2_overlap & olaps_sig_df$th17_areg_overlap],
           '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap/tisTreg_th17_areg_ilc2_te_overlap.txt')

# single te homer results dotplot -----------------------------------------

homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
binding_domains = binding_domains[!binding_domains == '?']
binding_domains[binding_domains=='T-box'] = 'T.box'
binding_domains = binding_domains[!binding_domains %in% c("Paired,Homeobox", "POU,Homeobox")]

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-25-browser_tracks/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_df$padj = p.adjust(homer_df$`P-value`, method = 'BH')
  homer_list[[i]] = homer_df[homer_df$padj < 0.05, ] #homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
sapply(homer_list, nrow)
# homer_motif_list = lapply(homer_list, '[[', 1)
# jj_plot_upsetr(homer_motif_list)

ggl = list()
for(k in seq_along(homer_list)){
  homer_plot_df = homer_list[[k]]
  homer_plot_df$enrichment = pmin(homer_plot_df$enrichment, 5)
  homer_plot_df$motif = sapply(strsplit(homer_plot_df$motif, split = '/'), '[[', 1)
  homer_plot_df$motif_short = sapply(strsplit(homer_plot_df$motif, split = '\\('), '[[', 1)
  homer_plot_df$log_padj = -log10(homer_plot_df$padj) #pmin(-log10(homer_plot_df$`q-value (Benjamini)`), 5)
  homer_plot_df$bindingDomain = 'Other'
  for(i in binding_domains){
    homer_plot_df$bindingDomain[grepl(i, homer_plot_df$motif)] = i  
  # jj_plot_volcano(homer_list[[1]], logfc_column = 'enrichment',
  #                 pval_column = 'padj', symbol_column = 'motif', labs_range = c(1, 5, 0, 5), marker_thres = c(1.5))
  text_df = head(homer_plot_df, 10) #homer_plot_df[homer_plot_df$`q-value (Benjamini)` < 0.05 & homer_plot_df$enrichment > 1.5, ]
  if( names(homer_list)[k] %in% c('mouse_tnkilc_RSINE1_736','mouse_tnkilc_RLTR48A_628', 'mouse_tnkilc_MERV1_LTR_614')){
    text_df = rbind(head(homer_plot_df,10), homer_plot_df[homer_plot_df$bindingDomain %in% c('bZIP') & homer_plot_df$enrichment > 4, ])
    text_df = text_df[!duplicated(text_df), ]
  }
  if( names(homer_list)[k] %in% c('mouse_tnkilc_RSINE1_736')){
    text_df = rbind(text_df, homer_plot_df[homer_plot_df$enrichment > 4 & homer_plot_df$log_padj > 30, ])
    text_df = text_df[!duplicated(text_df), ]
  }
  
  ggl[[ names(homer_list)[k] ]] = ggplot() + geom_point(data = homer_plot_df, aes(x=enrichment, y=log_padj, colour=bindingDomain), size = 3, alpha = 0.8) + 
    geom_text_repel(data = text_df, aes(x=enrichment, y=log_padj, label=motif_short), max.overlaps = 100) + 
    theme_minimal() + 
    scale_colour_manual(values = jj_get_colours(homer_plot_df$bindingDomain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
    labs(x = 'Enrichment', y = 'FDR (-log10)', colour = 'Binding domain', title = names(homer_list)[k] ) + 
    expand_limits(x=1, y=1) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 1)
  }
}

#7D,F
pdf(paste0(storeFigPath,'single_homer_known_motifs.pdf'),  width = 8, height =5)
ggl
dev.off()

# #TFs that are significant only in one subset
# diff_list = list()
# for(i in seq_along(homer_motif_list)){
#   diff_list[[i]] = setdiff(homer_motif_list[[i]], unlist(homer_motif_list[-i]))
# }
# names(diff_list) = names(homer_motif_list)

homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]

gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
#for subsets use
gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'mouse_tnkilc_ERV3-16A3_LTR_279'= 'ERV3-16A3_LTR (279)',
  'mouse_tnkilc_MERV1_LTR_614'= 'MERV1_LTR (614)',
  'mouse_tnkilc_MIR3_173'= 'MIR3 (173)',
  'mouse_tnkilc_MLT2B2_1351'= 'MLT2B2 (1351)',
  'mouse_tnkilc_RLTR48A_628'= 'RLTR48A (628)',
  'mouse_tnkilc_RMER19C_2722'= 'RMER19C (2722)',
  'mouse_tnkilc_RSINE1_736'= 'RSINE1 (736)'
))
gsea_res_comb$enrichment = jj_cap_vals(gsea_res_comb$enrichment, cap_top = 20)

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_single_te_homer_dotplot.pdf'), width =14, height=3)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'TE', y = '', fill = 'Enrichment') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()



# fig 7a treg te pseudotime trajectory ------------------------------------

library(ArchR)
library(tidyverse)

getTrajectory = function(mat, trajectory, matName = 'Matrix', groupEvery = 1,
                         log2Norm = TRUE, scaleTo = 10000, smoothWindow = 11, return_trajectory = FALSE) 
{
  #trajectory: DF with one column 'tisTreg_development' and barcodes as rownames
  trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
  breaks <- seq(0, 100, groupEvery)
  if (!all(is.numeric(trajectory[, 1]))) {
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if (!all(trajectory[, 1] >= 0 & trajectory[, 1] <= 100)) {
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }
  groupList <- lapply(seq_along(breaks), function(x) {
    if (x == 1) {
      NULL
    }
    else {
      rownames(trajectory)[which(trajectory[, 1] > breaks[x - 
                                                            1] & trajectory[, 1] <= breaks[x])]
    }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], 
                             "_", breaks[-1])
  
  trajectory$Group = NA
  for(k in seq_along(groupList)){
    trajectory$Group[rownames(trajectory) %in% groupList[[k]] ] = names(groupList)[k]
  }
  if(return_trajectory){
    return(trajectory)
  }
  
  # featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  # matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], 
  #                                    paste0(useMatrix, "/Info/Class")))
  # message("Creating Trajectory Group Matrix..")
  # groupMat <- .getGroupMatrix(ArrowFiles = getArrowFiles(ArchRProj), 
  #                             featureDF = featureDF, groupList = groupList, threads = threads, 
  #                             verbose = FALSE, useMatrix = useMatrix)
  mat = mat[, match(rownames(trajectory), colnames(mat))]
  stopifnot(identical(rownames(trajectory), colnames(mat)))
  groupMat = jj_summarize_sparse_mat(mat, summarize_by_vec = trajectory$Group, method = 'sum')
  cnames = unlist(lapply(strsplit(colnames(groupMat), '_|\\.'), '[[', 2))
  groupMat =  groupMat[, gtools::mixedorder(cnames)]
  
  if (!is.null(scaleTo)) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }
    else {
      groupMat <- t(t(groupMat)/colSums(groupMat)) * scaleTo
    }
  }
  if (log2Norm) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }
    else {
      groupMat <- log2(groupMat + 1)
    }
  }
  if (!is.null(smoothWindow)) {
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) ArchR:::.centerRollMean(x, 
                                                                                         k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))
    rownames(groupMat) = as.character(rownames(groupMat))
    rownames(smoothGroupMat) = as.character(rownames(smoothGroupMat))
    seTrajectory <- SummarizedExperiment(assays = SimpleList(smoothMat = as.matrix(smoothGroupMat), 
                                                             mat = as.matrix(groupMat)))#, rowData = featureDF)
    # if ("name" %in% colnames(featureDF)) {
    #   rownames(seTrajectory) <- paste0(featureDF$seqnames, 
    #                                    ":", featureDF$name)
    # }
    # else {
    #   rownames(seTrajectory) <- paste0(featureDF$seqnames, 
    #                                    ":", featureDF$start, "_", featureDF$end)
    # }
  }
  else {
    stop('smoothWindow must be set')
    # colnames(groupMat) <- paste0(colnames(groupMat))
    # seTrajectory <- SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)), 
    #                                      rowData = featureDF)
    # if ("name" %in% colnames(featureDF)) {
    #   rownames(seTrajectory) <- paste0(featureDF$seqnames, 
    #                                    ":", featureDF$name)
    # }
    # else {
    #   rownames(seTrajectory) <- paste0(featureDF$seqnames, 
    #                                    ":", featureDF$start, "_", featureDF$end)
    # }
  }
  metadata(seTrajectory)$Params <- list(useMatrix = matName, 
                                        matrixClass = "Sparse.Assays.Matrix",
                                        scaleTo = scaleTo, log2Norm = log2Norm, 
                                        smoothWindow = smoothWindow, date = Sys.Date())
  seTrajectory
}


plotTrajectoryHeatmap = function(seTrajectory = NULL, varCutOff = 0.9, maxFeatures = 25000, 
                                 scaleRows = TRUE, limits = c(-1.5, 1.5),# grepExclude = NULL, 
                                 pal = paletteContinuous(set = "solarExtra"), labelMarkers = NULL, 
                                 labelTop = 50, labelRows = FALSE, 
                                 rowOrder = NULL, #useSeqnames = NULL, 
                                 returnMatrix = FALSE, 
                                 force = FALSE) 
{
  # if (metadata(seTrajectory)$Params$matrixClass == "Sparse.Assays.Matrix") {
  #   if (is.null(useSeqnames) || length(useSeqnames) > 1) {
  #     .logMessage("useSeqnames is NULL or greater than 1 with a Sparse.Assays.Matrix trajectory input.", 
  #                 verbose = TRUE, logFile = logFile)
  #     if (force) {
  #       .logMessage("force=TRUE thus continuing", verbose = verbose, 
  #                   logFile = logFile)
  #     }
  #     else {
  #       useSeqnames <- rev(unique(rowData(seTrajectory)$seqnames))[1]
  #       .logMessage(paste0("force=FALSE thus continuing with subsetting useSeqnames = ", 
  #                          useSeqnames), verbose = TRUE, logFile = logFile)
  #     }
  #   }
  # }
  # if (!is.null(useSeqnames)) {
  #   seTrajectory <- seTrajectory[paste0(rowData(seTrajectory)$seqnames) %in% 
  #                                  paste0(useSeqnames), ]
  # }
  # if (nrow(seTrajectory) == 0) {
  #   .logStop("No features left in seTrajectory, please check input!", 
  #            logFile = logFile)
  # }
  mat <- assay(seTrajectory)
  # if (!is.null(grepExclude)) {
  #   idxExclude <- grep(grepExclude, rownames(mat))
  #   if (length(idxExclude) > 0) {
  #     mat <- mat[-grep(grepExclude, rownames(mat)), , drop = FALSE]
  #   }
  # }
  rSNA <- rowSums(is.na(mat))
  if (sum(rSNA > 0) > 0) {
    .logMessage("Removing rows with NA values...", verbose = TRUE, 
                logFile = logFile)
    mat <- mat[rSNA == 0, ]
  }
  #.logThis(mat, "mat-pre", logFile = logFile)
  varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))
  #.logThis(varQ, "varQ", logFile = logFile)
  orderedVar <- FALSE
  if (is.null(rowOrder)) {
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if (is.null(varCutOff) & is.null(maxFeatures)) {
      n <- nrow(mat)
    }else if (is.null(varCutOff)) {
      n <- maxFeatures
    }else if (is.null(maxFeatures)) {
      n <- (1 - varCutOff) * nrow(mat)
    } else {
      n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n), ]
  }
  
  
  #.logThis(mat, "mat-post", logFile = logFile)
  if(!is.null(labelTop)) {
    if (orderedVar) {
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    }
    else {
      idxLabel <- rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
    }
  }else {
    idxLabel <- NULL
  }
  #.logThis(idxLabel, "idxLabel", logFile = logFile)
  if (!is.null(labelMarkers)) {
    idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)),
                       nomatch = 0)
    idxLabel2 <- idxLabel2[idxLabel2 > 0]
  }
  else {
    idxLabel2 <- NULL
  }
  #.logThis(idxLabel2, "idxLabel2", logFile = logFile)
  idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])
  #.logThis(idxLabel, "idxLabel", logFile = logFile)
  
  # if (is.null(pal)) {
  #   if (is.null(metadata(seTrajectory)$Params$useMatrix)) {
  #     pal <- paletteContinuous(set = "solarExtra", n = 100)
  #   }
  #   else if (tolower(metadata(seTrajectory)$Params$useMatrix) == 
  #            "genescorematrix") {
  #     pal <- paletteContinuous(set = "blueYellow", n = 100)
  #   }
  #   else {
  #     pal <- paletteContinuous(set = "solarExtra", n = 100)
  #   }
  #}
  
  if (scaleRows) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
                 `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    #.logThis(mat, "mat-zscores", logFile = logFile)
  }
  if (nrow(mat) == 0) {
    stop("No Features Remaining!")
  }
  
  if(!is.null(rowOrder)) {
    idx <- rowOrder
  }else {
    if(anyNA(mat)){
      warning('replacing NA values in matrix with 0 because otherwise order function fails')
      mat[is.na(mat)] = 0
    }
    max_vals = apply(mat, 1, which.max)
    idx <- order(max_vals)
  }
  errl = list(mat = mat[idx, ], scale = FALSE, limits = c(min(mat), 
                                                          max(mat)), color = pal, clusterCols = FALSE, clusterRows = FALSE, 
              labelRows = labelRows, labelCols = FALSE, customRowLabel = match(idxLabel, 
                                                                               rownames(mat[idx, ])), showColDendrogram = TRUE, 
              name = metadata(seTrajectory)$Params$useMatrix, draw = FALSE)
  #.logThis(idx, "idx", logFile = logFile)
  ht <- tryCatch({
    ArchR:::.ArchRHeatmap(mat = mat[idx, ], scale = FALSE, limits = c(min(mat), max(mat)), 
                          color = pal, clusterCols = FALSE, clusterRows = FALSE, 
                          labelRows = labelRows, labelCols = FALSE, 
                          customRowLabel = match(idxLabel, rownames(mat[idx, ])), 
                          showColDendrogram = TRUE, 
                          name = metadata(seTrajectory)$Params$useMatrix, draw = T)
  }, error = function(e) {
    errorList = errl
    #.logError(e, fn = ".ArchRHeatmap", info = "", errorList = errorList, 
    #          logFile = logFile)
  })
  #.endLogging(logFile = logFile)
  if (returnMatrix) {
    return(errl)
  }
  else {
    return(ht)
  }
}


source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')
config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal_cd4)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj_use = loadArchRProject('mouse_cd4_treg_subset_monocle3_proj.RDS')
comb_norm_df = Matrix::t(prepare_comb_norm_df(pconfig, proj_use))

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_df$odds_ratio = (olaps_df$both * olaps_df$none) / (olaps_df$te_only * olaps_df$tisTreg_sig_only)
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
olaps_sig_df$logfdr = -log10(olaps_sig_df$fdr)
tes_label = olaps_sig_df$name[olaps_sig_df$odds_ratio > 2 | olaps_sig_df$logfdr > 5]
tes_keep = olaps_df$name[olaps_df$fdr < 0.05]

#tisTreg th17_areg ilc2 overlap (fig 7a)
tes_label = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap/tisTreg_th17_areg_ilc2_te_overlap.txt')


# trajectory_use = c('Treg_naive', 'tisTregST2_prog', 'tisTregST2')
# trajectory_use = c('C2', 'C4', 'C7', 'C9')
# proj_use <- addTrajectory(
#   ArchRProj = proj_use,
#   name = "tisTreg_development",
#   groupBy = "annotation",
#   #preFilterQuantile = 0.75,
#   #postFilterQuantile = 0.75,
#   #reducedDims = "IterativeLSI",
#   trajectory = trajectory_use,
#   embedding = 'UMAP', 
#   force = TRUE
# )
# p <- plotTrajectory(proj_use, trajectory = "tisTreg_development", colorBy = "cellColData", name = "signature_pct_overlap")
# p
# trajGSM <- ArchR::getTrajectory(ArchRProj = proj_use, name = "tisTreg_development", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
# p2 <- ArchR::trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), labelTop = 10)
#trajectory = proj_use@cellColData[, 'tisTreg_development', drop = FALSE]

trajectory = proj_use@cellColData[, c('pseudotime', 'annotation', 'Tissue'), drop = FALSE]
trajectory$pseudotime = 100*min_max_normalize(as.numeric(trajectory$pseudotime))
seTrajectory = getTrajectory(comb_norm_df,#[rownames(comb_norm_df) %in% tes_keep, ],
                             trajectory = trajectory, matName = 'Scaled TE accessibility',
                             groupEvery = 1#, 
                             #smoothWindow = 5
)
trajectory_annot_raw = getTrajectory(comb_norm_df,#[rownames(comb_norm_df) %in% tes_keep, ],
                                     trajectory = trajectory, matName = 'Scaled TE accessibility',
                                     groupEvery = 1, return_trajectory = T
)

# pdf(paste0(storeFigPath, 'treg_te_pseudotime_heatmap.pdf'), width = 8, height = 10)
# plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0, 
#                       scaleRows = T,
#                       labelMarkers = tes_keep, 
#                       labelTop = 0)
# plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0, 
#                       scaleRows = T,
#                       labelMarkers = tes_label, #tes_label, #c('RLTR48A', 'MERV1_LTR', 'B4A','RSINE1','ID_B1'), 
#                       labelTop = 0)
# dev.off()

trajectory_annot_raw$annotation = as.character(trajectory_annot_raw$annotation)
trajectory_annot_raw$Group = as.character(trajectory_annot_raw$Group)
trajectory_annot = trajectory_annot_raw %>% as.data.frame %>% dplyr::group_by(Group, annotation) %>% dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% dplyr::group_by(Group) %>% dplyr::mutate(total_n = sum(n)) %>% 
  dplyr::mutate(fraction = round(n / total_n, 4)) %>% 
  dplyr::ungroup() %>%
  complete(Group, annotation,
           fill = list(n = 0, fraction = 0))
trajectory_annot = trajectory_annot[!is.na(trajectory_annot$Group), ]

groupnames = unlist(lapply(strsplit(trajectory_annot$Group, '_|\\.'), '[[', 2))
trajectory_annot = trajectory_annot[gtools::mixedorder(groupnames), ]

# find_mode <- function(x) {
#   u <- unique(x)
#   tab <- tabulate(match(x, u))
#   u[tab == max(tab)]
# }
#trajectory_annot_tissue = trajectory_annot_raw %>% as.data.frame %>% dplyr::group_by(Group) %>% dplyr::summarise(main_tissue = find_mode(Tissue))

htmat = plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0, 
                              scaleRows = T,
                              labelMarkers = tes_keep, 
                              labelTop = 0, returnMatrix = T)

htmat2 = plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0, 
                               scaleRows = T,
                               labelMarkers = tes_label, 
                               labelTop = 0, returnMatrix = T)



trajectory_annot_tisTreg = as.data.frame(trajectory_annot[trajectory_annot$annotation == 'tisTregST2', ])
trajectory_annot_naiveTreg = as.data.frame(trajectory_annot[trajectory_annot$annotation == 'Treg_naive', ])
trajectory_annot_prec = as.data.frame(trajectory_annot[trajectory_annot$annotation == 'tisTregST2_prog', ])

column_col = colorRamp2(c(0, 1), c("white", "#E41A1C"))
column_col2 = colorRamp2(c(0, 1), c("white", "#4DAF4A"))
column_col3 = colorRamp2(c(0, 1), c("white", "#377EB8"))

# names(column_col) = trajectory_annot_tisTreg$fraction
# names(column_col2) = trajectory_annot_naiveTreg$fraction
# names(column_col3) = trajectory_annot_prec$fraction
ht1 = Heatmap(htmat$mat, cluster_columns = F, show_row_names = F, show_column_names = F, show_row_dend = F,
              col =  paletteContinuous(set = "horizonExtra"), cluster_rows = F, name = 'scaled\naccessibility',
              top_annotation = columnAnnotation(tisTreg_fraction = trajectory_annot_tisTreg$fraction,
                                                tisTreg_prog_fraction = trajectory_annot_prec$fraction,
                                                Treg_naive_fraction = trajectory_annot_naiveTreg$fraction,
                                                col = list(tisTreg_fraction =  column_col,
                                                           tisTreg_prog_fraction = column_col2,
                                                           Treg_naive_fraction = column_col3)))

ht2 = Heatmap(htmat2$mat, cluster_columns = F, show_row_names = F, show_column_names = F, show_row_dend = F,
              col =  paletteContinuous(set = "horizonExtra"), cluster_rows = F, name = 'scaled\naccessibility',
              top_annotation = columnAnnotation(tisTreg_fraction = trajectory_annot_tisTreg$fraction,
                                                tisTreg_prog_fraction = trajectory_annot_prec$fraction,
                                                Treg_naive_fraction = trajectory_annot_naiveTreg$fraction,
                                                col = list(tisTreg_fraction =  column_col,
                                                           tisTreg_prog_fraction = column_col2,
                                                           Treg_naive_fraction = column_col3)))

pdf(paste0(storeFigPath, 'mouse_cd4_treg_te_pseudotime_heatmap.pdf'), width = 6, height = 8)
ht1 + rowAnnotation(link =  ComplexHeatmap::anno_mark(at = htmat$customRowLabel, 
                                                      labels = rownames(htmat$mat)[htmat$customRowLabel], 
                                                      labels_gp = gpar(fontsize = 8)), 
                    width = unit(0.75, "cm") + max_text_width( htmat$customRowLabel))

ht2 + rowAnnotation(link =  ComplexHeatmap::anno_mark(at = htmat2$customRowLabel, 
                                                      labels = rownames(htmat2$mat)[htmat2$customRowLabel], 
                                                      labels_gp = gpar(fontsize = 8)), 
                    width = unit(0.75, "cm") + max_text_width( htmat2$customRowLabel))
dev.off()


##mouse atlas treg subset
# pconfig = yaml::read_yaml(config_list$mouse_normal)
# addArchRGenome(pconfig$GENOME)
# bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
# setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
# proj_treg = loadArchRProject('ArchR_Project_Treg_subset')
# proj_treg = proj_treg[!is.na(proj_treg$pseudotime), ]
# comb_norm_df = Matrix::t(prepare_comb_norm_df(pconfig, proj_treg))
# 
# trajectory = proj_treg@cellColData[, 'pseudotime', drop = FALSE]
# trajectory$pseudotime = 100*min_max_normalize(as.numeric(trajectory$pseudotime))
# seTrajectory = getTrajectory(comb_norm_df,#[rownames(comb_norm_df) %in% tes_keep, ],
#                              trajectory = trajectory, matName = 'Scaled TE accessibility',
#                              groupEvery = 1#,
#                              #smoothWindow = 5
# )
# 
# pdf(paste0(storeFigPath, 'mouse_atlas_treg_te_pseudotime_heatmap.pdf'), width = 8, height = 10)
# plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0,
#                       scaleRows = T,
#                       labelMarkers = tes_keep,
#                       labelTop = 0)
# plotTrajectoryHeatmap(seTrajectory, pal = paletteContinuous(set = "horizonExtra"), limits = c(-3,3), varCutOff = 0,
#                       scaleRows = T,
#                       labelMarkers = tes_label, #tes_label, #c('RLTR48A', 'MERV1_LTR', 'B4A','RSINE1','ID_B1'),
#                       labelTop = 0)
# dev.off()


# revisions ---------------------------------------------------------------


# functions ---------------------------------------------------------------

run_fisher = function(universe_gr, universe_subset_gr, query_gr){
  fisher_single = function(universe_gr, query_gr){
    universe_gr$query_pos = granges_overlap(universe_gr, query_gr, return_type = 'logical', olap_direction = 'a', minOverlap = 0)
    if(is.null(universe_gr$query_pos)){
      universe_gr$query_pos = FALSE
    }
    present_both = sum(universe_gr$subset_pos & universe_gr$query_pos)
    present_query_only = sum(!universe_gr$subset_pos & universe_gr$query_pos)
    present_subset_only = sum(universe_gr$subset_pos & !universe_gr$query_pos)
    present_none = sum(!universe_gr$subset_pos & !universe_gr$query_pos)
    
    cont_df = data.frame('query_yes' = c(present_both, present_query_only),
                         'query_no' = c(present_subset_only, present_none),
                         row.names = c('subset_yes', 'subset_no'))
    test <- fisher.test(cont_df)
    return(c(present_both, present_query_only, present_subset_only, present_none, test$p.value))
  }
  
  universe_gr$subset_pos = convert_granges(universe_gr) %in% convert_granges(universe_subset_gr)
  
  if(is.null(query_gr$name)){
    query_gr$name = 'query'
  }
  stopifnot(!anyNA(query_gr$name))
  
  query_n_df = as.data.frame(table(query_gr$name))
  colnames(query_n_df) = c('name', 'n_total_query')
  
  query_names = unique(query_gr$name)
  olaps_df = jj_initialize_df(5, length(query_names), init = 0,
                              col.names = c('both','query_only','subset_only','none', 'pval'),
                              row.names=query_names)
  for(i in 1:length(query_names)){
    if(i %% 50 == 0) message(i, '/', length(query_names))
    query_subset_gr = query_gr[query_gr$name == query_names[i]]
    olaps_df[i,] = fisher_single(universe_gr, query_subset_gr)
  }
  
  olaps_df$fraction_subset_olap = olaps_df$both / length(universe_subset_gr)
  olaps_df$fraction_inv_subset_olap = olaps_df$query_only / (length(universe_gr) - length(universe_subset_gr))
  #olaps_df$relative_enrichment = olaps_df$fraction_cons_sig_olap / olaps_df$fraction_other_peaks_olap #same as odds ratio
  #olaps_df = olaps_df[order(olaps_df$relative_enrichment, decreasing = T), ]
  olaps_df$odds_ratio = with(olaps_df, (both * none) / (subset_only * query_only))
  olaps_df$fdr = p.adjust(olaps_df$pval, method = 'fdr')
  olaps_df = olaps_df[, c('fdr','pval','odds_ratio','fraction_subset_olap','fraction_inv_subset_olap', 'both','query_only','subset_only','none')]
  olaps_df$n_subset_total = length(universe_subset_gr)
  olaps_df = olaps_df %>% 
    rownames_to_column('name') %>% 
    dplyr::left_join(query_n_df, by = 'name') %>% 
    dplyr::arrange(fdr)
  
  olaps_df
}


# revision part -----------------------------------------------------------




#########################################################################################################
# subsetting for Th1 signature similar to tisTreg signature ---------------
#########################################################################################################

####
# signature heatmap
####

source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

mouse_cd4_th1_sig = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/mouse_normal_cd4_th1_sig_lr_11_17_vs_4_7_8_npeaks_regressed.csv')
mouse_cd4_th1_sig_only = mouse_cd4_th1_sig[mouse_cd4_th1_sig$comparison == '12_17_versus_4_7_8', ]
mouse_cd4_th1_sig_gr = convert_granges(mouse_cd4_th1_sig_only$feature)

mouse_cd4_naive_sig_only = mouse_cd4_th1_sig[mouse_cd4_th1_sig$comparison == '4_7_8_versus_12_17', ]
mouse_cd4_naive_sig_gr = convert_granges(mouse_cd4_naive_sig_only$feature)


config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
proj = archr_add_peak_signatures(proj, signature_list = list(th1_vs_naive_sig = mouse_cd4_th1_sig_gr, naive_vs_th1_sig = mouse_cd4_naive_sig_gr), 'th1')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_th1_signature.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = c('z_th1_vs_naive_sig', 'z_naive_vs_th1_sig'))
dev.off()

#saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_mouse_cd4_th1_signature_peakset')
#proj = addPeakSet(proj, peakSet = mouse_cd4_th1_sig_gr, force=T)
#proj <- addPeakMatrix(proj,binarize = T)
#saveArchRProject(proj, outputDirectory = '/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_mouse_cd4_th1_signature_peakset')
proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_mouse_cd4_th1_signature_peakset')
pmat = get_peak_mat(proj)
cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
cluster_annotation_fine_df$cluster_annotation_fine[cluster_annotation_fine_df$cluster_annotation_fine == 'Th_Il21'] = 'Tfh_like'
proj$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
summarize_vec = proj$cluster_annotation_fine
mean_mat = jj_summarize_sparse_mat(pmat, summarize_vec)
mean_mat = mean_mat[, !colnames(mean_mat) %in% c('undefined')]
scaled_mat = t(scale(t(mean_mat)))
htmat = t(scaled_mat)
library(ComplexHeatmap)

ht = Heatmap(htmat, show_row_names=T, show_column_dend = F, show_row_dend = F, gap = unit(3,'mm'),
             show_column_names = F, column_split = 5, 
             name = 'scaled\naccessibility')
ht = draw(ht)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_th1_signature_subset_heatmap.pdf'), width =8, height=5)
ht
dev.off()


corder = column_order(ht)
subset_list = list(
  c1_myeloid =  colnames(htmat)[corder[[1]]],
  c2_Treg_Tfh_Th2_ILC2 =  colnames(htmat)[corder[[2]]],
  c3_mixed =  colnames(htmat)[corder[[3]]],
  c4_ILC_Th17 =  colnames(htmat)[corder[[4]]],
  c5_CD8mem_eff_NK_ILC1 = colnames(htmat)[corder[[5]]]
)

######
## homer results dotplot
######

# for(i in names(subset_list)){
#    jj_save_bed(convert_granges(subset_list[[i]]), file_name = paste0(storeFigPath, i, '.bed'))
# }

#for the subsets use:
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-01-th1_subsets/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_list[[i]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
sapply(homer_list, nrow)
homer_motif_list = lapply(homer_list, '[[', 1)
jj_plot_upsetr(homer_motif_list)


homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]

gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)
#for subsets use
gsea_res_comb$name = from_to(vec= gsea_res_comb$name, old_new_map_vec=c(
  'c1_myeloid'= '1',
  'c2_Treg_Tfh_Th2_ILC2'= '2',
  'c3_mixed'= '3',
  'c4_ILC_Th17'= '4',
  'c5_CD8mem_eff_NK_ILC1'= '5'
))

# Mouse atlas tisTreg signature subset homer known motif heatmap (and S7H)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_th1_signature_complete_homer.pdf'), width =5, height=8)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'Signature subset', y = '', fill = 'Enrichment')
dev.off()

########
#signature percentage overlap
########

config_list = get_config()
dataset_use = 'mouse_normal'
pconfig = yaml::read_yaml(config_list[[dataset_use]])
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df,
                                 nFrags_vec = proj$nFrags, signature_gr = mouse_cd4_th1_sig_gr,
                                 verbose = T, k = 100, count_thres = 2e5)
#write_rds(olap_df, paste0(storeFigPath, 'mouse_atlas_scATAC_th1_sig_olap_df.RDS'))
olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-19-footprints/mouse_atlas_scATAC_th1_sig_olap_df.RDS')

#plot_signature_olap(dr_df, olap_df)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
cluster_annotation_fine_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-01-26-treg_subsets/treg_fine_annotation.csv')
cluster_annotation_fine_df$cluster_annotation_fine[cluster_annotation_fine_df$cluster_annotation_fine == 'Th_Il21'] = 'Tfh_like'
cluster_annotation_fine_df$cluster_annotation_fine[cluster_annotation_fine_df$cluster_annotation_fine == 'CD4_other'] = 'Th17_Areg'

dr_df$cluster_annotation_fine = cluster_annotation_fine_df$cluster_annotation_fine
dr_df$pct_overlap = olap_df$signature_pct_overlap

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_th1_signature_overlap_umap.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = 'pct_overlap', return_gg_object = T)[[1]] + labs(colour = '% overlap') 
dev.off()
# celltypes_keep = dr_df$cluster_annotation_fine %>% table %>% .[. > 200] %>% names
# dr_df = dr_df[dr_df$cluster_annotation_fine %in% celltypes_keep, ]

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_th1_signature_overlap_boxplot.pdf'),  width = 6, height = 5)
jj_plot_numeric_by_group(dr_df[!dr_df$cluster_annotation_fine == 'undefined', ], 'pct_overlap', group_column = 'cluster_annotation_fine', 
                         order = T, flip_coordinates = T, type = 'boxplot',
                         custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  theme(legend.position = 'none') + labs(y='% overlap', x='Cell type')
dev.off()


#########################################################################################################
#########################################################################################################
# common and distinct features of tissue adaptation -----------------------
#########################################################################################################
#########################################################################################################

#########################################################################################################
#get marker peaks peripheral vs spleen for other immune cell subsets than T cells
#########################################################################################################

library(ArchR)
library(tidyverse)
library(jj)
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

proj$peripheral_tissue = ifelse(proj$Tissue == 'Spleen', 'spleen',  'peripheral')
proj$cluster_annotation_level00 = from_to(proj$cluster_annotation_level0, old_new_map_vec = c(Macrophage = "Macrophage/Monocyte", Monocyte = "Macrophage/Monocyte"))
proj$peripheral_cell_type = paste(proj$peripheral_tissue, proj$cluster_annotation_level00, sep='__')

dr_df = jj_get_reduction_coords(proj, redname='UMAP')

comparisons_keep = dr_df %>% 
  dplyr::group_by(cluster_annotation_level00) %>%
  dplyr::count(peripheral_tissue) %>%
  tidyr::spread(peripheral_tissue, n) %>%
  dplyr::filter(peripheral > 200 & spleen > 200 & !cluster_annotation_level00 %in%  c('undefined')) %>%
  dplyr::pull(cluster_annotation_level00)


markers_list = list()
for(i in unique(proj$cluster_annotation_level00)){
  message(i)
  if(!i %in%  comparisons_keep) next
  markers_list[[i]] <- getMarkerFeatures(
    ArchRProj = proj[proj$cluster_annotation_level00 == i, ],
    useMatrix = 'PeakMatrix',
    useGroups = 'peripheral',
    groupBy = "peripheral_tissue",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
}
head(assays(markers_list$`T cell`)[[1]])

#combine the markers from the individual comparisons in one assay, which is put in a summarized experiment
new_ass_list = list()
for(i in seq_along(assays(markers_list[[1]]))){
  new_name = names(assays(markers_list[[1]]))[i]
  ass_list = list()
  for(j in seq_along(markers_list)){
    ass_list[[j]] = assays(markers_list[[j]])[[i]]
  }
  ass_df = Reduce(cbind, ass_list)
  colnames(ass_df) = names(markers_list)
  new_ass_list[[new_name]] = ass_df[, colnames(ass_df) %in% comparisons_keep]
}
own_se = SummarizedExperiment(assays = new_ass_list,
                              metadata = list(Params = list(useMatrix ='PeakMatrix')), #metadata(markers_use_se),
                              rowData = rowData(markers_list[[1]]))

#write_rds(own_se, paste0(bigFilesDir, 'mouse_atlas_dc_macmono_b_t_nk_peripheral_vs_spleen_marker_peaks_se.RDS'))
own_se = read_rds(paste0(bigFilesDir, 'mouse_atlas_dc_macmono_b_t_nk_peripheral_vs_spleen_marker_peaks_se.RDS'))

marker_df = archr_get_markers_as_df(own_se, proj,  cutOff = "FDR <= 0.01 & Log2FC >= 0.58")
table(marker_df$comparison)
# comparisons_keep = comparisons_keep[comparisons_keep %in% names(which(table(marker_df$comparison) > 100))]


marker_list = split(marker_df$feature, marker_df$comparison)
pdf(paste0(storeFigPath, 'immune_atlas_archr_peripheral_vs_spleen_per_celltype_upset.pdf'),  width = 8, height = 6)
jj_plot_upsetr(marker_list)
dev.off()

#########################################################################################################
# Run Nebula on broader cell type annotation of entire dataset and get peripheral specific peaks --------
#########################################################################################################

### R4.2
#pgr_use = getPeakSet(proj)
#proj = addPeakSet(proj, peakSet = pgr_use, force=T)
#check current peak set with proj@peakSet
#the variance is 2.5 * mean
proj = loadArchRProject('/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchR_proj_nonbinary_peak_mat')
pmat = get_peak_mat(proj, binarize = FALSE)
pmat = ceiling(pmat / 2) #fixes problem of equal counts being more abundant, but still var = 1.25*mean
max(pmat) #25

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

proj$peripheral_tissue = ifelse(proj$Tissue == 'Spleen', 'spleen',  'peripheral')

###
dr_df = jj_get_reduction_coords(proj)
#dr_df$peripheral_tissue = factor(ifelse(dr_df$Tissue == 'Spleen', 'spleen',  'peripheral'), levels = c('spleen', 'peripheral'))
ref_celltype = 'B cell'
ref_tissue = 'Spleen'
dr_df$cluster_annotation_level0 = relevel(factor(dr_df$cluster_annotation_level0), ref = ref_celltype)
dr_df$Tissue = relevel(factor(dr_df$Tissue), ref = ref_tissue)

#dr_df$peripheral_tissue = factor(dr_df$peripheral_tissue, levels = c('spleen', 'peripheral'))
#df = model.matrix(~peripheral_tissue+cluster_annotation_level0, data = dr_df) 
df = model.matrix(~Tissue+cluster_annotation_level0, data = dr_df) 

re = nebula(pmat, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005, ncore =1, verbose = T)  
#saveRDS(re, paste0(bigFilesDir, 'nebula_immune_atlas_cell_type_coarse_peripheral_spleen.RDS'))

re = read_rds(paste0(bigFilesDir, 'nebula_immune_atlas_cell_type_coarse_peripheral_spleen.RDS'))
res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_level0', '', colnames(res_df))
colnames(res_df) = gsub('Tissue', '', colnames(res_df))

sig_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0.5)
sapply(sig_list,nrow)
celltype_da_peaks = unique(unlist(sapply(sig_list[!names(sig_list) %in% c('Colon','Skin','VAT')], '[[', 'feature')))
for(i in c("Colon", "VAT",  "Skin")){
  sig_list[[i]]$tissue_exclusive = !sig_list[[i]]$feature %in% celltype_da_peaks
}
sig_list = lapply(sig_list, function(x) x[order(x$FDR), ])
#jj_save_excel(sig_list, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc05_fdr05.xlsx') 
sig_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc05_fdr05.xlsx')

peripheral_sigs = lapply(sig_list, '[[', 'feature')

pdf(paste0(storeFigPath, 'immune_atlas_nebula_upset_ref_b_spleen_0.5.pdf'),  width = 10, height = 8)
jj_plot_upsetr(peripheral_sigs)
dev.off()

pdf(paste0(storeFigPath, 'immune_atlas_nebula_upset_only_tissue_ref_b_spleen_0.5.pdf'),  width = 8, height = 6)
jj_plot_upsetr(peripheral_sigs[c('Colon', 'Skin','VAT')])
dev.off()


# peripheral_sigs = c(list(
#   nebula_peripheral = peripheral_conserved$feature,
#   nebula_peripheral_exclusive = peripheral_conserved$feature[peripheral_conserved$exclusive]
# ), marker_list)

library(genomic_region_tools)
peripheral_sigs = lapply(peripheral_sigs, convert_granges)

proj = archr_add_peak_signatures(proj, signature_list = peripheral_sigs[c('Colon', 'Skin','VAT')], signature_name = 'nebula_tissue_')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
colnames(dr_df) = make.names(colnames(dr_df))

pdf(paste0(storeFigPath, 'immune_atlas_nebula_signatures.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = paste0('z_', c('Colon', 'Skin','VAT')), cap_top = 'q99', cap_bottom = 'q01')
dev.off()

##############
### plot enriched TFs in the tissue marker peaks
##############

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-30-homer_nebula/'
#save subset regions as bed
#sapply(seq_along(peripheral_sigs), function(x) jj_save_bed(peripheral_sigs[[x]], paste0(homer_res, gsub('/| ', '_', names(peripheral_sigs)[x]), '.bed')))
#exclusive sigs
tissue_sigs = lapply(sig_list[c('Colon','Skin','VAT')], '[[', 'feature')
tissue_excl_sigs = list()
for(i in c('Colon','Skin','VAT')){
  tissue_sig = tissue_sigs[[i]]
  other_sigs = unlist(tissue_sigs[!names(tissue_sigs) %in% i])
  tissue_sig = tissue_sig[!tissue_sig %in% other_sigs]
  message(i, ': ', length(tissue_sig))
  #jj_save_bed(tissue_sig , paste0(homer_res, i, '_exclusive.bed')) 
  tissue_excl_sigs[[i]] = tissue_sig
}

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap/' #tissue-associated peaks, but with shared fractions (fdr 0.05, logfc 0.5)
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-30-homer_nebula/' #tissue-associated peaks, which are not associated with another tissue (fdr 0.05, logfc 0.5)
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_list[[i]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 1.5, ]
}
sapply(homer_list, nrow)
homer_motif_list = lapply(homer_list, '[[', 1)
jj_plot_upsetr(homer_motif_list)

# #TFs that are significant only in one subset
# diff_list = list()
# for(i in seq_along(homer_motif_list)){
#   diff_list[[i]] = setdiff(homer_motif_list[[i]], unlist(homer_motif_list[-i]))
# }
# names(diff_list) = names(homer_motif_list)

homer_list_gsea_plot = lapply(homer_list, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]

# res = as.matrix(pivot_wider(gsea_res_comb, id_cols = motif, names_from =  name, values_from = enrichment, values_fill = NA))
# rnames = res[, 1]
# res = res[, -1]
# res = apply(res, 2, as.numeric)
# rownames(res) = rnames
# Heatmap(res)

gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)

# Nebula tissue peaks homer known motif heatmap
pdf(paste0(storeFigPath, 'Peripheral_tissue_nebula_peaks_homer_enrichment.pdf'), width =5, height=8)
pdf(paste0(storeFigPath, 'Peripheral_tissue_exclusive_nebula_peaks_homer_enrichment.pdf'), width =5, height=8)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'Associated peaks (Nebula)', y = '', fill = 'Enrichment') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##############
### go enrichment (results are very broad)
##############

# library(clusterProfiler)
go_files = c(
  hallmark="/omics/groups/OE0436/internal/msimon/scATAC/mh.all.v2022.1.Mm.symbols.gmt",
  curated_sets='/omics/groups/OE0436/internal/msimon/scATAC/m2.all.v2022.1.Mm.symbols.gmt',
  go_sets='/omics/groups/OE0436/internal/msimon/scATAC/m5.all.v2022.1.Mm.symbols.gmt'
  #celltype_signatures='/omics/groups/OE0436/internal/msimon/scATAC/m8.all.v2022.1.Mm.symbols.gmt'
)
go_list <- lapply(go_files, clusterProfiler::read.gmt)

#complet immune atlas nebula result (naming of file wrong)
sig_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc05_fdr05.xlsx')

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))


##nebula
closest_genes = list()
pset = getPeakSet(proj)
pset_peaks = convert_granges(pset)

for(i in c('Colon', 'Skin','VAT')){
  temp_set = sig_list[[i]]
  temp_set = temp_set %>% dplyr::filter(FDR < 0.001, logFC > 1.5) %>% pull(feature)
  closest_genes[[i]] = unique(na.omit(pset$nearestGene[pset_peaks %in% temp_set]))
}

# ##archr
# #taken from next section "comparison between archr and nebula"
# closest_genes = list()
# for(i in names(marker_list)){
#   temp_set = marker_list[[i]]
#   closest_genes[[i]] = unique(na.omit(pset$nearestGene[pset_peaks %in% temp_set]))
# }

sapply(closest_genes, length)

go_res_list = list()
for(i in seq_along(closest_genes)){
  message(i, '/', length(closest_genes))
  go_res_list[[i]] = lapply(go_list, function(x) {
    enricher(
      gene = closest_genes[[i]],
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = x)
  }
  )
}
names(go_res_list) = names(closest_genes)

qvalue_thres = 0.05
#start with 0 and increase
min_gene_ratio = 0.05


dp_list = list()
for(j in names(go_list)){
  hallmark_list = lapply(go_res_list, function(x) x[[j]]@result)
  for(i in seq_along(hallmark_list)){
    halltemp = hallmark_list[[i]]
    halltemp$pathway = halltemp$ID
    halltemp$GeneRatio = sapply(strsplit(halltemp$GeneRatio, '/'), function(x) as.integer(x[1]) / as.integer(x[2]))
    halltemp = halltemp[halltemp$qvalue < qvalue_thres, ]
    halltemp = halltemp[halltemp$GeneRatio > min_gene_ratio, ]
    hallmark_list[[i]] = halltemp
  }
  hallmark_list = hallmark_list[sapply(hallmark_list, nrow) > 0]
  if(length(hallmark_list) > 0){
    
    gsea_res_w_id <- lapply(1:length(hallmark_list), FUN = function(x) {
      dataset <- hallmark_list[[x]]
      dataset$id <- names(hallmark_list)[x]
      dataset$pval_na <- sapply(dataset$qvalue, FUN = function(y) { is.na(y) })
      dataset$logq <- (-1) * log10(dataset$qvalue)
      return(dataset)
    }) 
    
    # Combine the cell-type-specific data sets.
    gsea_res_comb <- do.call(rbind, gsea_res_w_id)
    gsea_res_comb$id <- factor(gsea_res_comb$id, levels = names(hallmark_list))
    gsea_res_comb$pathway <- factor(gsea_res_comb$pathway, levels = rev(sort(unique(gsea_res_comb$pathway))))
    
    # Generate the dot plot.
    dot_plot <- ggplot(gsea_res_comb) +
      aes(x = id, y = pathway, colour = GeneRatio, size = logq, shape = pval_na) +
      viridis::scale_color_viridis() + 
      #scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(0, max(gsea_res_comb$GeneRatio)), midpoint = max(gsea_res_comb$GeneRatio)/2, na.value = "black") +
      scale_shape_manual(breaks = c(T, F), values = c(4, 16)) +
      scale_size_area(na.value = 3) +
      geom_point(stroke = 2) +
      guides(colour = guide_colourbar(order = 1),
             size = guide_legend(title = "-log10(q)", order = 2),
             shape = "none") +
      xlab("Group") +
      ylab("Gene set") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"))
    
    dp_list[[j]] = dot_plot
  }
}

#gene ratio 0.05, qvalue 0.05, nebula fdr 0.001, logFC > 1.5
pdf(paste0(storeFigPath, 'nebula_tissue_signatures_go_enrichment_dotplots.pdf'),width = 8, height = 8)
dp_list
dev.off()

##############
###compare with archr markers
##############

#nebula
sig_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc05_fdr05.xlsx')
tissue_sigs = lapply(sig_list[c('Colon','Skin','VAT')], '[[', 'feature')
nebula_peripheral_peaks = unique(unlist(tissue_sigs))

#archr
own_se = read_rds(paste0(bigFilesDir, 'mouse_atlas_dc_macmono_b_t_nk_peripheral_vs_spleen_marker_peaks_se.RDS'))
marker_df = archr_get_markers_as_df(own_se, proj,  cutOff = "FDR <= 0.01 & Log2FC >= 0.58")
marker_list = split(marker_df$feature, marker_df$comparison)
archr_peripheral_peaks = unique(unlist(marker_list))

pdf(paste0(storeFigPath, 'nebula_archr_peripheral_peaks_comparison_upset.pdf'),  width = 8, height = 5)
jj_plot_upsetr(list(nebula = nebula_peripheral_peaks, archr=archr_peripheral_peaks))
dev.off()
#length(pset) 254545

###############
###test TE enrichment within the tissue signatures
###############

te_gr = read_te('mm10', return_gr = T)
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)

te_gr$family = family_df$Family[match(te_gr$name, family_df$TE)]

##use this to compare on the family level
#te_gr$name = te_gr$family


pset = getPeakSet(proj)

tissue_sigs_gr = lapply(tissue_sigs, convert_granges)

olaps_df_list = list()
for(i in seq_along(tissue_sigs_gr)){
  message(i)
  olaps_df_list[[i]] = run_fisher(pset, tissue_sigs_gr[[i]], te_gr)
}
names(olaps_df_list) = names(tissue_sigs_gr)

olaps_df = do.call(rbind, olaps_df_list)
olaps_df$Tissue = rep(names(olaps_df_list), sapply(olaps_df_list, nrow))

#write_csv(olaps_df, paste0(storeFigPath, 'peripheral_tissue_te_family_enrichment.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'peripheral_tissue_te_enrichment.csv'))
# olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/peripheral_tissue_te_family_enrichment.csv')
# olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/peripheral_tissue_te_enrichment.csv')

olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
#olaps_plot_df = olaps_df

#olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
#olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
olaps_plot_df$logfdr[olaps_plot_df$logfdr > 15] = 15
olaps_plot_df$Class = family_df$Class[match(olaps_plot_df$name, family_df$Family)]

#tissue signatures TE Family enrichment volcano plot
pdf(paste0(storeFigPath, 'nebula_tissue_sig_te_family_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = Tissue, shape = Class), size = 4) + 
  geom_line(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, group = name), colour = 'grey') + 
  scale_colour_manual(values = jj_get_colours(olaps_plot_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  #viridis::scale_color_viridis() +  
  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)')

dev.off()


olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
#olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
olaps_plot_df$logfdr[olaps_plot_df$logfdr > 15] = 15

pdf(paste0(storeFigPath, 'nebula_tissue_sig_te_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = Tissue), size = 4) + 
  geom_line(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, group = name), colour = 'grey') + 
  scale_colour_manual(values = jj_get_colours(olaps_plot_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv')) + 
  #viridis::scale_color_viridis() +  
  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)')


dev.off()


############################################################################
# human vs mouse tisTreg sig ----------------------------------------------
############################################################################

library(jj)
library(GenomicRanges)
library(tidyverse)
library(genomic_region_tools)
library(Seurat)
#Compare 4,411 human tissue Treg peaks to 7,500 murine tissue Treg peaks with regard to TE subfamily enrichment .

### mouse (already included in fig S7f)

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
#seurat_atac = read_rds(pick_content('mouse_normal_CD4', 'seurat_file'))
#mouse_cd4_peakset = rownames(GetAssayData(seurat_atac, assay = 'scATAC_raw'))
mouse_cd4_peakset = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-09-peripheral_sig/mouse_cd4_peakset.txt')

mouse_cd4_universe_gr = convert_granges(mouse_cd4_peakset)

te_gr = read_te('mm10', return_gr = T)

mouse_cd4_universe_gr$tisTreg_sig_olap = granges_overlap(mouse_cd4_universe_gr, mouse_cd4_tisTreg_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)

### family level comparison
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)

te_gr$family = family_df$Family[match(te_gr$name, family_df$TE)]
te_gr$name = te_gr$family

olaps_df = run_fisher(mouse_cd4_universe_gr, mouse_cd4_tisTreg_gr, te_gr)

#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv'))
olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')

olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
olaps_plot_df$Class = family_df$Class[match(olaps_plot_df$name, family_df$Family)]

#S7G tisTreg TE Family enrichment volcano plot
pdf(paste0(storeFigPath, 'tisTreg_sig_te_family_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

ggplot() + geom_point(data = olaps_plot_df[olaps_plot_df$both > 10, ], aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$both > 10, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

dev.off()


### human


seurat_atac = read_rds(pick_content('human_CD4', 'seurat_file'))
human_cd4_universe_gr = convert_granges(rownames(GetAssayData(seurat_atac, assay = 'scATAC_raw')))
rm(seurat_atac)

human_tisTreg_sig = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-02-08-cluster_comparison_3_vs_7_human_cd4/diff_results_7_versus_3_seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.csv')
human_tisTreg_sig = human_tisTreg_sig[human_tisTreg_sig$comparison == '3_versus_7', ]
human_cd4_tisTreg_gr = convert_granges(human_tisTreg_sig$feature)

te_gr = read_te('hg19', return_gr = T)

human_cd4_universe_gr$tisTreg_sig_olap = granges_overlap(human_cd4_universe_gr, human_cd4_tisTreg_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)

### family level comparison
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/hg19.te.bed', col_names = F)
family_df = te_bed_annot[, c('X4', 'X5', 'X6')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)

te_gr$family = family_df$Family[match(te_gr$name, family_df$TE)]

#for family level comparison use:
te_gr$name = te_gr$family
#otherwise leave the name column as is

###

olaps_df = run_fisher(human_cd4_universe_gr, human_cd4_tisTreg_gr, te_gr)

#write_csv(olaps_df, paste0(storeFigPath, 'human_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'human_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv'))
olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/human_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')

olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
olaps_plot_df$Class = family_df$Class[match(olaps_plot_df$name, family_df$Family)]

#human tisTreg TE Family enrichment volcano plot
pdf(paste0(storeFigPath, 'human_tisTreg_sig_te_family_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

ggplot() + geom_point(data = olaps_plot_df[olaps_plot_df$both > 10, ], aes(x = odds_ratio, y = logfdr, colour = log10(both), shape = Class), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$both > 10, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

dev.off()

#######
## human individual TEs
#######

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-01-th1_subsets/human_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')

olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 5] = 5
olaps_plot_df = olaps_plot_df[olaps_plot_df$both >= 10, ]
olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)

te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/hg19.te.bed', col_names = F)
family_df = te_bed_annot[, c('X4', 'X5', 'X6')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)
family_col_df = family_df %>% dplyr::select(Class, Family) %>% .[!duplicated(.), ] %>% dplyr::arrange(Class, Family)
cols_use_family = jj_get_colours(family_col_df$Family, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

olaps_plot_df = olaps_plot_df %>% dplyr::left_join(family_df, by = c('name' = 'TE'))

#human tisTreg TE enrichment volcano plot
pdf(paste0(storeFigPath, 'human_tisTreg_sig_te_enrichment_fisher_volcano_plot.pdf'), width = 10, height = 8)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both)), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df, 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both)), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$logfdr > 10 | olaps_plot_df$odds_ratio > 3, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)')

ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = Family, shape = Class), size = 3) + 
  scale_colour_manual(values = cols_use_family) + theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$odds_ratio > 3 | olaps_plot_df$logfdr > 10, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + labs(x = 'Odds ratio', y = '-log10(FDR)')

dev.off()


#mouse individual TEs
olaps_df_mouse = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_df_mouse$odds_ratio = (olaps_df_mouse$both * olaps_df_mouse$none) / (olaps_df_mouse$te_only * olaps_df_mouse$tisTreg_sig_only)
olaps_sig_df = olaps_df_mouse[olaps_df_mouse$fdr < 0.05, ] 
olaps_sig_df$logfdr = -log10(olaps_sig_df$fdr)
mouse_family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
olaps_sig_df = olaps_sig_df %>% dplyr::left_join(mouse_family_df, by = c('name' = 'TE'))

family_col_df = family_df %>% dplyr::select(Class, Family) %>% .[!duplicated(.), ] %>% dplyr::arrange(Class, Family)
cols_use_family = jj_get_colours(family_col_df$Family, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

#6G Fisher test enrichment of TE insertions in the tisTreg signature
pdf(paste0(storeFigPath, 'tisTreg_sig_te_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)

ggplot() + geom_point(data = olaps_sig_df, aes(x = odds_ratio, y = logfdr, colour = Family, shape = Class), size = 3) + 
  scale_colour_manual(values = cols_use_family) + theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_sig_df[olaps_sig_df$odds_ratio > 2 | olaps_sig_df$logfdr > 5, ], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + labs(x = 'Odds ratio', y = '-log10(FDR)')

dev.off()



##comparison
# olaps_df_mouse = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')
# olaps_df_mouse = olaps_df_mouse[olaps_df_mouse$fdr < 0.05, ]
# olaps_df_mouse$odds_ratio[olaps_df_mouse$odds_ratio > 5] = 5
# olaps_df_mouse = olaps_df_mouse[olaps_df_mouse$both >= 10, ]
# olaps_df_mouse$logfdr = -log10(olaps_df_mouse$fdr)
# olaps_df_mouse$Class = family_df$Class[match(olaps_df_mouse$name, family_df$Family)]
# 
# olaps_df_human = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/human_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')
# olaps_df_human = olaps_df_human[olaps_df_human$fdr < 0.05, ]
# olaps_df_human$odds_ratio[olaps_df_human$odds_ratio > 5] = 5
# olaps_df_human = olaps_df_human[olaps_df_human$both >= 10, ]
# olaps_df_human$logfdr = -log10(olaps_df_human$fdr)
# olaps_df_human$Class = family_df$Class[match(olaps_df_human$name, family_df$Family)]
# 
# olaps_plot_df = olaps_df_mouse %>% dplyr::inner_join(olaps_df_human, by = 'name', suffix = c('.mouse', '.human'))
#
# pdf(paste0(storeFigPath, 'mouse_human_tisTreg_sig_te_comparison.pdf'), width = 8, height = 6)
# 
# ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio.mouse, y = odds_ratio.human, shape = Class.mouse), size = 4) + 
#   #viridis::scale_color_viridis() +  
#   theme_minimal() + 
#   #scale_colour_manual(values = cols_use_family) +
#   geom_hline(yintercept = 1) + geom_vline(xintercept = 1) + 
#   ggrepel::geom_text_repel(data=olaps_plot_df, 
#                            aes(x=odds_ratio.mouse, y=odds_ratio.human, label=name), max.overlaps = 500) + 
#   labs(x = 'Odds ratio mouse', y = 'Odds ratio human')
# 
# ggplot() + geom_point(data = olaps_plot_df, aes(x = logfdr.mouse, y = logfdr.human, shape = Class.mouse), size = 4) + 
#   #viridis::scale_color_viridis() +  
#   theme_minimal() + 
#   #scale_colour_manual(values = cols_use_family) +
#   ggrepel::geom_text_repel(data=olaps_plot_df, 
#                            aes(x=logfdr.mouse, y=logfdr.human, label=name), max.overlaps = 500) + 
#   labs(x = 'FDR mouse', y = 'FDR human')
# dev.off()

####
# Mouse vs human TE enrichment heatmap
####
olaps_df_mouse = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/mouse_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')
olaps_df_mouse = olaps_df_mouse[olaps_df_mouse$both >= 10, ]
olaps_df_human = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-12-16-mouse_human_tisTreg_te_comparison/human_cd4_tisTreg_sig_vs_universe_te_family_overlap_df.csv')
olaps_df_human = olaps_df_human[olaps_df_human$both >= 10, ]
olaps_plot_df = olaps_df_mouse %>% dplyr::full_join(olaps_df_human, by = 'name', suffix = c('.mouse', '.human'))

#fdr_mat = olaps_plot_df %>% dplyr::select(name, fdr.mouse, fdr.human) %>% as.data.frame
# library(ggrepel)
# fdr_mat[is.na(fdr_mat)] = 1
# fdr_mat$te_occurence = 'both'
# fdr_mat$te_occurence[fdr_mat$fdr.mouse == 1] = 'human'
# fdr_mat$te_occurence[fdr_mat$fdr.human == 1] = 'mouse'
# ggplot(fdr_mat, aes(x=-log10(fdr.mouse), y = -log10(fdr.human))) + geom_point(aes(colour=te_occurence)) + geom_text_repel(aes(label = name))

odds_mat = olaps_plot_df %>% dplyr::select(name, odds_ratio.mouse, odds_ratio.human) %>%
  dplyr::arrange(name) %>% as.data.frame %>% column_to_rownames('name') %>% as.matrix
fdr_mat = olaps_plot_df %>% dplyr::select(name, fdr.mouse, fdr.human) %>%
  dplyr::arrange(name) %>% as.data.frame %>% column_to_rownames('name') %>% as.matrix

# ComplexHeatmap::Heatmap(odds_mat, cluster_rows = F, cluster_columns = F, name = 'Odds ratio', 
#                         cell_fun = function(j, i, x, y, w, h, fill) {grid.text(round(fdr_mat[i,j],5), x, y)})

fdr_mat[is.na(fdr_mat)] = 2
colnames(odds_mat) = gsub('odds_ratio.', '', colnames(odds_mat))
ComplexHeatmap::Heatmap(odds_mat, cluster_rows = F, cluster_columns = F, name = 'Odds ratio', 
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(fdr_mat[i, j] < 0.001) {
                            grid.text("***", x, y)
                          }else if(fdr_mat[i, j] < 0.01) {
                            grid.text("**", x, y)
                          }else if(fdr_mat[i, j] < 0.05) {
                            grid.text("*", x, y)
                          }})


###########################################################################
# signature percentage overlap sampling -----------------------------------
###########################################################################

library(ArchR)
library(tidyverse)
library(jj)
library(chromVAR)
library(SummarizedExperiment)

source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')


mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
seurat_atac = read_rds(pick_content('mouse_normal_CD4', 'seurat_file'))
mouse_cd4_universe_gr = convert_granges(rownames(GetAssayData(seurat_atac, assay = 'scATAC_raw')))

rData = as(seurat_atac@assays$scATAC_raw@meta.features, 'DataFrame')
peak_mat = GetAssayData(seurat_atac, assay = 'scATAC_raw', slot = 'counts')
mouse_cd4_se = SummarizedExperiment::SummarizedExperiment(
  #rowData = rData, 
  rowRanges = mouse_cd4_universe_gr,
  assays = SimpleList(counts = peak_mat)
  #colData <- as(seurat_atac@meta.data, 'DataFrame')
)

library(BSgenome.Mmusculus.UCSC.mm10)
mouse_cd4_se <- addGCBias(mouse_cd4_se, 
                          genome = BSgenome.Mmusculus.UCSC.mm10)
bgpeaks <- getBackgroundPeaks(mouse_cd4_se)
tisTreg_index = rownames(peak_mat) %in% mouse_cd4_tisTreg_df$feature
bgpeaks_use = bgpeaks[tisTreg_index, ]
bgpeaks_mat = apply(bgpeaks_use, 2, function(x) rownames(peak_mat)[x])

config_list = get_config()
dataset_use = 'mouse_normal'
pconfig = yaml::read_yaml(config_list[[dataset_use]])
setwd(pconfig$ARCHR_DIR)
proj = loadArchRProject(pconfig$ARCHR_PROJECT)

emb_df = getReducedDims(proj, reducedDims = 'IterativeLSI')
pmat = get_peak_mat(proj)
olap_df_list = list()
for(i in 1:ncol(bgpeaks_mat)){
  message(i)
  olap_df_list[[i]] = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = emb_df,
                                             nFrags_vec = proj$nFrags, signature_gr = convert_granges(bgpeaks_mat[, i]),
                                             verbose = T, k = 100, count_thres = 2e5)
}

#write_rds(olap_df, paste0(storeFigPath, 'mouse_atlas_scATAC_tisTreg_sig_olap_df.RDS'))
#olap_df = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2022-11-16-percentage_signature/mouse_atlas_scATAC_tisTreg_sig_olap_df.RDS')

dr_df = jj_get_reduction_coords(proj, 'UMAP')
olap_df = do.call(cbind, lapply(olap_df_list,  '[[', 'signature_pct_overlap'))
olap_df$cluster_annotation_fine = dr_df$cluster_annotation_fine
olap_df_long = pivot_longer(olap_df, cols = -cluster_annotation_fine, names_to = 'set', values_to = 'pct_overlap')

ggplot(olap_df_long, aes(x = cluster_annotation_fine,  y = pct_overlap)) + geom_boxplot()

#Mouse atlas tisTreg signature percentage overlap boxplot
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tisTreg_signature_overlap_boxplot.pdf'),  width = 6, height = 5)
jj_plot_numeric_by_group(dr_df[!dr_df$cluster_annotation_fine == 'undefined', ], 'pct_overlap', group_column = 'cluster_annotation_fine', 
                         order = T, flip_coordinates = T, type = 'boxplot',
                         custom_colors = jj_get_colours(dr_df$cluster_annotation_fine, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')) + 
  theme(legend.position = 'none') + labs(y='% overlap', x='Cell type')
dev.off()


#######################################
## nebula on tisTreg signature
#######################################

pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
#contains the scATAC tisTreg regions as peakset
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets', pconfig$DONOR)) #use uncorrected version, because peakset is overwritten...
proj = addPeakSet(proj, peakSet = mouse_cd4_tisTreg_gr, force=T)
proj <- addPeakMatrix(proj,binarize = F)


proj2 = loadArchRProject('ArchRProject_t_nk_ilc_subset_corrected')

proj2 = proj2[!proj2$cluster_annotation_fine == 'undefined', ]
dr_df = jj_get_reduction_coords(proj2, 'UMAP')
rm(proj2)

proj = proj[proj$cellNames %in% rownames(dr_df), ]
dr_df = dr_df[match(proj$cellNames, rownames(dr_df)), ]
identical(rownames(dr_df), proj$cellNames)
dr_df$cluster_annotation_fine = factor(dr_df$cluster_annotation_fine, levels = c('CD4_Tnaive', unique(dr_df$cluster_annotation_fine)[!unique(dr_df$cluster_annotation_fine) == 'CD4_Tnaive']))
dr_df$Tissue = factor(proj$Tissue, levels = c('Spleen', 'Colon', 'VAT', 'Skin'))

df = model.matrix(~Tissue+cluster_annotation_fine, data = dr_df)

pmat = get_peak_mat(proj, binarize = FALSE)
pmat = ceiling(pmat / 2) #fixes problem of equal counts being more abundant, but still var = 1.25*mean
max(pmat) #25

re = nebula(pmat, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005, ncore =1, verbose = T)  
#saveRDS(re, paste0(bigFilesDir, 'nebula_immune_atlas_tissue_treg_signature.RDS'))

re = read_rds(paste0(bigFilesDir, 'nebula_immune_atlas_tissue_treg_signature.RDS'))
res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_fine', '', colnames(res_df))
colnames(res_df) = gsub('Tissue', '', colnames(res_df))

sig_list = convert_nebula_df_to_list(res_df, fdr_thres = 0.05, logFC_thres = 0.5)
sapply(sig_list,nrow)

sig_peak_list = lapply(sig_list, '[[', 'feature')

pdf(paste0(storeFigPath, 'immune_atlas_tissueTreg_signature_fdr05_fc0.5.pdf'),  width = 8, height = 6)
jj_plot_upsetr(sig_peak_list)
dev.off()


simplified_sig_peak_list = list(
  tisTreg_pTreg = unique(c(sig_peak_list$pTreg, sig_peak_list$tisTregST2, sig_peak_list$tisTregST2_prog)),
  Colon = sig_peak_list$Colon,
  Skin = sig_peak_list$Skin,
  VAT = sig_peak_list$VAT,
  other_cell_types = unlist(sig_peak_list[!names(sig_peak_list) %in% c('tisTregST2', 'tisTregST2_prog', 'pTreg', 'Colon', 'VAT', 'Skin')])
)
jj_plot_upsetr(simplified_sig_peak_list)



###############################################################################
# TE comparison for th17_areg vs Th17 and ILC2_Tbet vs ILC2
###############################################################################

tistreg_tes = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-19-footprints/tisTreg_sig_61_enriched_tes.txt')

library(ArchR)
library(tidyverse)
library(jj)
library(nebula)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')

dr_df_tnkilc = jj_get_reduction_coords(loadArchRProject(paste0(pconfig$ARCHR_DIR, 'ArchRProject_t_nk_ilc_subset_corrected')))
c1_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C1']
c2_bc = dr_df_tnkilc$cellNames[dr_df_tnkilc$Clusters_1.2 == 'C2']

###
#ilc2 tbet vs ilc2
###

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject("/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_t_nk_ilc_subset_nonbinary_peakset/")
proj = proj[proj$Tissue == 'Colon', ]

proj$tissue_cell_type = proj$cluster_annotation_fine
#ilc2 tbet vs ilc2
proj$tissue_cell_type[proj$cluster_annotation_fine=='ILC2' & proj$cellNames %in% c1_bc] = 'ILC2__C1'
proj$tissue_cell_type[proj$cluster_annotation_fine=='ILC2' & proj$cellNames %in% c2_bc] = 'ILC2__C2'
proj = proj[proj$tissue_cell_type %in% c('ILC2_C1', 'ILC2_C2'), ]

count_df = t(prepare_comb_norm_df(pconfig, proj, normalize = F))

dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$cluster_annotation_fine = relevel(factor(dr_df$tissue_cell_type), ref = 'ILC2__C1')

df = model.matrix(~cluster_annotation_fine, data = dr_df)

#cpc (default = 0.005)
#A non-negative threshold for filtering low-expressed genes. Genes with counts per cell smaller than the specified value will not be analyzed.
re = nebula(count_df, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005)

#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/nebula_TE_ILC2_C1_vs_C2.RDS')
re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/nebula_TE_ILC2_C1_vs_C2.RDS')

res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_fine', '', colnames(res_df))
res_df$name = rownames(res_df)

pdf(paste0(storeFigPath, 'ilc2_tbet_vs_ilc2_nebula_te_volcano.pdf'), width = 10, height = 8)
jj_plot_volcano(res_df[res_df$ILC2__C2_FDR < 0.05, ], 'ILC2__C2_logFC', 'ILC2__C2_FDR', symbol_column = 'name',labs_range = c(-2,2), marker_thres = c(1, 10)) #+ ylim(c(0,20))
jj_plot_volcano(res_df[res_df$ILC2__C2_FDR < 0.05, ], 'ILC2__C2_logFC', 'ILC2__C2_FDR', symbol_column = 'name',labs_range = c(-2,2),marker_thres = Inf, markers_highlight = tistreg_tes) #+ ylim(c(0,20))
dev.off()

###
#th17 areg vs th17
###

config_list = get_config()
pconfig = yaml::read_yaml(config_list$mouse_normal)
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject("/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/mouse_normal/ArchRProject_t_nk_ilc_subset_nonbinary_peakset/")
proj = proj[proj$Tissue == 'Colon', ]
proj = proj[proj$cluster_annotation_fine %in% c('Th17', 'Th17_Areg'), ]

dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$cluster_annotation_fine = relevel(dr_df$cluster_annotation_fine, ref = 'Th17')

count_df = t(prepare_comb_norm_df(pconfig, proj, normalize = F))

#cpc (default = 0.005)
#A non-negative threshold for filtering low-expressed genes. Genes with counts per cell smaller than the specified value will not be analyzed.
re = nebula(count_df, id = dr_df$Sample, pred=df, offset = dr_df$nFrags, cpc = 0.005)

#write_rds(re, '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/nebula_TE_Th17_vs_Th17_Areg.RDS')
re = read_rds('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-02-th1_subsets/nebula_TE_Th17_vs_Th17_Areg.RDS')

res_df = get_nebula_df(nebula_res_df = re$summary)
colnames(res_df) = gsub('cluster_annotation_fine', '', colnames(res_df))
res_df$name = rownames(res_df)

pdf(paste0(storeFigPath, 'th17_areg_vs_th17_nebula_te_volcano.pdf'), width = 10, height = 8)
jj_plot_volcano(res_df[res_df$Th17_Areg_FDR < 0.05, ], 'Th17_Areg_logFC', 'Th17_Areg_FDR', symbol_column = 'name',labs_range = c(-2,2), marker_thres = c(1, 10)) #+ ylim(c(0,20))
jj_plot_volcano(res_df[res_df$Th17_Areg_FDR < 0.05, ], 'Th17_Areg_logFC', 'Th17_Areg_FDR', symbol_column = 'name',labs_range = c(-2,2), marker_thres = Inf, markers_highlight = tistreg_tes) #+ ylim(c(0,20))
dev.off()



##############################################################################
# BATF Chip TE comparison
##############################################################################

#test tisTreg-associated TEs for enrichment of transcription factor chip seq peaks

tisTreg_ilc2_th17_areg_tes = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-06-te_pseudotime_heatmap/tisTreg_th17_areg_ilc2_te_overlap.txt')
tistreg_tes = readLines('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-09-19-footprints/tisTreg_sig_61_enriched_tes.txt')

te_gr = read_te('mm10', return_gr = T)
te_gr_use = te_gr[te_gr$name %in% tistreg_tes]

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)

pset = getPeakSet(loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected')))

#######
##GSE40918
#######

##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40918
#mm9 is used
# pCRMs_5TFs_th17.txt: 
#   pCRMs computed over ChiP-Seq datasets for BATF, IRF4, STAT3, c-MAF, RORC. Columns are: 
#   l - pCRM number, sorted based on chr and position
# chr - chromosome of pCRM
# expt - ChIP peaks comprising this pCRM, ordered by bp position
# expt.alphanum.sorted - same as expt column, but sorted alphabetically on TF name
# start,end - the start and end of the pCRM; s - median summit of peaks within pCRM
# pval - p-value of peaks, sorted as in expt column
# pval.mean - mean p-value
# span.tfs - span of each peak in pCRM
# span.l - span of entire pCRM
# peak.ids - peak IDs from MACS file
# trgt.prox - genes associated with pCRM (pCRM within +/- 5kb of TSS)
# trgt.dist - genes associated with pCRM (pCRM within gene or +/- 10kb of gene)
# dtss.prox - distances from pCRM to TSS of genes in column trgt.prox
# dtss.dist - distances from pCRM to TSS of genes in column trgt.dist.

gse40918 = read_delim('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/GSE40918_pCRMs_5TFs_th17.txt', delim = '\t')

chain = rtracklayer::import.chain('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm9ToMm10.over.chain')

olaps_df_list = list()
for(i in c('Batf', 'IRF4', 'Stat3', 'Maf', 'RORC')){
  #either use: also common binding sites
  batf_chip_peaks = gse40918[grepl(i, gse40918$expt, ignore.case = T), ]
  #or: exclusive binding sites
  #batf_chip_peaks = gse40918[gse40918$expt %in% i, ]
  
  batf_chip_peaks_gr = makeGRangesFromDataFrame(batf_chip_peaks, keep.extra.columns = T)
  print(length(batf_chip_peaks_gr))
  lift_res = rtracklayer::liftOver(batf_chip_peaks_gr, chain = chain)  
  #lift_res = lift_res[sapply(lift_res, length) == 1]
  lift_gr = unlist(lift_res)#do.call(c, lift_res)
  print(length(lift_gr))
  #get peaks overlapping with batf_chip_peaks_gr
  pset$batf_olap = granges_overlap(pset, lift_gr, return_type = 'logical', minOverlap = 0.5)
  pset_subset = pset[pset$batf_olap]
  pset_subset = pset_subset[granges_overlap(pset_subset, mouse_cd4_tisTreg_gr, return_type = 'logical')]
  print(length(pset_subset))
  olaps_df_list[[i]] = run_fisher(pset, pset_subset, te_gr_use)#te_gr_use)
}

lapply(olaps_df_list, function(x) x[x$fdr < 0.05, ])
#write_rds(olaps_df_list, paste0(storeFigPath, 'chip_gse40918_te_olaps_df_list.rds'))

pdf(paste0(storeFigPath, 'GSE40918_tf_exclusive_chip_tisTreg_TE_enrichments.pdf'), width = 10, height = 8)
for(i in names(olaps_df_list)){
  olaps_plot_df = olaps_df_list[[i]]
  #olaps_plot_df = olaps_df[olaps_df$fdr < 0.05, ]
  olaps_plot_df$logfdr = -log10(olaps_plot_df$fdr)
  
  print(ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both)), size = 4) + 
          viridis::scale_color_viridis() +  theme_minimal() + 
          #scale_colour_manual(values = cols_use_family) +
          geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + ggtitle(i) + 
          ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$logfdr > 0.5,], 
                                   aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
          labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)') + geom_hline(yintercept = -log10(0.05), colour='red', linetype = 'dashed'))
}
dev.off()



#####
#GSE121295
#####

GSE121295_batf = readr::read_tsv(file = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/GSM3430973_Treg_WT_BATF.bed', skip = 34)
colnames(GSE121295_batf) = c('chr', 'start', 'end', 'ID', 'dummy', 'strand')
GSE121295_batf_gr = convert_granges(GSE121295_batf)

chain = rtracklayer::import.chain('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm9ToMm10.over.chain')
lift_res = rtracklayer::liftOver(GSE121295_batf_gr, chain = chain)  
#lift_res = lift_res[sapply(lift_res, length) == 1]
lift_gr = unlist(lift_res)#do.call(c, lift_res)


pset = getPeakSet(loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected')))
pset$batf_olap = granges_overlap(pset, lift_gr, return_type = 'logical', minOverlap = 0)
table(pset$batf_olap)
pset_subset = pset[pset$batf_olap]
pset_subset = pset_subset[granges_overlap(pset_subset, mouse_cd4_tisTreg_gr, return_type = 'logical')]
print(length(pset_subset))
olaps_df = run_fisher(pset, pset_subset, te_gr_use)
#write_csv(olaps_df, paste0(storeFigPath, 'GSE121295_batf_te_overlap.csv'))

olaps_df[olaps_df$fdr < 0.05, ]

olaps_df$logfdr = -log10(olaps_df$fdr)
olaps_plot_df = olaps_df
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 10] = 10

pdf(paste0(storeFigPath, 'GSE121295_batf_chip_tisTreg_TE_enrichments.pdf'), width = 10, height = 8)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both)), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$logfdr > 0.5,], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)') + geom_hline(yintercept = -log10(0.05), colour='red', linetype = 'dashed')
dev.off()

######
###JunB
######
GSE121295_junb = readr::read_tsv(file = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/GSM3430972_Treg_WT_JunB.bed', skip = 34)
colnames(GSE121295_junb) = c('chr', 'start', 'end', 'ID', 'dummy', 'strand')
GSE121295_junb_gr = convert_granges(GSE121295_junb)

lift_res = rtracklayer::liftOver(GSE121295_junb_gr, chain = chain)  
lift_gr = unlist(lift_res)#do.call(c, lift_res)


pset = getPeakSet(loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected')))
pset$batf_olap = granges_overlap(pset, lift_gr, return_type = 'logical', minOverlap = 0)
table(pset$batf_olap)
pset_subset = pset[pset$batf_olap]
pset_subset = pset_subset[granges_overlap(pset_subset, mouse_cd4_tisTreg_gr, return_type = 'logical')]
print(length(pset_subset))
olaps_df = run_fisher(pset, pset_subset, te_gr_use)
#write_csv(olaps_df, paste0(storeFigPath, 'GSE121295_junb_te_overlap.csv'))

olaps_df$logfdr = -log10(olaps_df$fdr)
olaps_plot_df = olaps_df
olaps_plot_df$odds_ratio[olaps_plot_df$odds_ratio > 10] = 10

pdf(paste0(storeFigPath, 'GSE121295_junb_chip_tisTreg_TE_enrichments.pdf'), width = 10, height = 8)
ggplot() + geom_point(data = olaps_plot_df, aes(x = odds_ratio, y = logfdr, colour = log10(both)), size = 4) + 
  viridis::scale_color_viridis() +  theme_minimal() + 
  #scale_colour_manual(values = cols_use_family) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data=olaps_plot_df[olaps_plot_df$logfdr > 0.5,], 
                           aes(x=odds_ratio, y=logfdr,label=name), max.overlaps = 500) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'log10(sig. peaks with TE olap)') + geom_hline(yintercept = -log10(0.05), colour='red', linetype = 'dashed')
dev.off()



#########################################################
# TF enrichment in TE regions
#########################################################

te_gr = read_te('mm10')
pgr = getPeakSet(proj)
#te_gr$peak_olap = granges_overlap(te_gr, pgr, return_type = 'logical')
#table(te_gr$peak_olap)

sig_list = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-06-05-tissue_peaks_nebula/tnkilc_nebula_TE_lfc05_fdr05.xlsx')
peripheral_sigs = lapply(sig_list, '[[', 'feature')
common_tissue = convert_granges(Reduce(intersect, peripheral_sigs[c('Colon', 'VAT', 'Skin')] ))
any_tissue = convert_granges(unique(Reduce(c,  peripheral_sigs[c('Colon', 'VAT', 'Skin')])))

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_corrected'))
proj = archr_add_peak_signatures(proj, signature_list = list(all_peripheral = common_tissue, any_peripheral = any_tissue), 'peripheral_sig_')
dr_df = jj_get_reduction_coords(proj, 'UMAP')
#pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_nebula_peripheral_signatures.pdf'),  width = 10, height = 8)
jj_plot_features(dr_df, features = c('z_all_peripheral', 'z_any_peripheral'))
#dev.off()

te_gr$peripheral_tissue_olap = granges_overlap(te_gr, any_tissue, return_type = 'logical')
table(te_gr$peripheral_tissue_olap)
te_gr_use = te_gr[te_gr$peripheral_tissue_olap]
n_df = te_gr_use$name %>% table %>% as.data.frame
colnames(n_df) = c('name', 'n')
n_df = n_df %>% dplyr::filter(n > 10) %>% dplyr::arrange(desc(n))
for(i in n_df$name){
  jj_save_bed(te_gr_use[te_gr_use$name == i], paste0(storeFigPath, i, '.bed'))
}
paste(gsub('\\.bed', '', list.files('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-06-homer_peripheral_tes')), collapse = ',')

tes_plot = n_df %>% dplyr::filter(n >= 10) %>% dplyr::pull(name) %>% as.character

homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2024-01-06-homer_peripheral_tes/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = F)
homer_subsets = homer_subsets[homer_subsets %in% tes_plot]
homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(homer_res, i, '/knownResults.txt'))
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = i
  homer_df = as.data.frame(homer_df)
  homer_list[[i]] = homer_df
}

homer_list_use = lapply(homer_list, function(x) x[x$`q-value (Benjamini)` < 0.001, ] )# & x$enrichment > 2, ])

sapply(homer_list_use, nrow)
homer_list_use = homer_list_use[sapply(homer_list_use, nrow) > 0]
#homer_motif_list = lapply(homer_list_use, '[[', 1)
#jj_plot_upsetr(homer_motif_list)

homer_list_gsea_plot = lapply(homer_list_use, function(x) x[, c('name','motif', 'enrichment', 'rank')])
# Combine the cell-type-specific data sets.
gsea_res_comb <- do.call(rbind, homer_list_gsea_plot)
rownames(gsea_res_comb) = NULL
gsea_res_comb$rank = NULL
gsea_res_comb = gsea_res_comb[!duplicated(gsea_res_comb), ]
gsea_res_comb = gsea_res_comb %>% distinct() #%>% dplyr::group_by(name, motif) %>% count() %>% arrange(desc(freq))
gsea_res_comb = gsea_res_comb %>% dplyr::rename(TE = name)

#gsea_res_comb$motif_clean = make.names(replace_special_chars(sapply(strsplit(as.character(gsea_res_comb$motif), split = '/'), '[[', 1), pattern = '[^[:alnum:] ]'))
#gsea_res_comb$motif_short = sapply(strsplit(gsea_res_comb$motif_clean, split = '\\('), '[[', 1)

res = gsea_res_comb %>% 
  # .[complete.cases(.), ] %>% 
  # .[!duplicated(.), ] %>% 
  pivot_wider(., id_cols = motif, names_from =  TE, values_from = enrichment, values_fill = 0) %>% 
  as.data.frame %>% set_rownames(.$motif) %>% dplyr::select(-motif) %>% as.matrix


# res = gsea_res_comb %>%
#   group_by(TE, motif) %>%
#   dplyr::summarise(enrichment = mean(enrichment, na.rm = T)) %>%
#   tidyr::pivot_wider(names_from = TE, values_from = enrichment, 
#                      values_fill = 0) %>% 
#   as.data.frame %>% set_rownames(.$motif) %>% dplyr::select(-motif) %>% as.matrix

res[is.na(res)] = 0
rownames(res) = sapply(strsplit(rownames(res), split = '/'), '[[', 1)

homer_headers = read_homer('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/custom.motifs', headers_only = T)
homer_headers_names = sapply(strsplit(homer_headers, '\t'),'[[',2)
binding_domains = gsub('.*\\((.*)\\)', '\\1', sapply(strsplit(homer_headers_names, '/'), '[[',1))
binding_domains = names(table(binding_domains))[as.vector(table(binding_domains) > 2)]
binding_domains = binding_domains[!binding_domains == '?']
binding_domains[binding_domains=='T-box'] = 'T.box'
binding_domains = binding_domains[!binding_domains %in% c("Paired,Homeobox", "POU,Homeobox")]
name_df = data.frame(name = rownames(res))
name_df$bindingDomain = 'Other'
for(i in binding_domains){
  name_df$bindingDomain[grepl(i, name_df$name)] = i  
}
ha = HeatmapAnnotation(domain = name_df$bindingDomain, col = list(domain =  jj_get_colours(name_df$bindingDomain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')))

family_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-07-03-finalizing_te_analysis/mm10_te_family.csv')
family_annot = family_df$Family[match(colnames(res), family_df$TE)]

colcsv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv'
ra = rowAnnotation(family = family_annot, col = list(family =  jj_get_colours(family_annot, colour_csv = colcsv, comment_char = '$')))
# rnames = res[, 1]
# res = res[, -1]
# res = apply(res, 2, as.numeric)
# rownames(res) = rnames
library(circlize)
pdf(paste0(storeFigPath, 'te_homer_heatmap_tissue.pdf'), width = 40, height = 20)
Heatmap(t(res), 
        top_annotation = ha, 
        right_annotation = ra,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8), 
        name = 'Enrichment',
        col = colorRamp2(c(0, 0.0001, 5, 10),c("grey95",'#440154' , '#21918c', "#fde725")))
dev.off()


gsea_res_comb$motif = sapply(strsplit(gsea_res_comb$motif, split = '/'), '[[', 1) #gsub('(.*)/Homer', '\\1', gsub('(.*)\\(GSE.*', '\\1', gsub('(.*)-ChIP.*','\\1', gsea_res_comb$pathway)))

levels_use = gsea_res_comb %>% 
  dplyr::group_by(motif) %>% 
  dplyr::summarise(enr_sum = sum(enrichment)) %>% 
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = levels_use)

#pdf(paste0(storeFigPath, ''), width =5, height=8)
ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + 
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) + 
  labs(x = 'Associated peaks (Nebula)', y = '', fill = 'Enrichment') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()

