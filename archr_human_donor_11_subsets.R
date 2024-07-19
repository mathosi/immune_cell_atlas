library(jj)
library(ArchR)
#mouse subset analysis\
#addArchRThreads(threads = 1) 
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor11.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
plotEmbedding(proj, colorBy='GeneScoreMatrix', name=c('BCL6'), imputeWeights = getImputeWeights(proj), plotAs = 'points')
jj_plot_features(reduction = dr_df, meta_features = 'Clusters_0.6', pt.size = 0.5, label = T)
jj_plot_features(reduction=dr_df, meta_features='Tissue', label = F, pt.size = 0.5, custom_colors = jj_get_colours(dr_df$Tissue, '/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/colour_map.csv'))
jj_plot_features(reduction = dr_df, meta_features = 'singler_label', pt.size = 0.5, label = T)


proj$cluster_annotation = from_to(annot_vec = proj$Clusters_0.6, 
                                  old_new_map_vec = c(
                                    C1='Monocyte/Neutrophil', 
                                    C2='Monocyte/Neutrophil',
                                    C3='Macrophage/DC/Monocyte',
                                    C4='Macrophage/DC/Monocyte', 
                                    C5='Macrophage/DC/Monocyte', 
                                    C6='B cell', 
                                    C7='B cell', 
                                    C8='undefined', 
                                    C9='CD4 T cell', 
                                    C10='CD4 T cell/NK cell', 
                                    C11='CD4 T cell', 
                                    C12='CD8 T cell', 
                                    C13= 'CD4 T cell', 
                                    C14='CD4 T cell',
                                    C15='NK cell',
                                    C16='CD8 T cell',
                                    C17='MAIT cell', 
                                    C18='NK cell', 
                                    C19='NK cell'
                                  ))

proj$cluster_annotation_level2 = from_to(annot_vec = proj$Clusters_0.6, 
                                  old_new_map_vec = c(
                                    C1='Monocyte/Neutrophil',
                                    C2='CD16 Monocyte', 
                                    C3='', 
                                    C4='', 
                                    C5='', 
                                    C6='B memory', 
                                    C7='B naive', 
                                    C8='', 
                                    C9='',
                                    C10='', 
                                    C11='CD4 naive',
                                    C12='CD8 naive', 
                                    C13= 'CD4 naive/memory', 
                                    C14='CD4/Treg', 
                                    C15='NK cell',
                                    C16='CD8 memory',
                                    C17='MAIT cell', 
                                    C18='NK cell', 
                                    C19='NK cell'
                                  ))


#blueprint reference (major cell types such as fibroblasts included)
# singler_pred = read_csv('/omics/groups/OE0436/internal/msimon/scATAC/imm_atlas_human_donor11_normal_gene_activity_mat_singler_ref_blueprintencodedata_label.fine_predictions.csv')
# singler_pred = singler_pred[match(rownames(dr_df), singler_pred$X1), ]
# stopifnot(identical(singler_pred$X1, rownames(dr_df)))
# dr_df$singler_blueprint = singler_pred$labels
# jj_plot_features(reduction = dr_df, meta_features = 'singler_blueprint', pt.size = 0.5, label = T)

singler_pred = read_csv('/omics/groups/OE0436/internal/msimon/scATAC/imm_atlas_human_donor11_normal_gene_activity_mat_singler_ref_monacoimmunedata_label.fine_predictions.csv')
singler_pred = singler_pred[match(rownames(dr_df), singler_pred$X1), ]
stopifnot(identical(singler_pred$X1, rownames(dr_df)))
proj$singler_label_fine = singler_pred$labels
dr_df = get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'singler_label_fine', pt.size = 0.5, label = T)

##singler peaks pred
pmat = get_peak_mat(proj)
library(genomic_region_tools)
ass_df = GetAssayData(pbmc_multi, assay = 'ATAC')
unified_list = unify_peak_matrices(pmat_list=list(donor11=pmat, pbmc=ass_df))
summ_mat = jj_summarize_sparse_mat(unified_list$pbmc, pbmc_multi$annot)
library(SingleR)
preds <- SingleR(test = unified_list$donor11, ref = summ_mat, assay.type.test=1, labels = colnames(summ_mat))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
dr_df$singler_label_peaks = preds$labels
jj_plot_features(reduction=dr_df, meta_features = 'singler_label_peaks')

###

#b cells: cluster 2 could also be ILC, cluster 4 is also suspicious
proj_b = proj[proj$Clusters_0.6 %in% paste0('C',6:7), ]
proj_b = archr_dim_reduction(proj_b, cluster_res = NULL)
dr_df = jj_get_reduction_coords(proj_b, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 1)
#saveArchRProject(proj_b, outputDirectory='ArchRProject_donor11_B_subset', load=F)
proj =  loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_B_subset'))
proj = archr_clustering(proj, cluster_res = 1.1)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Clusters_1.1', pt.size = 0.5, label = T)
proj$cluster_annotation = from_to(vec= proj$Clusters_1.1, old_new_map_vec=c(
  'C1'= 'other',
  'C2'= 'Naive B cell',
  'C3'= 'Naive B cell',
  'C4'= 'Plasmablast',
  'C5'= 'Memory B cell',
  'C6'= 'Memory B cell'
))


#innate immune cells: 
proj_i = proj[proj$Clusters_0.6 %in% paste0('C',c(1:5, 8)), ]
proj_i = archr_dim_reduction(proj_i, cluster_res = NULL)
dr_df = jj_get_reduction_coords(proj_i, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
#saveArchRProject(proj_i, outputDirectory='ArchRProject_donor11_macrophage_dc_subset', load=F)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_macrophage_dc_subset'))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Clusters_2', pt.size = 0.5, label = T)

proj$cluster_annotation = from_to(vec= proj$Clusters_2, old_new_map_vec=c(
  'C1'= 'Neutrophil/Basophil',
  'C2'= 'other', 
  'C3'= 'Monocyte', 
  'C4'= 'Monocyte', 
  'C5'= 'Monocyte', 
  'C6'= 'Monocyte', 
  'C7'= 'Monocyte',
  'C8'= 'DC',
  'C9'= 'DC',
  'C10'= 'DC',
  'C11'= 'DC', 
  'C12'= 'pDC'
))
#saveArchRProject(proj)

marker_fine = unlist(marker_list_fine)
names(marker_fine) = gsub('^(.*)[0-9]$', '\\1', names(marker_fine))
jj_plot_heatmap(gmat, features_use = marker_fine, group_vec = proj$Clusters_2)

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation', pt.size = 0.5, label = T)

celltypist_markers = read_delim('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/encyclopedia_table.csv', delim = '\t')
marker_list_coarse = marker_list_fine = list()
for(i in 1:nrow(celltypist_markers)){
  marker_list_coarse[[celltypist_markers$`High-hierarchy cell types`[i] ]] = c(marker_list_coarse[[celltypist_markers$`High-hierarchy cell types`[i] ]], unlist(strsplit(celltypist_markers$`Curated markers`[i], ', ')))
}
for(i in 1:nrow(celltypist_markers)){
  marker_list_fine[[celltypist_markers$`High-hierarchy cell types`[i] ]] = unlist(strsplit(celltypist_markers$`Curated markers`[i], ', '))
}

gmat = get_gene_mat(proj)
lapply(marker_list_fine, function(x) x[!x %in% rownames(gmat)])
marker_list_fine = lapply(marker_list_fine, function(x) x[x %in% rownames(gmat)])
proj = addModuleScore(proj, features = marker_list_fine, useMatrix = 'GeneScoreMatrix')
dr_df = jj_get_reduction_coords(proj, 'UMAP')

colnames(dr_df) = replace_special_chars(colnames(dr_df))
gg = jj_plot_features(reduction = dr_df, meta_features = grep('Module', colnames(dr_df), value = T), cap_top = 'auto', return_gg_object = T)
pdf(paste0(storeFigPath, 'human_donor_11_myeloid_marker_modules.pdf'), width = 8, height=7)
gg
dev.off()

mlist = archr_plot_markers(proj, marker_genes = markers_plot, impute_weights = T, reduction = 'UMAP')
pdf(paste0(storeFigPath, 'human_donor_11_myeloid_marker_genes.pdf'), width = 8, height=7)
lapply(seq_along(mlist), function(x) try(mlist[[x]] + ggtitle(paste(gsub('^(.*?)[0-9]+$','\\1',names(markers_plot[x])), markers_plot[x], sep = ': '))))
dev.off()
gmat = get_gene_mat(proj)
markers_plot = unlist(marker_list)[unlist(marker_list) %in% rownames(gmat)]
tannot = HeatmapAnnotation(region = gsub('^(.*?)[0-9]+$','\\1',names(markers_plot)))
jj_plot_heatmap(gmat, features_use = markers_plot, group_vec = proj$Clusters_2,
                cluster_rows = T, cluster_columns = T, top_annot = tannot)

#tissue t, nk, ilcs
proj_t = proj[proj$Clusters_0.6 %in% paste0('C',c(9,10,18)), ]
proj_t = archr_dim_reduction(proj_t, cluster_res = NULL)
dr_df = jj_get_reduction_coords(proj_t, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
#saveArchRProject(proj_t, outputDirectory='ArchRProject_donor11_tissue_t_nk_subset', load=F)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_tissue_t_nk_subset'))
proj = archr_clustering(proj, cluster_res = 1)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Clusters_1', pt.size = 0.5, label = T)
proj$cluster_annotation = from_to(vec= proj$Clusters_1, old_new_map_vec=c(
  'C1'= 'other', 
  'C2'= 'NK cell', 
  'C3'= 'NK cell', 
  'C4'= 'CD8 memory T cell',
  'C5'= 'CD4 memory T cell',
  'C6'= 'CD8 memory T cell', 
  'C7'= 'tisTreg', 
  'C8'='CD4 memory T cell',
  'C9'='CD4 memory T cell', 
  'C10'='CD4 memory T cell' 
))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation', pt.size = 0.5, label = T)
#saveArchRProject(proj)

#blood t, nk, ilcs
proj_t2 = proj[proj$Clusters_0.6 %in% paste0('C',c(11:17,19)), ]
proj_t2 = archr_dim_reduction(proj_t2, cluster_res = NULL)
dr_df = jj_get_reduction_coords(proj_t2, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
jj_plot_features(reduction = dr_df, meta_features = 'singler_label', pt.size = 0.5, label=T)
#saveArchRProject(proj_t2, outputDirectory='ArchRProject_donor11_blood_t_nk_subset', load=F)
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_blood_t_nk_subset'))
proj = archr_clustering(proj, cluster_res = 1)
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'Clusters_1', pt.size = 0.5, label = T)
proj$cluster_annotation = from_to(vec= proj$Clusters_1, old_new_map_vec=c(
  'C1'= 'CD8 naive T cell',
  'C2'= 'CD4 naive T cell',
  'C3'= 'CD4 naive T cell', 
  'C4'= 'Treg cell',
  'C5'= 'Tfh/Th1', 
  'C6'= 'Tfh/Th1/Th2',
  'C7'= 'Th1/Th17', 
  'C8'= 'Th17 cell', 
  'C9'= 'MAIT cell',
  'C10'= 'CD8 memory T cell', 
  'C11'= 'CD8 memory T cell', 
  'C12'= 'CD8 memory T cell',
  'C13'= 'NK cell', 
  'C14'= 'NK cell',
  'C15'= 'NK cell' 
))
dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'cluster_annotation', pt.size = 0.5, label = T)
#saveArchRProject(proj)


# put annotations back into main object -----------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_macrophage_dc_subset'))
myeloid_df = data.frame(cellNames = proj$cellNames, cluster_annotation_level3 = proj$cluster_annotation)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_B_subset'))
b_df = data.frame(cellNames = proj$cellNames, cluster_annotation_level3 = proj$cluster_annotation)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_tissue_t_nk_subset'))
tissue_lymphoid_df = data.frame(cellNames = proj$cellNames, cluster_annotation_level3 = proj$cluster_annotation)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_donor11_blood_t_nk_subset'))
blood_lymphoid_df = data.frame(cellNames = proj$cellNames, cluster_annotation_level3 = proj$cluster_annotation)
cluster_annotation_level3_df = rbind(myeloid_df, b_df, tissue_lymphoid_df, blood_lymphoid_df)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
dr_df = as.data.frame(proj@cellColData)
dr_df = dr_df %>% dplyr::left_join(cluster_annotation_level3_df, by = 'cellNames')
rownames(dr_df) = dr_df$cellNames
#dr_df$cluster_annotation_level3[is.na(dr_df$cluster_annotation_level3)] = dr_df$cluster_annotation[is.na(dr_df$cluster_annotation_level3)]
stopifnot(!anyNA(dr_df$cluster_annotation_level3))
proj@cellColData = as(dr_df, 'DataFrame')

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_plot_features(reduction=dr_df, meta_features = 'cluster_annotation_level3')

proj$cluster_annotation_level2 = from_to(vec= proj$cluster_annotation_level3, old_new_map_vec=c(
  'CD4 memory T cell'= 'CD4_Tmem',
  'CD4 naive T cell'= 'CD4_Tnaive',
  'CD8 memory T cell'= 'CD8_Tmem',
  'CD8 naive T cell'= 'CD8_Tnaive',
  'DC'= 'DC',
  'MAIT cell'= 'MAIT',
  'Memory B cell'= 'B_memory',
  'Monocyte'= 'Monocyte',
  'Naive B cell'= 'B_naive',
  'Neutrophil/Basophil'= 'Neutrophil/Basophil',
  'NK cell'= 'NK',
  'other'= 'other',
  'pDC'= 'other',
  'Plasmablast'= 'Plasmablast',
  'T naive'= 'Tnaive',
  'Tfh/Th1'= 'CD4_Tmem',
  'Tfh/Th1/Th2'= 'CD4_Tmem',
  'Th1/Th17'= 'CD4_Tmem',
  'Th17 cell'= 'CD4_Tmem',
  'tisTreg'= 'tisTreg',
  'Treg cell'= 'Treg'
))

proj$cluster_annotation_level1 = from_to(vec= proj$cluster_annotation_level2, old_new_map_vec=c(
  'B_memory'= 'B',
  'B_naive'= 'B',
  'CD4_Tmem'= 'CD4_T',
  'CD4_Tnaive'= 'CD4_T',
  'CD8_Tmem'= 'CD8_T',
  'CD8_Tnaive'= 'CD8_T',
  'DC'= 'DC',
  'MAIT'= 'T',
  'Monocyte'= 'Monocyte',
  'Neutrophil/Basophil'= 'Neutrophil/Basophil',
  'NK'= 'NK',
  'other'= 'other',
  'Plasmablast'= 'B',
  'tisTreg'= 'CD4_T',
  'Tnaive'= 'T',
  'Treg'= 'CD4_T'
))
#saveArchRProject(proj)

# chromvar ----------------------------------------------------------------
proj = loadArchRProject(paste0(pconfig$ARCHR_DIR,'ArchRProject_donor11_tissue_t_nk_subset'))
dr_df = get_reduction_coords(proj, redname='UMAP')
#transfer from full project
#dr_df$z_nr1 = dr_df_full$z_nr1[match(dr_df$cellNames,dr_df_full$cellNames)]

if(pconfig$GENOME=='mm10'){
  signature_list = read_rds('/abi/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
  signature_list = signature_list[-c(4)]
  signature_list = lapply(signature_list, makeGRangesFromDataFrame)
  #early and late precursor sig cause error 
  # Error in Matrix::sparseMatrix(i = queryHits(overlapRegions), j = match(names(allPositions),  : 
  #                                                                          NA's in (i,j) are not allowed
  proj <- addPeakAnnotations(ArchRProj = proj, 
                             regions = signature_list[1:5],
                             name = "tisTreg",
                             force=T)
  #getPeakAnnotation(proj, 'tisTreg')
  #stores list of 3 anntoations: 
  #$Name: name of the signatures
  #$Positions: GRangesList with Granges for each signature
  #$Matches: SummarizedExperiment with logical sparse assay 'matches' that contains cells in rows and signatures in columns (True = Match of peak and signature)
  proj <- addBgdPeaks(proj, force = T, method = 'ArchR')
  
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
  # dev_score_mat = t(assays(sig_se)[['deviations']])

  #dr_df$cell_id = dr_df$cellNames
  dr_df = dr_df %>% dplyr::select(-ends_with('sig')) %>% 
    dplyr::left_join(rownames_to_column(as.data.frame(z_score_mat), 'cell_id'), by = 'cell_id')
  dr_df = dr_df %>% dplyr::select(-ends_with('sig')) %>% 
    dplyr::left_join(rownames_to_column(as.data.frame(z_score_mat), 'cellNames'), by = 'cellNames')
  
  
  #plot(dr_df$colon_treg_sig.x, dr_df$colon_treg_sig.y)
  features_plot = colnames(z_score_mat)
  gg = plotFeatures(seurat_obj = NULL, 
                    reduction = dr_df,
                    pt.size = 0.5,
                    meta_features = features_plot,
                    cont_or_disc = 'c', 
                    colorScale = 'viridis', 
                    topCap = 'q95', bottomCap = 'q05',
                    custom_theme = theme_minimal(), return_gg_object = T)
  
  # p <- plotEmbedding(
  #   ArchRProj = proj, 
  #   colorBy = "tisTregMatrix", 
  #   name = signatures, 
  #   embedding = "UMAP",
  #   imputeWeights = getImputeWeights(proj)
  # )
  
  pdf(paste0(storeFigPath, 't_subset_tisTreg_signatures.pdf'), width=12, height=8)
  #p
  gg
  dev.off()
}else{

  signature_list = read_rds('/abi/data2/simonma/projects/scATAC/analysis/2020-04-30-tconv_signature_substraction/tisTreg_vs_Tconv_filtering_plot_signatures_hg38.RDS')
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = signature_list, name = "tisTreg", force=T)
  proj <- addBgdPeaks(proj, force = T, method = 'ArchR')
  
  proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "tisTreg",
    force = TRUE
  )
  
  #plotVarDev <- getVarDeviations(proj, plot = TRUE, name = "tisTreg2Matrix")
  #plotVarDev
  
  # signatures <- getFeatures(proj, useMatrix = "tisTregMatrix")
  # signatures <- sort(grep("z:", signatures, value = TRUE))
  # signatures
  
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


z_score_mat = t(assays(sig_se)[['z']])
dev_score_mat = t(assays(sig_se)[['deviations']])
z_score_mat = z_score_mat[match(proj$cellNames, rownames(z_score_mat)), ]
identical(proj$cellNames, rownames(z_score_mat))
colnames(z_score_mat) = paste0('z_', colnames(z_score_mat))
dev_score_mat = dev_score_mat[match(proj$cellNames, rownames(dev_score_mat)), ]
identical(proj$cellNames, rownames(dev_score_mat))
colnames(dev_score_mat) = paste0('dev_', colnames(dev_score_mat))

proj = add_df_to_cellcoldata(proj, z_score_mat)
proj = add_df_to_cellcoldata(proj, dev_score_mat)
dr_df_full = get_reduction_coords(proj, 'UMAP')

# plot markers ------------------------------------------------------------
proj = addImputeWeights(proj)
plotEmbedding(proj, colorBy='GeneScoreMatrix', embedding='UMAP',name = 'FOXP3', size = 1, rastr = F, plotAs = 'points')
plotEmbedding(proj, colorBy='GeneScoreMatrix', embedding='UMAP',name = 'GZMA', size = 1, rastr = F, plotAs = 'points')

gene_mat = get_gene_mat(proj)
dr_df <- get_reduction_coords(proj, 'UMAP')
identical(colnames(gene_mat), rownames(dr_df))

marker_genes = c('Foxp3', 'Il1rl1', 'Tigit', 'Pdcd1', 'Il2ra', 'Cd4', 'Cd8a')
marker_genes = c('Nt5e', 'Cd40') #c('Cd19','Sdc1','Prdm1','Cd1d1','Havcr2','Pdcd1')
marker_genes = 'GZMA'

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
                 pt.size = 1,
                 meta_features = marker_genes, #'Foxp3',
                 cont_or_disc = 'c', colorScale = 'wbr', cap_top = 'auto', order=F,
                 custom_theme = theme_minimal())


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


sig_se = getPeakAnnotation(proj, 'tisTreg')
match_se = read_rds("/omics/groups/OE0436/internal/msimon/scATAC/archr/imm_cell_atlas/human_normal/ArchRProject_donor11_tissue_t_nk_subset/Annotations/tisTreg-Matches-In-Peaks.rds")
match_vec = assays(match_se)[[1]][, 1]
# FALSE   TRUE 
# 212025   3993
peak_mat = getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
peak_mat_sp =  assays(peak_mat)[[1]]
peak_mat_sp_keep = peak_mat_sp[match_vec, ]
Heatmap(as.matrix(peak_mat_sp_keep[1:100, 1:100]))


# misc --------------------------------------------------------------------

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters_0.5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

proj <- addArchRAnnotations(ArchRProj = proj, collection = "ATAC")
enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 15, transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "ATAC",
  force = TRUE
)

sig_se = getMatrixFromProject(proj, 'ATACMatrix')
z_score_mat = t(assays(sig_se)[['z']])
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))
proj = add_df_to_cellcoldata(proj, z_score_mat)
dr_df = jj_get_reduction_coords(proj, redname = 'UMAP')
jj_plot_features(reduction = dr_df, meta_features = 'IAtlas_B_Bulk')

plotVarDev <- getVarDeviations(proj, plot = TRUE, name = "ATACMatrix")
plotVarDev
