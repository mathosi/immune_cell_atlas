library(jj)
library(ArchR)
pconfig = yaml::read_yaml('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor13.yaml')
addArchRGenome(pconfig$GENOME)
bigFilesDir <- "/omics/groups/OE0436/internal/msimon/scATAC/"
setwd(pconfig$ARCHR_DIR) #important since archr works with relative paths
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))

library(clustree)
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df = dr_df[, colnames(dr_df) %in% paste0('Clusters_', seq(0,2, 0.1))]
clustree(dr_df[,1:20], prefix = "Clusters_", show_axis=T)
dr_df = get_reduction_coords(proj, 'UMAP')
dr_df$Clusters_1.2 = factor(dr_df$Clusters_1.2, levels=gtools::mixedsort(unique(dr_df$Clusters_1.2)))
jj_plot_features(reduction=dr_df, meta_features = 'Clusters_1.2', label_type='geom_label', label_col = 'white')

dr_df = jj_get_reduction_coords(proj, 'UMAP')
jj_arrange_ggplots(cols = 3, jj_plot_features(reduction = dr_df, meta_features = c('Tissue','Clusters_1.2','singler_label'), return_gg_object = T))
plotEmbedding(proj, colorBy='GeneScoreMatrix', name=c('BCL6'), imputeWeights = getImputeWeights(proj), plotAs = 'points')

jj_plot_features(reduction = dr_df, meta_features = 'Clusters_0.6', pt.size = 0.5, label_type = 'geom_label')
jj_plot_features(reduction = dr_df, meta_features = 'Tissue', pt.size = 0.5)
jj_plot_features(reduction = dr_df, meta_features = 'singler_label', pt.size = 0.5, label = T)

features_plot = list(
  Tcells= c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b'), 
  Bcells=c('Ms4a1', 'Cd79a','Cd79b'),
  Monocyte = c('Csf1r','Cx3cr1','Cd14'),
  DCs=c('LAMP3','HLA-DRA','CXCR3','FLT3','CD74','SIGLEC1'),
  NK=c('Gzma','Prf1', "Klrb1",'Gzmb'),
  neutrophil = c('Itgam','S100a9'), #'Ly6g' not available
  Mastcells_Eo_Basophil=c('Itgam',"Il4"),
  t_nk_ilc = c('Klrg1','Id2','Il7r')
) 
features_plot = lapply(features_plot, toupper)

gmat = get_gene_mat(proj)
seurat_dotplot(gmat=gmat, metadf = dr_df, features = unique(unname(unlist(features_plot))), group_column = 'Clusters_0.5')

proj$cluster_annotation_level3 = from_to(vec= dr_df$Clusters_0.5, old_new_map_vec=c(
  'C1'= 'other',
  'C2'= 'other', 
  'C3'= 'NK cell',
  'C4'= 'effector CD8 T cell',
  'C5'= 'naive CD4 T cell',
  'C6'= 'CD4 T cell',
  'C7'= 'Treg',
  'C8'= 'B cell', 
  'C9'= 'naive B cell',
  'C10'= 'Monocytes', 
  'C11'= 'Monocytes',  
  'C12'= 'Monocytes',
  'C13'= 'DC', 
  'C14'= 'DC', 
  'C15'= 'Myeloid progenitor/Neutrophil (netosis?)',
  'C16'= 'Neutrophil'
))

dr_df = as.data.frame(proj@cellColData)
dr_df$singler_monaco = dr_df$singler_label
dr_df$singler_nover = preds$labels
dr_df$singler_dice = preds2$labels
dr_df$singler_blueprint = preds3$labels
dr_df$singler_hpca = preds4$labels
#saveArchRProject(proj)
