
# functions ---------------------------------------------------------------

# set parameters ----------------------------------------------------------

dataset_use = 'human_normal_donor11'
dataset_use = 'mouse_normal'
#Sys.setenv(R_SESSION_TMPDIR='/abi/data2/simonma/temp')


# load libraries and set options ------------------------------------------
#run always

library(ArchR)
library(tidyverse)
source('/omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R')
#this will change to /omics/groups/OE0436/data2/simonma/projects/scATAC/scripts/useful_functions_sc.R

addArchRThreads(threads = 1) 

config_list = list(
  mouse_normal = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal.yaml',
  mouse_nfil3_gfp = '/omics/groups/OE0436/data2/simonma/projectsimm_cell_atlas/scripts/config_mouse_nfil3_gfp.yaml',
  mouse_areg_gfp = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_areg_gfp.yaml',
  human_normal = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal.yaml',
  human_normal_donor11 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor11.yaml', #use for now
  human_normal_donor13 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_donor13.yaml',
  mouse_normal_cd4 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal_cd4.yaml',
  human_normal_cd4 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_cd4.yaml',
  human_normal_cd8 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_normal_cd8.yaml',
  mouse_normal_cd8 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_normal_cd8.yaml',
  mouse_tumor_t = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_mouse_tumor_t.yaml',
  human_hcc = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/config_human_hcc.yaml'
)
pconfig = yaml::read_yaml(config_list[[dataset_use]])

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

# create archr project ----------------------------------------------------


#get fragment files
inputFiles = Sys.glob(sprintf('%s*/*/outs/fragments.tsv.gz', pconfig$CELLRANGER_OUTPUT_DIR))
names(inputFiles) = basename(dirname(dirname(inputFiles)))

#mouse normal cd4:
pheno_df = read_rds(pick_content('mouse_normal_CD4', 'seurat_pheno'))
inputFiles = inputFiles[names(inputFiles) %in% unique(pheno_df$sample_name)]
#human normal cd4:
pheno_df = read_rds(pick_content('human_CD4', 'seurat_pheno'))
inputFiles = inputFiles[names(inputFiles) %in% unique(pheno_df$sample_name)]
#mouse normal cd8
#mouse normal cd4:
pheno_df = read_rds(pick_content('mouse_normal_CD8', 'seurat_pheno'))
inputFiles = inputFiles[names(inputFiles) %in% unique(pheno_df$sample_name)]
#mouse her2 t cells
inputFiles = inputFiles[grepl('tumor', names(inputFiles),ignore.case = T)]
#human_hcc
pheno_df = read_rds(pick_content('human_hcc', 'seurat_pheno'))
inputFiles = inputFiles[names(inputFiles) %in% unique(pheno_df$sample_name)]

# # # #do this in clusterjob
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 1000,
  maxFrags = 1e+05,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force=T
)
# 
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "ArchRProject_raw", #relative to working directory
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#saveArchRProject(proj, outputDirectory='ArchRProject_raw', load=F)

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
#proj = loadArchRProject(paste0(pconfig$ARCHR_DIR, 'ArchRProject_mouse_normal_CD4'))

#mouse/human cd4
pheno_df = read_rds(pick_content('human_CD8', 'seurat_pheno'))
pheno_df$cellNames = paste(pheno_df$sample_name, pheno_df$cells, sep='#')
pheno_df$cellNames = gsub('^(.*)-[0-9]+$', '\\1-1', pheno_df$cellNames)
compare_vectors(pheno_df$cellNames, proj$cellNames)
proj = proj[proj$cellNames %in% pheno_df$cellNames, ]
dr_df = as.data.frame(getCellColData(proj))
dr_df = dr_df %>% rownames_to_column('cellNames') %>% dplyr::left_join(pheno_df, by = 'cellNames')
rownames(dr_df) = dr_df$cellNames #important to set the rownames!
proj@cellColData = as(dr_df, 'DataFrame')


proj$Tissue = sapply(strsplit(proj$Sample, '_'), '[[', 3)
proj$Donor = sapply(strsplit(proj$Sample, '_'), '[[', 4)
proj$log10_nFrags = log10(proj$nFrags)


# add in metadata from sample sheet ---------------------------------------

dr_df = as.data.frame(getCellColData(proj))
single_cell_csvs <- Sys.glob(sprintf('%sMD*/MD*/outs/singlecell.csv', pconfig$CELLRANGER_OUTPUT_DIR))
names(single_cell_csvs) = basename(dirname(dirname(single_cell_csvs)))
stopifnot(identical(names(inputFiles), names(single_cell_csvs)))
cells_metrics_df <- load_single_cell_csvs(single_cell_csvs, remove_no_barcode = T, return_df = T, full_names = F)
#make barcode+sample combination according to archr cellNames
cells_metrics_df$cellNames = paste(cells_metrics_df$Short.name, cells_metrics_df$barcode, sep='#')
dr_df$cellNames = rownames(dr_df)
stopifnot(all(dr_df$cellNames %in% cells_metrics_df$cellNames))
dr_df = dr_df %>% dplyr::left_join(cells_metrics_df, by = 'cellNames')
dr_df$Short.name = NULL
rownames(dr_df) = dr_df$cellNames
proj@cellColData = as(dr_df, "DataFrame")

#saveArchRProject(proj, outputDirectory='ArchRProject_raw', load=F) #overwrites previous 

# sample number of fragments and barcodes ---------------------------------

dr_df = as.data.frame(getCellColData(proj))
res = dr_df %>% 
   dplyr::select(Sample, nFrags) %>% 
   dplyr::group_by(Sample) %>% 
   dplyr::summarise(n_barcodes=n(), median_n_fragments=as.integer(median(nFrags)))
res = res %>% add_row(Sample = 'total', 
                      n_barcodes = nrow(dr_df), 
                      median_n_fragments = as.integer(median(dr_df$nFrags)))
res

# nr of cells and median fragments for different filters ------------------

ggl = cutoff_heatmap(dr_df_keep, 'Sample', 'nFrags', 'TSSEnrichment', return_df = F, 
                     nFrags_cutoff_seq = seq(1000, 10000, 1000),
                      tss_cutoff_seq = seq(3, 12, 1))
ggl[[1]] #n barcodes
ggl[[2]] #n fragments


#plot results for specific tss cutoff
res_df = cutoff_heatmap(dr_df_keep, 'Sample', 'nFrags', 'TSSEnrichment', return_df = T)
res_plot_df = res_df[res_df$tss_cutoff == 6, ]
res_plot_df = res_plot_df[res_plot_df$Sample != 'total', ]
p1 = ggplot(res_plot_df, aes(x=fragment_cutoff, y=n_barcodes, color=Sample)) + geom_line(size = 1.5) +
  theme_minimal() + scale_color_manual(values = msPickSampleColors(res_plot_df$Sample)) +
  scale_x_continuous(breaks = seq(min(res_plot_df$fragment_cutoff), max(res_plot_df$fragment_cutoff), by = 1000)) + 
  NoLegend()
p2 = ggplot(res_plot_df, aes(x=fragment_cutoff, y= median_n_fragments, color=Sample)) + geom_line(size = 1.5) +
  theme_minimal() + scale_color_manual(values = msPickSampleColors(res_plot_df$Sample)) + 
  scale_x_continuous(breaks = seq(min(res_plot_df$fragment_cutoff), max(res_plot_df$fragment_cutoff), by = 1000))
p3 = cowplot::get_legend(p2)
p2 = p2 + NoLegend()
grid.arrange(p1, p2, p3, ncol=2)

#plot result for specific fragment cutoff
res_plot_df = res_df[res_df$fragment_cutoff == 3000, ]
res_plot_df = res_plot_df[res_plot_df$Sample != 'total', ]
p1 = ggplot(res_plot_df, aes(x=tss_cutoff, y=n_barcodes, color=Sample)) + geom_line(size = 1.5) +
  theme_minimal() + scale_color_manual(values = msPickSampleColors(res_plot_df$Sample)) +
  scale_x_continuous(breaks = seq(min(res_plot_df$tss_cutoff), max(res_plot_df$tss_cutoff), by = 1)) + 
  NoLegend()
p2 = ggplot(res_plot_df, aes(x=tss_cutoff, y= median_n_fragments, color=Sample)) + geom_line(size = 1.5) +
  theme_minimal() + scale_color_manual(values = msPickSampleColors(res_plot_df$Sample)) + 
  scale_x_continuous(breaks = seq(min(res_plot_df$tss_cutoff), max(res_plot_df$tss_cutoff), by = 1))
p3 = cowplot::get_legend(p2)
p2 = p2 + NoLegend()
grid.arrange(p1, p2, p3, ncol=2)


# distribution of fragments -----------------------------------------------

dr_df_temp = dr_df
dr_df_temp$nFrags[dr_df_temp$nFrags > 30000] = 30000
ggplot(dr_df_temp, aes(x=nFrags)) + geom_histogram(bins=300) + facet_wrap(.~Sample) + theme_minimal() 

### fragment size distribution
fsize_df <- as.data.frame(plotFragmentSizes(ArchRProj = proj, returnDF = T))
fsize_df$Tissue = sapply(strsplit(fsize_df$group, '_'), '[[', 3)
p1 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% fragments', color='Sample') + 
  scale_color_manual(values = msPickSampleColors(fsize_df$group))

p2 = ggplot(fsize_df, aes(x = fragmentSize, y = fragmentPercent, color=group)) + 
  geom_line() + theme_minimal() + labs(x='Fragment length (bp)', y='% fragments', color='Sample') + 
  scale_color_manual(values = msPickSampleColors(fsize_df$group)) + facet_wrap(Tissue~., ncol = 1)

pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_fragment_size_distribution.pdf'), width = 8, height = 3)
  p1
dev.off()

pdf(paste0(storeFigPath, 'human_normal_donor11_fragment_size_distribution2.pdf'), width = 7, height = 7)
p2
dev.off()

#tss enrichment
#TSS enrichment profiles should show a clear peak in the center
#and a smaller shoulder peak right-of-center which is caused by the well-positioned +1 nucleosome
tss_df <- as.data.frame(plotTSSEnrichment(ArchRProj = proj, groupBy = 'Sample', returnDF = T))
p2 = ggplot(tss_df, aes(x = x, y = smoothValue, color=group)) + 
  geom_line() + theme_minimal() + 
  scale_color_manual(values = msPickSampleColors(tss_df$group)) + 
  labs(x='Distance from TSS (bp)', y = 'Normalized Insertion Profile')
p2
p2 + facet_wrap(.~group)


# plot QC stats -----------------------------------------------------------

cor(dr_df$nFrags, dr_df$passed_filters) #almost identical, passed_filters is 10-100 fragments higher
table(dr_df$nFrags == dr_df$passed_filters)

fragments_cutoff = pconfig$FRAGMENT_CUTOFF
tss_cutoff = pconfig$TSS_CUTOFF

nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', nFrags_cutoff = 1000, tss_cutoff = 8)

dr_df %>% dplyr::mutate(own_qc_pass = nFrags > fragments_cutoff & TSSEnrichment > tss_cutoff ) %>% 
  dplyr::mutate(cellranger_qc_pass = as.logical(as.integer(as.character(.$is__cell_barcode)))) %>% 
  dplyr::count(own_qc_pass, cellranger_qc_pass)

res_ = dr_df %>% 
  dplyr::select(Sample, nFrags, TSSEnrichment) %>% 
  dplyr::filter(nFrags > fragments_cutoff & TSSEnrichment > tss_cutoff ) 
res = res_ %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(n_barcodes=n(), median_n_fragments=as.integer(median(nFrags)))
res = res %>% add_row(Sample = 'total', 
                      n_barcodes = nrow(res_), 
                      median_n_fragments = as.integer(median(res_$nFrags)))
res

gr_thres_df = dr_df
gr_thres_df$gr_frags_thres = dr_df$nFrags > fragments_cutoff
gr_thres_df$gr_tss_thres = dr_df$TSSEnrichment > tss_cutoff
padd = 0.3
gr_thres_summary_df = gr_thres_df %>% 
  dplyr::select(gr_frags_thres, gr_tss_thres) %>% 
  dplyr::group_by_all() %>%
  dplyr::summarise(n_cells = n())

gr_thres_summary_df$x_pos = c((1-padd) * fragments_cutoff, (1-padd) * fragments_cutoff, (1 + padd) * fragments_cutoff, (1 + padd) * fragments_cutoff ) 
gr_thres_summary_df$y_pos = c((1-padd) * tss_cutoff, (1 + padd) * tss_cutoff, (1-padd) * tss_cutoff, (1 + padd) * tss_cutoff)

#TSS vs nr fragments
dr_df = as.data.frame(getCellColData(proj))
dr_df$log10_nFrags = log10(dr_df$nFrags)
ggplot(dr_df, aes(x=log10_nFrags, y = TSSEnrichment)) + 
  ggpointdensity::geom_pointdensity(size = 0.3) + 
  viridis::scale_color_viridis() + 
  theme_minimal() + 
  geom_vline(xintercept = log10(fragments_cutoff), col = 'red') + 
  geom_hline(yintercept = tss_cutoff, col = 'red') + 
  labs(title=sprintf('Cutoffs at %i fragments and %i TSSEnrichment', fragments_cutoff, tss_cutoff))+
  geom_label(data = gr_thres_summary_df, mapping = aes(x = log10(x_pos), y =  y_pos, label = n_cells), color='red', size = 3, alpha=0.7) 


gr_thres_summary_sample_df = gr_thres_df %>% 
  dplyr::select(Sample, gr_frags_thres, gr_tss_thres) %>%
  dplyr::group_by_all() %>% 
  dplyr::summarise(n_cells = n())
#gr_thres_summary_sample_df$x_pos = rep(c((1-padd) * fragments_cutoff, (1-padd) * fragments_cutoff, (1 + padd) * fragments_cutoff, (1 + padd) * fragments_cutoff ), length(unique(dr_df$Sample)))
#gr_thres_summary_sample_df$y_pos = rep(c((1-padd) * tss_cutoff, (1 + padd) * tss_cutoff, (1-padd) * tss_cutoff, (1 + padd) * tss_cutoff), length(unique(dr_df$Sample)))
gr_thres_summary_sample_df$x_pos = fragments_cutoff
gr_thres_summary_sample_df$y_pos = tss_cutoff

alpha_use=0.7;text_size=3
pdf(paste0(storeFigPath, 'human_normal_donor11_tss_vs_nFrags.pdf'), width = 12, height=9)
pdf(paste0(storeFigPath, 'mouse_imm_cell_atlas_tss_vs_nFrags.pdf'), width = 12, height=9)
#without log, bad visibility of low fragment numbers
# ggplot(dr_df, aes(x=nFrags, y = TSSEnrichment)) +
#   ggpointdensity::geom_pointdensity(size = 0.3) + 
#   viridis::scale_color_viridis() + 
#   facet_wrap(.~Sample) +
#   theme_minimal() + labs(x='fragment count', y = 'TSS score') + 
#   geom_vline(xintercept = fragments_cutoff, col = 'red') + 
#   geom_hline(yintercept = tss_cutoff, col = 'red') + 
#   geom_label(data = gr_thres_summary_sample_df, mapping = aes(x = x_pos, y =  y_pos, label = n_cells), color='red', size = 3, alpha=0.7) 

ggplot(dr_df, aes(x=log10_nFrags, y = TSSEnrichment)) +
  ggpointdensity::geom_pointdensity(size = 0.3) + 
  viridis::scale_color_viridis() + 
  facet_wrap(.~Sample) +
  theme_minimal() + labs(x='log10(fragment count)', y = 'TSS score', color='n neighbours') + 
  geom_vline(xintercept = log10(fragments_cutoff), col = 'red') + 
  geom_hline(yintercept = tss_cutoff, col = 'red') + 
  geom_label(data = gr_thres_summary_sample_df[!gr_thres_summary_sample_df$gr_frags_thres & !gr_thres_summary_sample_df$gr_tss_thres, ], hjust = 1, vjust=1, mapping = aes(x = log10(x_pos), y =  y_pos, label = n_cells), color='red', size = text_size, alpha=alpha_use) + 
  geom_label(data = gr_thres_summary_sample_df[!gr_thres_summary_sample_df$gr_frags_thres & gr_thres_summary_sample_df$gr_tss_thres, ], hjust = 1, vjust = 0, mapping = aes(x = log10(x_pos), y =  y_pos, label = n_cells), color='red', size = text_size, alpha=alpha_use) + 
  geom_label(data = gr_thres_summary_sample_df[gr_thres_summary_sample_df$gr_frags_thres & !gr_thres_summary_sample_df$gr_tss_thres, ], hjust = 0, vjust = 1,mapping = aes(x = log10(x_pos), y =  y_pos, label = n_cells), color='red', size = text_size, alpha=alpha_use) +
  geom_label(data = gr_thres_summary_sample_df[gr_thres_summary_sample_df$gr_frags_thres & gr_thres_summary_sample_df$gr_tss_thres, ], hjust = 0, vjust = 0, mapping = aes(x = log10(x_pos), y =  y_pos, label = n_cells), color='red', size = text_size, alpha=alpha_use) 

dev.off()


nfrags_tss_scatter(dr_df, 'Tissue', 'nFrags', 'TSSEnrichment', 4000, 6)
nfrags_tss_scatter(dr_df, 'Sample', 'nFrags', 'TSSEnrichment', 4000, 8)

### nr fragments

#plot without log for better understanding
p1 = plotGroups(
  ArchRProj = proj, 
  baseSize = 10, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "ridges"
)

# p1 <- plotGroups(
#   ArchRProj = proj, 
#   baseSize = 10, 
#   groupBy = "Sample", 
#   colorBy = "cellColData", 
#   name = "log10(nFrags)",
#   plotAs = "ridges"
# )
p1 = p1 + geom_vline(xintercept = fragments_cutoff, col = 'red') +
  scale_fill_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + NoLegend()
p2 <- plotGroups(
  ArchRProj = proj, 
  baseSize = 10, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2 = p2 + geom_hline(yintercept = log10(fragments_cutoff), col = 'red') + 
  scale_fill_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + 
  scale_color_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + NoLegend()
grid.arrange(p1, p2, ncol=2)

### tss enrichment
p1 <- plotGroups(
  ArchRProj = proj, 
  baseSize = 10, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1 = p1 + geom_vline(xintercept = tss_cutoff, col = 'red') + 
  scale_fill_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + NoLegend()
p2 <- plotGroups(
  ArchRProj = proj, 
  baseSize = 10, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2 = p2 + geom_hline(yintercept = tss_cutoff, col = 'red') + 
  scale_fill_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + 
  scale_color_manual(values = msPickSampleColors(gtools::mixedsort(unique(dr_df$Sample)))) + NoLegend()

grid.arrange(p1, p2, ncol=2)



# comparison to old filter steps ------------------------------------------

#define filter thresholds
unmapped_lowmapq_chimeric_threshold <- 0.2
duplicate_threshold <- 0.9
mito_fraction_threshold <- 0.1

#do the archr filtered cells accord with the previous qc filter steps?
gr_thres_df = gr_thres_df %>% dplyr::filter(gr_frags_thres, gr_tss_thres) %>% 
  dplyr::mutate(poor_mapping = ((unmapped+lowmapq+chimeric) / total) > unmapped_lowmapq_chimeric_threshold,
                high_duplicate_fraction = (duplicate / total) > duplicate_threshold,
                high_mito_fraction = (mitochondrial / total) > mito_fraction_threshold)
gr_thres_df %>% 
  dplyr::count(poor_mapping, high_duplicate_fraction, high_mito_fraction) %>% 
  dplyr::mutate(fraction = round(n / sum(n), 3))

upsetr_tb_list_new = list()
for(i in seq_along(unique(dr_df$Sample))){
  #exclude first row which is a summary across all barcodes
  upsetr_tb <- dr_df[dr_df$Sample == unique(dr_df$Sample)[i], ]
  upsetr_tb$low_fragment_number <- upsetr_tb$nFrags < fragments_cutoff
  upsetr_tb$low_tss <- upsetr_tb$TSSEnrichment < tss_cutoff
  upsetr_tb <- upsetr_tb %>% 
    dplyr::select(low_fragment_number:low_tss) %>% 
    apply(., 2, as.integer) %>% 
    as.data.frame
  upsetr_tb$initial_barcodes <- 1
  upsetr_tb$remaining_barcodes <- 0
  upsetr_tb$remaining_barcodes[apply(upsetr_tb,1,sum)==1] <- 1
  upsetr_tb_list_new[[i]] <- upsetr_tb %>% dplyr::select(initial_barcodes, everything())
}
#pdf(paste0(storeFigPath, fragments_cutoff, 'mouse_normal_fragments_upsetr_barcode_filtering_new.pdf'))
# for(i in seq_along(upsetr_tb_list_new)){
#   message(i)
#   print(UpSetR::upset(upsetr_tb_list_new[[i]], nsets=100,  mb.ratio = c(0.55, 0.45),
#                       sets = rev(colnames(upsetr_tb_list_new[[i]])), keep.order = T, 
#                       mainbar.y.label =  paste0(names(upsetr_tb_list_new)[i], '\nIntersection Size'), 
#                       order.by='freq', set_size.show = TRUE,   set_size.scale_max= nrow(upsetr_tb_list_new[[i]]) +1.5e5 ))
# }
# dev.off()



#apply basic filtering on the raw barcode matrix
upsetr_tb_list_old = list()
for(i in seq_along(unique(dr_df$Sample))){
  #exclude first row which is a summary across all barcodes
  upsetr_tb <- dr_df[dr_df$Sample == unique(dr_df$Sample)[i], ]
  upsetr_tb$low_fragment_number <- upsetr_tb$passed_filters < fragments_cutoff
  upsetr_tb$poor_mapping <- with(upsetr_tb, (unmapped+lowmapq+chimeric) / total > unmapped_lowmapq_chimeric_threshold)
  upsetr_tb$high_duplicate_fraction <- with(upsetr_tb, duplicate / total > duplicate_threshold)
  upsetr_tb$high_mito_fraction <- with(upsetr_tb, mitochondrial / total > mito_fraction_threshold)
  upsetr_tb$low_tss = upsetr_tb$TSSEnrichment < tss_cutoff
  upsetr_tb <- upsetr_tb %>% 
    dplyr::select(low_fragment_number:low_tss) %>% 
    apply(., 2, as.integer) %>% 
    as.data.frame
  upsetr_tb$initial_barcodes <- 1
  upsetr_tb$remaining_barcodes <- 0
  upsetr_tb$remaining_barcodes[apply(upsetr_tb,1,sum)==1] <- 1
  upsetr_tb_list_old[[i]] <- upsetr_tb %>% dplyr::select(initial_barcodes, everything())
}

barcodes_passed_qc_new = sapply(upsetr_tb_list_new, function(x) sum(x$remaining_barcodes))
barcodes_passed_qc_old = sapply(upsetr_tb_list_old, function(x) sum(x$remaining_barcodes))
gg_df = data.frame(Sample=unique(dr_df$Sample), n_cells_new_qc=barcodes_passed_qc_new, n_cells_old_qc = barcodes_passed_qc_old)
gg_df = gg_df %>% dplyr::arrange(Sample)
ggplot(gg_df, aes(x=n_cells_new_qc, y = n_cells_old_qc, color=Sample)) + theme_minimal() + 
  geom_point(size = 3) + geom_abline(intercept = 0, slope=1) + scale_color_manual(values = msPickSampleColors(gg_df$Sample)) 


# downfilter arrow files --------------------------------------------------

#human normal filtered only for frags > 1000 and excluded_reason == 0
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw_excluded_reason_filtered'))

#saveArchRProject(proj, outputDirectory='ArchRProject_raw', load=F) #overwrites previous 
proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_raw'))
#mouse normal: 104552 cells
#human normal: 484783 cells

table(proj$nFrags > pconfig$FRAGMENT_CUTOFF & proj$TSSEnrichment > pconfig$TSS_CUTOFF)
idxPass <- which(proj$nFrags > pconfig$FRAGMENT_CUTOFF & proj$TSSEnrichment > pconfig$TSS_CUTOFF)
cellsPass <- proj$cellNames[idxPass]
proj = proj[cellsPass, ]

#save under new directory and reload the filtered arrow files
#saveArchRProject(proj, outputDirectory='ArchRProject_filtered', load=T, dropCells = T, overwrite = T) #overwrites previous 
#mouse normal: 72766 cells
#mouse nfil3 gfp: 32742 cells
#mouse areg gfp: 22625 cells
#human normal: 79470 cells
#mouse her2 t: 14702 cells

# add doublet scores ------------------------------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered'))
table(proj$Sample)
getAvailableMatrices(proj)
ncol(getMatrixFromArrow(getArrowFiles(proj)[4]))

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
jj_summarize_dataframe(dr_df[, c('snn_harmony_res.1', 'DoubletEnrichment')], summarize_by_vec = dr_df$snn_harmony_res.1) %>% arrange(desc(DoubletEnrichment))

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

proj = archr_dim_reduction(proj, cluster_res = NULL)
#overwrites previous archr project, but does not copy all arrows if already present
#saveArchRProject(proj, outputDirectory='ArchRProject_filtered', overwrite = F) 

# dataset specific filtering ----------------------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered'))
#mouse normal: 72766  cells
#mouse areg gfp: 22625 cells
#mouse nfil3 gfp: 32742 cells
#human normal: 79470 cells

###general
dr_df = jj_get_reduction_coords(proj, 'UMAP')
plotFeatures(reduction=dr_df, meta_features = 'Sample', use_facets = 1, pt.size = 0.5)
plotFeatures(reduction=dr_df, meta_features = 'TSSEnrichment', colorScale = 'viridis', pt.size = 0.5, topCap = 'q99')
###doublet scores
#on complete umap
dr_df$doubletEnrich_above_cutoff = dr_df$DoubletEnrichment > pconfig$DOUBLET_ENRICH_CUTOFF
jj_plot_features(reduction=dr_df, meta_features = 'DoubletEnrichment', pt.size = 0.5, colorScale = 'viridis', cap_top =  'q95')
plotFeatures(reduction=dr_df, meta_features = 'doubletEnrich_above_cutoff', pt.size = 0.5)
dr_df$log10_nFrags = log10(dr_df$nFrags)
plotFeatures(reduction=dr_df, meta_features = 'log10_nFrags', colorScale = 'viridis', pt.size = 0.5)

#on per sample umaps
dr_df[, 1] = dr_df$UMAP_1_ss
dr_df[, 2] = dr_df$UMAP_2_ss
plotFeatures(reduction=dr_df, meta_features = 'DoubletEnrichment', pt.size = 0.5, colorScale = 'viridis', topCap = 'q95', use_facets = 'Sample')
plotFeatures(reduction=dr_df, meta_features = 'log10_nFrags', pt.size = 0.5, colorScale = 'viridis', topCap = 'q95', use_facets = 'Sample')
plotFeatures(reduction=dr_df, meta_features = 'doubletEnrich_above_cutoff', pt.size = 0.5, use_facets = 'Sample')

dr_df %>% dplyr::mutate(doubletEnrich_grThres = DoubletEnrichment > pconfig$DOUBLET_ENRICH_CUTOFF) %>% 
  dplyr::group_by(Sample, doubletEnrich_grThres) %>% 
  dplyr::summarise(ncells = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Sample) %>%  
  dplyr::mutate(fraction = ncells/sum(ncells)) %>% 
  dplyr::arrange(desc(doubletEnrich_grThres)) %>% 
  dplyr::filter(doubletEnrich_grThres)

###excluded reason
table(proj$excluded_reason)
#meaning:
#•	0: barcode is not excluded
#•	1: barcode is excluded because it is a gel bead multiplet
#•	2: barcode is excluded because it is low-targeting
#•	3: barcode is excluded because it is a barcode multiplet
dr_df = get_reduction_coords(proj, 'UMAP')
plotFeatures(reduction=dr_df, meta_features = 'excluded_reason', pt.size = 0.5)
plotFeatures(reduction=dr_df, meta_features = 'excluded_reason', pt.size = 0.5, use_facets = T, pointdensity = T)


###filter cells by excluded reason
idxPass <- which(proj$excluded_reason == 0)
cellsPass <- proj$cellNames[idxPass]
proj = proj[cellsPass, ]
#proj = archr_dim_reduction(proj)

#how many cells are retained after doublet filtering with different thresholds
count_retained_cells = function(drdf, logical_col, sample_col){
  drdf %>%
    dplyr::group_by(!!sym(sample_col), !!sym(logical_col)) %>% 
    dplyr::summarise(ncells = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(!!sym(sample_col)) %>%  
    dplyr::mutate(fraction = ncells/sum(ncells)) %>% 
    dplyr::filter( !!sym(logical_col)) %>% 
    dplyr::select(-!!sym(logical_col))
}

dr_df = get_reduction_coords(proj, 'UMAP')
cell_summary = list()
summary_type = c('ncells', 'fraction')[2]
for(thres in 1:7){
  thres_col = paste0('DE_',thres)
  dr_df[, thres_col] = dr_df$DoubletEnrichment < thres
  cell_summary[[thres_col]]  = pull(count_retained_cells(dr_df, logical_col = thres_col, sample_col = 'Sample'), summary_type)
}
cell_summary_df = data.frame(Sample = pull(count_retained_cells(dr_df, logical_col = thres_col, sample_col = 'Sample'), Sample),
                             do.call(cbind, cell_summary))

#cell_summary_df[nrow(cell_summary_df)+1, ] = c('total', apply(cell_summary_df[, 2:ncol(cell_summary_df)], 2, sum))
cell_summary_df = cell_summary_df %>% mutate_at(vars(-Sample), function(x) round(as.numeric(x),2))
cell_summary_df = tidyr::pivot_longer(cell_summary_df, !Sample, names_to = 'cutoff', values_to = 'ncells')
ggplot(cell_summary_df[cell_summary_df$Sample !='total', ], aes(x=cutoff, y = Sample, fill=ncells)) + 
  geom_bin2d(stat='identity') + geom_text(aes(label=ncells), size=3, color='black') + 
  viridis::scale_fill_viridis() + theme_minimal()  #scale_colour_gradient2(midpoint = get_scale_midpoint(res_df[res_df$Sample !='total', 'n_barcodes']), low = 'white', mid='orange', high = 'black')+




###filter cells by doublet score
proj = proj[proj$DoubletEnrichment < pconfig$DOUBLET_ENRICH_CUTOFF, ]
dr_df = proj@embeddings$UMAP@listData$df
dr_df = dr_df[match(rownames(proj), rownames(dr_df)), ]
dr_df <- cbind(dr_df, proj@cellColData)
dr_df[, 1] = dr_df$UMAP_1_ss
dr_df[, 2] = dr_df$UMAP_2_ss
plotFeatures(reduction=dr_df, 
             meta_features = 'DoubletEnrichment',
             pt.size = 1, colorScale = 'viridis', 
             topCap = 'q95')
plotFeatures(reduction=dr_df, 
             meta_features = 'DoubletEnrichment',
             pt.size = 1, colorScale = 'viridis', 
             topCap = 'q95', use_facets = 'Sample')
plotFeatures(reduction=dr_df, 
             meta_features = 'nFrags',
             pt.size = 1, colorScale = 'viridis', 
             topCap = 'q95')
plotFeatures(reduction=dr_df, 
             meta_features = 'Sample',
             pt.size = 1)
plotFeatures(reduction=dr_df, 
             meta_features = 'Sample',
             pt.size = 1, use_facets = T, pointdensity = T)



#rerun dimensionality reduction, this time with clustering
proj = archr_dim_reduction(proj, cluster_res = NULL)
#proj = archr_clustering(proj)

#save with copying arrow files since cells were excluded
#saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets', load=T, dropCells = T, overwrite = T) 

library(clustree)
dr_df = dr_df[, colnames(dr_df) %in% paste0('Clusters_', seq(0,2, 0.1))]
clustree(dr_df[,1:15], prefix = "Clusters_", show_axis=T)

# # add harmony -------------------------------------------------------------
# 
# proj <- addHarmony(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample",
#   theta=2,
#   lambda=1,
#   sigma=1,
#   do_pca=FALSE,
#   force=TRUE
# )
# 
# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "Harmony",
#   name = "UMAP_Harmony",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine",
#   force=TRUE
# )
# 
# dr_df = get_reduction_coords(proj, 'UMAP_Harmony')
# plotFeatures(reduction=dr_df, 
#              meta_features = c('Donor'),
#              pt.size = 0.5, use_facets = 'Donor', background_cells = T)

# plot metadata, filter out doublets? -----------------------------------------------

(proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets')))

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Tissue", embedding = "UMAP")

dr_df = get_reduction_coords(proj, 'UMAP')

storeFigPath = msMakeAnDir('mouse_areg_clustering', outputDir = '/abi/data2/simonma/projects/imm_cell_atlas/analysis/')
pdf(paste0(storeFigPath, dataset_use, '_clustering_tnkilc_subset.pdf'), width = 12, height = 10)
for(i in gtools::mixedsort(grep('Clusters_[0-9]', colnames(dr_df), value=T))){
  #dr_df[, i] = factor(dr_df[, i], levels = gtools::mixedsort(unique(dr_df[, i])))
  print(plotFeatures(reduction=dr_df, meta_features = i, pt.size = 0.5, cont_or_disc = 'd', custom_colors = msPickSampleColors(gtools::mixedsort(unique(dr_df[, i])))))
}
dev.off()
lowerRes(paste0(storeFigPath, dataset_use, '_clustering.pdf'), remove_original=T)

cluster_df = dr_df[, grepl('Clusters_[0-9]', colnames(dr_df))][, 1:8]
plot_sankey(cluster_df)

# run singler  -------------------------------------------------

#use script run_singler.R to save gene_activity matrix, run singler with singler_classify_cells_bsub.R

#load singler results
singler_df = read_csv(pconfig$SINGLER_RESULTS)
dr_df = get_reduction_coords(proj, 'UMAP')
identical(singler_df$X1, rownames(dr_df))
singler_df = singler_df[match(rownames(dr_df), singler_df$X1), ]
stopifnot(identical(singler_df$X1, rownames(dr_df)))
dr_df$singler_label = singler_df$labels


dr_df$singler_label[dr_df$singler_label %in% names(table(dr_df$singler_label))[table(dr_df$singler_label) < 20]] = 'other'
dr_df$singler_label[dr_df$singler_label %in% c('Eosinophils', 'Basophils', 'Mast cells', 'Tgd', 'Neutrophils')] = 'other'

dr_df$tuning_scores_second = singler_df$tuning.scores.second
#proj@cellColData  = as(dr_df, 'DataFrame')
gg = plotFeatures(reduction=dr_df, meta_features = 'singler_label', pt.size = 0.5, cont_or_disc = 'd' , return_gg_object = T)
gg[[1]]$data$singler_label = as.factor(gg[[1]]$data$singler_label)
pdf(paste0(storeFigPath, dataset_use, 'immgen_label_main_simplified_singler_predictions.pdf'), width = 14, height = 12)
  LabelClusters(gg[[1]], id='singler_label', size = 3, box=T) + NoLegend()
  plotFeatures(reduction=dr_df, meta_features = 'singler_label', pt.size = 0.5, cont_or_disc = 'd', 
               use_facets = 'singler_label', pointdensity = T, background_cells = T)
  plotFeatures(reduction=dr_df, meta_features = c('tuning_scores_second'),  pt.size = 1, colorScale = 'viridis')
dev.off()

proj$tuning_scores_second = singler_df$tuning.scores.second
proj$singler_label = singler_df$labels


# dataset specific filtering ----------------------------------------------

if(dataset_use=='mouse_normal'){
  #cluster 15 only occurs in Skin_78, but not in Skin_77, it is a subset of the big barcode multiplet cluster, which however was 
  #not annotated as multiplet by the cellranger-atac pipeline
  #to compare, load the proj before filtering out the barcode multiplets and get the cells within the barcode multiplet cluster
  #barcodes_bad = fractionsByCoordinates(seurat_obj = dr_df, x.coord = c(-12, 0), y.coord = c(-14, -3.5))
  #then load after filtering and check whether the strange cells belong to this cluster
  #dr_df$barcode_multiplet_cluster = rownames(dr_df) %in% rownames(barcodes_bad)
  #plotFeatures(reduction=dr_df, meta_features = c('barcode_multiplet_cluster'), pt.size = 1)
  proj = proj[!proj$Clusters_0.5 == 'C15', ]
  #saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets', load=T, dropCells = T, overwrite = T) 
}
if(dataset_use=='human'){
  #split by donor
  proj_donor11 = archr_dim_reduction(proj[proj$Donor=='Donor11', ]) #39171 cells
  #found that c2 is most likely still doublets
  proj_donor11 = archr_dim_reduction(proj[proj$Clusters_1.1 != 'C2', ]) #38783 cells
  #saveArchRProject(proj_donor11, outputDirectory='ArchRProject_filtered_no_doublets_donor11', load=T, dropCells = T, overwrite = T) 
  
  proj_donor13 = archr_dim_reduction(proj[proj$Donor=='Donor13', ]) #18674 cells
  #saveArchRProject(proj_donor13, outputDirectory='ArchRProject_filtered_no_doublets_donor13', load=T, dropCells = T, overwrite = T) 
}
if(dataset_use=='mouse_areg_gfp'){
  dr_df = get_reduction_coords(proj, redname = "UMAP")
  cells_exclude = rownames(fractionsByCoordinates(reduction = dr_df, x.coord = c(-2,-1), y.coord = c(-11, -10))) #13 cells
  dr_df$exclude = rownames(dr_df) %in% cells_exclude
  jj_plot_features(reduction = dr_df, meta_features = 'exclude')
  proj = proj[!proj$cellNames %in% cells_exclude, ]
  #saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets_peaks', load=T, dropCells = T, overwrite = T) 
}


# plot qc stats on umap ---------------------------------------------------


dr_df = get_reduction_coords(proj, 'UMAP')
features_plot = c('nMultiFrags','nMonoFrags','nFrags','nDiFrags','TSSEnrichment','ReadsInTSS','ReadsInPromoter','ReadsInBlacklist','PromoterRatio','PassQC','NucleosomeRatio','BlacklistRatio','total','duplicate','chimeric','unmapped','lowmapq','mitochondrial','nonprimary','passed_filters','is__cell_barcode','excluded_reason','TSS_fragments','DNase_sensitive_region_fragments','enhancer_region_fragments','promoter_region_fragments','on_target_fragments','blacklist_region_fragments','peak_region_fragments','peak_region_cutsites','fraction_fragments_overlapping_peaks')
gg = jj_plot_features(reduction=dr_df, meta_features=features_plot,
                      colorScale = 'viridis', cap_top = 'auto', return_gg_object = T)
pdf(paste0(storeFigPath, 'mouse_normal_qc.pdf'), width = 16, height= 10)
print(plot_grobs(gg, grobs=6, cols=3))
dev.off()


# summarise fractions per cluster -----------------------------------------

dr_df = get_reduction_coords(proj)
summary_df = summarise_fractions(dplyr::select(dr_df, 
                                               singler_label, 
                                               Sample 
                                               ),
                                 summarise_by = 'Sample')

ggplot(summary_df, aes(x=Sample, fill = singler_label, y = fraction)) +
  geom_bar(stat='identity', position="fill") + 
  #facet_wrap(.~clinical) +
  coord_flip()  + theme_minimal() + #NoLegend() +
  scale_fill_manual(values=msPickSampleColors(summary_df$singler_label))

summary_table = summary_df %>% dplyr::select(-c(total, ncells)) %>%  pivot_wider(names_from = 'singler_label', values_from = 'fraction')
summary_table[is.na(summary_table)] = 0
summary_table = summary_table %>% dplyr::mutate_at(vars(-Sample), round, 3)
summary_table = summary_table[, c(T, apply(summary_table[, 2:ncol(summary_table)],2, max) >= 0.001)]
summary_table %>% View


# add peaks with macs2 ----------------------------------------------------

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
#254545 peaks

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI_peaks", 
  force = T,
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI_peaks", 
  name = "UMAP_peaks", 
  force = T,
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

cluster_res = seq(0.1, 1, 0.1)
for(i in cluster_res){
    proj <- addClusters(
      input = proj,
      force = T,
      reducedDims = "IterativeLSI_peaks",
      method = "Seurat",
      maxClusters = Inf,
      name = paste0("Clusters_peaks_", i),
      resolution = i
    )
}

dr_df <- get_reduction_coords(proj, 'UMAP_peaks')
plotFeatures(reduction=dr_df, meta_features = 'Sample')

#saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets_peaks', load=F)
#saveArchRProject(proj, outputDirectory = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))

# plot tisTreg signatures (archr conversion to seurat) ---------------------------------------------------------

peak_se = getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
ass_df = assays(peak_se)[[1]]
rownames(ass_df) = paste(seqnames(rowData(peak_se)), ranges(rowData(peak_se)), sep= '-')
#in R 3.6
seurat_atac = Seurat::CreateSeuratObject(counts = ass_df, assay = 'scATAC_raw')
#write_rds(seurat_atac, paste0(bigFilesDir, 'imm_atlas_mouse_normal_peak_seurat.RDS'))

sig_df = read.csv('/abi/data2/simonma/projects/scATAC/analysis/2021-05-10-file_overview_cd8/mouse_imm_atlas_tisTreg_signatures_ATAC.csv')
sig_df$cell_id = readLines('/abi/data2/simonma/projects/scATAC/rownames_temp.txt')

dr_df <- get_reduction_coords(proj, 'UMAP')
dr_df = dr_df %>% rownames_to_column('cell_id') %>% dplyr::select(starts_with('UMAP'), everything()) %>%  dplyr::left_join(sig_df, by = 'cell_id')
rownames(dr_df) = dr_df$cell_id
#proj@cellColData = dr_df %>% dplyr::select(!contains('Iter')) %>% as(., 'DataFrame')
pdf(paste0(storeFigPath, dataset_use, '_tisTreg_signatures.pdf'), width=12, height = 10)
for(i in grep('_sig$', colnames(dr_df), value = T)){
  print(plotFeatures(reduction=dr_df, meta_features = i, pt.size = 1, cont_or_disc = 'c', colorScale = 'viridis', topCap = 'auto'))
}
dev.off()

dr_df$cluster = gsub('C', '', dr_df$Clusters_0.5)
signature_piechart(dr_df[, c('cluster', grep('_sig$', colnames(dr_df), value=T))])

passay = assays(pmat)[[1]]
feat_avail = ArchR::getPeakSet(proj)
#feat_avail = getSeqnames(proj, 'PeakMatrix')
peaknames = paste(seqnames(feat_avail), ranges(feat_avail), sep='-')
test_module = 
proj = addModuleScore(
  ArchRProj = proj,
  useMatrix = 'TileMatrix',
  name = "Module",
  features = list(a=c('chr1-1:1000', 'chr2-1:1000'), b=c('chr1-1:1000', 'chr2-1:1000')),
  nBin = 25,
  nBgd = 100,
  seed = 1
)


# browser tracks ----------------------------------------------------------
#https://www.archrproject.com/bookdown/plotting-marker-peaks-in-archr.html
gene_plot = 'Sdc1'
geneAnnotation = getGeneAnnotation(proj)
region <- geneAnnotation$genes
region <- region[which(tolower(mcols(region)$symbol) %in% 
                         tolower(gene_plot))]
region <- region[order(match(tolower(mcols(region)$symbol), 
                             tolower(gene_plot)))]
region = region + 1000

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters_0.5", useGroups = c('C1', 'C3'),
  region = region,
  #geneSymbol = c("CD19"),
  upstream = 5000, #upstream of tss
  downstream = 10000,  #downstream of tss
  plotSummary = c("bulkTrack", 
                  "scTrack",
                  "featureTrack", 
                  "geneTrack"),
)
grid::grid.draw(p[[1]])





# chromvar ----------------------------------------------------------------

#adds annotation with name 'Motif'
proj = addMotifAnnotations(proj, motifSet = 'homer')
proj <- addBgdPeaks(proj)
#~20 min
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
ggplot() + geom_point(data = plotVarDev$data, aes(x=rank, y=combinedVars, colour=combinedVars)) + 
  geom_text_repel(data = plotVarDev$data[1:20, ], aes(x=rank, y=combinedVars, label=name)) + 
  viridis::scale_colour_viridis() + theme_minimal()

#sig_se = get_archr_mat(proj, 'MotifMatrix')
sig_se = getMatrixFromProject(proj, 'MotifMatrix')
z_score_mat = t(assays(sig_se)[['z']])
dr_df = jj_get_reduction_coords(proj, redname='UMAP')
z_score_mat = z_score_mat[match(rownames(dr_df), rownames(z_score_mat)), ]
stopifnot(identical(rownames(dr_df), rownames(z_score_mat)))
dr_df = as.data.frame(cbind(dr_df, z_score_mat))

gg = jj_plot_features(reduction=dr_df, return_gg_object = T,
                 meta_features = c('BATF.bZIP_20','IRF.BATF.IRF.bZIP_19'),#,'Bach1.bZIP_16', 'Bach2.bZIP_17'),
                 colorScale = 'viridis', cap_top = 'q95', cap_bottom = 'q05')
plot_grobs(gg)

chromvar_summarized <- z_score_mat %>%
  as.data.frame %>% dplyr::mutate(cluster=dr_df$Clusters_0.6) %>% 
  dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean)) %>% 
  as.data.frame %>% column_to_rownames('cluster') %>% t

library(ComplexHeatmap)
h1 <- Heatmap(chromvar_summarized,
              name = 'chromvar deviation', 
              row_names_gp = gpar(fontsize = 3),
              top_annotation = HeatmapAnnotation(cluster = anno_text(x = colnames(chromvar_summarized))))
pdf(paste0(storeFigPath, dataset_use, 'chromvar_heatmap.pdf'), width = 8, height=25)
h1
#+ rowAnnotation(link = anno_mark(at = 1:nrow(chromvar_mat),
#labels = rownames(chromvar_summarized), 
#labels_gp = gpar(fontsize = 10))) 
dev.off()


#proj_t <- addMotifAnnotations(ArchRProj = proj_t, motifSet = "cisbp", name = "Motif")

#https://www.archrproject.com/bookdown/archr-and-custom-deviations.html

if(pconfig$GENOME=='mm10'){
  signature_list = read_rds('/abi/data2/simonma/projects/imm_cell_atlas/analysis/1_res/tisTreg_signature_gr_list.RDS')
  signature_list = signature_list[-c(4)]
  signature_list = lapply(signature_list, makeGRangesFromDataFrame)
  proj <- addPeakAnnotations(ArchRProj = proj, regions = signature_list[1:5], name = "tisTreg", force=T)
  #getPeakAnnotation(proj, 'tisTreg')
  #stores list of 3 anntoations: 
  #$Name: name of the signatures
  #$Positions: GRangesList with Granges for each signature
  #$Matches: SummarizedExperiment with logical sparse assay 'matches' that contains cells in rows and signatures in columns (True = Match of peak and signature)
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
    dplyr::left_join(rownames_to_column(as.data.frame(z_score_mat), 'cell_id'), by = 'cell_id')
  
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
  
  # p <- plotEmbedding(
  #   ArchRProj = proj, 
  #   colorBy = "tisTregMatrix", 
  #   name = signatures, 
  #   embedding = "UMAP",
  #   imputeWeights = getImputeWeights(proj)
  # )
  
  pdf(paste0('/abi/data2/simonma/projects/imm_cell_atlas/analysis/2021-10-11-mouse_normal_tisTreg_signatures/tisTreg_signatures.pdf'), width=12, height=8)
  #p
  gg
  dev.off()
}else{
  signature_list_hg19 = read_csv('/abi/data2/simonma/projects/scATAC/analysis/2020-04-30-tconv_signature_substraction/tisTreg_vs_Tconv_filtering_plot_signatures.csv')
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


#proj <- addArchRAnnotations(ArchRProj = proj, collection = "ATAC")

#saveArchRProject(proj, outputDirectory='ArchRProject_filtered_no_doublets_peaks', load=F)

# seurat module score -----------------------------------------------------

feat_avail = ArchR::getPeakSet(proj)
assay_sp_mat = assays(getMatrixFromProject(proj, 'PeakMatrix', binarize = T))[[1]]
rownames(assay_sp_mat) = paste(seqnames(feat_avail), ranges(feat_avail), sep='-')
z_score_df = AddChromatinModuleOwn(sparse_peak_mat = assay_sp_mat, features_gr_list = signature_list, genome = BSgenome.Mmusculus.UCSC.mm10)
dr_df = dr_df %>% dplyr::select(!ends_with('sig')) %>% dplyr::left_join(rownames_to_column(z_score_df, 'cell_id'), by = 'cell_id')
plotFeatures(reduction=dr_df, meta_features = 'core_tisTreg_sig', colorScale = 'viridis', pt.size = 0.5, topCap = 'auto')


# create bw tracks --------------------------------------------------------

get_cells_per_group <- function(dr_df, cell_id_col, sample_col, meta_col, idents_keep, drop_if_less_than=0.00,
                                cell_id_match = '^.*#([ACGT]+)-.*$', saveDir=NULL){
  #cell_id_match: extract substring from cell_id_col by defining a match group, which should result in strings present in the bam file
  dr_df = as.data.frame(dr_df)
  stopifnot(all(c(cell_id_col, sample_col, meta_col) %in% colnames(dr_df)))
  dr_df_keep <- dr_df[dr_df[, meta_col] %in% idents_keep, c(cell_id_col, sample_col, meta_col), drop=F]
  dr_df_keep$cells <- gsub(cell_id_match, '\\1', dr_df_keep[, cell_id_col])
  print(head(dr_df_keep))
  stopifnot(!anyNA(dr_df_keep$cells))
  if(!is.null(saveDir)){
    saveDir <- gsub('/$', '', saveDir)
    system(sprintf("mkdir -p %s", saveDir))
    dr_df_spl <- split(dr_df_keep, dr_df_keep[, sample_col])
    dr_df_spl <- dr_df_spl[sapply(dr_df_spl, nrow) >= drop_if_less_than * nrow(dr_df_keep)]
    ret_0 <- sapply(seq_along(dr_df_spl), function(x) {
      savename <- sprintf('%s/%s__barcodes.txt', saveDir, names(dr_df_spl)[x])
      message(savename)
      writeLines(dr_df_spl[[x]]$cells, savename)
    })
    writeLines(sprintf("meta_col: %s\nidents_keep: %s\ndrop_if_less_than: %f", meta_col, idents_keep, drop_if_less_than), 
               sprintf("%s/params_used_in_seurat_get_group.txt", saveDir))
    return('All files saved.')
  }

  dr_df_keep %>% 
    dplyr::group_by_(meta_col) %>% 
    dplyr::summarise(n_cells=n()) %>% 
    dplyr::mutate(n_total=sum(n_cells)) %>% 
    dplyr::mutate(fraction=round(n_cells/n_total, 3)) %>% 
    arrange(desc(fraction)) %>% 
    dplyr::filter(fraction >= drop_if_less_than)
}

#store cells from each cluster in a separate textfile by sample (used for Snakefile_sinto_filterbarcodes.py)
for(i in unique(dr_df[, pconfig$CLUSTERS])){
  get_cells_per_group(dr_df, cell_id_col='cellNames', #'cell_id', 
                      sample_col='Sample', meta_col=pconfig$CLUSTERS,
                      idents_keep=i, cell_id_match = '^.*#([ACGT]+)-.*$', drop_if_less_than=0.00, 
                      saveDir = paste0("/abi/data2/simonma/projects/imm_cell_atlas/analysis/2021-11-05-human_normal_donor11/", i))
}
#/abi/data2/simonma/projects/imm_cell_atlas/analysis/2021-10-14-mouse_normal_bigwig/

#store peakset as bedfile
feat_avail = ArchR::getPeakSet(proj)
feat_avail$name = paste(seqnames(feat_avail), ranges(feat_avail), sep='-')
cluster_peaks = feat_avail$name #[names(feat_avail) == 'C1'] #peaks for one cluster only
write_bed_from_string(cluster_peaks, paste0(storeFigPath, 'mouse_normal_peakset.bed'))

proj = loadArchRProject(path = paste0(pconfig$ARCHR_DIR, 'ArchRProject_filtered_no_doublets_peaks', pconfig$DONOR))
