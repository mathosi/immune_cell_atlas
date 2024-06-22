library(openxlsx)
library(Signac)
library(GenomicRanges)
library(signify)
library(jj)
library(tidyverse)

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


# te enrichment in mouse cd4 tisTreg sig ----------------------------------

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
seurat_atac = read_rds(pick_content('mouse_normal_CD4', 'seurat_file'))
mouse_cd4_universe_gr = convert_granges(rownames(GetAssayData(seurat_atac, assay = 'scATAC_raw')))

te_gr = read_te('mm10', return_gr = T)

mouse_cd4_universe_gr$tisTreg_sig_olap = granges_overlap(mouse_cd4_universe_gr, mouse_cd4_tisTreg_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)


### family level comparison
te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')
family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)
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
View(olaps_df)

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

#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv'))
#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_naiveTreg_sig_vs_universe_te_overlap_df.csv'))
olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
ggplot(olaps_sig_df, aes(x = odds_ratio, y = -log10(fdr), colour = log(n_te))) + geom_point() + viridis::scale_color_viridis()

pdf(paste0(storeFigPath, 'tisTreg_sig_te_enrichment_fisher_volcano_plot.pdf'), width = 8, height = 6)
jj_plot_volcano(olaps_sig_df, logfc_column = 'odds_ratio',
                pval_column = 'fdr', symbol_column = 'name',  marker_thres = c(2, Inf), labs_range = c(0,10, 0, 10)) + 
  labs(x = 'Odds ratio', y = '-log10(FDR)', colour = 'Odds ratio') + 
  scale_colour_manual(
    values = c(">-2 & <2" = 'black', ">= 2" = 'blue'),
    breaks = c(">= 2")
  )
dev.off()

olaps_sig_df = olaps_sig_df[olaps_sig_df$odds_ratio > 2, ]

### single te
te_olap_list = list()
for(i in seq_along(olaps_sig_df$name)){
  te_use = olaps_sig_df$name[i]
  te_test_gr = te_gr[te_gr$name == te_use]
  te_olap_res = granges_overlap(mouse_cd4_tisTreg_gr, te_test_gr, olap_direction = 'both',  minOverlap = 0)
  te_olap_list[[te_use]] = mouse_cd4_tisTreg_df[te_olap_res$a_index, ]
  te_olap_list[[te_use]]$pct_peak_olap_with_te = te_olap_res$pct_a_olap_b
  te_olap_list[[te_use]]$pct_te_olap_with_peak = te_olap_res$pct_b_olap_a
}
jj_save_excel(te_olap_list, paste0(storeFigPath, 'te_scATAC_tisTreg_signature_overlaps.xlsx'))

# conserved tisTreg sig te enrichment ------------------------------------------------------

library(openxlsx)
library(Signac)
library(GenomicRanges)
library(signify)
library(jj)
library(tidyverse)
library(Seurat)

cons_df = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-04-16-liftover_minmatch02_2000bp_upstream/liftover_tisTreg_annotation_distance_to_gene_body_and_2000bp_upstream.xlsx')
cons_df = cons_list$peaks_lifted_homologs_dist0_sam
cons_gr = convert_granges(cons_df$mm10_original_peak)
mm_cons_peaks = unique(cons_df$mm10_original_peak[cons_df$comparison_mm10_original == '16_23_versus_0_3_14' ]) #273 unique peaks
all_peaks = rownames(GetAssayData(read_rds('/omics/groups/OE0436/internal/msimon/scATAC/final/seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.RDS'), assay = 'scATAC_raw'))
all_peaks_gr = convert_granges(all_peaks)
mm_cons_gr = all_peaks_gr[ all_peaks %in% mm_cons_peaks ]

te_gr = read_te('mm10', return_gr = T)
te_gr = te_gr[seqnames(te_gr) %in% seqnames(all_peaks_gr)] 
mm_cons_gr = mm_cons_gr[seqnames(mm_cons_gr) %in% seqnames(te_gr)]

#te_names = unique(te_gr$name)
olaps_df = run_fisher(all_peaks_gr, mm_cons_gr, te_gr)

#write_csv(olaps_df, paste0(storeFigPath, 'mouse_conserved_tisTreg_sig_vs_universe_te_overlap_df.csv'))

olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-04-14-fig1to4_additional_plots/mouse_conserved_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
ggplot(olaps_sig_df, aes(x = odds_ratio, y = -log10(fdr), colour = log(n_total_query))) + 
  geom_point() + viridis::scale_color_viridis()
jj_volcano_plot(olaps_sig_df, logfc_column = 'odds_ratio',
                pval_column = 'fdr', symbol_column = 'name', 
                labs_range = c(0,10),
                marker_thres = c(5,10)) + geom_vline(xintercept = 1, colour = 'red') + labs(colour = 'Odds ratio')

### single te
cons_df$TE_olap = ''
for(i in seq(cons_gr)){
  cons_sig_gr = cons_gr[i]
  cons_df[i, 'TE_olap']  = paste(te_gr$name[granges_overlap(te_gr, cons_sig_gr, minOverlap = 0, olap_direction = 'a', return_type = 'logical')], collapse = ', ')
}
write_csv(cons_df, paste0(storeFigPath, 'mouse_conserved_tisTreg_sig_te_annotation.csv'))


# batf te enrichment ------------------------------------------------------

library(openxlsx)
library(Signac)
library(GenomicRanges)
library(signify)
library(jj)
library(tidyverse)

tab7h_peaks = read.xlsx('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/1-s2.0-S1074761319304984-mmc8.xlsx', 
                        sheet = 'Figure 7H', startRow = 3,colNames = T, rowNames = F)
batf_pos_peaks = tab7h_peaks[tab7h_peaks$padj < 0.01 & tab7h_peaks$log2FoldChange < -2, ] #10529 regions that are only present in batf+ cells
batf_neg_peaks = tab7h_peaks[tab7h_peaks$padj < 0.01 & tab7h_peaks$log2FoldChange > 2, ] #1260 regions that are only present in batf- cells
tab7h_peaks$da_batf_pos = tab7h_peaks$index %in% batf_pos_peaks$index
tab7h_peaks_gr = convert_granges(tab7h_peaks$index)
tab7h_peaks_gr$index = tab7h_peaks$index
tab7h_peaks_gr$da_batf_pos = tab7h_peaks$da_batf_pos
#tab7h_peaks_gr$da_batf_pos = factor(tab7h_peaks$da_batf_pos, levels = c('FALSE','TRUE'))

batf_pos_gr = tab7h_peaks_gr[tab7h_peaks_gr$da_batf_pos]
#batf_pos_gr = StringToGRanges(batf_pos_peaks$index, sep = c(':','-'))

te_gr = read_te('mm10', return_gr = T)
te_gr = te_gr[seqnames(te_gr) %in% seqnames(batf_pos_gr)] 
batf_pos_gr = batf_pos_gr[seqnames(batf_pos_gr) %in% seqnames(te_gr)] #nothing filtered out

te_names = unique(te_gr$name)


olaps_df = jj_initialize_df(5, length(te_names), init = 0,
                            col.names = c('both','te_only','batf_sig_only','none', 'pval'), row.names=te_names)
for(i in 1:length(te_names)){
  if(i %% 25 == 0) message(i, '/', length(te_names))
  te_test_gr = te_gr[te_gr$name == te_names[i]]
  tab7h_peaks_gr$te_olap = granges_overlap(tab7h_peaks_gr, te_test_gr, return_type = 'logical', olap_direction = 'a',  minOverlap = 0)
  if(is.null(tab7h_peaks_gr$te_olap)){
    tab7h_peaks_gr$te_olap = FALSE
  }
  present_both = sum(tab7h_peaks_gr$da_batf_pos & tab7h_peaks_gr$te_olap)
  present_te_only = sum(!tab7h_peaks_gr$da_batf_pos & tab7h_peaks_gr$te_olap)
  present_batf_only = sum(tab7h_peaks_gr$da_batf_pos & !tab7h_peaks_gr$te_olap)
  present_none = sum(!tab7h_peaks_gr$da_batf_pos & !tab7h_peaks_gr$te_olap)
  
  cont_df = data.frame('te_yes' = c(present_both, present_te_only),
                       'te_no' = c(present_batf_only, present_none),
                       row.names = c('batf_sig_yes', 'batf_sig_no'))
  
  test <- fisher.test(cont_df)
  olaps_df[i,] = c(present_both, present_te_only, present_batf_only, present_none, test$p.value)
}

olaps_df$fraction_batf_sig_olap = olaps_df$both / length(batf_pos_gr)
olaps_df$fraction_other_peaks_olap = olaps_df$te_only / (length(tab7h_peaks_gr) - length(batf_pos_gr))
olaps_df$relative_enrichment = olaps_df$fraction_batf_sig_olap / olaps_df$fraction_other_peaks_olap
olaps_df = olaps_df[order(olaps_df$relative_enrichment, decreasing = T), ]
olaps_df$odds_ratio = with(olaps_df, (both * none) / (batf_sig_only * te_only))
View(olaps_df)

te_n_df_mouse = as.data.frame(table(te_gr$name))
colnames(te_n_df_mouse) = c('name', 'n_te')
olaps_df = olaps_df %>% 
  rownames_to_column('name') %>% 
  dplyr::left_join(te_n_df_mouse, by = 'name') %>% 
  dplyr::mutate(fraction_olap_both = both / n_te)
olaps_df$fdr = p.adjust(olaps_df$pval, method = 'fdr')
olaps_df = olaps_df %>% 
  relocate(name, n_te, both, te_only, batf_sig_only, none, pval, fdr, relative_enrichment) %>% 
  arrange(fdr)

#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv'))
#/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-29-te_analysis
olaps_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
olaps_sig_df = olaps_df[olaps_df$fdr < 0.05, ] 
ggplot(olaps_sig_df, aes(x = relative_enrichment, y = -log10(fdr), colour = log(n_te))) + geom_point() + viridis::scale_color_viridis()
jj_volcano_plot(olaps_sig_df, logfc_column = 'relative_enrichment',
                pval_column = 'fdr', symbol_column = 'name', labs_range = c(1,5, 0, 10), marker_thres = c(2, 5))

olaps_sig_df = olaps_sig_df[olaps_sig_df$relative_enrichment > 2, ]

olaps_tisTreg_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-03-06-fig1_to_4_additional_plots/mouse_cd4_tisTreg_sig_vs_universe_te_overlap_df.csv')
plot_df = olaps_df %>% dplyr::full_join(olaps_tisTreg_df, by = 'name', suffix = c('_batf', '_tisTreg'))
plot_df[is.na(plot_df)] = 0
jj_fc_fc_plot(plot_df, logfc_column1 = 'odds_ratio_tisTreg', logfc_column2 = 'odds_ratio_batf',
                 symbol_column = 'name', marker_thres = 2, labs_range = c(0,8), use_text = T)
plot_df$fdr_batf = -log10(plot_df$fdr_batf)
plot_df$fdr_tisTreg = -log10(plot_df$fdr_tisTreg)
jj_fc_fc_plot(plot_df, logfc_column1 = 'fdr_tisTreg', logfc_column2 = 'fdr_batf',
              symbol_column = 'name', marker_thres = 2, labs_range = c(0,10), use_text = T) + labs(x = '-log10(fdr tisTreg)', y = '-log10(fdr BATF)')


olaps_sig_df = olaps_df[olaps_df$fdr < 0.05 & olaps_df$both >= 10, ] 
ggplot(mapping = aes(x = odds_ratio_batf, y = odds_ratio_tisTreg)) + 
  geom_point(data = plot_df, aes(size=log10(both_tisTreg), colour = -log10(fdr_tisTreg))) +
  geom_text_repel(data = plot_df[(plot_df$odds_ratio_batf > 2.5 | plot_df$odds_ratio_tisTreg > 2.5) & plot_df$both_batf >=5 & plot_df$both_tisTreg >=5, ],
                  aes(label = name)) + 
  geom_hline(yintercept = 1, colour='red2') + geom_vline(xintercept = 1, colour='red2') + 
  viridis::scale_colour_viridis() + theme_minimal()

# te list for the mouse scATAC tisTreg signature -------------------

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
te_gr = read_te('mm10', return_gr = T)
names(te_gr) = te_gr$name
te_gr$te_width = width(te_gr)

olap_df = granges_overlap(mouse_cd4_tisTreg_gr, te_gr, minOverlap = 0)
olap_df$pct_b_olap_a = round(olap_df$pct_b_olap_a, 3)
olap_df$te_width = te_gr$te_width[as.integer(olap_df$b_index)]
olap_summary_df = olap_df %>% dplyr::group_by(a_name) %>%
  dplyr::summarise(TE_overlap = paste(b_name, collapse = ', '),
                   TE_fraction_overlap = paste(pct_b_olap_a, collapse = ', '),
                   TE_region_size = paste(te_width, collapse = ', '))
mouse_cd4_tisTreg_df = mouse_cd4_tisTreg_df %>% dplyr::left_join(olap_summary_df, by = c('feature' = 'a_name'))
mouse_cd4_tisTreg_df$peak_size = width(convert_granges(mouse_cd4_tisTreg_df$feature))
#write_csv(mouse_cd4_tisTreg_df, paste0(storeFigPath, 'mouse_cd4_scATAC_tisTreg_16_23_versus_0_3_14_TE_olap.csv'))


# get sequences -----------------------------------------------------------

library(Biostrings)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)

target_regions = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/3_mail/Target regions.xlsx')[[1]]
target_regions = jj_load_excel('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/3_mail/Target regions_3.xlsx')[[1]]
crispr_gr = convert_granges(gsub(':','-', gsub(',','', target_regions$CRISPRa.Target.Region)))

crispr_gr = corrected_regions_gr
target_regions = crispr_use

crispr_seqs <- Views(BSgenome.Mmusculus.UCSC.mm10, crispr_gr)
show(crispr_seqs)
crispr_seqs = as.data.frame(crispr_seqs)
target_regions$crispr_dna_seq = crispr_seqs$dna
#write_csv(target_regions, paste0(storeFigPath, 'crispr_dna_sequences.csv'))

myseq=as.data.frame(Views(BSgenome.Mmusculus.UCSC.mm10, convert_granges('chr19-40599442-40754525')))
myseq$dna
writeLines(myseq$dna, paste0(storeFigPath, 'batf_sequence.txt'))


# homer on TE part of mouse tisTreg signature -----------------------------
library(signify)

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
mouse_cd4_tisTreg_gr$width = width(mouse_cd4_tisTreg_gr)
sum(mouse_cd4_tisTreg_gr$width) #11446376
summary(mouse_cd4_tisTreg_gr$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 230.0   751.5  1145.0  1495.3  1796.0 31678.0 

length(mouse_cd4_tisTreg_gr)
te_gr = read_te('mm10', return_gr = T)

te_tisTreg_overlap = granges_overlap(mouse_cd4_tisTreg_gr, te_gr, minOverlap = 0)
length(unique(te_tisTreg_overlap$a_index)) #5701/7655 peaks overlap with a TE

#te_tisTreg_overlap$a_width = width(convert_granges(te_tisTreg_overlap$a_name))
#te_tisTreg_overlap$b_width = width(convert_granges(te_tisTreg_overlap$b_name))

#require 50% overlap in at least one direction (14534 TEs)
te_olaps = unique(te_tisTreg_overlap$b_name[te_tisTreg_overlap$pct_a_olap_b >= 0.5 | te_tisTreg_overlap$pct_b_olap_a >= 0.5])
sort(table(te_gr$name[as.integer(te_tisTreg_overlap$b_index)]), decreasing = T)
te_olap_gr = convert_granges(te_olaps)
te_olap_gr = GenomicRanges::reduce(te_olap_gr) #make nonoverlapping granges
jj_save_bed(te_olap_gr, file_name = paste0(storeFigPath, 'mouse_tisTreg_overlapping_tes.bed'))

te_olap_gr$width = width(te_olap_gr)
sum(te_olap_gr$width) #2529154
summary(te_olap_gr$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.0   105.0   147.0   193.5   215.0  7086.0 

#2529154 / 11446376 = 0.2209567 -> 22% of the original width make up the te width of the tisTreg signature

#get the nonoverlapping part
overlap_tisTreg_sig_gr = convert_granges(te_tisTreg_overlap$a_name)
overlap_te_gr = convert_granges(te_tisTreg_overlap$b_name)
nonolap_region_list = list()
for(i in seq_along(overlap_tisTreg_sig_gr)){
  if(i %% 50 == 0) message(i)
  d_gr = disjoin(c(overlap_tisTreg_sig_gr[i], overlap_te_gr[i]))
  nonolap_region_list[[i]] = d_gr[!granges_overlap(d_gr, overlap_te_gr[i], minOverlap = 0, return_type = 'logical')]
}
gr_list = list(tisTreg = overlap_tisTreg_sig_gr[i], te = overlap_te_gr[i], disjoin = d_gr, disjoin_keep = nonolap_region_list[[i]])
plot_granges(gr_list)
nonolap_gr = Reduce(c, nonolap_region_list)
jj_save_bed(nonolap_gr, file_name = paste0(storeFigPath, 'mouse_tisTreg_no_overlapping_tes.bed'))


### naive Treg signature
mouse_cd4_naive_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '0_3_14_versus_16_23', ]
mouse_cd4_naive_gr = convert_granges(mouse_cd4_naive_df$feature)
length(mouse_cd4_naive_gr)

te_naive_treg_overlap = granges_overlap(mouse_cd4_naive_gr, te_gr, minOverlap = 0)
length(unique(te_naive_treg_overlap$a_index)) #5649/7299 peaks overlap with a TE

#require 50% overlap in at least one direction (14534 TEs)
te_naive_olaps = unique(te_naive_treg_overlap$b_name[te_naive_treg_overlap$pct_a_olap_b >= 0.5 | te_naive_treg_overlap$pct_b_olap_a >= 0.5])
sort(table(te_gr$name[as.integer(te_naive_treg_overlap$b_index)]), decreasing = T)
te_olap_gr = convert_granges(te_naive_olaps)
te_olap_gr = GenomicRanges::reduce(te_olap_gr) #make nonoverlapping granges
jj_save_bed(te_olap_gr, file_name = paste0(storeFigPath, 'mouse_naive_Treg_overlapping_tes.bed'))


#get the nonoverlapping part
overlap_naiveTreg_sig_gr = convert_granges(te_naive_treg_overlap$a_name)
overlap_naive_te_gr = convert_granges(te_naive_treg_overlap$b_name)
nonolap_region_list = list()
for(i in seq_along(overlap_naiveTreg_sig_gr)){
  if(i %% 50 == 0) message(i)
  d_gr = disjoin(c(overlap_naiveTreg_sig_gr[i], overlap_naive_te_gr[i]))
  nonolap_region_list[[i]] = d_gr[!granges_overlap(d_gr, overlap_naive_te_gr[i], minOverlap = 0, return_type = 'logical')]
}
gr_list = list(tisTreg = overlap_tisTreg_sig_gr[i], te = overlap_naive_te_gr[i], disjoin = d_gr, disjoin_keep = nonolap_region_list[[i]])
plot_granges(gr_list)
nonolap_naive_gr = Reduce(c, nonolap_region_list)
jj_save_bed(nonolap_naive_gr, file_name = paste0(storeFigPath, 'mouse_naiveTreg_no_overlapping_tes.bed'))

##
olaps_df = run_fisher(universe_gr = convert_granges(mouse_cd4_scATAC_sig_df$feature), 
                      universe_subset_gr =  mouse_cd4_tisTreg_gr, 
                      query_gr = te_gr)
#write_csv(olaps_df, paste0(storeFigPath, 'mouse_cd4_tisTreg_sig_vs_naive_treg_sig_te_fisher_results.csv'))

View(olaps_df[olaps_df$fdr < 0.3, ])

#fig 6a
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-03-peripheral_markers/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = T)
homer_res2 = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-04-peripheral_markers/'
homer_subsets2 = list.dirs(homer_res2, recursive = F, full.names = T)
homer_subsets = c(homer_subsets, homer_subsets2)

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
pdf(paste0(storeFigPath, 'mouse_tisTreg_sig_te_olap_nonolap_subsets_homer_heatmap.pdf'),  width = 5, height = 10)
gg + domain_legend #+ plot_layout(widths = c(3, 1))
dev.off()





joined_df = dplyr::full_join(homer_list$mouse_tisTreg_overlapping_tes, homer_list$mouse_tisTreg_no_overlapping_tes, by = 'Motif Name', suffix= c('_te', '_no_te'))
joined_df$name_use = sapply(base::strsplit(joined_df$`Motif Name`, split = '/'), '[[', 1)
joined_df$enrichment_te[is.na(joined_df$enrichment_te)] = 0.9
joined_df$enrichment_no_te[is.na(joined_df$enrichment_no_te)] = 0.9
ggplot(joined_df, aes(x = enrichment_te, y = enrichment_no_te)) + geom_point() + geom_text_repel(aes(label = name_use)) + 
  coord_fixed() + geom_abline(intercept = 0, slope = 1)

joined_df$rank_te[is.na(joined_df$rank_te)] = nrow(joined_df)+1
joined_df$rank_no_te[is.na(joined_df$rank_no_te)] = nrow(joined_df)+1
ggplot(joined_df, aes(x = rank_te, y = rank_no_te)) + geom_point() + geom_text_repel(aes(label = name_use), max.overlaps = 1000) + 
  coord_fixed() + geom_abline(intercept = 0, slope = 1)

joined_df$`Log P-value_te`[is.na(joined_df$`Log P-value_te`)] = 0
joined_df$`Log P-value_no_te`[is.na(joined_df$`Log P-value_no_te`)] = 0
joined_df$`Log P-value_te`[joined_df$`Log P-value_te` < -1000] = -1000
joined_df$`Log P-value_no_te`[joined_df$`Log P-value_no_te` < -1000] = -1000
ggplot(joined_df, aes(x = `Log P-value_te`, y = `Log P-value_no_te`)) + geom_point() + geom_text_repel(aes(label = name_use)) + 
  coord_fixed() + geom_abline(intercept = 0, slope = 1)


# parse homer de novo motif results for tisTreg signature subsets ---------

#use R 4.2.0
library(ggplot2)
library(ggseqlogo)
library(flextable)
library(patchwork)
set_flextable_defaults(font.family = '') #standard font not found when saving as pdf


for(i in homer_subsets){
  pwm_list = pwm_gg_list = match_df_list = match_df_gg = combined_gg = list()
  for(j in seq(10)){
    similar_motifs = readLines(sprintf('%s/homerResults/motif%s.info.html',i, j))
    #get the statistics for the de novo motif
    motif_pval = gsub('.*(1e-[0-9]+).*', '\\1', grep('<TD>p-value:</TD>', similar_motifs, value = T))
    motif_target = gsub('.*?([0-9]+\\.[0-9]+\\%).*', '\\1', grep('Percentage of Target Sequences with motif', similar_motifs, value = T))
    motif_bg = gsub('.*?([0-9]+\\.[0-9]+\\%).*', '\\1', grep('Percentage of Background Sequences with motif', similar_motifs, value = T))
    title_use = sprintf('Rank: %i, P-value: %s, Target: %s, Background: %s', j, motif_pval, motif_target, motif_bg)
    #rows are the positions specific probabilities for each nucleotide (A/C/G/T)
    #read the pwm
    pwm_list[[j]] = t(read.delim(sprintf('%s/homerResults/motif%s.motif',i, j),
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
  pdf(paste0(storeFigPath, basename(i), '_top5_de_novo.pdf'), width = 9, height = 9)
  print(all_gg)
  dev.off()
} 

# homer on TE part of mouse tisTreg signature per family analysis -----------------------------
library(signify)

mouse_cd4_scATAC_sig_df = read_csv('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-03-16-1/diff_results_16_23_versus_0_3_14_seurat_mouse_normal_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar_clustered.csv')
mouse_cd4_tisTreg_df = mouse_cd4_scATAC_sig_df[mouse_cd4_scATAC_sig_df$comparison == '16_23_versus_0_3_14', ]
mouse_cd4_tisTreg_gr = convert_granges(mouse_cd4_tisTreg_df$feature)
mouse_cd4_tisTreg_gr$width = width(mouse_cd4_tisTreg_gr)
sum(mouse_cd4_tisTreg_gr$width) #11446376
summary(mouse_cd4_tisTreg_gr$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 230.0   751.5  1145.0  1495.3  1796.0 31678.0 

length(mouse_cd4_tisTreg_gr)
te_gr = read_te('mm10', return_gr = T)


te_bed_annot = read_tsv('/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/data/2_download/mm10.te.bed', col_names = F)
family_df = te_bed_annot[, c('X5', 'X6', 'X7')] %>% .[!duplicated(.),]
colnames(family_df)  = c('TE', 'Class', 'Family')

family_df$Class = gsub('\\?', '', family_df$Class)
family_df$Family = gsub('\\?', '', family_df$Family)

mouse_cd4_tisTreg_gr@elementMetadata$width = NULL
te_tisTreg_overlap = granges_overlap(mouse_cd4_tisTreg_gr, te_gr, minOverlap = 0)
length(unique(te_tisTreg_overlap$a_index)) #5701/7655 peaks overlap with a TE

#te_tisTreg_overlap$a_width = width(convert_granges(te_tisTreg_overlap$a_name))
#te_tisTreg_overlap$b_width = width(convert_granges(te_tisTreg_overlap$b_name))

#require 50% overlap in at least one direction (14534 TEs)
te_olaps = unique(te_tisTreg_overlap$b_name[te_tisTreg_overlap$pct_a_olap_b >= 0.5 | te_tisTreg_overlap$pct_b_olap_a >= 0.5])
sort(table(te_gr$name[as.integer(te_tisTreg_overlap$b_index)]), decreasing = T)
te_olap_gr = convert_granges(te_olaps)
te_olap_gr$name = te_gr$name[match(te_olaps, convert_granges(te_gr))]
te_olap_gr$family = family_df$Family[match(te_olap_gr$name, family_df$TE)]
sort(table(te_olap_gr$family))
te_families_use =  names(table(te_olap_gr$family))[table(te_olap_gr$family) > 100]
#REMEMBER: -mask needs to be removed from homer call to not filter out the majority of these regions!
for(i in te_families_use){
  te_olap_gr_use = GenomicRanges::reduce(te_olap_gr[te_olap_gr$family == i]) #make nonoverlapping granges
  storeFigPath = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-30-tisTreg_TE_overlap_per_family_subsets'
  jj_save_bed(te_olap_gr_use, file_name = sprintf('%s/mouse_tisTreg_overlapping_tes_%s_%i.bed', storeFigPath, i, length(te_olap_gr_use)))
}

## plot results
#enrichment cutoff increased from 1.5 to 2
homer_res = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/analysis/2023-05-30-tisTreg_TE_overlap_per_family_subsets/'
homer_subsets = list.dirs(homer_res, recursive = F, full.names = T)

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
domain_cols = jj_get_colours(mapping_df$domain, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/imm_cell_atlas/scripts/colour_map.csv', comment_char = '$')

homer_list = list()
for(i in homer_subsets){
  homer_df  = read_tsv(paste0(i, '/knownResults.txt'))
  homer_df$domain = mapping_df$domain[match(homer_df$`Motif Name`, mapping_df$homer_name)]
  homer_df$domain_col = domain_cols[match(homer_df$domain, names(domain_cols))]
  homer_df$motif = homer_df$`Motif Name`
  homer_df$enrichment = as.numeric(gsub('%','', homer_df$`% of Target Sequences with Motif`)) /  as.numeric(gsub('%','', homer_df$`% of Background Sequences with Motif`))
  homer_df$rank = 1:nrow(homer_df)
  if(nrow(homer_df)>20){
    homer_df$norm_rank = c(seq(from = 2, to = .2, length.out = 20), rep(0.2, nrow(homer_df) -20))
  }else{
    homer_df$norm_rank = head(seq(from = 2, to = .2, length.out = 20), nrow(homer_df))
  }
  #1= size 2, 20+ = size 0.2
  homer_df$name = basename(i)
  homer_df = as.data.frame(homer_df)
  homer_list[[basename(i)]] = homer_df[homer_df$`q-value (Benjamini)` < 0.001 & homer_df$enrichment > 3, ]
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
levels_use = gsea_res_comb %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(enr_sum = sum(enrichment)) %>%
  dplyr::arrange( enr_sum) %>%
  pull(motif)
gsea_res_comb$name = gsub('mouse_tisTreg_overlapping_tes_', '', gsea_res_comb$name)

gsea_res_comb_wide = pivot_wider(gsea_res_comb, id_cols = c('name', 'motif'), values_from = 'enrichment',  values_fill = 1, values_fn = mean) %>% 
  as.data.frame %>% column_to_rownames('motif') %>% as.matrix
gsea_res_comb_wide = jj_plot_heatmap(gsea_res_comb_wide, features_use = rownames(gsea_res_comb_wide), group_vec = colnames(gsea_res_comb_wide), return_matrix = T)
gsea_res_comb$name = factor(gsea_res_comb$name, levels = rownames(gsea_res_comb_wide))
gsea_res_comb$motif = factor(gsea_res_comb$motif, levels = colnames(gsea_res_comb_wide))

gg = ggplot(gsea_res_comb, aes(x = name, y = motif, fill = enrichment)) +
  geom_tile(aes(colour = "black")) +
  scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + theme_minimal() + coord_flip() +
  scale_fill_gradientn(colours = paletteContinuous(set = 'comet', n=100)) +
  labs(x = 'Signature subset', y = '', fill = 'Enrichment') +  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1))  + 
  annotate(
    xmin = 12,
    xmax = 12.5,
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
#S7I TE family olapping tisTreg signautre TF homer heatmap
pdf(paste0(storeFigPath, 'mouse_tisTreg_sig_te_family_olap_qval0001_enr3_homer_heatmap.pdf'),  width = 13, height = 5)
gg + domain_legend + plot_layout(widths = c(12, 1))
dev.off()

