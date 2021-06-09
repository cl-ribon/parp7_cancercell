## Cancer Cell Manuscript Comp Bio Data, Code, Figures

### Analysis to define the PARP7 inhibition response groups from large cell panel screen using multiple PARP7 inhibitors and PARP1 inhibitor, niraparib.
*  Large cell panel screen using PARP7 inhibitors (3 PARP7 inhibitors used in paper), where ~50% showed growth response and the other half did not show a response.
*  Source data: data/parp7.horizon_cell_line_response_metrics_wdf.publication.tsv
*  Code: parp7_cellpanel_screen_analysis.R


### PARP7, Niraparib Waterfall plot (Fig5A)
*  figs/Fig5A.parp7_compound_screened_cell_lines.horizon.niraparib_pub.waterfall.pdf
*  Source data: data/parp7.horizon_cell_line_response_metrics_wdf_response_groups.publication.tsv
*  Code: parp7_cellpanel_screen_analysis.R


### PARP7 inhibitor basal gene expression predictive features analysis (Analysis for Fig5B)
*  Description: Using the PARP7 inhibition response groups
*  Source data: parp7_cell_line_sensitivity_basal_dge.deseq.RData (Note: requires CCLE expression data)
*  Code: parp7_cellscreen_basalexp_analysis.R


### PARP7 inhibitor basal gene expression predictive features analysis plot (Fig 5B)
*  Source data: parp7_cell_line_sensitivity_basal_dge.deseq.RData
*  Code: parp7_cellscreen_basalexp_analysis.R


### PARP7 expression vs PARP7i response group colored by published ISG expression score (Fig 5C)
*  Source data: parp7_cell_line_biomarker_data.tsv
*  Code: parp7_biomarker_analysis.R


### Fig 5F
1. RNA-seq volcano: figs/Fig5F_FigS7.scatter.parp7_trt_ctrl.merged_volcano.pub.pdf
2. GSEA barplot: figs/Fig5F_FigS7.gse_barplot.parp7_trt_ctrl.NCIH1373-NCIH647.24hr.gse.barplot.pdf - geneset enrichment barplot with select genesets highlighted. Contains both NCI-H1373 and NCI-H647, only NCI-H647 barplot is used in this figure.
3. Predictive vs. Post-dose gene expression venn
*  figs/FigS7.raw_venn.parp7_cell_line_pred_pd_exp.ncih647.venn.pdf - raw venn diagram generated from R script ()
*  figs/FigS7.formatted_venn.parp7_cell_line_pred_pd_exp.ncih647.venn.png - formatted venn using raw venn diagram and additional labels and graphics in Powerpoint.

Source data:
*  data/NCIH1373_NCIH647.6_24hr.dge_list.RData
*  data/NCIH1373_NCIH647.6_24hr.gse_list.RData
*  data/parp7_cell_line_sensitivity_basal_dge.deseq.RData
Code: parp7_ncih647_ncih1373_rnaesq.R

### Fig S7
1. RNA-seq volcano: figs/Fig5F_FigS7.scatter.parp7_trt_ctrl.merged_volcano.pub.pdf
2. GSEA barplot:
*  figs/Fig5F_FigS7.gse_barplot.parp7_trt_ctrl.NCIH1373-NCIH647.24hr.gse.barplot.pdf
*  Geneset enrichment barplot with select genesets highlighted. Contains both NCI-H1373 and NCI-H647, only NCI-H647 barplot is used in this figure.
3. Predictive vs. Post-dose gene expression venn
*  figs/FigS7.raw_venn.parp7_cell_line_pred_pd_exp.ncih647.venn.pdf - raw venn diagram generated from R script (parp7_ncih647_ncih1373_rnaseq.R)
*  figs/FigS7.formatted_venn.parp7_cell_line_pred_pd_exp.ncih647.venn.png - formatted venn using raw venn diagram and additional labels and graphics in Powerpoint.

Source data:
*  data/NCIH1373_NCIH647.6_24hr.dge_list.RData
*  data/NCIH1373_NCIH647.6_24hr.gse_list.RData
*  data/parp7_cell_line_sensitivity_basal_dge.deseq.RData
Code: parp7_ncih647_ncih1373_rnaesq.R

### Fig 6A NCI-H1373 CRISPRi/a scatter
*  See pptx with formatted figured to match the additions to the figure that are needed (Fig6A.raw_plot.crispria_phenotype_scores.innate_immune_response_scatter.pdf)
*  Source data: parp7_ncih1373_crispria_data.RData
*  Code: parp7_ncih1373_crispria_analysis.R

### Fig 6B NCI-H1373 CRISPRi/a geneset enrichment barplot
*  See pptx with formatted figured to match the additions to the figure that are needed (Fig6B.raw_plot.crispria_fgsea_all.bar.pdf)
*  Source data: parp7_ncih1373_crispria_data.RData
*  Code: parp7_ncih1373_crispria_analysis.R

### Fig 6C CRISPR array screen waterfall plot of RRA scores by gene
*  figs/Fig6C.ctg_aggregate_intx_stats_v2_v3.rra_bar.pdf
*  Source data: data/ctg_aggregate_intx_stats_v2_v3.Rdata
*  Code: parp7_ctg_crispr_array_analysis.R

### Fig 6D Heatmap of merged CRISPRi, CRISPR array, Substrate ID, Synthetic lethality direction
*  figs/Fig6D.parp7_crispr_subid_merged.chm.pdf
*  Source data: data/invitro_substrate_id.RData, data/parp7_crispr_array_pstat1_aggregate_intx_stats_v2_v3.Rdata, data/parp7_crispr_array_ctg_aggregate_intx_stats_v2_v3.Rdata, data/parp7_ncih1373_crispria_data.RData
*  Code: parp7_crispr_subid_merged_analysis.R

### Fig 7A Heatmap of Preclinical PARP7i gene expression measured by nanoString
*  figs/Fig7A.parp7_nanostring_100gp.foldchanges.RBN2397.20191023.pub.chm.pdf
*  Source data: data/parp7_nanostring_100gp.RData
*  Code: parp7_preclinical_nanostring_analysis.R

### Fig S1 TCGA and CCLE plots for PARP7 biomarkers in select indications
*  figs/FigS1A.parp7_tcga_cna_freq.bar.pdf
*  figs/FigS1B.parp7_tcga_cna_exp.pt.pdf
*  figs/FigS1C.parp7_tcga_exp_isg_cor.pdf
*  figs/FigS1D.parp7_parp1_cell_line_crispr_dep.scatter.pdf
*  figs/FigS1E.parp7_cell_line_crispr_dep_exp.scatter.pdf
*  Source data: parp_cnv.tcga_pancancer.cnv_summary.tsv
*  Code: parp7_ccle_tcga_analysis.R
