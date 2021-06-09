
deseqDGE <- function(param.list, response.metrics.df) {
  # Run DESeq2 on cell line basal expression between the PARP7i responder and nonresponder groups

  source("/home/rabo/github_repos/scripts/cell_line_utils.R")

  require("DESeq2")
  require("data.table")
  require("fgsea")
  require("GSVA")
  require("GSEABase")
  require("BiocParallel", lib.loc="/home/rabo/R/x86_64-redhat-linux-gnu-library/3.4")
  register(MulticoreParam())
  require("AnnotationDbi")
  require("org.Hs.eg.db")

  ccle.counts.exp <- loadCCLERnaseq(exp.format = "counts", gene.identifier = "ensids", save.file = TRUE, return.value = "dataframe")
  ccle.meta.df <- getCCLEMeta()
  exp.cell.names <- colnames(ccle.counts.exp)
  exp.cell.names <- sapply(exp.cell.names, function(x) unlist(strsplit(x, "_"))[1])
  cell.match.idxs <- match(response.metrics.df$norm.name, exp.cell.names)
  resp.cell.idxs <- which(!is.na(cell.match.idxs))
  cell.match.idxs <- na.omit(cell.match.idxs)
  ccle.counts.subset.df <- ccle.counts.exp[, cell.match.idxs]
  colnames(ccle.counts.subset.df) <- exp.cell.names[cell.match.idxs]
  colData.df <- response.metrics.df[resp.cell.idxs,]

  rownames(ccle.counts.subset.df) <- unlist(lapply(ccle.counts.exp[,1], function(x) unlist(strsplit(x, '.', fixed = TRUE))[1]))

  dds <- DESeqDataSetFromMatrix(countData = ccle.counts.subset.df,
                                colData = colData.df,
                                design = ~ Response.group)
  dds$Response.group <- factor(dds$Response.group, levels = c("Nonresponder", "Responder"))
  dds$Response.group <- relevel(dds$Response.group, ref = "Nonresponder")

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  deseq.res <- results(dds, name = "Response.group_Responder_vs_Nonresponder")
  deseq.res <- deseq.res[order(deseq.res$pvalue),]

  # res.lfc <- lfcShrink(dds, coef = "Response.group_Responder_vs_Nonresponder", type = "apeglm")

  deseq.res$ensid <- rownames(deseq.res)
  deseq.res$symbol <- mapIds(org.Hs.eg.db, 
                       keys = deseq.res$ensid, 
                       column = c("SYMBOL"), 
                       keytype = "ENSEMBL",
                       multiVals = "first")
  deseq.res$entrez <- mapIds(org.Hs.eg.db, 
               keys = deseq.res$ensid, 
               column = "ENTREZID", 
               keytype = "ENSEMBL",
               multiVals = "first")

  out.fn <- file.path(param.list$processed.data.dir, "parp7_cell_line_sensitivity_dge.deseq2.tsv")
  write.table(deseq.res, file = out.fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  gmt <- "/home/rabo/ref_data/genesets/c2_c5.all.v6.2.entrez.gmt"
  gs <- getGmt(gmt)

  res <- data.frame(deseq.res)
  z <- na.omit(data.frame(res[,c("entrez", "log2FoldChange")]))
  values <- z$log2FoldChange
  values <- values + runif(length(values), 0, 1e-6)
  names(values) <- z$entrez
  deseq.fgsea.res <- fgseaMultilevel(pathways=geneIds(gs), 
                     stats = sort(values, dec = FALSE),
                     minSize = 5,
                     maxSize = 500)
  deseq.fgsea.res <- deseq.fgsea.res[sort.list(deseq.fgsea.res$padj, dec = FALSE), ]

  deseq.list <- list(dge = deseq.res,
                     fgsea = deseq.fgsea.res)

  save.fn <- file.path(param.list$output.data.dir, "interim", "parp7_cell_line_sensitivity_basal_dge.deseq.RData")
  save(deseq.list, file = save.fn)
  return(deseq.list)
}

# Function to generate Fig5B.
publicationBasalExpResponseAnalysis <- function(param.list) {
  # Plot the differentially expressed genes and highlight genes in key enriched pathways

  # Load PARP7 inhibitor response metrics and clustered response groups based on publication compounds (RBN011364, RBN011595, RBN011628)
  # Load this file locally or from Quilt package
  response.metrics.df <- read.table("data/parp7.horizon_cell_line_response_metrics_wdf_response_groups.publication.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Differential expressed genes in enriched pathways plotted
  deseq.list <- deseqDGE(param.list, response.metrics.df)
  deseq.res.df <- data.frame(deseq.list$dge)
  deseq.fgsea.df <- deseq.list$fgsea
  lfc.df <- deseq.res.df[which(abs(deseq.res.df$log2FoldChange) > 0.8 & deseq.res.df$padj < 0.01 & !is.na(deseq.res.df$entrez)),]
  lfc.df <- lfc.df[sort.list(lfc.df$log2FoldChange, dec=FALSE),]
  lfc.df$ensid <- factor(lfc.df$ensid, levels = lfc.df$ensid)
  terms <- list("TYPE_I_INTERFERON" = "Type I Interferon Signaling Pathway", 
                 "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION" = "Antigen Processing and Presentation", 
                 "GO_EXTRACELLULAR_MATRIX" = "Extracellular Matrix")
  gs.terms <- rep(NA, nrow(lfc.df))
  sig.fgsea.df <- deseq.fgsea.df[which(deseq.fgsea.df$padj < 0.05 & deseq.fgsea.df$NES > 0),]
  for ( term in names(terms) ) {
    term.genes <- unlist(sig.fgsea.df$leadingEdge[grep(term, sig.fgsea.df$pathway)])
    gene.idxs <- match(term.genes, lfc.df$entrez)
    gs.terms[gene.idxs] <- terms[[term]]
  }
  lfc.df$gs.terms <- gs.terms
  gp.lfc <- ggplot(lfc.df, aes(x = ensid, y = log2FoldChange, fill = gs.terms)) +
            geom_bar(stat="identity") +
            scale_fill_manual(values=c("blue", "red", "purple"), na.value="lightgrey") +
            xlab("Differentially Expressed Genes") +
            ylab("Expression Fold-change (log2)") +
            ylim(-6, 6) +
            theme_minimal() + 
            coord_flip() +
            theme(axis.text.y=element_blank(), text=element_text(size=22), legend.title=element_blank(), legend.text=element_text(size=20), panel.grid.major=element_blank())

  plot.fn <- file.path(param.list$report.dir, "parp7_cell_line_sensitivity_deseq_dge_pathway.bar.pdf")
  pdf(plot.fn, width=10, height = 6)
  plot(gp.lfc)
  dev.off()

  # Differential expressed genes in enriched pathways plotted
  edger.list <- edgerDGE(param.list, response.metrics.df)
  edger.res.df <- data.frame(edger.list$dge)
  edger.fgsea.df <- edger.list$fgsea
  lfc.df <- edger.res.df[which(abs(edger.res.df$logFC) > 0.5 & edger.res.df$FDR < 0.05 & !is.na(edger.res.df$entrez)),]
  lfc.df <- lfc.df[sort.list(lfc.df$logFC, dec=FALSE),]
  lfc.df$ensid <- factor(lfc.df$ensid, levels = lfc.df$ensid)
  terms <- list("TYPE_I_INTERFERON" = "Type I Interferon Signaling Pathway", 
                 "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION" = "Antigen Processing and Presentation", 
                 "GO_EXTRACELLULAR_MATRIX" = "Extracellular Matrix")
  gs.terms <- rep(NA, nrow(lfc.df))
  sig.fgsea.df <- edger.fgsea.df[which(edger.fgsea.df$padj < 0.05 & edger.fgsea.df$NES > 0),]
  for ( term in names(terms) ) {
     term.genes <- unlist(sig.fgsea.df$leadingEdge[grep(term, sig.fgsea.df$pathway)])
     gene.idxs <- match(term.genes, lfc.df$entrez)
     gs.terms[gene.idxs] <- terms[[term]]
  }
  lfc.df$gs.terms <- gs.terms
  gp.lfc <- ggplot(lfc.df, aes(x = ensid, y = logFC, fill = gs.terms)) +
            geom_bar(stat="identity") +
            scale_fill_manual(values=c("blue", "red", "purple"), na.value="lightgrey") +
            xlab("Differentially Expressed Genes") +
            ylab("Expression Fold-change (log2)") +
            ylim(-6, 6) +
            theme_minimal() + 
            coord_flip() +
            theme(axis.text.y=element_blank(), text=element_text(size=22), legend.title=element_blank(), legend.text=element_text(size=20), panel.grid.major=element_blank())

  plot.fn <- file.path(param.list$report.dir, "parp7_cell_line_sensitivity_edger_dge_pathway.bar.pdf")
  pdf(plot.fn, width = 10, height = 6)
  plot(gp.lfc)
  dev.off()
}