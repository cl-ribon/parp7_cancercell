#! /usr/bin/Rscript

# FigS1A
plotTCGAPARP7CNVFreq <- function() {
  #

  require("ggplot2")

  # Load CNV data and plot bar plots for freq across cancers
  cnv.fn <- "data/parp_cnv.tcga_pancancer.cnv_summary.tsv"
  cnv.df <- read.table(cnv.fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  cnv.df <- cnv.df[which(cnv.df$gene.symbol == "TIPARP"),]

  cna.df <- cnv.df[,c("cancer.code", "cancer.name", "all.amp.high.freq")]
  cna.df <- cna.df[-which(cna.df$all.amp.high == 0),]
  cna.df <- cna.df[sort.list(cna.df$all.amp.high, dec=TRUE),]
  cna.df$cancer.code <- toupper(cna.df$cancer.code)
  cna.df$cancer.code <- factor(cna.df$cancer.code, levels = cna.df$cancer.code)
  gp <- ggplot(cna.df, aes(x = cancer.code, y = all.amp.high.freq*100)) +
        geom_bar(stat = "identity") +
        ylim(0, 40) +
        xlab("") +
        ylab("Amplification Frequency (%)") +
        # ggtitle("PARP7 Copy-number Amplification Frequency in TCGA") +
        theme_minimal() +
        theme(text = element_text(size = 18), axis.text.x = element_text(angle =45, hjust = 1))

  fn <- "figs/FigS1A.parp7_tcga_cna_freq.bar.pdf"
  pdf(fn, width = 8, height = 4)
  plot(gp)
  dev.off()
}

# FigS1B
plotTCGAPARP7CNVExp <- function() {

  source("/home/rabo/github_repos/scripts/utils.R")
  source("/home/rabo/github_repos/scripts/cnv_utils.R")
  require("ggplot2")

  getAmpCancers <- function(threshold = 0.1) {
    #

    cnv.fn <- "data/parp_cnv.tcga_pancancer.cnv_summary.tsv"
    cnv.df <- read.table(cnv.fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cnv.df <- cnv.df[which(cnv.df$gene.symbol == "TIPARP"),]
    cna.df <- cnv.df[,c("cancer.code", "cancer.name", "all.amp.high.freq")]
    cna.df <- cna.df[-which(cna.df$all.amp.high == 0),]
    cna.df <- cna.df[sort.list(cna.df$all.amp.high, dec=TRUE),]
    amp.cancer.codes <- cna.df$cancer.code[which(cna.df$all.amp.high.freq >= threshold)]
    return(amp.cancer.codes)
  }

  getCNVExpValues <- function(cnv.df, cancer.code, gene.values) {
    # 

    merged.df <- NA
    exp.data <- getNormalizedTCGAexp(cancer.code)
    exp.df <- exp.data$exp.data$nofiltered.normalized.logcpm
    if ( gene.values$ensembl.id %in% rownames(exp.df) ) {
      exp.sample.ids <- colnames(exp.df)
      exp.sample.subids <- unlist(lapply(exp.sample.ids, function(x) substr(x, start=1, stop=16)))
      gene.idx <- which(rownames(exp.df) == gene.values$ensembl.id)
      exp.df <- data.frame(sample.ids = exp.sample.ids,
                           exp.values = as.numeric(exp.df[gene.idx, ]),
                           sample.subids = exp.sample.subids)

      cnv.df$sample.ids <- gsub(".", "-", cnv.df$sample.ids, fixed = TRUE)
      cnv.df$sample.subids <- unlist(lapply(cnv.df$sample.ids, function(x) substr(x, start=1, stop=16)))

      merged.df <- merge(cnv.df, exp.df, by = "sample.subids", all = TRUE)
      # merged.df <- merge(merged.df, exp.data$exp.manifest.data$df[,c("cases.0.samples.0.submitter.id", "cases.0.samples.0.sample.type")], by.x = "sample.subids", by.y = "cases.0.samples.0.submitter.id")
    } else {
      print("Gene", gene.values$symbol, gene.values$ensembl.id, "missing from expression matrix!")
    }

    return(merged.df)
  }

  getOverallCNV <- function(gistic.data, gene.values) {
    # Extract the overall CNV values and calls for a given gene and cancer 
    # type. Note that this differs from the focal CNV values because these 
    # include segments that would be considered "broad" or arm-level. So 
    # these calls are a combination of broad and focal CNV values.
    # 
    # Args:
    #   gistic.data:
    #   gene.values:
    #
    # Returns:
    # 

    sample.cutoffs.df <- gistic.data$gistic.sample.cutoffs
    geneCNV.values.df <- gistic.data$gistic.geneCNV$values
    geneCNV.calls.df <- gistic.data$gistic.geneCNV$calls
    gene.idx <- which(geneCNV.values.df[,c("Gene Symbol")] == gene.values$symbol)

    all.summary.df <- data.frame(cancer.code = gistic.data$cancer,
                                 gene.symbol = gene.values$symbol,
                                 all.amp.high = 0, 
                                 all.amp.low = 0,
                                 all.amp.high.freq = 0,
                                 all.amp.low.freq = 0,
                                 all.amp.max = 0,
                                 all.del.high = 0, 
                                 all.del.low = 0,
                                 all.del.high.freq = 0,
                                 all.del.low.freq = 0,
                                 all.del.min = 0,
                                 all.cnv.exp.cor = 0,
                                 all.cnv.exp.cor.pvalue = 0,
                                 all.amphigh.exp.lfc = 0,
                                 all.amplow.exp.lfc = 0,
                                 all.amp.exp.lfc = 0,
                                 nsamples = 0)

    gene.cnv.df <- data.frame(cancer.code = gistic.data$cancer,
                              sample.ids = NA, 
                              all.cnv.values = NA, 
                              all.cnv.calls = NA)

    if ( length(gene.idx) > 0 ) {
      sample.ids <- colnames(geneCNV.values.df)[-c(1:3)]
      gene.cnv.values <- as.numeric(geneCNV.values.df[gene.idx, -c(1:3)])
      gene.cnv.calls <- as.numeric(geneCNV.calls.df[gene.idx, -c(1:3)])
      gene.cnv.df <- data.frame(cancer.code = gistic.data$cancer,
                                sample.ids = sample.ids, 
                                all.cnv.values = gene.cnv.values, 
                                all.cnv.calls = gene.cnv.calls)

      # Amplifications
      all.summary.df$all.amp.high <- sum(gene.cnv.df$all.cnv.calls == 2)
      all.summary.df$all.amp.low <- sum(gene.cnv.df$all.cnv.calls == 1)
      all.summary.df$nsamples <- nrow(gene.cnv.df)
      all.summary.df$all.amp.high.freq <- all.summary.df$all.amp.high / all.summary.df$nsamples
      all.summary.df$all.amp.low.freq <- all.summary.df$all.amp.low / all.summary.df$nsamples
      all.summary.df$all.amp.max <- max(gene.cnv.values)

      # Deletions
      all.summary.df$all.del.high <- sum(gene.cnv.df$all.cnv.calls == -2)
      all.summary.df$all.del.low <- sum(gene.cnv.df$all.cnv.calls == -1)
      all.summary.df$all.del.high.freq <- all.summary.df$all.del.high / all.summary.df$nsamples
      all.summary.df$all.del.low.freq <- all.summary.df$all.del.low / all.summary.df$nsamples
      all.summary.df$all.del.min <- min(gene.cnv.values)

      # Merge CNV and Expression data sets
      cnv.exp.df <- getCNVExpValues(gene.cnv.df, cancer.code, gene.values)
      cor.res <- suppressWarnings(cor.test(cnv.exp.df$all.cnv.values, cnv.exp.df$exp.values, method="spearman", use="complete.obs"))
      all.summary.df$all.cnv.exp.cor <- cor.res$estimate
      all.summary.df$all.cnv.exp.cor.pvalue <- cor.res$p.value

      if ( all.summary.df$all.amp.high >= 10 ) {
        highamp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls == 2)], na.rm = TRUE)
        noamp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls == 0)], na.rm = TRUE)
        all.summary.df$all.amphigh.exp.lfc <- highamp.meanexp - noamp.meanexp
      }
      if ( all.summary.df$all.amp.low >= 10 ) {
        lowamp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls == 1)], na.rm = TRUE)
        noamp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls == 0)], na.rm = TRUE)
        all.summary.df$all.amplow.exp.lfc <- lowamp.meanexp - noamp.meanexp
      }
      if ( (all.summary.df$all.amp.high + all.summary.df$all.amp.low) >= 10 ) {
        amp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls > 0)], na.rm = TRUE)
        noamp.meanexp <- mean(cnv.exp.df$exp.values[which(cnv.exp.df$all.cnv.calls == 0)], na.rm = TRUE)
        all.summary.df$all.amp.exp.lfc <- amp.meanexp - noamp.meanexp
      }
    }

    overall.list <- list(summary = all.summary.df, 
                         values = gene.cnv.df,
                         cnv.exp.df = cnv.exp.df)
    return(overall.list)
  }

  cancer.codes <- getAmpCancers()
  gene.values <- getGeneList("parp")
  gene.values <- gene.values[which(gene.values$symbol == "TIPARP"),]

  cnv.exp.df <- c()
  cor.df <- c()
  for ( cancer.code in cancer.codes ) {
    print(cancer.code)

    gistic.data <- getGisticData(cancer.code, "amp")
    overall.values <- getOverallCNV(gistic.data, gene.values)
    df <- overall.values$cnv.exp.df[,c("all.cnv.values", "exp.values")]
    df$cancer.code <- toupper(cancer.code)
    corv <- cor.test(df$all.cnv.values, df$exp.values, method = "spearman", use ="complete.obs")
    cor.str <- paste("Spearman rho = ", format(corv$estimate, digits = 2), "\n", "P = ", format(corv$p.value, digits = 2), sep = "")
    cdf <- data.frame(cancer.code = toupper(cancer.code), cor.label = cor.str)
    cor.df <- rbind(cor.df, cdf)
    cnv.exp.df <- rbind(cnv.exp.df, df)
  }

  gp <- ggplot(cnv.exp.df, aes(x = exp.values, y = all.cnv.values)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method="lm", se=FALSE, color = "black", alpha = 0.5) +
        geom_text(data = cor.df, aes(x = 7, y = 3, label = cor.label, group = cancer.code), hjust = 0, size = 4) +
        xlab("Expression Level (log2TPM)") +
        xlim(0, 14) +
        ylab("CNV Level") +
        facet_wrap(~ cancer.code) +
        theme_minimal() +
        theme(text = element_text(size= 18))

  fn <- "/home/rabo/github_repos/p7/analysis/pub_analysis/p7_tcga_cna_exp.pt.pdf"
  pdf(fn, width = 10, height = 6)
  plot(gp)
  dev.off()
}

# FigS1C
plotTCGAPARP7ExpISGExpCorr <- function() {
  #

  source("/home/cbio/tap/src/tcga_utils.R")

  param.list <- getParams()
  # Load processed data
  fn <- file.path(param.list$data.dir, "interim", "tcga_L1.RData")
  load(fn)
  # Write merged molecular data for PARP7 to file
  exportMergedSummaryData(param.list, NA, "TIPARP")
  # Read in merged data file
  merged.fn <- file.path(param.list$data.dir, "processed", "tcga", "tcga_merged_moldata.TIPARP.tsv")
  parp7.ldf <- unique(read.table(merged.fn, header = TRUE, sep="\t", stringsAsFactors = FALSE))

  isg.list <- getISGSignature(param.list)
  metadata.df <- param.list$data.list$metadata

  isg.exp.df <- isg.list$score.df
  isg.exp.df$TCGA.ID <- gsub(".", "-", isg.exp.df$TCGA.ID, fixed = TRUE)
  parp7.ldf <- merge(parp7.ldf, isg.exp.df, by = "TCGA.ID", all.x = TRUE)

  p7_isg_cor.tumor.df <- c()
  for ( study in na.omit(unique(parp7.ldf$Study.Abbreviation)) ) {
    sdf.tumor <- parp7.ldf[which(parp7.ldf$Study.Abbreviation == study & parp7.ldf$Sample.type != "Solid Tissue Normal"),]
    cor.res <- cor.test(sdf.tumor$exp.tpm.value, sdf.tumor$score, method = "spearman", use = "complete.obs")
    nsamples <- nrow(na.omit(sdf.tumor[, c("exp.tpm.value", "score")]))
    sig.label <- ifelse(cor.res$p.value < 0.05, "Significant (P < 0.05)", "NS")
    cdf <- data.frame(Study.Abbreviation = study, Spearman.rho = cor.res$estimate, Spearman.Pval = cor.res$p.value, Nsamples = nsamples, Sig = sig.label)
    p7_isg_cor.tumor.df <- rbind(p7_isg_cor.tumor.df, cdf)
    print(paste(study, cor.res$estimate, cor.res$p.value, nsamples))
  }

  p7_isg_cor.tumor.df <- p7_isg_cor.tumor.df[sort.list(p7_isg_cor.tumor.df$Spearman.rho, dec=TRUE),]
  p7_isg_cor.tumor.df$Study.Abbreviation <- factor(p7_isg_cor.tumor.df$Study.Abbreviation, levels = p7_isg_cor.tumor.df$Study.Abbreviation)
  gp2 <- ggplot(p7_isg_cor.tumor.df, aes(x = Study.Abbreviation, y = Spearman.rho, fill = Sig)) +
        geom_bar(stat = "identity") +
        ylab("Correlation (Spearman Rho)") +
        xlab("") +
        ggtitle("PARP7 Expression - ISG Expression Score Correlations") +
        scale_fill_manual(values = c("Significant (P < 0.05)" = "blue", "NS" = "darkgrey")) +
        theme_minimal() +
        theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) 

  plot.fn <- "figs/FigS1C.parp7_tcga_exp_isg_cor.pdf"
  pdf(plot.fn, width = 9, height = 5)
  plot(gp2)
  dev.off()
}

# FigS1D/FigS1E
plotCCLEPARP7PARP1Dependency <- function() {
  #

  require("ggplot2")
  require("gridExtra")
  require("ggrepel")
  require("writexl")

  source("/home/cbio/tap/src/cell_line_utils.R")

  param.list <- getParams()
  # Load processed data
  fn <- file.path(param.list$data.dir, "interim", "cell_line_L1.RData")
  load(fn)
  # Write merged molecular data for PARP7 to file
  exportGeneGroupMerged(param.list, NA, "TIPARP")
  # Read in merged data file
  merged.fn <- file.path(param.list$data.dir, "processed", "cell_line", "cellline_merged_moldata.TIPARP.tsv")
  parp7.ldf <- read.table(merged.fn, header = TRUE, sep="\t", stringsAsFactors = FALSE)

  exportGeneGroupMerged(param.list, NA, "PARP1")
  # Read in merged data file
  merged.fn <- file.path(param.list$data.dir, "processed", "cell_line", "cellline_merged_moldata.PARP1.tsv")
  parp1.ldf <- read.table(merged.fn, header = TRUE, sep="\t", stringsAsFactors = FALSE)

  cols <- c("DepMap.ID", "CCLE.Name", "crispr.score.value", "crispr.prob")
  p7.dep.ldf <- parp7.ldf[, cols]
  colnames(p7.dep.ldf) <- c("DepMap.ID", "CCLE.Name", "PARP7.crispr.score.value", "PARP7.crispr.prob")
  p1.dep.ldf <- parp1.ldf[, cols]
  colnames(p1.dep.ldf) <- c("DepMap.ID", "CCLE.Name", "PARP1.crispr.score.value", "PARP1.crispr.prob")
  p7_p1.dep.ldf <- merge(p7.dep.ldf, p1.dep.ldf, all = TRUE, by = c("DepMap.ID", "CCLE.Name"))
  p7_p1.dep.ldf$Dep.label <- "No Dependency"
  p7_p1.dep.ldf$Dep.label[which(p7_p1.dep.ldf$PARP7.crispr.prob > 0.5 & p7_p1.dep.ldf$PARP1.crispr.prob > 0.5)] <- "PARP7 & PARP1 Dependent"
  p7_p1.dep.ldf$Dep.label[which(p7_p1.dep.ldf$PARP7.crispr.prob > 0.5 & p7_p1.dep.ldf$PARP1.crispr.prob < 0.5)] <- "PARP7 Dependent"
  p7_p1.dep.ldf$Dep.label[which(p7_p1.dep.ldf$PARP7.crispr.prob < 0.5 & p7_p1.dep.ldf$PARP1.crispr.prob > 0.5)] <- "PARP1 Dependent"

  cor.dep <- cor.test(p7_p1.dep.ldf$PARP7.crispr.score.value, p7_p1.dep.ldf$PARP1.crispr.score.value, method = "spearman")
  gp <- ggplot(p7_p1.dep.ldf, aes(x = PARP7.crispr.score.value, y = PARP1.crispr.score.value)) +
        geom_point(alpha = 0.6, size = 3, aes(color = Dep.label)) +
        geom_smooth(method = "lm", se = FALSE, color = "darkgrey") +
        annotate("text", x = -1, y = 0.25, label = paste("Spearman rho = ", round(cor.dep$estimate, 2), " (P = ", round(cor.dep$p.value, 2), ")", sep = ""), size = 5) +
        xlab("PARP7 Dependency Score (CRISPR)") +
        ylab("PARP1 Dependency Score (CRISPR)") +
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        scale_color_manual(values = c("PARP1 Dependent" = "blue", "PARP7 Dependent" = "red", "PARP7 & PARP1 Dependent" = "purple", "No Dependency" = "grey")) +
        theme_minimal() +
        theme(text = element_text(size = 16), legend.position = c(0.8, 0.8), legend.title = element_blank(), axis.line = element_line(size = 1))

  plot.fn <- "figs/FigS1D.parp7_parp1_cell_line_crispr_dep.scatter.pdf"
  pdf(plot.fn, width = 10, height = 7)
  plot(gp)
  dev.off()

  # PARP7 dependency vs PARP7 mRNA
  cor.val <- cor.test(parp7.ldf$crispr.score.value, parp7.ldf$exp.value, method = "spearman")
  cor.sig <- round(cor.val$p.value, 2)
  if ( cor.sig < 1e-15) {
    cor.sig <- "< 2e-16"
  } else {
    cor.sig <- paste(" = ", cor.sig, sep ="")
  }
  gp <- ggplot(parp7.ldf, aes(x = crispr.score.value, y = exp.value)) +
        geom_point(alpha = 0.6, size = 3) +
        geom_smooth(method = "lm", se = FALSE, color = "darkgrey") +
        annotate("text", x = -1.5, y = 2, label = paste("Spearman rho = ", round(cor.val$estimate, 2), " (P", cor.sig, ")", sep = ""), size = 5.5, hjust = 0) +
        xlab("PARP7 Dependency Score (CRISPR)") +
        ylab("PARP7 Expression (log2)") +
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        theme_minimal() +
        theme(text = element_text(size = 16), legend.position = c(0.8, 0.8), legend.title = element_blank(), axis.line = element_line(size = 1))

  plot.fn <- "figs/FigS1E.parp7_cell_line_crispr_dep_exp.scatter.pdf"
  pdf(plot.fn, width = 10, height = 7)
  plot(gp)
  dev.off()
}

plotTCGAPARP7CNVFreq()
plotTCGAPARP7CNVExp()
plotTCGAPARP7ExpISGExpCorr()
plotCCLEPARP7PARP1Dependency()