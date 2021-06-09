#! /usr/bin/Rscript

# Plot Fig5F and FigS7 volcanoes
plotVolcanos <- function() {
  #

  # Load differential expression list with all results
  # Loads variable dge.list
  load("data/NCIH1373_NCIH647.6_24hr.dge_list.RData")

  df <- c()
  parp.df <- c()
  for ( dataset.name in names(dge.list) ) {
    df <- rbind(df, dge.list[[dataset.name]]$df)
    parp.df <- rbind(parp.df, dge.list[[dataset.name]]$parp.res)
  }

  df$label2 <- gsub("parp7_trt_ctrl.", "", df$label)
  df$label2 <- gsub(".", " ", df$label2, fixed=TRUE)
  df$label2 <- factor(df$label2, levels=c("NCI-H1373 24hr", "NCI-H647 24hr", "NCI-H1373 6hr", "NCI-H647 6hr"))

  parp.df$label2 <- gsub("parp7_trt_ctrl.", "", parp.df$label)
  parp.df$label2 <- gsub(".", " ", parp.df$label2, fixed=TRUE)
  parp.df$label2 <- factor(parp.df$label2, levels=c("NCI-H1373 24hr", "NCI-H647 24hr", "NCI-H1373 6hr", "NCI-H647 6hr"))

  volcano.df <- df[-which(is.na(df$padj) | is.na(df$log2FoldChange)),]
  volcano.df$pt.size <- ifelse(volcano.df$threshold == "TRUE", volcano.df$log2FoldChange + -log10(volcano.df$padj), 0.5)

  # Publication plot
  h1373.df <- volcano.df[grep("H1373", volcano.df$label2),]
  h1373.df$label2 <- factor(h1373.df$label2, levels = c("NCI-H1373 6hr", "NCI-H1373 24hr"))
  h1373.df$threshold.label <- h1373.df$threshold
  h1373.df$threshold.label[which(h1373.df$threshold)] <- "Padj<0.05, |log2FC|>1"
  h1373.df$threshold.label[-which(h1373.df$threshold)] <- "NS"
  gp.h1373 <- ggplot(h1373.df, aes(x=log2FoldChange, y=-log10(padj), colour=threshold.label, size = pt.size)) +
              geom_point(alpha=0.5) +
              scale_size(guide = "none") +
              xlab("Expression Fold-change (log2)") + 
              ylab("Significance") +
              facet_wrap(~ label2, scale="free") +
              theme_minimal() +
              scale_color_manual(values = c("Padj<0.05, |log2FC|>1" = "#FF5733", "NS" = "#bcbfc4")) +
              guides(colour = guide_legend(override.aes = list(size=10))) +
              theme(legend.title = element_blank(), text=element_text(size=20), strip.text = element_text(size = 24), axis.line = element_line(size = 0.25))

  h647.df <- volcano.df[grep("H647", volcano.df$label2),]
  h647.df$label2 <- factor(h647.df$label2, levels = c("NCI-H647 6hr", "NCI-H647 24hr"))
  h647.df$threshold.label <- h647.df$threshold
  h647.df$threshold.label[which(h647.df$threshold)] <- "Padj<0.05, |log2FC|>1"
  h647.df$threshold.label[-which(h647.df$threshold)] <- "NS"
  gp.h647 <- ggplot(h647.df, aes(x=log2FoldChange, y=-log10(padj), colour=threshold.label, size = pt.size)) +
              geom_point(alpha=0.5) +
              scale_size(guide = "none") +
              xlab("Expression Fold-change (log2)") + 
              ylab("Significance") +
              facet_wrap(~ label2, scale="free") +
              theme_minimal() +
              scale_color_manual(values = c("Padj<0.05, |log2FC|>1" = "#FF5733", "NS" = "#bcbfc4")) +
              guides(colour = guide_legend(override.aes = list(size=10))) +
              theme(legend.title = element_blank(), text=element_text(size=20), strip.text = element_text(size = 24), axis.line = element_line(size = 0.25))

  plot.fn <- "figs/Fig5F_FigS7.parp7_trt_ctrl.merged_volcano.pub.pdf")
  pdf(plot.fn, width=12, height=10)
  grid.arrange(gp.h1373, gp.h647, nrow = 2)
  dev.off()
}

# Plot Fig 5F and FigS7 geneset enrichment barplots
plotGSEs <- function(param.list, gse.list) {
  # Plot the significant genesets from up-regulated genes from treatment
  # as a barplot with size of bar reflecting significance value.
  # 
  # Highlight specific genesets related to biology themes that were 
  # consistent between cell lines or within a cell line.

  require("RColorBrewer")

  jaccard <- function(a, b) {
    c = intersect(a, b)
    jc <- length(c) / (length(a) + length(b) - length(c))
    return(jc)
  }

  deprecatedCode <- function() {
    # 
    # Deprecated code
    # gps <- list()
    # select.clusters.list <- list(NCIH1373.24hr = list(clusters = list("1" = "type I interferon signaling pathway", "2" = "DNA replication", "3" = "Influenza A", "4" = "One carbon pool by folate", "5" = "Negative regulation of gene expression"), cols = c()),
    #                              NCIH647.24hr = list(clusters = list("1" = "type I interferon signaling pathway", "2" = "regulation of viral genome replication", "3" = "extracellular matrix organization", "4" = "cytokine activity", "5" = "Cytosolic DNA-sensing pathway"), cols = c("type I interferon signaling pathway" = "#e2361f", "regulation of viral genome replication" = "#e2761e", "extracellular matrix organization" = "#aae21e", "cytokine activity" = "#1e93e2")))

    # bar.df <- c()
    # for ( dataset in c("NCIH1373.24hr", "NCIH647.24hr") ) {
    #   clusters <- list()
    #   df <- gse.list[[dataset]]$up$sig.enrichr
    #   for ( i in 1:nrow(df) ) {
    #     if ( df$Adjusted.P.value[i] < 0.05 ) {
    #       genes <- unlist(strsplit(df$Genes[i], ";"))
    #       print(paste(i, df$Term[i]))
    #       print(genes)
    #       hit <- FALSE
    #       if ( length(clusters) > 0 ) {
    #         for ( cluster in 1:length(clusters) ) {
    #           cluster.jaccard <- jaccard(genes, clusters[[cluster]]$genes)
    #           # print(paste("Checking cluster", cluster, cluster.jaccard, length(intersect(genes, clusters[[cluster]]$genes))))
    #           if ( cluster.jaccard > 0.35 & !hit ) {
    #             # Merge
    #             print("Merging")
    #             clusters[[cluster]]$genes <- unique(c(clusters[[cluster]]$genes, genes))
    #             print(clusters[[cluster]]$genes)
    #             clusters[[cluster]]$terms <- c(clusters[[cluster]]$terms, df$Term[i])
    #             bar.df <- rbind(bar.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = cluster, dataset = dataset))
    #             hit <- TRUE
    #           }
    #         }
    #       }
    #       if ( ! hit ) {
    #         clusters[[length(clusters)+1]] <- list(genes = genes, terms = df$Term[i])
    #         bar.df <- rbind(bar.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = length(clusters), dataset = dataset))
    #       }
    #     }
    #   }
      # cols <- c("type I interferon" = "red", "other" = "lightgrey")
      # bar.df$cluster.name <- bar.df$cluster
      # bar.df$cluster.name[which(bar.df$cluster == "1")] <- "type I interferon"
      # bar.df$cluster.name[which(bar.df$cluster != "1")] <- "other"
      # gp.bar <- ggplot(bar.df, aes(x = reorder(term, -log10(pval), FUN=max), y = -log10(pval), fill = cluster.name)) +
      #           geom_bar(stat="identity") +
      #           coord_flip() +
      #           xlab("") +
      #           ggtitle(dataset) +
      #           scale_fill_manual(values = cols) +
      #           theme_minimal() +
      #           theme(axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), legend.position = "none")
      # gps[[dataset]] <- gp.bar
    # }

    # gs.freq.df <- data.frame(table(bar.df$term))
    # gs.freq.df <- gs.freq.df[sort.list(gs.freq.df$Freq, dec=T),]
    # sig.terms <- as.character(gs.freq.df$Var1[which(gs.freq.df$Freq == 2)])

    # clusters <- list()
    # for ( dataset in c("NCIH1373.24hr", "NCIH647.24hr") ) {
    #   df <- gse.list[[dataset]]$up$sig.enrichr
    #   df <- df[which(df$Term %in% sig.terms),]
    #   for ( i in 1:nrow(df) ) {
    #     if ( df$Adjusted.P.value[i] < 0.05 ) {
    #       genes <- unlist(strsplit(df$Genes[i], ";"))
    #       # print(paste(i, df$Term[i]))
    #       # print(genes)
    #       hit <- FALSE
    #       if ( length(clusters) > 0 ) {
    #         for ( cluster in 1:length(clusters) ) {
    #           cluster.jaccard <- jaccard(genes, clusters[[cluster]]$genes)
    #           # print(paste("Checking cluster", cluster, cluster.jaccard, length(intersect(genes, clusters[[cluster]]$genes))))
    #           if ( cluster.jaccard > 0.25 & !hit ) {
    #             # Merge
    #             print("Merging")
    #             clusters[[cluster]]$genes <- unique(c(clusters[[cluster]]$genes, genes))
    #             print(clusters[[cluster]]$genes)
    #             clusters[[cluster]]$terms <- c(clusters[[cluster]]$terms, paste(dataset, df$Term[i], sep="."))
    #             # bar.df <- rbind(bar.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = cluster, dataset = dataset))
    #             hit <- TRUE
    #           }
    #         }
    #       }
    #       if ( ! hit ) {
    #         clusters[[length(clusters)+1]] <- list(genes = genes, terms = paste(dataset, df$Term[i], sep="."))
    #         # bar.df <- rbind(bar.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = length(clusters), dataset = dataset))
    #       }
    #     }
    #   }
    # }

    # cols <- c("type I interferon" = "red", "other" = "lightgrey")
    # bar.df$cluster.name <- bar.df$cluster
    # bar.df$cluster.name[which(bar.df$cluster == "1")] <- "type I interferon"
    # bar.df$cluster.name[which(bar.df$cluster != "1")] <- "other"
    # gp.bar <- ggplot(bar.df, aes(x = reorder(term, -log10(pval), FUN=max), y = -log10(pval), fill = cluster.name)) +
    #           geom_bar(stat="identity") +
    #           coord_flip() +
    #           xlab("") +
    #           scale_fill_manual(values = cols) +
    #           facet_wrap(~ dataset, scales = "free", ncol = 2) +
    #           theme_minimal() +
    #           theme(axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank())
    # merged.dir <- checkDir(file.path(param.list$analysis.dir, param.list$analysis.label, "merged"))
    # plot.fn <- file.path(merged.dir, "parp7_trt_ctrl.NCIH1373-NCIH647.24hr.gse.barplot.pdf")
    # pdf(plot.fn, width=8, height=6)
    # grid.arrange(gps[[1]], gps[[2]], ncol=2)
    # dev.off()
  }

  clusterGenesetsByGenes <- function(gse.list) {
    #

    clusters <- list()
    cluster.df <- c()
    for ( dataset in c("NCIH1373.24hr", "NCIH647.24hr") ) {
      clusters[[dataset]] <- list()
      df <- gse.list[[dataset]]$up$sig.enrichr
      df <- df[which(df$Adjusted.P.value[i] < 0.05),]
      for ( i in 1:nrow(df) ) {
        genes <- unlist(strsplit(df$Genes[i], ";"))
        hit <- FALSE
        if ( length(clusters[[dataset]]) > 0 ) {
          for ( cluster in 1:length(clusters[[dataset]]) ) {
            cluster.jaccard <- jaccard(genes, clusters[[dataset]][[cluster]]$genes)
            if ( cluster.jaccard > 0.35 & !hit ) {
              clusters[[dataset]][[cluster]]$genes <- unique(c(clusters[[dataset]][[cluster]]$genes, genes))
              clusters[[dataset]][[cluster]]$terms <- c(clusters[[dataset]][[cluster]]$terms, df$Term[i])
              cluster.df <- rbind(cluster.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = cluster, dataset = dataset))
              hit <- TRUE
            }
          }
          if ( ! hit ) {
            clusters[[dataset]][[length(clusters)+1]] <- list(genes = genes, terms = df$Term[i])
            cluster.df <- rbind(cluster.df, data.frame(term = df$Term[i], pval = df$Adjusted.P.value[i], ngenes = length(genes), cluster = length(clusters[[dataset]]), dataset = dataset))
          }
        }
      }
    }

          # Determine the significant genesets from up-regulated genes 
    gs.freq.df <- data.frame(table(bar.df$term))
    gs.freq.df <- gs.freq.df[sort.list(gs.freq.df$Freq, dec=T),]

    res.list <- list(clusters = clusters, overlap.df = gs.freq.df)
    return(res.list)
  }

  # Load geneset enrichment list object
  # Loads variable gse.list
  load("data/NCIH1373_NCIH647.6_24hr.gse_list.RData")
  gse.clusters <- clusterGenesetsByGenes(gse.list)

  select.clusters <- list("dsRNA response/binding" = c("dsRNA response/binding", "nuclease activity", "double-stranded RNA"), 
                          "Type I interferon" = c("type I interferon"),
                          "Chemotaxis" = c("chemotaxis"),
                          "Cell proliferation" = c("cell proliferation"),
                          "Cytokine/Chemokine" = c("cytokine", "chemokine"),
                          "Viral response" = c("viral", "Influenza", "Herpes", "Hepatitis", "Measles"),
                          "Extracellular matrix" = c("extracellular", "ECM", "adhesion"))

  darkcols <- brewer.pal(7, "Dark2")
  select.cols <- c("dsRNA response/binding" = darkcols[1], 
                   "Other" = "lightgrey",
                   "Type I interferon" = darkcols[2],
                   "Chemotaxis" = darkcols[3],
                   "Cell proliferation" = darkcols[4],
                   "Viral response" = darkcols[5],
                   "Cytokine/Chemokine" = darkcols[6],
                   "Extracellular matrix" = darkcols[7])

  gps <- list()
  for ( dataset in c("NCIH1373.24hr", "NCIH647.24hr") ) {
    df <- gse.list[[dataset]]$up$sig.enrichr
    df <- df[which(df$Adjusted.P.value < 0.05),]
    df$dataset <- dataset
    df$group <- "Other"

    for ( select.cluster in names(select.clusters)) {
      print(select.cluster)
      search.terms <- select.clusters[[select.cluster]]
      select.idxs <- c()
      for ( search.term in search.terms ) {
        select.idxs <- c(select.idxs, grep(search.term, df$Term))
      }
      select.idxs <- unique(select.idxs)
      if ( length(select.idxs) > 0 ) {
        df$group[select.idxs] <- select.cluster
      }
    }

    df$group <- factor(df$group, levels = c("Type I interferon", "Viral response", "Cytokine/Chemokine", "dsRNA response/binding", "Chemotaxis", "Cell proliferation", "Extracellular matrix", "Other"))
    lp <- ifelse(dataset == "NCIH1373.24hr", "bottom", "bottom")
    gp.bar <- ggplot(df, aes(x = reorder(Term, log10(Adjusted.P.value), FUN=max), y = -log10(Adjusted.P.value), fill = group)) +
              geom_bar(stat = "identity") +
              xlab("") +
              ylab("-log10(Pvalue)") +
              ggtitle(dataset) +
              scale_fill_manual(values = select.cols) +
              theme_minimal() +
              theme(text = element_text(size=16), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), legend.position = lp, legend.title = element_blank(), legend.text = element_text(size = 10))
    gps[[dataset]] <- gp.bar
  }

  plot.fn <- "figs/Fig5F_FigS7.gse_barplot.parp7_trt_ctrl.NCIH1373-NCIH647.24hr.gse.barplot.pdf"
  pdf(plot.fn, width=16, height=4)
  grid.arrange(gps[[1]], gps[[2]], ncol=2)
  dev.off()
}

# Plot Fig5F and FigS7 overlap between cell panel screen and cell line deg
publicationComparePredictivePostDoseExp <- function(param.list) {
  # Requires use of Rv3.4.1 (/usr/bin/R) for proper font size in Venn diagram

  library("venn")

  # Load RNA-seq differential gene expression results
  load("data/NCIH1373_NCIH647.6_24hr.dge_list.RData")
  invitro_trx_dge.list <- dge.list

  # NCI-H1373 24hr DEG results
  ncih1373.24hr.df <- invitro_trx_dge.list$NCIH1373.24hr$df
  ncih1373.24hr.df$hits <- 0
  hit.filter <- which(ncih1373.24hr.df$padj < 0.05 & ncih1373.24hr.df$log2FoldChange > 1)
  ncih1373.24hr.df$hits[hit.filter] <- 1
  ncih1373.24hr.df <- ncih1373.24hr.df[,c("entrez", "symbol", "log2FoldChange", "padj", "hits")]

  # Cell line proliferation screen
  # Loads variable deseq.list
  load("data/parp7_cell_line_sensitivity_basal_dge.deseq.RData")
  # parp7i.dge.df <- read.table(parp7i.dge.fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Filter columns
  parp7i.dge.df <- data.frame(deseq.list$dge)
  parp7i.dge.df$hits <- 0
  hit.filter <- which(parp7i.dge.df$padj < 0.05 & parp7i.dge.df$log2FoldChange > 1)
  parp7i.dge.df$hits[hit.filter] <- 1
  parp7i.dge.df <- parp7i.dge.df[, c("entrez", "symbol", "log2FoldChange", "padj", "hits")]
  colnames(parp7i.dge.df) <- c("Entrez.GeneID", "Gene.symbol", "PARP7i.DGE.log2FC", "PARP7i.DGE.padj", "PARP7i.DGE.hits")

  x <- ncih1373.24hr.df[which(ncih1373.24hr.df$hits == 1 & !is.na(ncih1373.24hr.df$entrez)),]
  y <- parp7i.dge.df[which(parp7i.dge.df$PARP7i.DGE.hits == 1 & !is.na(parp7i.dge.df$Entrez.GeneID)),]
  m <- merge(x, y, by.x = "entrez", by.y = "Entrez.GeneID", all = TRUE)

  # xx <- ncih1373.24hr.df[which(!is.na(ncih1373.24hr.df$entrez)),]
  # yy <- parp7i.dge.df[which(!is.na(parp7i.dge.df$Entrez.GeneID)),]
  # no <- c()
  # for ( i in 1:10000 ) {
  #   xx$shits <- sample(xx$hits)
  #   yy$shits <- sample(yy$PARP7i.DGE.hits)
  #   xx.temp <- xx[which(xx$shits == 1),]
  #   yy.temp <- yy[which(yy$shits == 1),]
  #   mm <- merge(xx.temp, yy.temp, by.x = "entrez", by.y = "Entrez.GeneID")
  #   no <- c(no, nrow(mm))
  # }

  mv <- m[,c("hits", "PARP7i.DGE.hits")]
  colnames(mv) <- c("PD", "Pred")
  mv[is.na(mv)] <- 0

  pdf("figs/Fig5F.raw_venn.parp7_cell_line_pred_pd_exp.ncih1373.venn.pdf")
  venn(mv, ilab=TRUE, zcolor="style", cexil=2, cexsn=1)
  dev.off()

  # NCI-H647
  ncih647.24hr.df <- invitro_trx_dge.list$NCIH647.24hr$df
  ncih647.24hr.df$hits <- 0
  hit.filter <- which(ncih647.24hr.df$padj < 0.05 & ncih647.24hr.df$log2FoldChange > 1)
  ncih647.24hr.df$hits[hit.filter] <- 1
  ncih647.24hr.df <- ncih647.24hr.df[,c("entrez", "symbol", "log2FoldChange", "padj", "hits")]

  x2 <- ncih1373.24hr.df[which(ncih647.24hr.df$hits == 1 & !is.na(ncih647.24hr.df$entrez)),]
  m2 <- merge(x2, y, by.x = "entrez", by.y = "Entrez.GeneID", all = TRUE)

  mv2 <- m2[,c("hits", "PARP7i.DGE.hits")]
  colnames(mv2) <- c("PD", "Pred")
  mv2[is.na(mv2)] <- 0

  pdf("figs/FigS7.raw_venn.parp7_cell_line_pred_pd_exp.ncih647.venn.pdf")
  venn(mv2, ilab=TRUE, zcolor="style", cexil=2, cexsn=1)
  dev.off()
}

plotVolcanos()
plotGSEs()
publicationComparePredictivePostDoseExp()