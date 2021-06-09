#! /usr/bin/Rscript

# Fig7A heatmap of expression fold changes due to PARP7 inhibition
plotNormalizedValues <- function() {
  #

  getIFNGenes <- function() {
    #

    source("/home/rabo/github_repos/scripts/utils.R")

    require("ontologyIndex")
    require("data.table")

    entrez.gene2go.fn <- "/home/rabo/ref_data/ncbi_gene/gene2go.gz"
    entrez.gene2go.df <- fread(entrez.gene2go.fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Filter down to human taxon 9606
    entrez.gene2go.df <- entrez.gene2go.df[which(entrez.gene2go.df[,1] == "9606"),]

    obo.fn <- "/home/rabo/ref_data/gene_ontology/go-basic.obo"
    obo <- get_ontology(obo.fn, extract_tags = "everything", propagate_relationships = c("is_a", "part_of", "regulates"))

    iter <- 0
    # Innate immune response
    new.gos <- c("GO:0045087")
    all.go <- new.gos
    while ( length(new.gos) > 0 ) {
      buffer <- c()
      for ( go.id in new.gos ) {
        idx <- which(obo[["id"]] == go.id)
        go.term <- obo[["name"]][idx]
        desc.go <- get_descendants(obo, go.id, exclude_roots = TRUE)
        buffer <- unique(c(buffer, desc.go))
        print(go.term)
      }
      if ( length(buffer) > 0 ) {
        if ( length(which(buffer %in% all.go)) > 0 ) {
          buffer <- buffer[-which(buffer %in% all.go)]
        }
        new.gos <- buffer
        all.go <- unique(c(all.go, buffer))
      } else {
        new.gos <- c()
      }
      iter <- iter + 1
    }
    go.df <- entrez.gene2go.df[which(entrez.gene2go.df$GO_ID %in% all.go),]
    go.df$GeneID <- as.character(go.df$GeneID)
    go.genevalues.df <- updateGeneValuesBatch(go.df$GeneID, "entrezgene")
    go.dff <- merge(go.df, unique(go.genevalues.df[,c("query", "entrezgene", "symbol")]), by.x = "GeneID", by.y = "query")
    go.df <- go.dff[, c("entrezgene", "symbol", "GO_term")]
    colnames(go.df) <- c("Entrez.GeneID", "Gene.symbol", "Source")
    go.df$Source <- paste("GO", go.df$Source, sep=".")
    ifn.df <- go.df[grep("interf", go.df$Source),]
    genes <- c(ifn.df$Gene.symbol, "CXCL10")
    return(genes)
  }

  require("ComplexHeatmap")
  require("RColorBrewer")
  require("circlize")

  # Load normalized nanostring expression data
  # Loads variable normalized.ldf
  load("data/parp7_nanostring_100gp.RData")

  ifn.genes <- getIFNGenes()
  nano.ifn.ldf <- normalized.ldf[which(normalized.ldf$Probe.name %in% ifn.genes),]
  nano.ifn.ldf <- nano.ifn.ldf[-which(nano.ifn.ldf$Probe.name == "EGR1"),]
  nano.ifn.ldf <- nano.ifn.ldf[-which(nano.ifn.ldf$Cell.line == "NCIH596"),]

  baseline.ldf <-  nano.ifn.ldf[which(nano.ifn.ldf$Concentration.nm == 0),]
  baseline.ldf$GI50.bin <- ifelse(baseline.ldf$GI50.nM < 100, "responder", "nonresponder")
  exp.df <- aggregate(value ~ GI50.bin, data = baseline.ldf, mean)
  exp.res.fc <- exp.df[2,2] / exp.df[1,2]
  res.pval <- wilcox.test(value ~ GI50.bin, baseline.ldf)$p.value
  print(paste("Baseline ISG expression GI50 Responder Groups", format(exp.res.fc, digits=3), format(res.pval, digits = 2)))

  # Heatmap
  fc.df <- c()
  for ( cell.line in unique(nano.ifn.ldf$Cell.line) ) {
    cell.df <- nano.ifn.ldf[which(nano.ifn.ldf$Cell.line == cell.line),]
    for ( gene in unique(cell.df$Probe.name) ) {
      gene.df <- cell.df[which(cell.df$Probe.name == gene),]
      ctrl.value <- gene.df$value[which(gene.df$Concentration.nm == 0)]
      gene.df$fc <- gene.df$value / ctrl.value
      fc.df <- rbind(fc.df, gene.df)
    }
  }
  fc.df <- fc.df[-which(fc.df$Concentration.nm == 0),]

  sample.df <- nano.ifn.ldf[which(nano.ifn.ldf$Concentration.nm == 0 & nano.ifn.ldf$Probe.name == "STAT1"), c("Cell.line", "Cancer.type", "GI50.nM")]
  sample.df$GI50.bin <- ifelse(sample.df$GI50.nM < 100, "GI50 < 100 nM", "GI50 >= 100 nM")
  gp.cell.wdf <- dcast(data = fc.df, formula = Probe.name ~ Cell.line, value.var = "fc")
  rownames(gp.cell.wdf) <- gp.cell.wdf$Probe.name
  gp.cell.wdf <- gp.cell.wdf[,-1]

  # PARP7 publication
  sample.df <- sample.df[match(colnames(gp.cell.wdf), sample.df$Cell.line), ]
  nanostring.lfc <- as.matrix(log2(gp.cell.wdf))
  nanostring.lfc <- nanostring.lfc[,-which(colnames(nanostring.lfc) == "CFPAC1")]
  sample.df <- sample.df[-which(sample.df$Cell.line == "CFPAC1"),]

  ha.nano = HeatmapAnnotation(df = data.frame(GI50 = sample.df$GI50.bin), col = list(GI50 = c("GI50 < 100 nM" = "darkolivegreen", "GI50 >= 100 nM" = "lightgrey")))
  hm.cols <- colorRamp2(c(-2, 0, 7), c("blue", "white", "red"))
  hm.nano <- Heatmap(nanostring.lfc, 
                     name="Nanostring log2 FC", 
                     km=1, 
                     col=hm.cols,
                     cluster_rows = TRUE,
                     cluster_columns = TRUE,
                     show_row_dend = TRUE,
                     show_column_names=TRUE,
                     row_names_side = "right",
                     column_title="Nanostring Treatment Fold-Changes",
                     column_title_gp = gpar(fontsize = 18),
                     row_names_gp = gpar(fontsize = 10),
                     show_row_names=TRUE,
                     top_annotation = ha.nano,
                     heatmap_legend_param=list(legend_direction="vertical", color_bar="continuous"))

  plot.fn <- "figs/Fig7A.parp7_nanostring_100gp.foldchanges.RBN2397.20191023.pub.chm.pdf"
  pdf(plot.fn, width=9, height=6)
  draw(hm.nano)
  dev.off()
}

plotNormalizedValues()