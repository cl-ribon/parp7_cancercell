#! /usr/bin/Rscript

updateGeneValuesBatch.mu <- function(genes, type = "symbol") {
  # 

  require("jsonlite")
  require("httr")

  all.df <- c()
  start.idxs <- seq(1, length(genes), 1000)
  for ( start.idx in start.idxs ) {
    end.idx <- min(length(genes), start.idx + 1000 - 1)
    print(paste("Querying genes: ", start.idx, "-", end.idx, sep=""))
    gene.query <- paste(genes[start.idx:end.idx], collapse = ",")
    post.params <- list(q = gene.query, scopes = "symbol", fields = "name,symbol,taxid,entrezgene", species = "mouse", size = 1)
    res <- POST("http://mygene.info/v3/query", body = post.params, content.type = ("application/json"))
    res.df <- fromJSON(rawToChar(res$content))

    if ( !"notfound" %in% colnames(res.df) ) {
      res.df$notfound <- NA
    }
    all.df <- rbind(all.df, res.df)
  }

  missing.df <- all.df[which(all.df$notfound),]
  missing.genes <- missing.df$query
  alias.df <- c()
  start.idxs <- seq(1, length(missing.genes), 1000)
  for ( start.idx in start.idxs ) {
    end.idx <- min(length(missing.genes), start.idx + 1000 - 1)
    print(paste("Querying genes: ", start.idx, "-", end.idx, sep=""))
    gene.query <- paste(missing.genes[start.idx:end.idx], collapse = ",")
    post.params <- list(q = gene.query, scopes = "alias", fields = "name,symbol,taxid,entrezgene", species = "mouse")
    res <- POST("http://mygene.info/v3/query", body = post.params, content.type = ("application/json"))
    res.df <- fromJSON(rawToChar(res$content))

    if ( !"notfound" %in% colnames(res.df) ) {
      res.df$notfound <- NA
    }
    alias.df <- rbind(alias.df, res.df)
  }

  if ( length(alias.df) > 0 ) {
    alias.df <- alias.df[, colnames(all.df)]
  }

  dup.idxs <- which(duplicated(alias.df$query))
  dup.genes <- unique(alias.df$query[dup.idxs])
  for(dup.query.name in dup.genes) {
    dup.df <- alias.df[which(alias.df$query == dup.query.name),]
    keep.idxs <- which(! dup.df$entrezgene %in% all.df$entrezgene)
    if ( length(keep.idxs) > 0 ) {
      # print(dup.df[keep.idxs[1],])
      alias.df <- alias.df[-which(alias.df$query == dup.query.name),]
      alias.df <- rbind(alias.df, dup.df[keep.idxs[1],])
    } else {
      alias.df <- alias.df[-which(alias.df$query == dup.query.name),]
      print("NA")
    }
  }

  midx <- match(alias.df$query, all.df$query)
  all.df[midx,] <- alias.df

  return(all.df)
}

updateGeneValuesBatch <- function(genes, type = "symbol") {
  # 

  require("jsonlite")
  require("httr")

  all.df <- c()
  start.idxs <- seq(1, length(genes), 1000)
  for ( start.idx in start.idxs ) {
    end.idx <- min(length(genes), start.idx + 1000 - 1)
    print(paste("Querying genes: ", start.idx, "-", end.idx, sep=""))
    gene.query <- paste(genes[start.idx:end.idx], collapse = ",")
    post.params <- list(q = gene.query, scopes = "symbol", fields = "name,symbol,taxid,entrezgene", species = "human", size = 1)
    res <- POST("http://mygene.info/v3/query", body = post.params, content.type = ("application/json"))
    res.df <- fromJSON(rawToChar(res$content))

    if ( !"notfound" %in% colnames(res.df) ) {
      res.df$notfound <- NA
    }
    all.df <- rbind(all.df, res.df)
  }

  missing.df <- all.df[which(all.df$notfound),]
  missing.genes <- missing.df$query
  alias.df <- c()
  start.idxs <- seq(1, length(missing.genes), 1000)
  for ( start.idx in start.idxs ) {
    end.idx <- min(length(missing.genes), start.idx + 1000 - 1)
    print(paste("Querying genes: ", start.idx, "-", end.idx, sep=""))
    gene.query <- paste(missing.genes[start.idx:end.idx], collapse = ",")
    post.params <- list(q = gene.query, scopes = "alias", fields = "name,symbol,taxid,entrezgene", species = "human")
    res <- POST("http://mygene.info/v3/query", body = post.params, content.type = ("application/json"))
    res.df <- fromJSON(rawToChar(res$content))

    if ( !"notfound" %in% colnames(res.df) ) {
      res.df$notfound <- NA
    }
    alias.df <- rbind(alias.df, res.df)
  }

  dup.idxs <- which(duplicated(alias.df$query))
  dup.genes <- unique(alias.df$query[dup.idxs])
  for(dup.query.name in dup.genes) {
    dup.df <- alias.df[which(alias.df$query == dup.query.name),]
    keep.idxs <- which(! dup.df$entrezgene %in% all.df$entrezgene)
    if ( length(keep.idxs) > 0 ) {
      # print(dup.df[keep.idxs[1],])
      alias.df <- alias.df[-which(alias.df$query == dup.query.name),]
      alias.df <- rbind(alias.df, dup.df[keep.idxs[1],])
    } else {
      alias.df <- alias.df[-which(alias.df$query == dup.query.name),]
      print("NA")
    }
  }

  midx <- match(alias.df$query, all.df$query)
  all.df[midx,] <- alias.df

  return(all.df)
}

getSubIDMerged <- function(metric) {
  # 
  
  require("RobustRankAggreg")

  # Load substrated ID object
  # Loads variable param.list
  load("data/invitro_substrate_id.RData")
  invitro_substrate_id.list <- param.list$data
  invitro_substrate_id.formatted.list <- list()
  pub_table.ldf <- c()
  merged.wdf <- c()

  p7.experiments <- c("skmes", "skmes.mg132.af1521", "ncih1373", "skmes.parp14m3", "skmes.dox")
  for ( exp in p7.experiments ) {
    df <- invitro_substrate_id.list[[exp]]$df

    # Deal with duplicate genes
    dup.idxs <- which(duplicated(df$genename))
    if ( length(dup.idxs) > 0 ) {
      dup.genes <- unique(df$genename[dup.idxs])
      for ( dup.gene in dup.genes ) {
        dup.df <- df[which(df$genename == dup.gene),]
        new.df <- collapseDups(dup.df)
        df <- df[-which(df$genename == dup.gene),]
        df <- rbind(df, new.df)
      }
    }

    # Calculate raw fold-change
    fc.constant <- 0.01
    trt.fc.col <- which(colnames(df) == "trt.fc")
    colnames(df)[trt.fc.col] <- "trt.fc.binding"
    spc.cols <- grep("_SpC", colnames(df))
    df$trt.fc <- (df[,spc.cols[2]] + fc.constant) / (df[,spc.cols[1]] + fc.constant)

    # Map gene
    if ( exp == "unscc" ) {
      df.genevals <- updateGeneValuesBatch.mu(df$genename)
      df <- merge(df, df.genevals[, c("query", "symbol", "entrezgene")], by.x = "genename", by.y = "query", all.x = TRUE)
      colnames(df)[which(colnames(df) == "symbol")] <- "gene.symbol.xref"
    } else {
      df.genevals <- updateGeneValuesBatch(toupper(df$genename))
      df <- merge(df, df.genevals[, c("query", "symbol")], by.x = "genename", by.y = "query", all.x = TRUE)
      colnames(df)[which(colnames(df) == "symbol")] <- "gene.symbol.xref"
    }
    # substrate.list[[exp]] <- df
    invitro_substrate_id.formatted.list[[exp]] <- df[,c("entrez.id", "genename", metric, "gene.symbol.xref")]

    cnames <- colnames(df)[grep("_SpC", colnames(df))]
    pt.ldf <- df[, c("entrez.id", "genename", "protein_desc", "Molecular_Weight")]
    pt.ldf$protein_desc <- gsub('"', "", pt.ldf$protein_desc)

    if ( exp == "skmes" ) {
      merged.wdf <- df[, c("entrez.id", cnames, metric)]
      colnames(merged.wdf) <- c("entrez.id", cnames, paste(exp, metric, sep="."))
    } else {
      metric.df <- df[, c("entrez.id", cnames, metric)]
      colnames(metric.df) <- c("entrez.id", cnames, paste(exp, metric, sep="."))
      merged.wdf <- merge(merged.wdf, metric.df, by = "entrez.id", all = TRUE)
    }
    pub_table.ldf <- rbind(pub_table.ldf, pt.ldf)

    fn <- file.path("/home/rabo/github_repos/p7/analysis/pub_analysis", paste(exp, "subid_formatted.tsv", sep = "."))
    write.table(df, file = fn, quote = FALSE, row.names = FALSE, sep = "\t")
  }

  wdf <- c()
  for ( exp in p7.experiments ) {
    df <- invitro_substrate_id.formatted.list[[exp]]
    colnames(df)[3:4] <- paste(colnames(df)[3:4], exp, sep=".")
    missing.symbols.idx <- which(is.na(df$gene.symbol.xref))
    df$gene.symbol.xref[missing.symbols.idx] <- df$genename[missing.symbols.idx]
    missing.entrez.idx <- which(is.na(df$entrez.id))
    df$entrez.id[missing.entrez.idx] <- df$gene.symbol.xref[missing.entrez.idx]
    df <- df[,-which(colnames(df) == "genename")]
    df <- df[,-which(colnames(df) == paste("gene.symbol.xref", exp, sep="."))]
    if ( length(wdf) == 0 ) {
      wdf <- df
    } else {
      if ( exp == "unscc") {
        df$gene.symbol.xref <- toupper(df$gene.symbol.xref)
        wdf <- merge(wdf, df, by = "gene.symbol.xref", all = TRUE)
        wdf$entrez.id <- unlist(apply(wdf[, grep("entrez", colnames(wdf)),], 1, function(x) unique(x[which(!is.na(x))])[1]))
        idxs <- grep("entrez", colnames(wdf))
        wdf <- wdf[,-c(idxs[1:2])]
      } else {
        df <- df[-which(colnames(df) == "gene.symbol.xref"),]
        wdf <- merge(wdf, df, by = "entrez.id", all = TRUE)
        wdf$gene.symbol.xref <- unlist(apply(wdf[, grep("xref", colnames(wdf)),], 1, function(x) unique(x[which(!is.na(x))])[1]))
        idxs <- grep("xref", colnames(wdf))
        wdf <- wdf[,-c(idxs[1:2])]
      }
    }
  }

  # Flip SKMES dox to match others
  skmes.dox.idx <- grep("skmes.dox", colnames(wdf))
  flip.idx <- which(!is.na(wdf[,skmes.dox.idx]))
  if ( metric == "spc.trt_ctrl.delta" | metric == "spc.ctrl_fraction_delta" ) {
    wdf[flip.idx, skmes.dox.idx] <- -1*wdf[flip.idx, skmes.dox.idx]
  } else {
    wdf[flip.idx, skmes.dox.idx] <- 1/wdf[flip.idx, skmes.dox.idx]
  }

  xdf <- data.frame(entrez.id = wdf$entrez.id, gene.symbol.xref = as.character(wdf$gene.symbol.xref), apply(wdf[,grep(metric, colnames(wdf))], 2, function(x) as.numeric(x)), stringsAsFactors = FALSE)

  glist <- list()
  for ( i in 3:ncol(xdf) ) {
    rdf <- xdf[, c(2,i)]
    rdf <- rdf[sort.list(rdf[,2], dec=FALSE),]
    glist[[colnames(xdf)[i]]] <- rdf$gene.symbol.xref[which(!is.na(rdf[,2]))]
  }
  rra <- aggregateRanks(glist=glist, method="RRA")
  rra$fdr <- p.adjust(rra$Score, method = "fdr")
  colnames(rra) <- c("gene.symbol.xref", "SpC.Delta.RRA.Score", "SpC.Delta.RRA.Score.FDR")
  rra.df <- merge(xdf, rra, by = "gene.symbol.xref")

  return(rra.df)
}

getArrayMetrics <- function() {
  #

  pstat.rdata.fn <- "data/pstat1_aggregate_intx_stats_v2_v3.Rdata"
  ctg.rdata.fn <- "data/ctg_aggregate_intx_stats_v2_v3.Rdata"

  load(pstat.rdata.fn)
  load(ctg.rdata.fn)

  ctg.rra.df <- ctg.list$rra
  ctg.stats.df <- ctg.list$intx.wdf
  ctg.stats.df <- merge(ctg.stats.df, ctg.rra.df[,c(1:3)], by.x = "gene", by.y = "Gene.Symbol")
  ctg.stats.df[, grep("FDR", colnames(ctg.stats.df))] <- 10^(-1*abs(ctg.stats.df[, grep("FDR", colnames(ctg.stats.df))]))

  pstat.rra.df <- pstat.list$rra
  pstat.stats.df <- pstat.list$intx.wdf
  pstat.stats.df <- merge(pstat.stats.df, pstat.rra.df[,c(1:3)], by.x = "gene", by.y = "Gene.Symbol")
  pstat.stats.df[, grep("FDR", colnames(pstat.stats.df))] <- 10^(-1*abs(pstat.stats.df[, grep("FDR", colnames(pstat.stats.df))]))

  array.genevalues <- updateGeneValuesBatch(ctg.stats.df$gene, type = "symbol")
  ctg.stats.df$Entrez.GeneID <- array.genevalues$entrezgene[match(array.genevalues$query, ctg.stats.df$gene)]
  pstat.stats.df$Entrez.GeneID <- array.genevalues$entrezgene[match(array.genevalues$query, pstat.stats.df$gene)]

  m <- merge(ctg.stats.df, pstat.stats.df[, -which(colnames(pstat.stats.df) == "Entrez.GeneID")], by = "gene")
  array.list <- list(pstat = pstat.stats.df, ctg = ctg.stats.df, all = m)
  return(array.list)
}

getCRISPRia <- function() {
  # Load CRISPRia data object
  # Loads variable param.list
  load("data/parp7_ncih1373_crispria_data.RData")
  crispr.data <- param.list
  crispri.sl.df <- crispr.data$crispri.df
  crispri.sl.df$crispri.lfc.score <- crispri.sl.df$pos.score
  crispri.sl.df$crispri.lfc.score[which(crispri.sl.df$neg.lfc < 0)] <- crispri.sl.df$neg.score[which(crispri.sl.df$neg.lfc < 0)]
  crispri.sl.df$crispri.lfc <- crispri.sl.df$pos.lfc
  crispri.sl.df$crispri.lfc[which(crispri.sl.df$neg.lfc < 0)] <- crispri.sl.df$neg.lfc[which(crispri.sl.df$neg.lfc < 0)]
  # Filter
  crispri.sl.df$hits <- 0
  hit.filter <- which(crispri.sl.df$crispri.lfc.score < 0.01)
  crispri.sl.df$hits[hit.filter] <- 1
  crispri.sl.df <- crispri.sl.df[, c("Entrez.GeneID", "id", "crispri.lfc", "Phenotype.Score", "crispri.lfc.score", "hits")]
  colnames(crispri.sl.df) <- c("Entrez.GeneID", "Gene.symbol", "PARP7i.CRISPRi.LFC", "PARP7i.CRISPRi.Phenotype.Score", "PARP7i.CRISPRi.RRA.Score", "PARP7.CRISPRi.hits")

  crispri.ko.df <- crispr.data$crispri.qc.df
  crispri.ko.df$crispri.lfc.score <- crispri.ko.df$pos.score
  crispri.ko.df$crispri.lfc.score[which(crispri.ko.df$neg.lfc < 0)] <- crispri.ko.df$neg.score[which(crispri.ko.df$neg.lfc < 0)]
  crispri.ko.df$crispri.lfc <- crispri.ko.df$pos.lfc
  crispri.ko.df$crispri.lfc[which(crispri.ko.df$neg.lfc < 0)] <- crispri.ko.df$neg.lfc[which(crispri.ko.df$neg.lfc < 0)]
  # Filter
  crispri.ko.df$hits <- 0
  hit.filter <- which(crispri.ko.df$crispri.lfc.score < 0.01)
  crispri.ko.df$hits[hit.filter] <- 1
  crispri.ko.df <- crispri.ko.df[, c("Entrez.GeneID", "id", "crispri.lfc", "Phenotype.Score", "crispri.lfc.score", "hits")]
  colnames(crispri.ko.df) <- c("Entrez.GeneID", "Gene.symbol", "CRISPRi.LFC", "CRISPRi.Phenotype.Score", "CRISPRi.RRA.Score", "CRISPRi.hits")

  crispra.sl.df <- crispr.data$crispra.df
  crispra.sl.df$crispra.lfc.score <- crispra.sl.df$pos.score
  crispra.sl.df$crispra.lfc.score[which(crispra.sl.df$neg.lfc < 0)] <- crispra.sl.df$neg.score[which(crispra.sl.df$neg.lfc < 0)]
  crispra.sl.df$crispra.lfc <- crispra.sl.df$pos.lfc
  crispra.sl.df$crispra.lfc[which(crispra.sl.df$neg.lfc < 0)] <- crispra.sl.df$neg.lfc[which(crispra.sl.df$neg.lfc < 0)]
  # Filter
  crispra.sl.df$hits <- 0
  hit.filter <- which(crispra.sl.df$crispra.lfc.score < 0.01)
  crispra.sl.df$hits[hit.filter] <- 1
  crispra.sl.df <- crispra.sl.df[, c("Entrez.GeneID", "id", "crispra.lfc", "Phenotype.Score", "crispra.lfc.score", "hits")]
  colnames(crispra.sl.df) <- c("Entrez.GeneID", "Gene.symbol", "PARP7i.CRISPRa.LFC", "PARP7i.CRISPRa.Phenotype.Score", "PARP7i.CRISPRa.RRA.Score", "PARP7.CRISPRa.hits")

  crispra.ko.df <- crispr.data$crispra.qc.df
  crispra.ko.df$crispra.lfc.score <- crispra.ko.df$pos.score
  crispra.ko.df$crispra.lfc.score[which(crispra.ko.df$neg.lfc < 0)] <- crispra.ko.df$neg.score[which(crispra.ko.df$neg.lfc < 0)]
  crispra.ko.df$crispra.lfc <- crispra.ko.df$pos.lfc
  crispra.ko.df$crispra.lfc[which(crispra.ko.df$neg.lfc < 0)] <- crispra.ko.df$neg.lfc[which(crispra.ko.df$neg.lfc < 0)]
  # Filter
  crispra.ko.df$hits <- 0
  hit.filter <- which(crispra.ko.df$crispra.lfc.score < 0.01)
  crispra.ko.df$hits[hit.filter] <- 1
  crispra.ko.df <- crispra.ko.df[, c("Entrez.GeneID", "id", "crispra.lfc", "Phenotype.Score", "crispra.lfc.score", "hits")]
  colnames(crispra.ko.df) <- c("Entrez.GeneID", "Gene.symbol", "CRISPRa.LFC", "CRISPRa.Phenotype.Score", "CRISPRa.RRA.score", "CRISPRa.hits")

  all.sl.df <- unique(merge(crispri.sl.df[,c(1:5)], crispra.sl.df[,c(1:5)], by = c("Entrez.GeneID", "Gene.symbol"), all = TRUE))
  all.ko.df <- unique(merge(crispri.ko.df[,c(1:5)], crispra.ko.df[,c(1:5)], by = c("Entrez.GeneID", "Gene.symbol"), all = TRUE))
  all.df <- merge(all.sl.df, all.ko.df, by = c("Entrez.GeneID", "Gene.symbol"), all = TRUE)
  crispria.cellline.list <- list(crispri.sl.df = crispri.sl.df, 
                                 crispra.sl.df = crispra.sl.df,
                                 crispri.df = crispri.ko.df,
                                 crispra.df = crispra.ko.df,
                                 all.df = all.df)
  return(crispria.cellline.list)
}

collapseDups <- function(dup.df) {
  #

  collapse.cols <- c("protein_desc", "Accession_Number", "Molecular_Weight", "label", "uniprotswissprot", "entrez.id", "genename")
  sum.cols <- colnames(dup.df)[grep("_SpC", colnames(dup.df))]
  avg.cols <- c("ctrl.binding", "trt.binding", "trt.fc", "molweight")
  recalc.cols <- c("spc.ctrl_fraction_delta", "spc.trt_ctrl.delta")
  new.values <- list()
  for ( i in 1:ncol(dup.df) ) {
    cname <- colnames(dup.df)[i]
    if ( cname %in% collapse.cols ) {
      new.value <- as.character(paste(unique(gsub('"', "", dup.df[,i])), collapse = ","))
    } else if ( cname %in% sum.cols ) {
      new.value <- sum(dup.df[,i])
    } else if ( cname %in% avg.cols ) {
      new.value <- mean(dup.df[,i])
    } else if ( cname %in% recalc.cols ) {
      new.value <- NA
    }
    new.values[[cname]] <- new.value
  }
  new.df <- data.frame(lapply(new.values, function(x) t(data.frame(x))))
  spc.idxs <- grep("_SpC", colnames(dup.df))
  new.df$spc.trt_ctrl.delta <- new.df[,spc.idxs[2]] - new.df[,spc.idxs[1]]
  new.df$spc.ctrl_fraction_delta <- (new.df[,spc.idxs[2]] - new.df[,spc.idxs[1]]) / sum(new.df[,spc.idxs])
  return(new.df)
}

pubHeatmap <- function() {
  #

  require("ComplexHeatmap")
  require("RColorBrewer")
  require("circlize")

  subid.df <- getSubIDMerged("spc.trt_ctrl.delta")
  array.list <- getArrayMetrics()
  crispria.list <- getCRISPRia()

  x <- merge(array.list$all, crispria.list$all.df, by= "Entrez.GeneID", all = TRUE)
  x <- merge(x, subid.df, by.x = "Entrez.GeneID", by.y = "entrez.id", all = TRUE)

  gs <- x[,c("gene", "Gene.symbol", "gene.symbol.xref")]
  nsymbols <- apply(gs, 1, function(x) length(unique(x[which(!is.na(x))])))
  gene.symbols <- apply(gs, 1, function(x) unique(x[which(!is.na(x))])[1])
  gene.symbols[which(nsymbols == 2)] <- ifelse(!is.na(gs[which(nsymbols == 2),3]), gs[which(nsymbols == 2),3], gs[which(nsymbols == 2),1])
  x <- x[,-which(colnames(x) %in% c("gene", "Gene.symbol", "gene.symbol.xref"))]
  x$Gene.Symbol <- gene.symbols

  subid.idxs <- which(!is.na(x$SpC.Delta.RRA.Score))
  array.idxs <- which(!is.na(x$CTG.SL.Sensitivity.RRA.FDR))
  ia.idxs <- which(!is.na(x$PARP7i.CRISPRa.RRA.Score) | !is.na(x$PARP7i.CRISPRi.RRA.Score))

  x$subid <- 0
  x$subid[subid.idxs] <- 1
  x$crispr.array <- 0
  x$crispr.array[array.idxs] <- 1
  x$crispria <- 0
  x$crispria[ia.idxs] <- 1

  # fn <- file.path("/home/rabo/github_repos/p7/analysis/integration/parp7_subid_crispr_merged.csv")
  # write.csv(x, file = fn, quote = FALSE, row.names = FALSE)

  # array.fn <- file.path("/home/rabo/github_repos/p7/analysis/integration/parp7_crispr_array.csv")
  # write.csv(array.list$all, file = array.fn, quote = FALSE, row.names = FALSE)

  # subid.fn <- file.path("/home/rabo/github_repos/p7/analysis/integration/parp7_subids.csv")
  # subid.df$NSubIDHits <- apply(subid.df[,grep("trt_ctrl.delta", colnames(subid.df))], 1, function(x) sum(x < 0, na.rm = TRUE))
  # write.csv(subid.df, file = subid.fn, quote = FALSE, row.names = FALSE)

  # Publication reduction
  # 1. Keep overlap among all datasets, remove pSTAT1 array results
  x <- x[,-c(12:17)]
  rx <- x[which(x$subid == 1 & x$crispr.array == 1 & x$crispria == 1),]

  # 1. Reduce CTG synthetic lethality 
  array.intx.col.idxs1 <- intersect(grep("mult", colnames(x)), grep("Synthetic", colnames(x)))
  ctg.mult.summary.values <- apply(rx[,array.intx.col.idxs1], 1, function(x) ifelse(sum(x < 0) > sum(x > 0), x[x<0][which.max(abs(x[x<0]))], x[x>0][which.max(abs(x[x>0]))]))
  array.intx.col.idxs2 <- intersect(grep("Synthetic.lethality", colnames(rx)), grep("add", colnames(x)))
  ctg.add.summary.values <- apply(rx[,array.intx.col.idxs2], 1, function(x) x[which.max(abs(x))])

  subid.hits <- apply(rx[,grep("trt_ctrl.delta", colnames(rx))], 1, function(x) sum(x < 0, na.rm = TRUE))

  df <- data.frame(Gene.Symbol = rx$Gene.Symbol, CRISPRi.SL.Score = rx$PARP7i.CRISPRi.Phenotype.Score, CRISPRa.SL.Score = rx$PARP7i.CRISPRa.Phenotype.Score, NSubExpHits = subid.hits, CRISPRarray.SL.mult.max = ctg.mult.summary.values)
  df <- df[which(df$NSubExpHits > 0),]

  fn <- file.path("/home/rabo/github_repos/p7/analysis/integration/parp7_subid_crispr_merged.all_overlap.annotated_reduced1.tsv")
  write.table(df, file = fn, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  df2 <- df[,c("Gene.Symbol", "CRISPRi.SL.Score", "CRISPRarray.SL.mult.max", "NSubExpHits")]
  df2$dir <- apply(df2[,c(2,3)], 1, function(x) ifelse(sum(x < 0) == 2, "Sensitization", ifelse(sum(x>0)==2, "Resistant", "None")))
  df2 <- df2[-which(df2$dir == "None"),]
  df2 <- df2[sort.list(df2$dir),]

  annotation.list <- list("AHR" = "Transcription factor", 
                          "TRIM32" = "Innate immune response", 
                          "ADAR" = "Innate immune response",
                          "ATRX" = "DNA helicase",
                          "DDX21" = "Innate immune response, helicase",
                          "DDX3X" = "RNA helicase",
                          "DDX41" = "Innate immune response, cytosolic sensor, helicase",
                          "DDX46" = "RNA helicase",
                          "DHX9" = "DNA helicase",
                          "EIF4A3" = "Innate immune response, RNA helicase",
                          "EP300" = "Innate immune response, RIG-I mediated",
                          "IFI16" = "Innate immune response, cytosolic sensor",
                          "JAK1" = "Innate immune response",
                          "PARP4" = "Protein ADP-ribosylation",
                          "PRKDC" = "Innate immune response, cytosolic sensor",
                          "TIPARP" = "Protein ADP-ribosylation",
                          "XRCC6" = "Innate immune response, cytosolic sensor",
                          "ZC3HAV1" = "Protein ADP-ribosylation")

  df2$annotation <- ""
  for ( gene in names(annotation.list) ) {
    df2$annotation[which(df2$Gene.Symbol == gene)] <- annotation.list[[gene]]
  }

  df.mat <- df2[,c(2,3)]
  rownames(df.mat) <- df2$Gene.Symbol
  colnames(df.mat) <- c("CRISPRi", "CRISPRarray")

  ha_anno = rowAnnotation(text = row_anno_text(df2$annotation, just = "left", gpar(fontsize = 8)), annotation_width = unit(0.1, "cm"), width = max_text_width(df2$annotation, gp = gpar(fontsize = 9)))
  hax <- Heatmap(df2$NSubExpHits, name = "RBN011364 SubID Exp Hits", col = colorRamp2(c(0, 5), c("white", "darkorange3")), width = unit(5, "mm"), show_column_names = FALSE)
  hay <- Heatmap(df2$dir, name = "RBN-2397 Synth.Lethal Direction", width = unit(5, "mm"), show_column_names = FALSE, col = c("goldenrod3", "darkorchid4"))
  # ha_meta = rowAnnotation(df = data.frame(Nhits = df2$NSubExpHits, Dir = df2$dir), width = unit(1.5, "cm"), show_annotation_name = TRUE)
  ha2 = HeatmapAnnotation(cn = anno_text(colnames(df.mat), rot = 0, gp = gpar(fontsize = 10)))
  # hm.cols <- colorRamp2(c(-2, 0, 1), c("royalblue4", "khaki1", "red4"))
  hm <- Heatmap(as.matrix(df.mat), 
                 name="RBN-2397 Synth.Lethal Score", 
                 km=1, 
                 # col=hm.cols,
                 cluster_rows = FALSE,
                 top_annotation = ha2,
                 cluster_columns = FALSE,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE,
                 row_names_side = "left",
                 column_title = "",
                 column_names_side = "top",
                 show_row_names = TRUE,
                 show_column_names = FALSE,
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 heatmap_legend_param=list(legend_direction="vertical", color_bar="continuous"))
  plot.fn <- "figs/Fig6D.parp7_crispr_subid_merged.chm.pdf")
  pdf(plot.fn, height = 4, width = 9)
  draw(hm + hax + hay + ha_anno, heatmap_legend_side = "right")
  dev.off()
}

pubHeatmap()
