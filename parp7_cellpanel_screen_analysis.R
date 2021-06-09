#! /usr/bin/Rscript


# Function contains code to define response groups
loadCellLineResponseData.pub <- function(param.list, data.list, plot.pca = FALSE) {
  #

  require("cluster")

  # Load response metric file
  # response.metrics.df <- load_response_metrics()
  response.metrics.df <- read.table("data/parp7.horizon_cell_line_response_metrics_wdf.publication.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Format the PARP7 specific metrics
  p7.metrics.df <- subset(response.metrics.df, grepl("PARP7", response.metrics.df$Compound.Label))
  p7.metrics.df$norm.name <- normalize_cellline_name(p7.metrics.df$cell.line)

  p7.select.cmpds <- c("RBN011628-001", "RBN011595-001", "RBN011364-002")
  select.metrics <- c("norm.name", "compound", "GI50_nm", "IC50_nm", "Delta.1um", "EMax.1um.fixed")
  p7.select.df <- subset(p7.metrics.df, compound %in% p7.select.cmpds, select = select.metrics)

  # Create a wide data.frame with columns for each PARP7 compound and corresponding response metric values.
  wdf.list <- list()
  for ( metric in select.metrics[-c(1,2)] ) {
    print(metric)
    wdf <- dcast(p7.select.df, norm.name ~ compound, value.var = metric)
    colnames(wdf)[-1] <- paste(colnames(wdf)[-1], "horizon", metric, sep=".")
    wdf.list[[metric]] <- wdf
  }
  # Merge all metrics into one data.frame
  p7.horizon.metrics.wdf <- Reduce(function(df1, df2) merge(df1, df2, by="norm.name", all.x=TRUE), wdf.list)

  # Cluster profiled cell lines by all response metrics
  cluster.df <- p7.horizon.metrics.wdf[,-1]
  rownames(cluster.df) <- p7.horizon.metrics.wdf$norm.name
  colnames(cluster.df) <- toupper(letters[1:ncol(cluster.df)])
  set.seed(123)
  pam.3cluster <- pam(cluster.df, 3, stand = TRUE)
  # m1$pam.cluster <- pam.df$clustering
  pam.2cluster <- pam(cluster.df, 2, stand=TRUE)
  # m1$pam.2cluster.cluster <- pam.2cluster.df$clustering

  # Merge the CCLE meta data with the PARP7 reponse metrics
  # m1 <- merge(p7.horizon.metrics.wdf, ccle.meta, by="norm.name", all.x=T)

  # Write metrics and response cluster to file
  metrics.df <- p7.horizon.metrics.wdf
  metrics.df$Response.group <- ifelse(pam.3cluster$clustering == 2, "Nonresponder", "Responder")
  write.csv(metrics.df, file = "/home/rabo/github_repos/p7/analysis/cell_line_isg/parp7.horizon_cell_line_response_metrics.csv", quote=FALSE)

  # Added for supplemental table for cancer cell resubmission
  niraparib.metrics.df <- subset(response.metrics.df, grepl("PARP1.RBN010061", response.metrics.df$Compound.Label))
  niraparib.metrics.df$norm.name <- normalize_cellline_name(niraparib.metrics.df$cell.line)
  select.metrics <- c("norm.name", "compound", "GI50_nm", "IC50_nm", "Delta.1um", "EMax.1um.fixed")
  niraparib.select.df <- subset(niraparib.metrics.df, select = select.metrics)
  colnames(niraparib.select.df)[-c(1,2)] <- paste("niraparib", "horizon", colnames(niraparib.select.df)[-c(1,2)], sep=".")
  niraparib.select.df <- niraparib.select.df[,-2]
  metrics.df <- merge(metrics.df, niraparib.select.df, by = "norm.name")
  write.csv(metrics.df, file = "/home/rabo/github_repos/p7/analysis/cell_line_isg/parp7.horizon_cell_line_response_metrics_niraparib.csv", quote=FALSE)

  if ( plot.pca ) {
    x.df <- prcomp(cluster.df, scale=T, center=T)
    pca.df <- data.frame(PC1=x.df$x[,1], PC2=x.df$x[,2], pam.cluster=pam.3cluster$clustering, response.group = ifelse(pam.3cluster$clustering == 2, "nonresponder", "responder"))
    gp.pca <- ggplot(pca.df, aes(x=PC1, y=PC2, color=as.factor(response.group), label=rownames(x.df$x))) +
              geom_point(size = 3, alpha = 0.8) +
              scale_color_manual(values = c("responder" = "#e838eb", "nonresponder" = "darkgrey")) +
              theme_minimal() +
              theme(legend.title = element_blank())
    plot.fn <- file.path("/home/rabo/github_repos/p7/analysis/cell_line_isg", "parp7.horizon_cell_line_responses.pca.pdf")
    pdf(plot.fn, width = 8, height = 6)
    plot(gp.pca)
    dev.off()
  }

  response.data.list <- list(p7.metrics.df = p7.horizon.metrics.wdf, 
                             p7.3cluster = pam.3cluster,
                             p7.2cluster = pam.2cluster)
}

# Function to generate waterfall plot Fig 5A
parp7_niraparib_waterfall_gi50 <- function() {
  #
  # Read in source data
  # response.metrics.df <- read.table("data/parp7_parp1_response_metrics.tsv", sep = "\t", header = TRUE)

  ldf <- melt(data.list$response.metrics.df[,c("cell.line", "compound", "GI50_nm")], id.vars = c("cell.line", "compound"))
  ldf <- ldf[which(ldf$compound %in% c("RBN011364-002", "RBN010061-002")),]
  ldf$compound <- gsub("RBN011364-002", "PARP7 inhibitor (RBN011364)", ldf$compound)
  ldf$compound <- gsub("RBN010061-002", "Niraparib", ldf$compound)
  ldf$compound <- factor(ldf$compound, levels=c("PARP7 inhibitor (RBN011364)", "Niraparib"))
  p7.order <- ldf[grep("RBN011364", ldf$compound),]
  p7.order <- p7.order[sort.list(p7.order$value, dec = FALSE),]
  ldf$cell.line <- factor(ldf$cell.line, levels = unique(p7.order$cell.line))
  gp.gi50_pt <- ggplot(ldf, aes(x = cell.line, y=value, color=compound)) +
                geom_hline(yintercept=c(1, 10, 100, 1000, 10000), col="lightgrey", alpha=0.6) +
                geom_point(alpha=0.75, size=3) + 
                scale_color_manual(breaks=c("Niraparib", "PARP7 inhibitor (RBN011364)"), values=c("#E84C22", "#011589")) +
                scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
                theme_minimal() +
                xlab("Cell Lines") +
                ylab(expression(paste("GI" ["50"], " (nM)"))) +
                theme(axis.ticks.x=element_blank(), 
                      axis.text.x = element_blank(), 
                      text=element_text(size = 18), 
                      axis.text.y=element_text(size = 20), 
                      legend.title=element_blank(), 
                      panel.grid.major=element_blank(), 
                      legend.text=element_text(size=22), 
                      legend.position=c(0.7, 0.25),
                      axis.line.x = element_line(size=1, color = "black"), 
                      axis.line.y = element_line(size=1, color = "black"))
  plot.fn <- file.path(param.list$report.dir, "parp7", "figs/Fig5A.parp7_compound_screened_cell_lines.horizon.niraparib_pub.waterfall.pdf")
  pdf(plot.fn, width = 8, height = 6)
  plot(gp.gi50_pt)
  dev.off()
}