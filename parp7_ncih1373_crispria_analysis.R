#! /usr/bin/Rscript

crisprScatter <- function(param.list, select.geneids, label) {
  # Function generates Fig6A scatter plot of CRISPRi/a results

  require("ggplot2")

  merged.df <- param.list$crispr.merged.df
  resistant.df <- merged.df[which(merged.df$crispri.Phenotype.Score > 0 & merged.df$crispra.Phenotype.Score < 0),]
  resistant.df$Group <- "Resistant"
  sens.df <- merged.df[which(merged.df$crispri.Phenotype.Score < 0 & merged.df$crispra.Phenotype.Score > 0),]
  sens.df$Group <- "Sensitive"
  color.df <- rbind(resistant.df, sens.df)
  select.df <- merged.df[match(select.geneids, merged.df$Entrez.GeneID),]
  gp.ia <- ggplot(merged.df, aes(x = crispra.Phenotype.Score, y = crispri.Phenotype.Score)) +
           geom_point(alpha = 0.5, color = "lightgrey") +
           geom_point(data = color.df, aes(x = crispra.Phenotype.Score, y = crispri.Phenotype.Score, color = Group), alpha = 0.5) +
           geom_vline(xintercept = 0, linetype = "dotted", size = 0.5, color = "darkgrey") +
           geom_hline(yintercept = 0, linetype = "dotted", size = 0.5, color = "darkgrey") +
           geom_point(data = select.df, aes(x = crispra.Phenotype.Score, y = crispri.Phenotype.Score), color = "orangered", alpha = 0.75, size = 2) +
           xlab("CRISPRa Phenotype Score") +
           ylab("CRISPRi Phenotype Score") +
           scale_color_manual(values = c("Resistant" = "dodgerblue3", "Sensitive" = "darkorange3")) +
           theme_minimal() +
           theme(text = element_text(size = 18), legend.position = "none", plot.title = element_text(size = 20), legend.title = element_blank())

  plot.fn <- "figs/Fig6A.crispria_phenotype_scores.", label, "_scatter.pdf", sep=""))
  pdf(plot.fn, width = 8, height = 8)
  plot(gp.ia)
  dev.off()
}

crispriaWaterfall <- function(param.list) {
  # Function to plot GSEA results as barplot - Fig6B

  require("ggplot2")

  fgsea.df <- param.list$ldiag.fgsea.df
  gs.df <- fgsea.df[which(fgsea.df$padj <= 0.05),]
  gs.df$Group <- "Res"
  gs.df$Group[which(gs.df$NES < 0)] <- "Sens"

  sens.df <- gs.df[which(gs.df$NES < 0),]
  res.df <- gs.df[which(gs.df$NES >= 0),]
  sens.df <- sens.df[order(sens.df$padj, decreasing=TRUE),]
  res.df <- res.df[order(res.df$padj),]
  plot.df <- rbind(res.df, sens.df)
  plot.df$pathway <- factor(plot.df$pathway, levels = rev(plot.df$pathway))
  gp.cia.fgsea <- ggplot(plot.df, aes(x = pathway, y = NES)) +
                  geom_bar(stat = "identity", fill = "grey") +
                  xlab("") +
                  ylab("Normalized Enrichment Score") +
                  ggtitle("CRISPRia Enriched Genesets") +
                  geom_bar(data = sens.df, aes(x = pathway, y = NES, alpha = as.numeric(-log10(padj))), stat = "identity", fill = "darkorange3") +
                  geom_bar(data = res.df, aes(x = pathway, y = NES, alpha = as.numeric(-log10(padj))), stat = "identity", fill = "dodgerblue3") +
                  scale_alpha_continuous() +
                  coord_flip() +
                  theme_minimal() +
                  theme(axis.text.y = element_blank(), text = element_text(size = 16), legend.title = element_blank(), legend.text=element_text(size=12), legend.position = "none", panel.grid = element_blank())

  plot.fn <- "figs/Fig6B.crispria_fgsea_all.bar.pdf"
  pdf(plot.fn, width = 5, height = 6)
  plot(gp.cia.fgsea)
  dev.off()
}

# Loads varable param.list object with all CRISPRi/a data
load("data/parp7_ncih1373_crispria_data.RData")
crisprScatter(param.list, param.list$ldiag.fgsea.df$leadingEdge[which(param.list$ldiag.fgsea.df$pathway == "GO_INNATE_IMMUNE_RESPONSE")][[1]], "innate_immune_response")
crispriaWaterfall(param.list)