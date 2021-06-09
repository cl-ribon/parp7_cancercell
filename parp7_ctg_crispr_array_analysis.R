#! /usr/bin/Rscript

require("ggplot2")

runRRA <- function(wdf, direction = "neg") {
  # 

  require("RobustRankAggreg")

  dec.val <- ifelse(direction == "neg", FALSE, TRUE)
  glist <- list()
  for ( i in 2:ncol(wdf) ) {
    rdf <- wdf[, c(1,i)]
    rdf <- rdf[sort.list(rdf[,2], dec = dec.val),]
    glist[[colnames(wdf)[i]]] <- rdf$gene[which(!is.na(rdf[,2]))]
  }
  rra <- aggregateRanks(glist=glist, method="RRA")
  rra$fdr <- p.adjust(rra$Score, method = "BH")
  return(rra)
}

# Load the CTG data list 
# Loads variable ctg.list <- list(rra = rra.intx.fdr, intx.wdf = intx.wdf, all.wdf = all.wdf)
load("data/ctg_aggregate_intx_stats_v2_v3.Rdata")

intx.wdf <- intx.wdf[-which(intx.wdf$gene == "Lethal"),]
trt.wdf <- trt.wdf[-which(trt.wdf$gene == "Lethal"),]

rra.intx.sens <- runRRA(intx.wdf, "neg")
colnames(rra.intx.sens) <- c("gene", "CTG.SL.Sensitivity.RRA.Score", "CTG.SL.Sensitivity.RRA.FDR")
rra.intx.res <- runRRA(intx.wdf, "pos")
colnames(rra.intx.res) <- c("gene", "CTG.SL.Resistance.RRA.Score", "CTG.SL.Resistance.RRA.FDR")
rra.intx <- merge(rra.intx.sens, rra.intx.res, by = "gene")
rra.intx <- merge(rra.intx, intx.wdf, by = "gene")

rra.intx.fdr <- -log10(rra.intx[,c(3,5)])
rownames(rra.intx.fdr) <- rra.intx$gene
rra.mult <- rra.intx[,grep("mult", colnames(rra.intx))]
rownames(rra.mult) <- rra.intx$gene
rra.add <- rra.intx[,grep("add", colnames(rra.intx))]
rownames(rra.add) <- rra.intx$gene

# Order of genes
mintx <- all.wdf[which(all.wdf$Stat.method %in% c("mult.intx") & all.wdf$Variable.type == "Synthetic.lethality"),]
aintx <- all.wdf[which(all.wdf$Stat.method %in% c("add.intx") & all.wdf$Variable.type == "Synthetic.lethality"),]
mintx.l2fc.mean <- aggregate(data = mintx, logFC ~ Gene.Symbol, function(x) log2(mean(exp(x))))
colnames(mintx.l2fc.mean) <- c("Gene.Symbol", "mean.log2FC")
aintx.delta.mean <- aggregate(data = aintx, logFC ~ Gene.Symbol, function(x) mean(x))
colnames(aintx.delta.mean) <- c("Gene.Symbol", "mean.Delta")
intx.mean.df <- merge(mintx.l2fc.mean, aintx.delta.mean, by = "Gene.Symbol")

intx.mean.df <- intx.mean.df[sort.list(intx.mean.df$mean.log2FC, dec=FALSE),]
# z$Gene.Symbol <- factor(z$Gene.Symbol, levels = zz$Gene.Symbol)
# zz <- aggregate(data = z, logFC ~ Gene.Symbol, mean)
# zz <- zz[sort.list(zz$logFC, dec=F),]
# z$Gene.Symbol <- factor(z$Gene.Symbol, levels = zz$Gene.Symbol)

all.wdf$Gene.Symbol <- factor(all.wdf$Gene.Symbol, levels = intx.mean.df$Gene.Symbol)
all.wdf$Variable.label <- paste(all.wdf$Variable.type, all.wdf$Stat.method, sep=".")

rra.intx.fdr[,grep("Sensitivity", colnames(rra.intx.fdr))] <- rra.intx.fdr[,grep("Sensitivity", colnames(rra.intx.fdr))]*-1
rra.intx.fdr$max.fdr <- apply(rra.intx.fdr, 1, function(x) x[which.max(abs(x))])
rra.intx.fdr$Gene.Symbol <- rownames(rra.intx.fdr)
rra.intx.fdr <- merge(rra.intx.fdr, intx.mean.df, by = "Gene.Symbol")
rra.intx.fdr <- rra.intx.fdr[sort.list(rra.intx.fdr$max.fdr, dec=FALSE),]
rra.intx.fdr$sig.cutoff <- ifelse(abs(rra.intx.fdr$max.fdr) > -log10(0.05), "Sig", "NS")
rra.intx.fdr$sens.dir <- apply(rra.intx.fdr[,c("mean.log2FC", "mean.Delta")], 1, function(x) sum(x < 0))
rra.intx.fdr$res.dir <- apply(rra.intx.fdr[,c("mean.log2FC", "mean.Delta")], 1, function(x) sum(x > 0))

# rra.intx.fdr$Gene.Symbol <- factor(rra.intx.fdr, levels = rra.intx.fdr$Gene.Symbol)

sig.sens.df <- rra.intx.fdr[which(rra.intx.fdr$sig.cutoff == "Sig" & rra.intx.fdr$sens.dir == 2 & rra.intx.fdr$max.fdr < 0),]
sig.res.df <- rra.intx.fdr[which(rra.intx.fdr$sig.cutoff == "Sig" & rra.intx.fdr$res.dir == 2 & rra.intx.fdr$max.fdr > 0),]

rra.intx.fdr$sig.label <- "NS.RRA"
rra.intx.fdr$sig.label[which(rra.intx.fdr$Gene.Symbol %in% sig.sens.df$Gene.Symbol)] <- "sig.sens.RRA"
rra.intx.fdr$sig.label[which(rra.intx.fdr$Gene.Symbol %in% sig.res.df$Gene.Symbol)] <- "sig.res.RRA"

require("ggrepel")
gp.rra.intx <- ggplot(rra.intx.fdr, aes(x = reorder(Gene.Symbol, max.fdr, FUN = min), y = max.fdr, fill = sig.label)) +
               geom_bar(stat = "identity") +
               # geom_text_repel(data = sig.sens.df, aes(x = Gene.Symbol, y = max.fdr, label = Gene.Symbol), direction = "y", nudge_x = 35, segment.alpha = 0.5, segment.size = 0.5, segment.color = "darkgrey", size= 3, force = 2) +
               # geom_text_repel(data = sig.res.df, aes(x = Gene.Symbol, y = max.fdr, label = Gene.Symbol), direction = "y", ylim = c(-5, 13), nudge_x = -35, segment.alpha = 0.5, segment.size = 0.5, segment.color = "darkgrey", size = 3) +
               scale_fill_manual(values = c("sig.res.RRA" = "dodgerblue2", "NS.RRA" = "darkgrey", "sig.sens.RRA" = "#FF8C00")) +
               xlab("Genes") +
               ylab("RRA Score Signifcance (-log10)") +
               theme_minimal() +
               theme(axis.text.x = element_blank(), text = element_text(size = 18), legend.position = "none", panel.grid = element_blank(), axis.line = element_line(size = 0.25))

fn <- "figs/ctg_aggregate_intx_stats_v2_v3.rra_bar.pdf"
pdf(fn, width = 6, height = 5)
plot(gp.rra.intx)
dev.off()
