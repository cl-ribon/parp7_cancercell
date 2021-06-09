
publicationPARP7ExpPARP7iResponseGroups <- function(param.list) {
  #

  require("ggplot2")
  require("ggbeeswarm")
  require("viridis")

  biomarker.df <- read.table("data/parp7_cell_line_biomarker_data.tsv")
  biomarker.df$Group.label <- "PARP7i Non-responder"
  biomarker.df$Group.label[which(biomarker.df$Response.group == "Responder")] <- "PARP7i Responder"
  gp2 <- ggplot(biomarker.df, aes(x = Group.label, y = TIPARP.log2TPM, color = ISG.score)) +
         geom_quasirandom(size = 2) +
         stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.5) +
         xlab("") +
         ylab("PARP7 Expression (log2 TPM)") +
         scale_color_viridis() +
         guides(size = FALSE) +
         labs(color = "ISG Expression Score") +
         theme_minimal() +
         theme(text = element_text(size = 16), axis.line = element_line(size = 1), legend.position = "bottom")

  # gp3 <- ggplot(biomarker.df, aes(x = horizon.PARP7.RBN011364.GI50_nm, y = TIPARP.log2TPM, color = ISG.score)) +
  #        geom_point(alpha = 0.8) +
  #        geom_smooth(method = "lm", se = FALSE, color = "darkgrey") +
  #        xlab(expression(paste("PARP7i (RBN-011364) GI" ["50"], " (nM)"))) +
  #        ylab("PARP7 Expression (log2 TPM)") +
  #        scale_color_viridis() +
  #        guides(size = FALSE) +
  #        labs(color = "ISG Expression Score") +
  #        theme_minimal() +
  #        theme(text = element_text(size = 16), axis.line = element_line(size = 1), legend.position = "bottom")

  plot.fn <- file.path("figs/Fig5C.parp7_response_group_vs_mrna_isg.pub.pdf")
  pdf(plot.fn, height = 6, width = 6)
  # grid.arrange(gp2, gp3, ncol = 2)
  plot(gp2)
  # plot(gp3)
  dev.off()
}