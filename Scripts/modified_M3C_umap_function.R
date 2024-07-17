modified_M3C_umap = function (mydata, labels = FALSE, printres = FALSE, seed = FALSE,
          axistextsize = 18, legendtextsize = 18, dotsize = 5, textlabelsize = 4,
          legendtitle = "Group", controlscale = FALSE, scale = 1,
          low = "grey", high = "red", colvec = c("skyblue", "gold",
                                                 "violet", "darkorchid", "slateblue", "forestgreen",
                                                 "violetred", "orange", "midnightblue", "grey31", "black"),
          printheight = 20, printwidth = 22, text = FALSE,
          umap_config)
{
  if (controlscale == TRUE && class(labels) %in% c("character",
                                                   "factor") && scale %in% c(1, 2)) {
    stop("when categorical labels, use scale=3")
  }
  if (controlscale == TRUE && class(labels) %in% c("numeric") &&
      scale %in% c(3)) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if (controlscale == FALSE && scale %in% c(2, 3)) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c("character",
                                                     "factor")) {
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- "Unknown"
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c("character",
                                                 "factor")) {
    warning("there is NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- "Unknown"
  }
  message("***UMAP wrapper function***")
  message("running...")
  if (seed != FALSE) {
    set.seed(seed)
  }
  if (labels[1] == FALSE && text[1] == FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = umap_config)
    scores <- data.frame(umap$layout)
    # colnames(scores)[which(colnames(scores) == "X1")] = "UMAP 1"
    # colnames(scores)[which(colnames(scores) == "X2")] = "UMAP 2"
    p <- ggplot(data = scores, aes(x = X1, y = X2)) +
      geom_point(colour = "skyblue",size = dotsize) +
      theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.text.y = element_text(size = axistextsize,colour = "black"),
                         axis.text.x = element_text(size = axistextsize,colour = "black"),
                         axis.title.x = element_text(size = axistextsize),
                         axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = colvec)
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAP.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] == FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = umap_config)
    scores <- data.frame(umap$layout)
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot(data = scores, aes(x = X1, y = X2)) +
          geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize),
                             axis.title.y = element_text(size = axistextsize),
                             legend.title = element_text(size = legendtextsize),
                             legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_color_viridis_c()
      }
      else if (scale == 2) {
        p <- ggplot(data = scores, aes(x = X1, y = X2)) +
          geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize),
                             axis.title.y = element_text(size = axistextsize),
                             legend.title = element_text(size = legendtextsize),
                             legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_colour_gradient(low = low,
                                                             high = high)
      }
      else if (scale == 3) {
        labels_df = data.frame(label = as.character(labels), size = 2, alpha = 0.8)
        labels_df$size = ifelse(labels_df$label == "Target", 3, 2)
        labels_df$alpha = ifelse(labels_df$label == "Target", 1, 0.3)
        p <- ggplot(data = scores, aes(x = X1, y = X2)) +
          geom_point(aes(colour = labels), size = labels_df$size, alpha = labels_df$alpha) +
          theme_bw() + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize),
                             axis.title.y = element_text(size = axistextsize),
                             legend.title = element_text(size = legendtextsize),
                             legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_colour_manual(values = colvec)
      }
    }
    else {
      p <- ggplot(data = scores, aes(x = X1, y = X2)) +
        geom_point(aes(colour = labels), size = dotsize) +
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                          colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                        colour = "black"), axis.title.x = element_text(size = axistextsize),
                           axis.title.y = element_text(size = axistextsize),
                           legend.title = element_text(size = legendtextsize),
                           legend.text = element_text(size = legendtextsize)) +
        labs(colour = legendtitle)
    }
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAPlabeled.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] != FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = umap_config)
    scores <- data.frame(umap$layout)
    scores$label <- text
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot(data = scores, aes(x = X1, y = X2,
                                       label = label)) + geom_point(aes(colour = labels),
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize),
                                                                                                         axis.title.y = element_text(size = axistextsize),
                                                                                                         legend.title = element_text(size = legendtextsize),
                                                                                                         legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_color_viridis_c() +
          geom_text(vjust = "inward", hjust = "inward",
                    size = textlabelsize)
      }
      else if (scale == 2) {
        p <- ggplot(data = scores, aes(x = X1, y = X2,
                                       label = label)) + geom_point(aes(colour = labels),
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize),
                                                                                                         axis.title.y = element_text(size = axistextsize),
                                                                                                         legend.title = element_text(size = legendtextsize),
                                                                                                         legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_colour_gradient(low = low,
                                                             high = high) + geom_text(vjust = "inward",
                                                                                      hjust = "inward", size = textlabelsize)
      }
      else if (scale == 3) {
        p <- ggplot(data = scores, aes(x = X1, y = X2,
                                       label = label)) + geom_point(aes(colour = labels),
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize),
                                                                                                         axis.title.y = element_text(size = axistextsize),
                                                                                                         legend.title = element_text(size = legendtextsize),
                                                                                                         legend.text = element_text(size = legendtextsize)) +
          labs(colour = legendtitle) + scale_colour_manual(values = colvec) +
          geom_text(vjust = "inward", hjust = "inward",
                    size = textlabelsize)
      }
    }
    else {
      p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) +
        geom_point(aes(colour = labels), size = dotsize) +
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                                                                                          colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                                        colour = "black"), axis.title.x = element_text(size = axistextsize),
                           axis.title.y = element_text(size = axistextsize),
                           legend.title = element_text(size = legendtextsize),
                           legend.text = element_text(size = legendtextsize)) +
        labs(colour = legendtitle) + geom_text(vjust = "inward",
                                               hjust = "inward", size = textlabelsize)
    }
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAPlabeled.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] == FALSE && text[1] != FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = umap_config)
    scores <- data.frame(umap$layout)
    scores$label <- text
    p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) +
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))),
                 size = dotsize) + theme_bw() + theme(legend.position = "none",
                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                      axis.text.y = element_text(size = axistextsize,
                                                                                 colour = "black"), axis.text.x = element_text(size = axistextsize,
                                                                                                                               colour = "black"), axis.title.x = element_text(size = axistextsize),
                                                      axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = colvec) + geom_text(vjust = "inward",
                                                       hjust = "inward", size = textlabelsize)
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAP.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  message("done.")
  p = p +
    xlab("UMAP 1") +
    ylab("UMAP 2")
  return(p)
}
