# Supplementary Figures ---------------------------------------------------
make_figure_S2 = function(){
  boxplot_data = list()
  for (i in 1:length(intens_distr_comp_all_fdr_10)){
    boxplot_data[[i]] = intens_distr_comp_all_fdr_10[[i]][["qq_test"]]
    boxplot_data[[i]][["ds_id"]] = names(intens_distr_comp_all_fdr_10)[i]
  }
  boxplot_data = dplyr::bind_rows(boxplot_data)
  boxplot_data$log_pval = -1 * log10(boxplot_data$p.value)
  boxplot_data$log_pval[is.infinite(boxplot_data$log_pval)] = 4
  boxplot_data$q = as.character(boxplot_data$q)

  boxplot_data = boxplot_data %>% dplyr::select(q, n1, n2, ds_id, p.value, p.crit, log_pval, `est1-est.2`)
  boxplot_data$delta_type = ifelse(boxplot_data$`est1-est.2` > 0, "Positive", "Negative")
  boxplot_data$sig_type = ifelse(boxplot_data$p.value < boxplot_data$p.crit, "Sig", "Non-Sig")
  boxplot_data = boxplot_data[!is.infinite(boxplot_data$`est1-est.2`),]

  beeswarm(boxplot_data$`est1-est.2` ~ boxplot_data$q,
           pch = 19,
           corral = "wrap", method = "swarm",
           pwcol = as.numeric(factor(boxplot_data$sig_type)),xlab = "Quantiles",
           ylab = "Delta", main = "Distribution of HD differences")
  # legend("bottomleft",legend = levels(factor(boxplot_data$sig_type)),
  #        col = 1:2, pch = 19, title = "",horiz = T, cex = 1,y.intersp = -0.7,
  #        adj = -1.2, bty = "n",title.adj = 0, x.intersp = -1.5)
} #beeswarm HD qq US legal landscape
make_figure_S3 = function(FDR_pct = 10,
                          plot_eval = c("MAP", "annot")){
  make_MAP_boxplot = function(df){
    df = df[which(df$metric_type == "MAP"),]
    p = ggboxplot(df, x = "features_type", y = "metric_value",
                  color = "score_type", palette = c("#A6CEE3", "#FB9A99"), add = "jitter",
                  facet.by = c("catmonotonic"), short.panel.labs = F, ylab = "MAP") +
      theme_pubr() +
      theme(strip.background = element_rect(fill = "white", linetype = "dashed"),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16)) +
      scale_y_continuous(n.breaks = 10)
    return(p)
  }
  make_annot_boxplot = function(annot_df,FDR_pct = 10){
    annot_df = annot_df[which(annot_df$FDR_pct == FDR_pct),]
    annot_df = annot_df %>% tidyr::separate(group_name, into = c("ds_id", "adduct"),
                                            sep = ",")
    annot_df = annot_df %>% group_by(features_type, catmonotonic, ds_id) %>%
      summarise(MSM = sum(msm_fdr_annots),
                `METASPACE-ML` = sum(pred_fdr_annots)) %>%
      as.data.frame()
    annot_df$Delta = annot_df$`METASPACE-ML` / annot_df$MSM
    annot_df$Delta[is.infinite(annot_df$Delta)] = max(annot_df$Delta[!is.infinite(annot_df$Delta)]) + 10
    annot_df = annot_df[order(annot_df$Delta, decreasing = T),]
    annot_df$LFC = log2(annot_df$Delta)
    annot_df$LFC[is.infinite(annot_df$LFC)] = 0

    p = ggboxplot(annot_df, x = "features_type", y = "LFC",
                  color = "catmonotonic", palette = c("#A6CEE3", "#FB9A99"), add = "jitter",
                  ylab = "MAP") +
      theme_pubr() +
      theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16))
    return(p)

  }

  combs = create_combs(all_param_optim_eval, c("includes_bad_rank", "ranking_type"))

  plot_list_MAP = list()
  plot_list_annot = list()
  for (i in 1:nrow(combs)){
    plt_name = paste0(combs[i,1], "_", combs[i,2])
    plot_list_MAP[[plt_name]] = make_MAP_boxplot(df = plyr::match_df(all_param_optim_eval,
                                                                     combs[i,]))
    plot_list_annot[[plt_name]] = make_annot_boxplot(annot_df = plyr::match_df(all_param_optim_annot,
                                                                               combs[i,]),
                                                     FDR_pct = FDR_pct)
  }

  if (plot_eval == "MAP"){
    pm = ggmatrix(
      plot_list_MAP,
      nrow = 2, ncol = 2,
      xAxisLabels = c("Exclude Low annot adducts: Yes", "Exclude Low annot adducts : No"),
      yAxisLabels = c("Ranking type : Separate", "Ranking type : Single"),
      showStrips = T,
      legend = c(1,2),
      showAxisPlotLabels = T
    ) +
      theme(strip.background = element_rect(fill = "white", linetype = "dashed"),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_blank(),
            legend.position = "bottom")
    return(pm)
  }
  else{
    pm = ggmatrix(
      plot_list_annot,
      nrow = 2, ncol = 2,
      xAxisLabels = c("Exclude Low annot adducts: Yes", "Exclude Low annot adducts : No"),
      yAxisLabels = c("Ranking type : Separate", "Ranking type : Single"),
      showStrips = T,
      legend = c(1,2),
      showAxisPlotLabels = T
    ) +
      theme(strip.background = element_rect(fill = "white", linetype = "dashed"),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_blank(),
            legend.position = "bottom")
    return(pm)
  }
} #Hyperparam optim US legal landscape

make_figure_S5 = function(plot_type = c("Sankey", "tissue_heatmap"),ds_ids = NULL){
  metaspace_meta = metaspace_metadata
  colnames(metaspace_meta)[1] = "ds_id"
  metaspace_meta = metaspace_meta[,c(1,3,4,5)]
  testing_datasets = dplyr::left_join(new_ds_testing, metaspace_meta)
  if (!is.null(ds_ids)){
    testing_datasets = testing_datasets[which(testing_datasets$ds_id %in% ds_ids),]
  }
  testing_datasets = testing_datasets %>% dplyr::left_join(testing_mz_range)
  testing_datasets$mz_range[is.na(testing_datasets$mz_range)] = "High"
  rNames = testing_datasets$ds_id
  testing_datasets = testing_datasets %>% dplyr::select(is_public, polarity, source, analyzer, rp_range, mz_range,
                                                        Sample_Information.Organism,
                                                        Sample_Information.Organism_Part) %>%
    as.data.frame()
  rownames(testing_datasets) = rNames
  colnames(testing_datasets) = c("Public", "Polarity", "Ionization Source", "Analyzer",
                                 "Resolving power range","mz range", "Organism", "Tissue")
  testing_datasets$`Ionization Source` = ifelse(testing_datasets$`Ionization Source` %in% c("MALDI","DESI"), testing_datasets$`Ionization Source`, "Other")
  testing_datasets$`Resolving power range` = str_to_title(testing_datasets$`Resolving power range`)
  for (i in 1:nrow(testing_datasets)){
    if (testing_datasets$Organism[i] %nin% c("Rattus norvegicus (rat)", "Mus musculus (mouse)","Homo sapiens (human)","Mouse ")){
      next()
    }
    testing_datasets$Organism[i] = switch(testing_datasets$Organism[i],"Rattus norvegicus (rat)" = "Rat",
                                          "Mus musculus (mouse)" = "Mouse",
                                          "Homo sapiens (human)" = "Human",
                                          "Mouse " = "Mouse")
  }
  testing_datasets$Organism = ifelse(testing_datasets$Organism %in% c("Mouse","Human","Rat"), testing_datasets$Organism, "Other")
  testing_datasets$Public = ifelse(testing_datasets$Public == "TRUE", "True", "False")

  testing_datasets_long = testing_datasets %>% ggsankey::make_long(Public, Polarity, `Ionization Source`,
                                                                   Analyzer, `Resolving power range`,`mz range`, Organism)

  cols = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))

  # for(col in colnames(testing_datasets_long)){
  #   testing_datasets_long[[col]] = as.factor(gsub(" ","\n",testing_datasets_long[[col]]))
  # }

  a = ggplot(testing_datasets_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_alluvial(flow.alpha = .6, space = 1) +
    geom_alluvial_text(size = 8, color = "Black", angle = 90) +
    scale_fill_discrete(type = cols) +
    #scale_fill_brewer(palette = "Paired", type = "qual") +
    theme_alluvial(base_size = 20) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5),
          axis.text.x = element_text(size = 30, color = "Black"),
          axis.title.y = element_text(size = 30, color = "Black"),
          axis.text.y = element_text(size = 30, color = "Black")) +
    scale_x_discrete(labels = addline_format(colnames(testing_datasets))) +
    ggtitle("Testing Datasets") +
    ylab("Number of Datasets")

  #testing_datasets = rename(testing_datasets, Tissue = Sample_Information.Organism_Part)
  testing_datasets$Tissue = testing_datasets$Tissue %>% str_to_title() %>% trimws()
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Brain")] = "Brain"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Whole")] = "Whole Organism"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Ovary With Endometriosis")] = "Ovary"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Artetry")] = "Artery"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "N/A")] = "None"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Stomach Mucus")] = "Stomach"
  testing_datasets$Tissue[str_which(testing_datasets$Tissue, "Palm Of Hand")] = "Skin"
  testing_datasets$Tissue[is.na(testing_datasets$Tissue)] = "None"



  wide_mat = testing_datasets %>% group_by(Organism, Tissue) %>% summarise(n = n()) %>%
    as.data.frame()
  # wide_mat$Tissue[is.na(wide_mat$Tissue)] = "None"
  wide_mat = wide_mat %>% pivot_wider(names_from = Organism, values_from = n)
  rNames = wide_mat$Tissue
  wide_mat = wide_mat[,-1] %>% as.matrix()
  rownames(wide_mat) = rNames
  b = pheatmap(wide_mat, cluster_rows = F, cluster_cols = F, na_col = "White",
               display_numbers = wide_mat, number_color = "white",
               color =  colorRampPalette(rev(brewer.pal(n = 11, name =
                                                          "RdBu")))(40),
               fontsize = 20, legend = F)

  if(plot_type == "Sankey"){
    return(a)
  }
  else{
    return(b)
  }
} #Sankey and tissue heatmap 15x20 landscape
make_figure_S6 = function(){
  plot_data = testing_param_optim_annot
  plot_data = plot_data %>%
    dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots,db) %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",") %>%
    group_by(ds_id,db, FDR_pct) %>%
    summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
              pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
    as.data.frame()
  plot_data$Delta = plot_data$pred_fdr_per_ds / plot_data$msm_fdr_per_ds
  plot_data$Delta[is.infinite(plot_data$Delta)] = max(plot_data$Delta[!is.infinite(plot_data$Delta)]) + 10
  plot_data = plot_data[order(plot_data$Delta, decreasing = T),]
  plot_data$FDR_pct = as.character(plot_data$FDR_pct)
  #plot_data = plot_data[which(plot_data$FDR_pct == "5" ),]
  plot_data$LFC = log2(plot_data$Delta)
  plot_data$Delta_class = ifelse(plot_data$Delta > 1, "Higher", "Lower")
  plot_data$Delta_class[which(plot_data$Delta == 1)] = "Equal"
  plot_data = plot_data[order(plot_data$LFC,decreasing = F),]

  new_data = plot_data %>% dplyr::select(ds_id, msm_fdr_per_ds, pred_fdr_per_ds, Delta,
                                         LFC, Delta_class, FDR_pct,db)
  new_data$diff = new_data$pred_fdr_per_ds - new_data$msm_fdr_per_ds
  #new_data$point_col = ifelse(new_data$diff > 0, "darkgreen", "darkred")
  #new_data$point_col[which(new_data$diff == 0)] = "blue"

  new_data = new_data[order(new_data$diff, decreasing = F),]


  p2 = ggplot(new_data[new_data$FDR_pct == "10",], aes(x=LFC, y=diff)) +
    geom_hex(binwidth = c(0.5,50)) +
    scale_y_continuous(breaks = seq(-50, 1000, 50),limits = c(-50, 700)) +
    scale_x_continuous(breaks = seq(-3, 7, 1)) +
    scale_fill_binned(type = "viridis", breaks = seq(0,40,5)) +
    theme_pubr(legend = "right") +
    theme(plot.title = element_text(size=18, hjust = 0.5),
          plot.subtitle = element_text(size=16, hjust = 0.5),
          axis.title.y = element_text(size=18, hjust = 0.5),
          axis.title.x = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=16),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(size = 18),
          strip.background = element_rect(colour = "white")) +
    ggtitle('Density of added annotations across testing datasets',
            subtitle = "FDR = 10%") +
    xlab("Log2(Fold change number of annotations)") +
    ylab("Difference in number of annotations") +
    facet_wrap(~ db, scales = "free")
  return(p2)
} #Hexplot Testing datasets 15x20 landscape
make_figure_S7 = function(perf_type = c("Duration", "memory"),
                          plot_type = c("line", "paired")){
  perf_data = performance_per_ds_pilot
  perf_data$memory_left = performance_per_ds_pilot$gb_secs_left / performance_per_ds_pilot$duration_left
  perf_data$memory_right = performance_per_ds_pilot$gb_secs_right / performance_per_ds_pilot$duration_right

  duration_data = perf_data %>% dplyr::select(ds_id, duration_left, duration_right)
  colnames(duration_data) = c("ds_id", "MSM", "METASPACE-ML")

  memory_data = perf_data %>% dplyr::select(ds_id, memory_left, memory_right)
  colnames(memory_data) = c("ds_id", "MSM", "METASPACE-ML")

  if(perf_type == "duration"){
    plot_data = duration_data
    y_lab = "Time/s"
  }
  else{
    plot_data = memory_data
    y_lab = "Memory/Gb"
  }
  if (plot_type == "paired"){
    p = plot_data %>% ggpaired(cond1 = "METASPACE-ML", cond2 = "MSM",
                               color = "condition", linetype = "solid",
                               line.size = 1, point.size = 3,
                               palette = c("#A6CEE3", "#FB9A99"),
                               ggtheme = theme_pubr(legend = "bottom"),
                               line.color = "lightgrey",
                               xlab = FALSE, ylab = y_lab, title = paste0(perf_type, "_Performance"),
                               short.panel.labs = T) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.title = element_blank(),
            strip.text = element_text(size = 20),
            strip.background = element_rect(fill = "white"))

    ggpar(p, legend.title = list(color = "Score_type"),
          x.text.angle = 0)
  }
  else{
    plot_data = plot_data %>% gather(-ds_id, key = "Approach", value = "Time")
    plot_data = plot_data[order(plot_data$Time, decreasing = T),]
    plot_data$ds_id = as.factor(plot_data$ds_id)
    dss_order = as.character(plot_data$ds_id)
    p =  plot_data %>%
      ggplot(aes(x=ds_id, y=Time, color = Approach, fill = Approach,
                 group = Approach)) +
      geom_line() +
      geom_point(shape=21, size=3) +
      theme_pubr() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.title = element_blank()) +
      xlab("Dataset") +
      ylab("Time/s") +
      scale_y_continuous(breaks = seq(0,30,2)) +
      scale_color_discrete(type = c("#A6CEE3", "#FB9A99")) +
      scale_fill_discrete(type = c("#A6CEE3", "#FB9A99"))

  }
  return(p)
} #Performance 10x8 portrait

#AP boxplot testing different databases
make_figure_S8 = function(eval_df,model_param_comb = NULL,
                          metric_type = "MAP"){
  if (!is.null(model_param_comb)){
    eval_df = eval_df %>% dplyr::filter(model_params == model_param_comb)
  }
  plot_data = eval_df
  plot_data = plot_data %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",") %>%
    dplyr::group_by(ds_id, model_params,score_type,db) %>%
    dplyr::summarise(metric_value = mean(metric_value)) %>%
    as.data.frame() %>%
    dplyr::select(db, score_type, metric_value)
  plot_data$db = sub("-.*", "", plot_data$db)
  p = ggplot(plot_data, aes(x=db, y=metric_value, color=score_type)) +
    geom_boxplot(outlier.shape = NA, size = 1.25) +
    scale_color_manual(values = c("#A6CEE3", "#FB9A99"))+
    geom_point(position=position_jitterdodge(dodge.width = 1)) +
    theme_pubr(legend = "bottom") +
    theme(axis.title.y = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=16),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
    scale_y_continuous(breaks = seq(0,1,0.1)) +
    xlab("") +
    ylab("MAP")
  return(p)
}

#Add probs = c(0.025, 0.25, 0.5, 0.75, 0.975) to quantile_list in #121 densityheatmap
