# Functions ---------------------------------------------------------------
plot_eval_metrics_boxplot = function(df, metric_type = c("MAP", "nDCG")){
  df = df %>% dplyr::select(metric_type, metric_value, score_type)
  df = df[which(df$metric_type == metric_type),]
  bp = ggplot(df, aes(x=score_type, y=metric_value,
                      color = score_type)) +
    geom_boxplot(width = 0.1, size = 1.5, position = "dodge",
                 notch = F, outlier.color = "Black") +
    scale_color_discrete(type = c("#A6CEE3", "#FB9A99")) +
    geom_jitter(size=4, alpha=0.6,width = 0.2) +
    theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size=18, hjust = 0.5),
          plot.subtitle = element_text(size=16, hjust = 0.5),
          axis.title.y = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=16),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
    ggtitle('Mean Average Precision',subtitle = "Cross validation")+
    xlab("") +
    ylab("MAP")
  return(bp)
}
plot_eval_metrics_testing = function(df, add_err_bar = T){
  df = df %>% dplyr::select(db, metric_value, score_type, metric_val_sd)
  bp = ggplot(df, aes(x=db, y=metric_value, fill=score_type)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 0.2, padding = 0),
             alpha = 1, width = 0.5) +
    scale_fill_discrete(type = c("#A6CEE3", "#FB9A99")) +
    theme_pubr(legend = "right") +
    theme(plot.title = element_text(size=18, hjust = 0.5),
          axis.title.y = element_text(size=20, hjust = 0.5),
          legend.text = element_text(size=20),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    ggtitle('Database comparison')+
    xlab("") +
    ylab("MAP") +
    rotate_x_text(angle = 45)

  if (add_err_bar){
    bp = bp +
      geom_errorbar(aes(x=db, ymin=metric_value-metric_val_sd, ymax=metric_value+metric_val_sd),
                    width=0.4, colour="black", alpha=0.9, size=1.3, position = position_dodge(0.9))
  }

  return(bp)
}

Get_intensity_dist_data = function(annot_w_ions_df, ds_id,
                                   intens_res_path,
                                   fdr_pct = 0.1,intens_type = c("median_nz", "q95_nz",
                                                                 "q99_nz", "total_intensity", "mean_nz")){
  annot_data = annot_w_ions_df %>%
    dplyr::select(ds_id, formula, modifier, msm_fdr, pred_score_fdr)
  ds_path = file.path(intens_res_path, paste0(ds_id, ".csv"))
  if (!file.exists(ds_path)){
    return(NULL)
  }
  intens_df = read.csv(ds_path)
  if (nrow(intens_df) == 0){
    return(NULL)
  }
  intens_df = intens_df[,c("formula", "adduct", intens_type)]
  colnames(intens_df)[3] = "total_intensity"
  intens_df = intens_df[!duplicated(intens_df[,c(1,2)]),]
  colnames(intens_df)[which(colnames(intens_df) == "adduct")] = "modifier"

  ds_fdr_ions = annot_data[which(annot_data$ds_id == ds_id),]
  common_ions = plyr::match_df(ds_fdr_ions, intens_df,
                               on = c("formula","modifier"))
  common_ions = dplyr::left_join(common_ions, intens_df, by = c("formula","modifier"))
  common_ions$ion = paste0(common_ions$formula, common_ions$modifier)
  common_ions = common_ions[!duplicated(common_ions),]

  common_ions$total_intensity = log10(common_ions$total_intensity)

  msm_ions = common_ions$ion[which(common_ions$msm_fdr <= fdr_pct)] %>% unique()
  ML_ions = common_ions$ion[which(common_ions$pred_score_fdr <= fdr_pct)] %>% unique()
  ML_only = ML_ions[which(ML_ions %nin% msm_ions)]
  msm_only = msm_ions[which(msm_ions %nin% ML_ions)]

  hist_data = list()
  for (i in c("MSM", "ML", "ML_only", "MSM_only")){
    sel_ions = switch(i, "MSM" = msm_ions, "ML" = ML_ions,
                      "ML_only" = ML_only, "MSM_only" = msm_only)
    if (length(sel_ions) == 0){
      hist_data[[i]] = NULL
    }
    else{
      df = common_ions[which(common_ions$ion %in% sel_ions),] %>%
        dplyr::select(ion, total_intensity)
      df$Group = paste0(i, " (", length(sel_ions), ")")
      hist_data[[i]] = df
    }
  }
  hist_data = dplyr::bind_rows(hist_data)
  return(hist_data)

}

prepare_data_umap = function(raw_res, ds_id = NULL, sf = NULL, fix_mz_err = T){
  if (!is.null(ds_id)){
    filtered = raw_res[which(raw_res$ds_id == ds_id),]
  }
  else{
    filtered = raw_res
  }
  if (!is.null(sf)){
    filtered = filtered[which(filtered$formula %in% sf),]
  }

  if (fix_mz_err){
     filtered$mz_err_abs = 1 - abs(filtered$mz_err_abs)
     filtered$mz_err_rel = 1 - abs(filtered$mz_err_rel)
  }

  filtered$ion_id = paste0(filtered$group_name, "_", filtered$formula, "_", filtered$modifier)
  labels = filtered %>% dplyr::select(ion_id, target)
  labels = labels[!duplicated(labels),]
  labels = ifelse(labels$target == 1, "Target", "Decoy")

  if (any(str_detect(colnames(filtered), "abserr"))){
    feature_mat = filtered %>% dplyr::select(ion_id, chaos, spatial, spectral,
                                             mz_err_abs_abserr,mz_err_rel_abserr)
  }
  else{
    feature_mat = filtered %>% dplyr::select(ion_id, chaos, spatial, spectral,
                                             mz_err_abs,mz_err_rel)
    if (fix_mz_err){
      feature_mat$mz_err_abs = 1 - abs(feature_mat$mz_err_abs)
      feature_mat$mz_err_rel = 1 - abs(feature_mat$mz_err_rel)
    }
  }
  feature_mat = feature_mat[!duplicated(feature_mat),]
  feature_mat = feature_mat %>% pivot_longer(!ion_id, names_to = "feature", values_to = "feat_value")
  feature_mat$feature[str_which(feature_mat$feature, "mz_err_abs")] = "mz_err_abs"
  feature_mat$feature[str_which(feature_mat$feature, "mz_err_rel")] = "mz_err_rel"

  feature_mat = feature_mat %>% pivot_wider(names_from = ion_id, values_from = feat_value)
  rNames = feature_mat$feature
  feature_mat = feature_mat[,-1] %>% as.data.frame()
  rownames(feature_mat) = rNames

  scores_mat = filtered %>% dplyr::select(ion_id, msm, pred_score)
  scores_mat = scores_mat[!duplicated(scores_mat),]
  scores_mat = scores_mat %>% pivot_longer(!ion_id, names_to = "score_type", values_to = "score_value") %>%
    as.data.frame()
  scores_mat = scores_mat[!duplicated(scores_mat),]
  scores_mat = scores_mat %>% pivot_wider(names_from = ion_id, values_from = score_value)
  rNames = scores_mat$score_type
  rNames = str_replace_all(rNames, c("msm" = "MSM", "pred_score" = "CatBoost"))
  scores_mat = scores_mat[,-1] %>% as.data.frame()
  rownames(scores_mat) = rNames

  return(list("feature_mat" = feature_mat, "labels" = labels, "score_mat" = scores_mat))

}
create_combs = function(df, comb_cols){
  if (all(comb_cols %nin% colnames(df))){
    stop("Error : some combination columns are missing from df")
  }
  else{
    combs = list()
    for (j in 1:length(comb_cols)){
      combs[[comb_cols[j]]] = unique(df[,comb_cols[j]])
    }
    combs = expand.grid(combs)
    return(combs)
  }
}
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

clean_up_eval_res = function(df, optim_type = c("old", "new"), dss = NULL){
  if(optim_type == "old"){
    if (length(intersect(c("includes_bad_rank", "ranking_type", "features_type", "catmonotonic"),
                         colnames(df))) == 0){
      if ("score_type" %in% colnames(df)){
        df$score_type[which(df$score_type == "Catboost")] =
          "METASPACE-ML"
      }
    }
    else{
      df$model_params = paste0("Incl_Bad_", df$includes_bad_rank,
                               "_rank_", df$ranking_type,
                               "_feat_type_", df$features_type,
                               "_catmonotonic_", df$catmonotonic)
      df = df %>% dplyr::select(-includes_bad_rank, -ranking_type,-features_type,
                                -catmonotonic)
      if ("score_type" %in% colnames(df)){
        df$score_type[which(df$score_type == "Catboost")] =
          "METASPACE-ML"
      }
    }

  }
  else{
    if ("model_params" %in% colnames(df)){
      if ("score_type" %in% colnames(df)){
        df$score_type[which(df$score_type == "Catboost")] =
          "METASPACE-ML"
      }
    }
    else{
      df$model_params = paste0("mzmedian_", sub(".*_","",df$mz_train_thresh),
                               "_dsz_", df$dsz_type,
                               "_chaos_", df$include_chaos,
                               "_OD_", df$include_OD,
                               "_objW_", df$is_objW)
      df = df %>% dplyr::select(-mz_train_thresh, -dsz_type,-include_chaos,
                                -include_OD, -is_objW)
      if ("score_type" %in% colnames(df)){
        df$score_type[which(df$score_type == "Catboost")] =
          "METASPACE-ML"
      }
    }
  }
  if (any(str_detect(colnames(df), "group"))){
    colnames(df)[str_which(colnames(df), "group")] = "group_name"
  }
  if ("ds_id" %in% colnames(df)){
    if (any(str_detect(df$ds_id, "ml_training")))
      df$ds_id = sub("_ml_training.*", "", df$ds_id)
  }
  if (!is.null(dss) & any(str_detect(colnames(df), "group"))){
    mapping_df = data.frame(grp = df$group_name,
                            ds_id = sub(",.*", "", df$group_name))
    if (any(str_detect(mapping_df$ds_id, "ml_training"))){
      mapping_df$ds_id = sub("_ml_training.*", "", mapping_df$ds_id)
    }
    df = df %>% dplyr::filter(group_name %in%
                                mapping_df$grp[which(mapping_df$ds_id %in% dss)])
  }
  return(df)
}
fisher_annot_ML = function(annot_w_ions_df,dataset_id, bg_list,
                           core_metabo_only = T, fdr_pct = 0.1,
                           pathway_rel_only = T){
  ds_df = annot_w_ions_df %>% filter(ds_id == dataset_id)
  if ("db" %in% colnames(ds_df)){
    ds_df = ds_df[which(ds_df$db == "CoreMetabolome-v3"),]
  }
  msm_ions = ds_df$formula[which(ds_df$msm_fdr <= fdr_pct)] %>% unique()
  ML_ions = ds_df$formula[which(ds_df$pred_score_fdr <= fdr_pct)] %>% unique()
  all_detectable_ions = ds_df$formula %>% unique()
  if (pathway_rel_only){
    all_detectable_ions = all_detectable_ions[which(all_detectable_ions %in% path_rel_sf)]
  }

  ML_ions_only = ML_ions[which(ML_ions %nin% msm_ions)]
  if(pathway_rel_only){
    ML_ions_only = ML_ions_only[which(ML_ions_only %in% path_rel_sf)]
  }
  #not_ML_ions = all_ions[all_ions %nin% ML_ions] %>% unique()

  final_res = list()
  for (i in 1:length(bg_list)){
    term_name = names(bg_list)[i]
    term_sf = bg_list[[i]]
    term_sf = term_sf[which(term_sf %in% ds_df$formula)]
    if (core_metabo_only){
      term_sf = term_sf[which(term_sf %in% core_metab$formula)]
    }
    if(pathway_rel_only){
      term_sf = term_sf[which(term_sf %in% path_rel_sf)]
    }
    if (length(term_sf) == 0){
      final_res[[i]] = NULL
    }
    else{
      TP = length(intersect(ML_ions_only, term_sf))
      TP_sf = intersect(ML_ions_only, term_sf)
      FP = length(which(ML_ions_only %nin% term_sf))
      FN = length(which(term_sf %nin% ML_ions_only))

      TN = all_detectable_ions[which(all_detectable_ions %nin% term_sf)]
      TN = length(TN[which(TN %nin% ML_ions_only)])

      fisher_mat = matrix(c(TP, FP, FN, TN), nrow = 2, ncol = 2, byrow = T)
      fisher_res = hypergea::hypergeom.test(fisher_mat)

      term_res = data.frame(ds_id = dataset_id, term = term_name,
                            TP = TP, FP = FP, FN = FN , TN = TN,
                            OR = fisher_res$estimate %>% as.numeric(),
                            pval = fisher_res$p.value,
                            TP_markers = paste(TP_sf, collapse = ","))
      term_res$FE = (term_res$TP / (term_res$TP + term_res$FP)) /
        ((term_res$TP + term_res$FN) / (term_res$TP + term_res$FP + term_res$FN + term_res$TN))
      final_res[[i]] = term_res
    }
  }
  final_res = dplyr::bind_rows(final_res)
  return(final_res)
}

Run_enrichment_pipeline = function(annot_w_ions_df,
                                   core_metabo_only = T,
                                   pathway_rel_only = T){
  all_ds = unique(annot_w_ions_df$ds_id)
  message("\n Running Enrichment ... \n")
  pb = txtProgressBar(min = 0, max = length(all_ds), style = 3)
  all_ds_enrich_res_subclass = list()
  for (i in 1:length(all_ds)){
    all_ds_enrich_res_subclass[[i]] = fisher_annot_ML(annot_w_ions_df = annot_w_ions_df,
                                                      dataset_id = all_ds[i],
                                                      bg_list = bg, core_metabo_only = core_metabo_only,
                                                      pathway_rel_only = pathway_rel_only)
    setTxtProgressBar(pb, i)
  }
  all_ds_enrich_res_subclass = dplyr::bind_rows(all_ds_enrich_res_subclass)
  return(all_ds_enrich_res_subclass)
}

select_eval_results = function(model_type = c("new", "old"),
                               ds_type = c("testing", "CV"),
                               abserr = T,
                               dsz_type = c("Fixed", "nonzero"),
                               mzmedian = c("50k", "10k")){
  eval_res = switch(model_type, "new" = all_new, "old" = all_old)

  eval_res = eval_res[[ds_type]]
  if (!dsz_type %in% names(eval_res)){
    stop("Parameter combination not found")
  }
  eval_res = eval_res[[dsz_type]]
  if (model_type == "old"){
    for (i in 1:length(eval_res)){
      eval_res[[i]] = clean_up_eval_res(df = eval_res[[i]], optim_type = "old",
                                        dss = NULL)
    }
    return(eval_res)
  }
  else{
    abserr_name = ifelse(abserr, "abserr", "no_abserr")
    eval_df = eval_res[["eval"]][[abserr_name]]
    annot_df = eval_res[["annot"]][[abserr_name]]
    annot_w_ions_df = eval_res[["annot_w_ions"]][[abserr_name]]

    final = list("eval" = eval_df[[str_which(names(eval_df), mzmedian)]],
                 "annot" = annot_df[[str_which(names(annot_df), mzmedian)]],
                 "annot_w_ions" = annot_w_ions_df[[str_which(names(annot_w_ions_df), mzmedian)]])

    sel_dss = mz_median_dss[[ds_type]]
    if (mzmedian == "10k"){
      sel_dss = sel_dss$ds_id[which(sel_dss$total_mz_peaks_median < 10000)]
    }
    else{
      sel_dss = sel_dss$ds_id[which(sel_dss$total_mz_peaks_median < 50000)]
    }

    for (i in 1:length(final)){
      if (!"model_params" %in% colnames(final)){
        final[[i]]$model_params = paste0("mzmedian_", mzmedian,
                                         "_dsz_", dsz_type,
                                         "_chaos_", "Yes",
                                         "_OD_", "No",
                                         "_objW_", "No")
      }
      final[[i]] = clean_up_eval_res(df = final[[i]], optim_type = "new",
                                     dss = sel_dss)
    }
    return(final)
  }

}
plot_intensity_bubble_plot = function(intens_df, x_colname, y_colname,
                                      wilcox_alternative = "less",
                                      return_pval = F){
  intens_data = intens_df %>% dplyr::select(-n_ions) %>%
    spread(key = "Group", value = "median_intens")
  intens_data = intens_data[,c("ds_id", x_colname, y_colname)]

  n_ions_data = intens_df %>%
    mutate(Group = paste0(Group, "_nions")) %>%
    dplyr::select(-median_intens) %>%
    spread(key = "Group", value = "n_ions")

  x_ion_colname = paste0(x_colname, "_nions")
  y_ion_colname = paste0(y_colname, "_nions")

  n_ions_data = n_ions_data[,c("ds_id", x_ion_colname, y_ion_colname)]


  plot_data = intens_data %>% left_join(n_ions_data, by  = "ds_id")
  plot_data$abline_pos = ifelse(plot_data[,y_colname] > plot_data[,x_colname],
                                "Higher", "Lower")
  plot_data$abline_pos[which(plot_data[,y_colname] == plot_data[,x_colname])] = "Equal"

  t_test_stats = t.test(plot_data[,x_colname], plot_data[,y_colname],
                        alternative = wilcox_alternative)
  wilcox_stats = wilcox.test(plot_data[,x_colname], plot_data[,y_colname],
                        alternative = wilcox_alternative)
  if (return_pval){
    return(list("t_test" = t_test_stats,
                "wilcox" = wilcox_stats))
  }

  abline_data = table(plot_data$abline_pos) %>% as.data.frame()
  colnames(abline_data) = c("abline_pos", "n_ds")
  abline_data$n_ds = paste0("N = ", abline_data$n_ds)
  abline_data = abline_data[which(abline_data$abline_pos != "Equal"),]
  abline_data$x = c(1,10)
  abline_data$y = c(10,1)

  p = ggplot() +
    geom_point(data = plot_data, aes(x = .data[[x_colname]],
                                     y = .data[[y_colname]],
                                     size = .data[[x_ion_colname]],
                                     colour = x_colname),
               alpha = 0.3) +
    geom_point(data = plot_data, aes(x = .data[[x_colname]],
                                     y = .data[[y_colname]],
                                     size = .data[[y_ion_colname]],
                                     colour = y_colname),
               alpha = 0.3) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(n.breaks = 10,limits = c(0,12)) +
    scale_y_continuous(n.breaks = 10, limits = c(0,12)) +
    scale_size_binned_area(breaks = c(1,10,50,100,300,500),
                           name = "Number of ions", max_size = 10) +
    scale_colour_manual("",
                        breaks = c(x_colname, y_colname),
                        values = c("#1F78B4","#FB9A99")) +
    xlab(paste0(x_colname, "\nMedian(Log10(Total Intensity))")) +
    ylab(paste0(y_colname, "\nMedian(Log10(Total Intensity))")) +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.position = "left") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    geom_text(data = abline_data, aes(x = x, y = y, label = n_ds), size = 5)

  p2 = ggExtra::ggMarginal(p,type = "histogram",
                           xparams = list("fill" = "#1F78B4"),
                           yparams = list("fill" = "#FB9A99"))
  return(p2)
}

plot_intensity_boxplot = function(intens_df){

  habal = intens_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(n=n_distinct(ds_id))

  intens_df = intens_df %>% dplyr::left_join(habal, by = "Group")
  intens_df$x_text = paste0(intens_df$Group, "\n", "(", intens_df$n, ")")

  unique_grps = unique(intens_df$x_text)

  sig_comparison_list = list(c(str_subset(unique_grps, "METASPACE-ML\n")
                               , str_subset(unique_grps, "METASPACE-ML_only\n")),
                             c(str_subset(unique_grps, "MSM\n")
                               , str_subset(unique_grps, "METASPACE-ML_only\n")),
                               c(str_subset(unique_grps, "MSM\n")
                                 , str_subset(unique_grps, "METASPACE-ML\n")),
                                 c(str_subset(unique_grps, "MSM\n")
                                   , str_subset(unique_grps, "MSM_only\n")))

  ggplot(data = intens_df, mapping = aes(x = x_text, y = median_intens,
                                                            color=log10(n_ions))) +
    geom_boxplot(outlier.shape = NA, size = 1) +
    geom_jitter(size = 3,alpha=0.8,
                position = position_jitterdodge(dodge.width = 0.9)) +
    scale_y_continuous(n.breaks = 10, limits = c(0,10)) +
    scale_color_continuous(type = "viridis", name = "Log10 (Number of ions)") +
    # scale_size_binned_area(breaks = c(1,10,50,100,300,500),
    #                        name = "Number of ions", max_size = 10) +
    theme_pubr() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))+
    ylab("Median(Log10(Intensity))") +
    xlab("") +
    geom_signif(comparisons = sig_comparison_list,
                map_signif_level = F, na.rm = T,textsize = 8,
                y_position = c(10,8, 9, 10))

}

make_density_hm = function(all_annot, all_eval, FDR = 10,
                           vals = c("Diff", "LogDiff", "LFC")){
  all_annot$diff_sign = ifelse(all_annot$Difference < 0,-1,1)

  all_annot = all_annot %>%
    dplyr::mutate(LFC = log2((pred_fdr_per_ds + 1) / (msm_fdr_per_ds + 1)),
                  pct_diff = (Difference + 1) / (msm_fdr_per_ds + 1),
                  log_diff = diff_sign * log10(abs(Difference) + 1))

  #all_annot$log_diff[is.infinite(all_annot$log_diff)] = 0


  val_col = switch(vals, "Diff" = "Difference", "LogDiff" = "log_diff",
                   "LFC" = "LFC")

  common_ds = all_annot %>%
    dplyr::select(ds_id, model_params) %>%
    dplyr::group_by(ds_id) %>%
    dplyr::summarise(n = n_distinct(model_params))
  common_ds = common_ds$ds_id[which(common_ds$n == length(unique(all_annot$model_params)))]

  d_hm_mat = all_annot %>%
    dplyr::filter(FDR_pct == FDR, ds_id %in% common_ds)
  d_hm_mat = d_hm_mat[,which(colnames(d_hm_mat) %in% c("ds_id", "model_params", val_col))]
  colnames(d_hm_mat)[which(colnames(d_hm_mat) == val_col)] = "values"
  d_hm_mat = d_hm_mat %>%
    pivot_wider(names_from = model_params, values_from = values)
  rNames = d_hm_mat$ds_id
  d_hm_mat = d_hm_mat[,-1] %>% as.matrix()

  x = colnames(d_hm_mat)
  anno_df = data.frame(mzmedian = sub("_.*","",sub(".*mzmedian_", "",x)),
                       dsz = sub("_.*","",sub(".*dsz_", "",x)),
                       chaos = sub("_.*","",sub(".*chaos_", "",x)),
                       OD = sub("_.*","",sub(".*OD_", "",x)),
                       objW = sub("_.*","",sub(".*objW_", "",x)))
  rownames(anno_df) = x

  y = all_eval
  y$ds_id = sub(",.*","",y$group_name)
  habal = y %>% dplyr::filter(ds_id %in% common_ds) %>%
    group_by(model_params, score_type) %>%
    summarise(MAP = mean(metric_value), AP_stdev = sd(metric_value)) %>%
    dplyr::filter(score_type %in% c("Catboost", "METASPACE-ML")) %>% as.data.frame()
  rownames(habal) = habal$model_params
  habal = habal[rownames(anno_df),] %>% dplyr::select(MAP)


  p = densityHeatmap(data = d_hm_mat, show_column_names = F, cluster_columns = T,
                     show_quantiles = T, title = paste0(FDR, " % FDR"),
                     ylab = vals,ylab_gp = gpar(fontsize = 16),
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 16),
                                                 labels_gp = gpar(fontsize = 14))) %v%
    HeatmapAnnotation(df = habal,
                      col = list(MAP = circlize::colorRamp2(
                        breaks = c(min(habal$MAP), max(habal$MAP)),
                        colors = c("white", "Blue"))),
                      annotation_name_gp = gpar(fontsize = 16),
                      annotation_legend_param = list(title_gp = gpar(fontsize = 16),
                                                     labels_gp = gpar(fontsize = 14))) %v%
    HeatmapAnnotation(df = anno_df,
                      col = list(mzmedian = c("10k" = "#66C2A5", "50k" = "#FC8D62",
                                              "all" = "#E5C494"),
                                 dsz = c("Fixed" = "#A6CEE3", "nonzero" = "#1F78B4"),
                                 chaos = c("1" = "#B2DF8A", "0" = "#33A02C"),
                                 OD = c("1" = "#CAB2D6", "0" = "#6A3D9A"),
                                 objW = c("1" = "#FB9A99", "0" = "#E31A1C")),
                      annotation_name_gp = gpar(fontsize = 16),
                      annotation_legend_param = list(title_gp = gpar(fontsize = 16),
                                          labels_gp = gpar(fontsize = 14)))

  return(p)


}

Fbeta_score <- function(tp, fp, fn, tn) {
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)

  beta <- 0.5  # Beta value for F0.5 score

  f_score <- ((1 + beta^2) * precision * recall) / ((beta^2 * precision) + recall)

  f <- matrix(f_score, ncol = 1)
  colnames(f) <- "F1_score"
  return(f)
}



calc_reliability_score = function(annot_ions_df){
  dataset_id = unique(annot_ions_df$ds_id)

  pb = txtProgressBar(min = 0, max = length(dataset_id), style = 3)

  counter = 0

  final = list()
  for (ds in dataset_id){
    counter = counter + 1
    habal = annot_ions_df[annot_ions_df$ds_id == ds,]
    y = cutpointr::cutpointr(data = habal,
                             x = pred_score_fdr,
                             class = target, metric = F1_score,
                             pos_class = 1, neg_class = 0) %>% suppressMessages()
    z = y$roc_curve[[1]]
    z$fbeta = Fbeta_score(z$tp, z$fp, z$fn, z$tn)

    res_fdr = list()
    for (f in c(0.05, 0.1, 0.2, 0.5)){
      closest = which.min(abs(z$x.sorted - f))

      fbeta_thresh = z$fbeta[closest]
      max_fbeta = max(z$fbeta, na.rm = T)
      optim_fdr = z$x.sorted[which.max(z$fbeta)]

      rel_score = (fbeta_thresh/max_fbeta) * (1 - optim_fdr)
      res = data.frame(ds_id = ds, rel_score = rel_score, FDR_pct = f*100,
                       fbeta_thresh = fbeta_thresh, max_fbeta = max_fbeta, optim_fdr = optim_fdr)
      res_fdr[[as.character(f)]] = res
    }
    res_fdr = res_fdr %>% dplyr::bind_rows()

    final[[ds]] = res_fdr
    setTxtProgressBar(pb, counter)
  }
  final = final %>% dplyr::bind_rows()
  return(final)
}

# Nat comms Data manipulation functions ----------------------------------------------
prepare_plot_df = function(eval_df, annot_df,
                           data_type = c("Training", "Testing"),
                           context_size = 50,
                           kingdom = c("Animal", "Plant"),
                           mean_eval_context = F,
                           metrics_per_grp = F,
                           FDR.pct = NULL){
  if (data_type == "Training"){
    ds_context_list = training_dss[[match.arg(kingdom)]][[as.character(context_size)]]
  } else if (data_type == "Testing"){
    ds_context_list = testing_dss[[match.arg(kingdom)]]
  }
  context_df = ds_context_list %>%
    as.data.frame(check.names = F) %>%
    gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F)

  eval_df$score_type = str_replace_all(eval_df$score_type,
                                       c("Catboost" = "METASPACE-ML"))
  eval_df = eval_df %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",",
                    remove = F) %>%
    dplyr::left_join(context_df, by = "ds_id")

  if(!metrics_per_grp){
    eval_df = eval_df %>%
      dplyr::group_by(context, score_type, metric_type, ds_id) %>%
      dplyr::summarise(metric_sd = sd(metric_value,
                                       na.rm = T),
                       metric_value = mean(metric_value,na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(context_df, by = c("ds_id", "context"))
  }

  if(mean_eval_context){
    eval_df = eval_df %>%
      dplyr::group_by(context, score_type, metric_type) %>%
      dplyr::summarise(metric_sd = sd(metric_value,na.rm = T),
                       metric_value = mean(metric_value, na.rm = T)) %>%
      dplyr::ungroup() %>%
      tidyr::separate(col = context,
                      into = c("polarity","mz_class",
                               "analyzer", "source",
                               "Species", "Sample_type"),
                      sep = ",", remove = F)
  }

  if("cv_split" %in% colnames(annot_df)){
    colnames(annot_df)[colnames(annot_df) == "cv_split"] = "db"
  }

  if (!metrics_per_grp){
    annot_df = annot_df %>%
      dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots, db) %>%
      tidyr::separate(group_name,
                      into = c("ds_id", "adduct"), sep = ",",
                      remove = F) %>%
      dplyr::group_by(ds_id, FDR_pct,db) %>%
      dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                       pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }
  else{
    colnames(annot_df)[colnames(annot_df) == "msm_fdr_annots"] = "msm_fdr_per_ds"
    colnames(annot_df)[colnames(annot_df) == "pred_fdr_annots"] = "pred_fdr_per_ds"
  }
  annot_df = annot_df %>%
    dplyr::mutate(Difference = pred_fdr_per_ds - msm_fdr_per_ds,
                  diff_sign = ifelse(Difference < 0,-1,1),
                  LogDiff = diff_sign * log10(abs(Difference) + 1),
                  FC = (pred_fdr_per_ds + 1) / (msm_fdr_per_ds + 1),
                  LFC = log2(FC))
  if(!is.null(FDR.pct)){
    annot_df = annot_df %>%
      dplyr::filter(FDR_pct == FDR.pct)
  }
  if(metrics_per_grp){
    all_res = left_join(eval_df, annot_df, by = c("group_name")) %>% suppressWarnings()
  }
  else{
    if(mean_eval_context){
      eval_df$source_analyzer = paste0(eval_df$source, "_", eval_df$analyzer)
      return(eval_df)
    }
    all_res = left_join(eval_df, annot_df, by = c("ds_id")) %>% suppressWarnings()
  }
  all_res$source_analyzer = paste0(all_res$source, "_", all_res$analyzer)
  all_res$FDR_pct = as.character(all_res$FDR_pct)
  return(all_res)
}

prepare_bulk_plot_df = function(eval_df, annot_df,
                           mean_eval_context = F,
                           metrics_per_grp = F,
                           FDR.pct = NULL){

  eval_df$score_type = str_replace_all(eval_df$score_type,
                                       c("Catboost" = "METASPACE-ML"))
  eval_df = eval_df %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",",
                    remove = F)

  if(!metrics_per_grp){
    eval_df = eval_df %>%
      dplyr::group_by(score_type, metric_type, ds_id) %>%
      dplyr::summarise(metric_sd = sd(metric_value,
                                      na.rm = T),
                       metric_value = mean(metric_value,na.rm = T)) %>%
      dplyr::ungroup()
  }

  if(mean_eval_context){
    eval_df = eval_df %>%
      dplyr::group_by(score_type, metric_type) %>%
      dplyr::summarise(metric_sd = sd(metric_value,na.rm = T),
                       metric_value = mean(metric_value, na.rm = T)) %>%
      dplyr::ungroup()
  }

  if("cv_split" %in% colnames(annot_df)){
    colnames(annot_df)[colnames(annot_df) == "cv_split"] = "db"
  }

  if (!metrics_per_grp){
    annot_df = annot_df %>%
      dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots, db) %>%
      tidyr::separate(group_name,
                      into = c("ds_id", "adduct"), sep = ",",
                      remove = F) %>%
      dplyr::group_by(ds_id, FDR_pct,db) %>%
      dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                       pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }
  else{
    colnames(annot_df)[colnames(annot_df) == "msm_fdr_annots"] = "msm_fdr_per_ds"
    colnames(annot_df)[colnames(annot_df) == "pred_fdr_annots"] = "pred_fdr_per_ds"
  }
  annot_df = annot_df %>%
    dplyr::mutate(Difference = pred_fdr_per_ds - msm_fdr_per_ds,
                  diff_sign = ifelse(Difference < 0,-1,1),
                  LogDiff = diff_sign * log10(abs(Difference) + 1),
                  FC = (pred_fdr_per_ds + 1) / (msm_fdr_per_ds + 1),
                  LFC = log2(FC))
  if(!is.null(FDR.pct)){
    annot_df = annot_df %>%
      dplyr::filter(FDR_pct == FDR.pct)
  }
  if(metrics_per_grp){
    all_res = left_join(eval_df, annot_df, by = c("group_name")) %>% suppressWarnings()
  }
  else{
    if(mean_eval_context){
      eval_df$source_analyzer = paste0(eval_df$source, "_", eval_df$analyzer)
      return(eval_df)
    }
    all_res = left_join(eval_df, annot_df, by = c("ds_id")) %>% suppressWarnings()
  }
  all_res$FDR_pct = as.character(all_res$FDR_pct)
  return(all_res)
}

prepare_data_umap <- function(raw_res, ds_id = NULL, sf = NULL, fix_mz_err = T) {
  if (!is.null(ds_id)) {
    filtered <- raw_res[which(raw_res$ds_id == ds_id), ]
  } else {
    filtered <- raw_res
  }
  if (!is.null(sf)) {
    filtered <- filtered[which(filtered$formula %in% sf), ]
  }

  if (fix_mz_err) {
    filtered$mz_err_abs <- 1 - abs(filtered$mz_err_abs)
    filtered$mz_err_rel <- 1 - abs(filtered$mz_err_rel)
  }

  filtered$ion_id <- paste0(filtered$group_name, "_", filtered$formula, "_", filtered$modifier)
  labels <- filtered %>% dplyr::select(ion_id, target)
  labels <- labels[!duplicated(labels), ]
  labels <- ifelse(labels$target == 1, "Target", "Decoy")

  if (any(str_detect(colnames(filtered), "abserr"))) {
    feature_mat <- filtered %>% dplyr::select(
      ion_id, chaos, spatial, spectral,
      mz_err_abs_abserr, mz_err_rel_abserr
    )
  } else {
    feature_mat <- filtered %>% dplyr::select(
      ion_id, chaos, spatial, spectral,
      mz_err_abs, mz_err_rel
    )
    if (fix_mz_err) {
      feature_mat$mz_err_abs <- 1 - abs(feature_mat$mz_err_abs)
      feature_mat$mz_err_rel <- 1 - abs(feature_mat$mz_err_rel)
    }
  }
  feature_mat <- feature_mat[!duplicated(feature_mat), ]
  feature_mat <- feature_mat %>% pivot_longer(!ion_id, names_to = "feature", values_to = "feat_value")
  feature_mat$feature[str_which(feature_mat$feature, "mz_err_abs")] <- "mz_err_abs"
  feature_mat$feature[str_which(feature_mat$feature, "mz_err_rel")] <- "mz_err_rel"

  feature_mat <- feature_mat %>% pivot_wider(names_from = ion_id, values_from = feat_value)
  rNames <- feature_mat$feature
  feature_mat <- feature_mat[, -1] %>% as.data.frame()
  rownames(feature_mat) <- rNames

  scores_mat <- filtered %>% dplyr::select(ion_id, msm, pred_score)
  scores_mat <- scores_mat[!duplicated(scores_mat), ]
  scores_mat <- scores_mat %>%
    pivot_longer(!ion_id, names_to = "score_type", values_to = "score_value") %>%
    as.data.frame()
  scores_mat <- scores_mat[!duplicated(scores_mat), ]
  scores_mat <- scores_mat %>% pivot_wider(names_from = ion_id, values_from = score_value)
  rNames <- scores_mat$score_type
  rNames <- str_replace_all(rNames, c("msm" = "MSM", "pred_score" = "CatBoost"))
  scores_mat <- scores_mat[, -1] %>% as.data.frame()
  rownames(scores_mat) <- rNames

  return(list("feature_mat" = feature_mat, "labels" = labels, "score_mat" = scores_mat))
}

parse_cv_train_error = function(out_file_path){
  out_data = read.delim(out_file_path, header = F,
                        na.strings = c("", "NA"))
  iter_time = out_data[,c(1,5)]
  colnames(iter_time) = c("iter", "time")

  iter_error = out_data[,c(1:3)]
  colnames(iter_error) = c("iter", "Training", "Testing")
  iter_error = iter_error[!is.na(iter_error$Testing),]

  iter_error$iter = sub(":.*", "", iter_error$iter) %>% as.numeric()
  iter_error$Training = sub(".*learn: ", "", iter_error$Training) %>%
    as.numeric() %>% suppressWarnings()
  iter_error$Testing = sub(".*test: ", "", iter_error$Testing) %>%
    as.numeric() %>% suppressWarnings()

  iter_error = iter_error %>%
    filter(!is.na(Training),
           !is.na(Testing)) %>%
    gather(-iter, key = "error_type", value = "error") %>%
    group_by(iter, error_type) %>%
    summarise(mean_error = mean(error),
              sd_error = sd(error))
  return(iter_error)

}

parse_final_train_time = function(out_file_path){
  out_data = read.delim(out_file_path, header = F,
                        na.strings = c("", "NA"))
  iter_time = out_data[,c(1,3)]
  colnames(iter_time) = c("iter", "time")

  iter_time = iter_time[!is.na(iter_time$time),]

  total_time = iter_time$time[str_which(iter_time$iter, "999")]
  total_time = sub(".*: ", "", total_time) %>%
    strsplit(" ") %>%
    unlist() %>%
    str_extract("[0-9]+") %>%
    as.numeric()

  if(length(total_time) == 3){
    total_time_mins = (total_time[1] * 60) + total_time[2] + (total_time[3] / 60)
  }
  else{
    total_time_mins = total_time[1] + (total_time[2] / 60)
  }



  return(total_time_mins)

}

get_arg_max_rel = function(rel_df){
  arg_max_rel = list()
  all_dss = unique(rel_df$ds_id)
  for (ds in all_dss){
    filt = rel_df[rel_df$ds_id == ds,]
    max_rel = max(filt$rel_score,na.rm = T)
    filt = filt[filt$rel_score == max_rel,]
    if(nrow(filt) > 1){
      filt = filt[which.min(filt$FDR_pct),]
    }
    arg_max_rel[[ds]] = filt
  }
  arg_max_rel = arg_max_rel %>% dplyr::bind_rows()
  return(arg_max_rel)

}

merge_rel_annot_FDR = function(max_rel_df, clean_annot_df,
                               kingdom = c("Animal", "Plant")){
  habal = max_rel_df[[match.arg(kingdom)]]
  habal$FDR_pct = as.character(habal$FDR_pct)
  habal = habal %>% left_join(clean_annot_df)
  habal$FDR_pct = "Optim FDR"
  habal = habal[,colnames(Annot_MAP_animal_nofdr)]

  final = rbind.data.frame(clean_annot_df, habal) %>%
    dplyr::filter(score_type == "METASPACE-ML")
  final = final[!is.na(final$FDR_pct),]
  return(final)
}

get_conting_bulk =  function(in_bulk, not_in_bulk, exclusive = T,
                             approach = c("METASPACE-ML", "MSM")){

  q_column = switch(approach,
                    "METASPACE-ML" = "passed_pred",
                    "MSM" = "passed_msm")
  r_column = switch(approach,
                    "METASPACE-ML" = "passed_msm",
                    "MSM" = "passed_pred")

  if(!exclusive){
    TP = in_bulk$ion[which(in_bulk[q_column] == 1)] %>%
      unique()

    FN = in_bulk$ion[in_bulk[q_column] == 0] %>%
      unique()
    FN = FN[FN %nin% TP]

    FP = not_in_bulk$ion[which(not_in_bulk[q_column] == 1)] %>%
      unique()
    FP = FP[FP %nin% c(TP, FN)]
  }
  else{
    TP = in_bulk$ion[which(in_bulk[q_column] == 1 &
                             in_bulk[r_column] == 0)] %>%
      unique()

    FN = in_bulk$ion[in_bulk[q_column] == 0] %>%
      unique()
    FN = FN[FN %nin% TP]

    FP = not_in_bulk$ion[which(not_in_bulk[q_column] == 1 &
                                 not_in_bulk[r_column] == 0)] %>%
      unique()
    FP = FP[FP %nin% c(TP, FN)]

  }

  TP = length(TP)
  FP = length(FP)
  FN = length(FN)

  TN = length(unique(c(not_in_bulk$ion, in_bulk$ion))) - (TP + FN + FP)

  # approach_label = ifelse(exclusive, paste(approach, "only", sep = "_"),
  #                         approach)

  exclusive_label = ifelse(exclusive, "Exclusive", "Inclusive")

  res = data.frame(approach = approach,incl_label = exclusive_label,
                   TP = TP, FN = FN, FP = FP, TN = TN)
  return(res)

}

preprocess_bulk_conting = function(eval_ions, bulk_data, fdr = 0.1,
                                   exclusive = T,
                                   approach = c("METASPACE-ML", "MSM"),
                                   ds_id = NULL){
  if (!is.null(ds_id)){
    eval_ions = eval_ions[eval_ions$ds_id == ds_id,]
  }
  eval_target = eval_ions %>%
    dplyr::filter(target == 1, modifier != "+K") %>%
    mutate(passed_msm = ifelse(msm_fdr < fdr, 1, 0),
           passed_pred = ifelse(pred_score_fdr < fdr, 1, 0))


  eval_bulk = eval_target
  eval_bulk$modifier[which(eval_bulk$modifier == "+H")] = "[M+H]+"
  eval_bulk$modifier[which(eval_bulk$modifier == "+Na")] = "[M+Na]+"
  colnames(eval_bulk)[which(colnames(eval_bulk) == "modifier")] = "Adduct.type"
  colnames(eval_bulk)[which(colnames(eval_bulk) == "formula")] = "Formula"

  not_in_bulk = mismatch_df(eval_bulk, bulk_data, on = c("Formula", "Adduct.type")) %>%
    dplyr::select(ion,ds_id,passed_msm, passed_pred) %>%
    mutate(ion = paste(ion, ds_id, sep = "_")) %>%
    dplyr::select(-ds_id) %>%
    distinct()

  eval_bulk = plyr::match_df(eval_bulk, bulk_data, on = c("Formula", "Adduct.type")) %>%
    dplyr::select(ion,ds_id,passed_msm, passed_pred) %>%
    mutate(ion = paste(ion, ds_id, sep = "_")) %>%
    dplyr::select(-ds_id) %>%
    distinct()

  conting_data = get_conting_bulk(in_bulk = eval_bulk,
                                  not_in_bulk = not_in_bulk,
                                  exclusive = exclusive,
                                  approach = approach)
  return(conting_data)
}

# Nat comms Plotting Functions ---------------------------------------------------------------
plot_context_dataexplorer <- function(meta_df, kingdom = c("Animal", "Plant")) {
  plot_df <- meta_df %>%
    dplyr::filter(Kingdom == kingdom) %>%
    dplyr::select(ds_id, polarity, source, analyzer, mz_class, Sample_type, Species) %>%
    DataExplorer::plot_bar(ggtheme = theme_pubr(), nrow = 2, ncol = 3)

  plot_df <- plot_df[[1]][["data"]]
  plot_df$value = str_replace_all(plot_df$value, c("_and_" = "\nand\n"))
  plot_df$value = str_replace_all(plot_df$value, c("Small_molecules" =
                                                         "Small molecules"))


  plot_df$h_just <- ifelse(plot_df$frequency > 0.25 * max(plot_df$frequency),
                           1, 0
  )

  p <- plot_df %>%
    ggplot(aes(x = reorder(value, frequency), y = frequency)) +
    coord_flip() +
    geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
    geom_text(aes(label = frequency,size = 3),
              position = position_dodge(width = 0.9),
              vjust = -0.25, hjust = plot_df$h_just,
              show.legend = F
    ) +
    facet_wrap(~variable, scales = "free") +
    theme_pubr() +
    # theme(axis.text = element_text(size = 16),
    #       strip.text = element_text(size = 14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 14),
          axis.title.x = element_text(size = 14)) +
    xlab("") +
    ylab("Number of datasets")
  return(p)
}
plot_sunburst_meta <- function(meta_df, kingdom = c("Animal", "Plant"),
                               pctile = 0.5) {
  trial <- meta_df %>%
    dplyr::filter(Kingdom == kingdom) %>%
    count(polarity, source, analyzer, mz_class, Sample_type, Species) %>%
    distinct() %>%
    plotme:::create_all_col_params(T, T)

  y <- str_count(trial$ids, "->")
  z <- trial$values[y == 5] %>% sort()
  diff_cum_sum <- cumsum(z) - sum(z)
  diff_cum_sum = abs(diff_cum_sum)

  idx <- which.min(diff_cum_sum[diff_cum_sum != 0] %% (pctile * sum(z)))

  thresh <- z[idx]

  # Make a list of 24 different colors in hex format
  # colors = c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))
  # names(colors) = unique(trial$labels)
  #
  # trial$colors = apply(trial, 1, function(x){
  #   col = colors[as.character(x["labels"])]
  #   col_list = list("colors" = as.character(col))
  #   as.character(col)
  # })

  trial <- trial[trial$values > thresh, ]

  trial$labels <- gsub("_", "<br>", trial$labels)
  trial$labels <- paste0(trial$labels, "\n", trial$values)



  x <- purrr::exec(plotly::plot_ly, !!!trial,
                   type = "sunburst",
                   branchvalues = "total"
  )
  return(x)
}
plot_new_sankey <- function(meta_df,ds_context_list = NULL,
                            kingdom, sel_dss, plot_title) {

  if(!is.null(ds_context_list)) {
    context_df = ds_context_list %>%
      unlist() %>% as.data.frame() %>%
      rownames_to_column("context")

    colnames(context_df)[2] = "ds_id"
    context_df$context = sub(".*[.]", "", context_df$context)
    context_df$context = sub("[0-9]+.*", "", context_df$context)

    context_df = context_df %>%
      distinct() %>%
      tidyr::separate(col = context,
                      into = c("polarity","mz_class",
                               "analyzer", "source",
                               "Species", "Sample_type"),
                      sep = ",", remove = T)
  }
  else{
    context_df = meta_df
  }

  fig <- sankey_plot_dss(
    meta_df = context_df,
    sel_dss = sel_dss,
    plot_title = plot_title,
    kingdom = kingdom
  )
  return(fig)
}

sankey_plot_dss = function(meta_df, sel_dss, plot_title, kingdom){

  if ("mz_class" %in% colnames(meta_df)) {
    cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = T) %>% unlist()
    context_cols = c("polarity","mz_class","source","analyzer","Species","Sample_type")
    meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
    meta_df_long = meta_df %>% ggsankey::make_long(polarity,mz_class, source,
                                                   analyzer,Species,Sample_type)
    # if (length(unique(meta_df$mz_class)) > 1){
    #
    # }
    # else{
    #   cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = F) %>% unlist()
    #   context_cols = c("polarity","source","analyzer","Species","Sample_type")
    #   meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
    #   meta_df_long = meta_df %>% ggsankey::make_long(polarity, source,
    #                                                  analyzer,Species,Sample_type)
    # }
  }
  else{
    cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = F) %>% unlist()
    context_cols = c("polarity","source","analyzer","Species","Sample_type")
    meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
    meta_df_long = meta_df %>% ggsankey::make_long(polarity, source,
                                                   analyzer,Species,Sample_type)
  }
  names(cols) = sub(".*[.]","",names(cols))

  p = ggplot(meta_df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                               fill = factor(node), label = node)) +
    geom_alluvial(flow.alpha = .6, space = 1) +
    geom_alluvial_text(color = "Black", angle = 0) +
    scale_fill_manual(values = cols) +
    #scale_fill_brewer(palette = "Paired", type = "qual") +
    theme_alluvial(base_size = 15) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5),
          axis.text.x = element_text(size = 15, color = "Black"),
          axis.title.y = element_text(size = 15, color = "Black"),
          axis.text.y = element_text(size = 15, color = "Black")) +
    # scale_x_discrete(labels = addline_format(colnames(train_ds))) +
    ggtitle(plot_title) +
    ylab("Number of Datasets")

  return(p)

}



get_ggmat_genbase <- function(df, plot_type = c(
  "box", "bar", "scatter",
  "intens_box", "2d_bin"
),
x_column = "score_type",
y_column = "metric_value",
color_column = "score_type", in_plot_text_size = 2) {

  if(x_column == "FDR_pct" & plot_type == "box") {
    df$FDR_pct = factor(df$FDR_pct, levels = c("5", "10", "20", "50"))
    if (y_column %in% c("LogDiff", "LFC")){
      df = df[df$score_type == "METASPACE-ML",]
    }
  }
  if(plot_type == "scatter" | plot_type == "2d_bin"){
    df = df[df$score_type == "METASPACE-ML",]
  }

  if (plot_type == "box") {
    counts = df %>%
      dplyr::group_by(.data[[x_column]]) %>%
      dplyr::summarise(n=n())
    n_labels = paste0(counts[[x_column]], "\n(n=", counts$n, ")")
    names(n_labels) = counts[[x_column]]

    gg_base <- ggplot(
      data = df,
      aes(
        x = .data[[x_column]], y = .data[[y_column]],
        color = .data[[color_column]]
      )
    ) +
      geom_boxplot(
        width = 0.5, position = "dodge",
        notch = F, outlier.color = "Black", outlier.shape = NA) +
      geom_jitter(size = 0.5, alpha = 0.3, width = 0.1) +
      stat_summary(fun=mean, geom="errorbar", aes(ymax = ..y..,
                                                  ymin = ..y..), size = 1,
                   linetype = "dashed", width = 0.5)+
      scale_x_discrete(labels = n_labels)
  } else if (plot_type == "bar") {
    if("mz_class" %in% colnames(df)){
      base = ggplot(
        data = df,
        aes(
          x = .data[[x_column]], y = .data[[y_column]],
          fill = .data[[color_column]],linetype = mz_class
        )
      )
    }
    else{
      base = ggplot(
        data = df,
        aes(
          x = .data[[x_column]], y = .data[[y_column]],
          fill = .data[[color_column]]
        )
      )
    }
    gg_base <- base +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = round(.data[[y_column]], 2)),
                position = position_dodge(width = 0.9),
                vjust = 1, hjust = 0, size = in_plot_text_size
      ) +
      geom_errorbar(aes(ymin=.data[[y_column]]-.data[["metric_sd"]], ymax=.data[[y_column]]+.data[["metric_sd"]]),
                    width=.2,
                    position=position_dodge(.9))
    if("mz_class" %in% colnames(df)){
      gg_base = gg_base +
        scale_linetype_manual(values = c("Lipids_and_small_moleclues" = "solid",
                                         "Lipids" = "dashed"))
    }
    gg_base = gg_base +
      scale_fill_manual(values = c("#A6CEE3", "#FB9A99")) +
      theme(legend.position = "none")
  } else if (plot_type == "scatter") {
    gg_base <- ggplot(
      data = df,
      aes(
        x = .data[[x_column]], y = .data[[y_column]],
        color = .data[[color_column]]
      )
    ) +
      geom_point(size = 2, alpha = 0.6) +
      # geom_line(aes(group=.data[[color_column]]), linetype=2) +
      # geom_smooth(method = "lm", se = F, size = 1) +
      # scale_color_manual(values = c("#A6CEE3", "#FB9A99")) +
      theme(legend.position = "none")
  } else if (plot_type == "intens_box") {
    gg_base <- ggplot(data = df, mapping = aes(
      x = .data[[x_column]], y = .data[[y_column]],
      color = log10(.data[[color_column]])
    )) +
      geom_boxplot(outlier.shape = NA, size = 1) +
      geom_jitter(
        size = 3, alpha = 0.8,
        position = position_jitterdodge(dodge.width = 0.9)
      ) +
      scale_y_continuous(n.breaks = 10, limits = c(0, 10)) +
      scale_color_continuous(type = "viridis", name = "Log10 (Number of ions)") +
      # scale_size_binned_area(breaks = c(1,10,50,100,300,500),
      #                        name = "Number of ions", max_size = 10) +
      theme_pubr() +
      # theme(legend.text = element_text(size = 14),
      #       legend.title = element_text(size = 16),
      #       axis.title = element_text(size = 16),
      #       axis.text = element_text(size = 14))+
      geom_signif(
        comparisons = list(c("MSM", "METASPACE-ML_only")),
        map_signif_level = sigFunc, na.rm = T, textsize = 8,
        y_position = c(10, 8, 9, 10)
      )
  } else if (plot_type == "2d_bin"){
    gg_base = ggplot(df, aes(x = .data[[x_column]], y = .data[[y_column]])) +
      geom_bin_2d(binwidth = c(0.5, 0.1)) +
      stat_bin2d(geom = "text", aes(label = ..count..),
                 binwidth = c(0.5, 0.1), size = in_plot_text_size) +
      scale_fill_gradient(low = "snow", high = "red") +
      xlim(-3, 3) +
      ylim(0, 1) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      scale_x_continuous(breaks = seq(-3, 3, 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted")
  }
  else {
    stop("Invalid plot type")
  }
  return(gg_base)
}
subPlot <- function(plot_df, pol,
                    species, plot_type,
                    x_column = "score_type",
                    y_column = "metric_value",
                    color_column = "score_type",
                    hide_x_axis = F, hide_x_axis_label = T,
                    in_plot_text_size = 2){
  plot_df <- plot_df %>%
    filter(polarity == pol, Species == species)

  plot_df$source_analyzer = gsub("_", "\n", plot_df$source_analyzer)
  plot_df$Sample_type = gsub(" ", "\n", plot_df$Sample_type)

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  y_label = switch(y_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change ")
  x_label = switch(x_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change", "score_type" = "Score type")

  if(plot_type %nin% c("paired_box", "2d_bin")){
    gg_base <- get_ggmat_genbase(df = plot_df, plot_type = plot_type,
                                 x_column = x_column,y_column = y_column,
                                 color_column = color_column,
                                 in_plot_text_size = in_plot_text_size)

    p = gg_base +
      scale_color_discrete(type = c("#A6CEE3", "#FB9A99")) +
      ylab(y_label) +
      xlab(x_label) +
      facet_grid(source_analyzer ~ Sample_type,
                 scales = "free"
      ) +
      theme_bw() +
      ggtitle(species) +
      theme(
        strip.background = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      )

    if(hide_x_axis){
      p = p + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    if(hide_x_axis_label){
      p = p + theme(axis.title.x = element_blank())
    }
    p
  }
  else if (plot_type == "2d_bin"){
    gg_base <- get_ggmat_genbase(df = plot_df, plot_type = plot_type,
                                 x_column = x_column,y_column = y_column,
                                 color_column = color_column,
                                 in_plot_text_size = in_plot_text_size)

    p = gg_base +
      ylab(y_label) +
      xlab(x_label) +
      facet_grid(source_analyzer ~ Sample_type,
                 scales = "free"
      ) +
      theme_bw() +
      ggtitle(species) +
      theme(
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )

    if(hide_x_axis){
      p = p + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    if(hide_x_axis_label){
      p = p + theme(axis.title.x = element_blank())
    }
    p
  }
  else{
    p = ggpubr::ggpaired(data = plot_df, x = x_column, y = y_column,
                         color = color_column, linetype = "dotted",
                         line.size = 0.5, point.size = 2,
                         palette = c("#A6CEE3", "#FB9A99"),
                         facet.by = c("source_analyzer","Sample_type"),
                         line.color = "#F4C08B",
                         xlab = FALSE, ylab = y_label,
                         title = species,
                         short.panel.labs = T, id = "ds_id",
                         ggtheme = theme_bw()) +
      theme(
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      stat_compare_means(paired = T, label = "p.format", label.x = 1.5, label.y = 0.65)
    p

  }
}




ggmat_wrapper <- function(plot_df, plot_type,
                          x_column = "score_type",
                          y_column = "metric_value",
                          color_column = "score_type",
                          hide_x_axis = F,hide_x_axis_label = T,
                          in_plot_text_size = 2) {
  all_species <- unique(plot_df$Species)

  plot_df = plot_df %>%
    dplyr::filter(!is.nan(.data[[x_column]]),
                  !is.na(.data[[y_column]]))

  plotList <- list()
  for (i in c("Negative", "Positive")) {
    counter = 0
    for (j in all_species) {
      counter = counter + 1
      p = subPlot(
        plot_df = plot_df,
        pol = i,
        species = j,
        plot_type = plot_type,
        x_column = x_column,
        y_column = y_column,
        color_column = color_column,
        hide_x_axis = hide_x_axis,
        hide_x_axis_label = hide_x_axis_label,
        in_plot_text_size = in_plot_text_size
      )
      p = p + theme(legend.position = "none")
      plotList[[paste0(i, "_", j)]] = p
    }
  }

  plotList <- plotList[!sapply(plotList, is.null)]

  p1 <- cowplot::plot_grid(plotlist = plotList[1:3], nrow = 1, ncol = 3)
  p2 <- cowplot::plot_grid(plotlist = plotList[4:6], nrow = 1, ncol = 3)

  ggmatrix(list(p1, p2),
           nrow = 2, ncol = 1,
           yAxisLabels = c("Negative", "Positive"),
           showStrips = T
  ) +
    theme_pubr(border = T)
}

context_upset <- function(annot_df, eval_df, dss_list, dss = NULL,FDR_pct = 10,
                          include_mz_class = F) {
  create_upset_plot <- function(train_ds_df, intersect_size = 10) {
    comb_mat_list <- list()
    for (context in names(train_ds_df)[-1]) {
      sets <- unique(train_ds_df[, context])
      for (set in sets) {
        comb_mat_list[[set]] <- train_ds_df$ds_id[which(train_ds_df[, context] == set)]
      }
    }
    m <- make_comb_mat(comb_mat_list, mode = "intersect")

    comb_sizes <- comb_size(m)[comb_size(m) >= intersect_size]
    comb_sets <- lapply(comb_name(m), function(nm) extract_comb(m, nm))
    comb_sets <- comb_sets[which(comb_degree(m) >= comb_degree_cutoff)]
    comb_sets <- comb_sets[lengths(comb_sets) >= intersect_size]
    comb_sets <- lapply(comb_sets, function(ds) {
      ds$LFC <- get_LFC_per_dataset(ds)
      ds$diff <- get_diff_per_dataset(ds)
      ds$MAP <- get_MAP_per_dataset(ds)
      ds
    })

    subgroup <- train_ds_df %>%
      dplyr::select(-ds_id) %>%
      distinct() %>%
      gather(key = "cat", value = "val") %>%
      distinct()

    subgroup <- setNames(subgroup$cat, subgroup$val)

    colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B9DC6E")
    colors = colors[1:(ncol(train_ds_df) - 1)]
    names(colors) <- colnames(train_ds_df)[-1]

    filt_m <- m[comb_size(m) >= intersect_size]
    filt_m <- filt_m[comb_degree(filt_m) >= comb_degree_cutoff]

    # highlight = comb_degree(filt_m)
    # highlight[which(substr(names(highlight),5,5) == "1")] = 1
    # highlight[which(substr(names(highlight),5,5) != "1")] = 2
    #
    # highlight_cols = c("red", "black")[highlight]

    upset_plot <- UpSet(filt_m,
                        comb_order = rev(order(comb_size(filt_m))),
                        # comb_col = highlight_cols,
                        top_annotation = upset_top_annotation(filt_m, add_numbers = TRUE),
                        left_annotation = rowAnnotation(
                          Context = subgroup[set_name(filt_m)],
                          show_annotation_name = FALSE, col = list("Context" = colors)
                        ),
                        right_annotation = upset_right_annotation(filt_m),
                        bottom_annotation = HeatmapAnnotation(
                          # LFC = anno_boxplot(lapply(comb_sets, function(ds) ds$LFC), outline = FALSE),
                          Log10_Difference = anno_boxplot(lapply(comb_sets, function(ds) ds$diff), outline = FALSE),
                          median_diff = sapply(comb_sets, function(ds) median(ds$diff)),
                          MAP = sapply(comb_sets, function(ds) unique(ds$MAP)),
                          annotation_name_side = "left", col = list(
                            "median_diff" = circlize::colorRamp2(
                              seq(0, 3, 1),
                              brewer.pal(4, "Blues")
                            ),
                            "MAP" = circlize::colorRamp2(
                              seq(0.1, 0.9, 0.1),
                              brewer.pal(9, "Reds")
                            )
                          ),
                          annotation_label = c(
                            "Log10\nDifference", "Median Difference",
                            "MAP"
                          )
                        )
    )
    return(upset_plot)
  }
  get_LFC_per_dataset <- function(ds_ids) {
    LFC <- c()
    for (i in ds_ids) {
      LFC <- c(LFC, LFC_annots$Delta[which(LFC_annots$ds_id == i)])
    }
    return(log2(LFC))
  }
  get_diff_per_dataset <- function(ds_ids) {
    diff <- c()
    for (i in ds_ids) {
      diff <- c(diff, LFC_annots$diff[which(LFC_annots$ds_id == i)])
    }
    return(diff)
  }
  get_MAP_per_dataset <- function(ds_ids) {
    filtered <- eval_df[which(eval_df$ds_id %in% ds_ids), ] %>%
      dplyr::filter(score_type == "METASPACE-ML")
    filtered <- filtered[!is.nan(filtered$metric_value), ]
    filtered <- filtered[!is.na(filtered$metric_value), ]
    MAP <- mean(filtered$metric_value)
    return(MAP)
  }
  # debug(get_LFC_per_dataset)

  testing_ds = dss_list %>%
    as.data.frame(check.names = F) %>%
    gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = T)
  if (!include_mz_class){
    testing_ds <- testing_ds %>%
      dplyr::select(-mz_class)
  }
  testing_ds = testing_ds %>%
    distinct() %>%
    relocate(ds_id)


  comb_degree_cutoff = ncol(testing_ds) - 1

  if (!is.null(dss)) {
    testing_ds <- testing_ds[which(testing_ds$ds_id %in% dss), ]
  }

  eval_df <- eval_df %>% tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",")
  eval_df$score_type = str_replace_all(eval_df$score_type,
                                       c("Catboost" = "METASPACE-ML"))

  LFC_annots <- annot_df[which(annot_df$FDR_pct == as.character(FDR_pct)), ] %>%
    dplyr::select(-FDR_pct) %>%
    dplyr::distinct() %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",") %>%
    group_by(ds_id) %>%
    summarise(
      msm_fdr_per_ds = sum(msm_fdr_annots),
      pred_fdr_per_ds = sum(pred_fdr_annots)
    ) %>%
    as.data.frame()

  LFC_annots$diff <- LFC_annots$pred_fdr_per_ds - LFC_annots$msm_fdr_per_ds
  LFC_annots$diff_sign <- ifelse(LFC_annots$diff < 0, -1, 1)
  LFC_annots$diff <- LFC_annots$diff_sign * log10(abs(LFC_annots$diff) + 1)

  LFC_annots$Delta <- LFC_annots$pred_fdr_per_ds / LFC_annots$msm_fdr_per_ds
  LFC_annots$Delta[is.infinite(LFC_annots$Delta)] <- max(LFC_annots$Delta[!is.infinite(LFC_annots$Delta)])
  LFC_annots <- LFC_annots[order(LFC_annots$Delta, decreasing = T), ]

  if (any(str_detect(LFC_annots$ds_id, "_ml_training"))) {
    LFC_annots$ds_id <- sub("_ml_training.*", "", LFC_annots$ds_id)
  }

  testing_ds <- testing_ds[which(testing_ds$ds_id %in% LFC_annots$ds_id), ]

  # stat_data = left_join(testing_ds, LFC_annots[,c("ds_id", "diff")])
  # stat_data$rp_grp = ifelse(stat_data$rp_range == "Low", "Low", "Medium-High")
  # stat_data = stat_data %>% dplyr::select(ds_id, rp_grp, diff)
  # stat_boxplot = ggbetweenstats(data = stat_data, x = rp_grp, y = diff,
  #                               type = "nonparametric", p.adjust.method = "BH",
  #                               package = "ggsci", palette = "default_jama", results.subtitle = F,
  #                               centrality.plotting = T) +
  #   theme(axis.text = element_text(size = 18),
  #         axis.title = element_text(size = 20)) +
  #   xlab("") +
  #   ylab("Log10 Difference") +
  #   geom_signif(comparisons = list(c("Low", "Medium-High")), map_signif_level = T, textsize = 8)
  #
  #
  # stat_data = stat_data %>% spread(key = rp_grp, value = diff)
  #
  # low_diff = stat_data$Low[!is.na(stat_data$Low)]
  # not_low_diff = stat_data$`Not low`[!is.na(stat_data$`Not low`)]
  #
  # stat_res = wilcox.test(low_diff, not_low_diff, alternative = "less", conf.int = T)
  #
  # if (return_stats_only){
  #   return(list("stat_res" = stat_res,
  #               "stat_boxplot" = stat_boxplot))
  # }

  p <- create_upset_plot(train_ds_df = testing_ds, intersect_size = 10)
  return(p)
} # upset plot landscape 10x15

plot_PR_per_ds <- function(raw_res, dataset_id, plot_type = c("PR", "ROC"),
                           feat_msm = NULL, plot = T) {
  ds_auc <- list()
  if (is.null(feat_msm)) {
    trial <- raw_res %>%
      dplyr::select(ds_id, msm, pred_score, target) %>%
      dplyr::filter(ds_id == dataset_id) %>%
      distinct()
  } else {
    trial <- raw_res[, c("ds_id", feat_msm, "pred_score", "target")]
    colnames(trial)[2] <- "msm"
    trial <- trial %>%
      dplyr::filter(ds_id == dataset_id) %>%
      distinct()
  }

  if(length(unique(trial$target)) == 1){
    return(NULL)
  }

  msm_pr <- PRROC::pr.curve(
    scores.class0 = trial$msm,
    weights.class0 = trial$target, curve = T,
    min.compute = T,
    max.compute = T, rand.compute = T
  )

  pred_pr <- PRROC::pr.curve(
    scores.class0 = trial$pred_score,
    weights.class0 = trial$target, curve = T,
    min.compute = T,
    max.compute = T, rand.compute = T
  )

  msm_roc <- roc.curve(
    scores.class0 = trial$msm,
    weights.class0 = trial$target,
    curve = T, min.compute = F,
    max.compute = F, rand.compute = T
  )

  pred_roc <- roc.curve(
    scores.class0 = trial$pred_score,
    weights.class0 = trial$target,
    curve = T, min.compute = T,
    max.compute = T, rand.compute = T
  )

  msm_feat = ifelse(is.null(feat_msm), "MSM", feat_msm)

  if (!plot){
    if (plot_type == "PR") {

      metric_res = data.frame(score_type = c(msm_feat, "METASPACE-ML"),
                              AUC = c(msm_pr$auc.integral,
                                      pred_pr$auc.integral))
      curve_data = data.frame(Recall = c(msm_pr$curve[, 1],
                                    pred_pr$curve[, 1]),
                              Precision = c(msm_pr$curve[, 2],
                                    pred_pr$curve[, 2]),
                              score_type = c(rep(msm_feat, nrow(msm_pr$curve)),
                                             rep("METASPACE-ML", nrow(pred_pr$curve))))
    }
    else {
      metric_res = data.frame(score_type = c(msm_feat, "METASPACE-ML"),
                              AUC = c(msm_roc$auc,
                                      pred_roc$auc))
      curve_data = data.frame(FPR = c(msm_roc$curve[, 1],
                                    pred_roc$curve[, 1]),
                              Sensitivity = c(msm_roc$curve[, 2],
                                    pred_roc$curve[, 2]),
                              score_type = c(rep(msm_feat, nrow(msm_roc$curve)),
                                             rep("METASPACE-ML", nrow(pred_roc$curve))))
    }
    final = left_join(curve_data, metric_res, by = "score_type")
    return(final)
  }

  if (plot_type == "PR") {
    msm_obj <- msm_pr
    pred_obj <- pred_pr

    plot(msm_obj,
         max.plot = F, min.plot = F,
         rand.plot = TRUE, fill.area = F,
         color = "#1c61b6", auc.main = FALSE,
         main = dataset_id
    )
    plot(pred_obj, add = TRUE, color = "#008600")
    text(0.03, 1, labels = paste("AUC = ", round(msm_obj$auc.integral, 4)), adj = c(0, 0.5), col = "#1c61b6")
    text(0.03, 0.9, labels = paste("AUC = ", round(pred_obj$auc.integral, 4)), adj = c(0, 0.5), col = "#008600")
    legend("topright", legend = c("MSM", "METASPACE-ML"), col = c("#1c61b6", "#008600"), lwd = 2)
  } else {
    msm_obj <- msm_roc
    pred_obj <- pred_roc

    plot(msm_obj,
         max.plot = F, min.plot = F,
         rand.plot = TRUE, fill.area = F,
         color = "#1c61b6", auc.main = FALSE,
         main = dataset_id
    )
    plot(pred_obj, add = TRUE, color = "#008600")
    text(0.03, 1, labels = paste("AUC = ", round(msm_obj$auc, 4)), adj = c(0, 0.5), col = "#1c61b6")
    text(0.03, 0.9, labels = paste("AUC = ", round(pred_obj$auc, 4)), adj = c(0, 0.5), col = "#008600")
    legend("bottomright", legend = c("MSM", "METASPACE-ML"), col = c("#1c61b6", "#008600"), lwd = 2)
  }
}


plot_PCA_per_ds <- function(pca_df_list, center_pca = T, scale_pca = F) {
  trial_pca <- prcomp(t(pca_df_list$feature_mat), scale. = scale_pca, center = center_pca)

  score_type_df <- pca_df_list$score_mat %>%
    t() %>%
    as.data.frame()

  decoy_target <- ggbiplot(trial_pca,
                           obs.scale = 0, var.scale = 1,
                           groups = as.factor(pca_df_list$labels), ellipse = T, circle = F,
                           choices = c(1, 2), alpha = 0.7, scale = 1, pc.biplot = T,
                           var.axes = T, varname.size = 5, varname.adjust = 2, varname.color = "black") +
    scale_color_discrete(name = "") +
    ggpubr::theme_pubr() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2)

  MSM_fill <- ggbiplot(trial_pca,
                       obs.scale = 0, var.scale = 1,
                       ellipse = T, circle = F,
                       choices = c(1, 2), alpha = 0.1, scale = 1, pc.biplot = T,
                       var.axes = T, varname.adjust = 2, varname.size = 5, varname.color = "black"
  ) +
    geom_point(aes(col = score_type_df$MSM, alpha = 0.3)) +
    scale_color_viridis_c() +
    ggpubr::theme_pubr() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2)

  ML_fill <- ggbiplot(trial_pca,
                      obs.scale = 0, var.scale = 1,
                      ellipse = T, circle = F,
                      choices = c(1, 2), alpha = 0.1, scale = 1, pc.biplot = T,
                      var.axes = T, varname.adjust = 2, varname.size = 5, varname.color = "black"
  ) +
    geom_point(aes(col = score_type_df$CatBoost, alpha = 0.3)) +
    scale_color_viridis_c() +
    ggpubr::theme_pubr() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2)

  return(list(
    "decoy_target" = decoy_target,
    "MSM_fill" = MSM_fill, "ML_fill" = ML_fill,
    "pca_res" = trial_pca,
    "orig_data" = t(pca_df_list$feature_mat)
  ))
}

PCA_per_ds <- function(raw_res, ds_id = NULL, sf = NULL, fix_mz_err = T,
                       center_pca = T, scale_pca = F) {
  pca_df_list <- prepare_data_umap(
    raw_res = raw_res, ds_id = ds_id, sf = sf,
    fix_mz_err = fix_mz_err
  )
  plot_res <- plot_PCA_per_ds(
    pca_df_list = pca_df_list, center_pca = center_pca,
    scale_pca = scale_pca
  )
  return(plot_res)
}


make_density_hm_context = function(df, FDR = 10,
                                   vals = c("Diff", "LogDiff", "LFC", "MAP"),
                                   kingdom = c("Animal", "Plant"), col_split = NULL,
                                   plot_title){

  val_col = switch(vals, "Diff" = "Difference", "LogDiff" = "LogDiff",
                   "LFC" = "LFC", "MAP" = "metric_value")

  y_label = switch(vals, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change ")

  anno_df = df %>%
    dplyr::select(context) %>%
    as.data.frame(check.names = F) %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    # dplyr::select(-mz_class) %>%
    distinct() %>%
    column_to_rownames(var = "context")

  ds_context = df %>%
    dplyr::select(ds_id, context) %>%
    distinct()

  d_hm_mat = df %>%
    dplyr::filter(FDR_pct == FDR,
                  ds_id %in% ds_context$ds_id,
                  score_type == "METASPACE-ML")
  d_hm_mat = d_hm_mat[,which(colnames(d_hm_mat) %in% c("ds_id", "context", val_col))]
  colnames(d_hm_mat)[which(colnames(d_hm_mat) == val_col)] = "values"
  d_hm_mat = d_hm_mat %>%
    pivot_wider(names_from = context, values_from = values) %>%
    column_to_rownames(var = "ds_id") %>%
    as.matrix()

  dm_cols = get_d_hm_colors_context(kingdom = kingdom, add_mz_class = T)


  p =
    densityHeatmap(data = d_hm_mat, show_column_names = F, cluster_columns = T,
                   show_quantiles = T, title = plot_title,
                   ylab = y_label, quantile_gp = gpar(fontsize = 10),
                   show_heatmap_legend = T,column_split = col_split
    ) %v%
    HeatmapAnnotation(df = anno_df,
                      col = dm_cols,show_legend = T,
                      annotation_legend_param = list("labels_gp" =
                                                       gpar(fontsize = 12),
                                                     "title_gp" =
                                                       gpar(fontsize = 14, fontface = "bold")))


  return(p)


}

make_density_hm_feat_imp = function(df, kingdom = c("Animal", "Plant"),
                                    shap_abs_prop = F, col_split = NULL){

  if(shap_abs_prop){
    df$abs_pred = rowSums(abs(df[,c("chaos", "spatial", "spectral", "mz_err_abs_abserr", "mz_err_rel_abserr")]))
    df = df %>%
      mutate(chaos = abs(chaos)/abs_pred,
             spatial = abs(spatial)/abs_pred,
             spectral = abs(spectral)/abs_pred,
             mz_err_abs_abserr = abs(mz_err_abs_abserr)/abs_pred,
             mz_err_rel_abserr = abs(mz_err_rel_abserr)/abs_pred)
  }

  df = df %>%
    distinct() %>%
    dplyr::group_by(ds_id, context) %>%
    dplyr::summarize(chaos = median(chaos),
              spatial = median(spatial),
              spectral = median(spectral),
              mz_err_abs_abserr = median(mz_err_abs_abserr),
              mz_err_rel_abserr = median(mz_err_rel_abserr),
              .groups = "drop")

  if(kingdom == "Animal"){
    df$context = fix_species_context(df$context)
  }

  ds_context = df %>%
    dplyr::select(ds_id, context) %>%
    distinct()

  anno_df = df %>%
    dplyr::select(context) %>%
    dplyr::distinct() %>%
    as.data.frame(check.names = F) %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    # dplyr::select(-mz_class) %>%
    dplyr::distinct() %>%
    column_to_rownames(var = "context")

  dm_cols = get_d_hm_colors_context(kingdom = kingdom, add_mz_class = T)


  d_hm_plots = list()
  for (feat in c("spectral", "spatial", "chaos", "mz_err_abs_abserr",
                 "mz_err_rel_abserr")){
    d_hm_mat = df[,which(colnames(df) %in% c("ds_id", "context", feat))] %>%
      distinct()
    colnames(d_hm_mat)[which(colnames(d_hm_mat) == feat)] = "values"
    d_hm_mat = d_hm_mat %>%
      pivot_wider(names_from = context, values_from = values) %>%
      column_to_rownames(var = "ds_id") %>%
      as.matrix()

    p =
      densityHeatmap(data = d_hm_mat, show_column_names = F, cluster_columns = T,
                     show_quantiles = F,
                     ylab = paste0(gsub("_abserr", "", feat),
                                   "\n", "SHAP"),
                     show_heatmap_legend = T,
                     title = "Feature importance\n Per context",column_split = col_split
      )
    d_hm_plots[[feat]] = p
  }

  final_plot = d_hm_plots[["spectral"]] %v%
    d_hm_plots[["spatial"]] %v%
    d_hm_plots[["mz_err_rel_abserr"]] %v%
    d_hm_plots[["mz_err_abs_abserr"]] %v%
    d_hm_plots[["chaos"]] %v%

    HeatmapAnnotation(df = anno_df,
                      col = dm_cols,show_legend = T,
                      annotation_legend_param = list("labels_gp" =
                                                       gpar(fontsize = 12),
                                                     "title_gp" =
                                                       gpar(fontsize = 14, fontface = "bold")))


  return(final_plot)

}

make_density_hm_param_optim = function(df, exclude_feat_type = T){

  if(exclude_feat_type){
    df = df %>%
      filter(features_type == "nonFDR") %>%
      select(-features_type)
  }

  anno_df = df[,colnames(df) %in%
                 c("includes_bad_rank", "ranking_type",
                   "catmonotonic", "features_type")] %>%
    mutate(optim_comb = do.call(paste, c(., sep = ","))) %>%
    distinct()

  d_hm_mat = df %>%
    left_join(anno_df) %>%
    dplyr::filter(score_type == "Catboost") %>%
    dplyr::select(optim_comb, metric_value, group_name) %>%
    distinct() %>%
    pivot_wider(names_from = optim_comb, values_from = metric_value) %>%
    column_to_rownames(var = "group_name") %>%
    as.matrix()

  anno_df = anno_df %>%
    column_to_rownames(var = "optim_comb")


  dm_cols = list("includes_bad_rank" = c("Yes" = "#A6CEE3", "No" = "#1F78B4"),
                 "ranking_type" = c("separate" = "#B2DF8A", "single" = "#33A02C"),
                 "catmonotonic" = c("Yes" = "#FB9A99", "No" = "#E31A1C"),
                 "features_type" = c("Both" = "#8DD3C7", "nonFDR" = "#FFFFB3", "FDRonly" = "#BC80BD"))


  p =
    densityHeatmap(data = d_hm_mat, show_column_names = F, cluster_columns = T,
                   show_quantiles = T, title = "Parameter optimization",
                   ylab = "MAP", quantile_gp = gpar(fontsize = 10),
                   show_heatmap_legend = T
    ) %v%
    HeatmapAnnotation(df = anno_df,
                      col = dm_cols,show_legend = T)

  return(p)


}


plot_context_hm = function(df,kingdom = "Animal",
                           cols_column = "context_size",
                           vals_column = "AUC"){

  anno_df = df[,colnames(df) %in% c("context","polarity","mz_class",
                                    "analyzer", "source",
                                    "Species", "Sample_type")] %>%
    distinct() %>%
    column_to_rownames(var = "context")

  mat_df = df %>% dplyr::select(ds_id,context, {{cols_column}}, {{vals_column}}) %>%
    group_by(context, .data[[cols_column]]) %>%
    summarise(vals = mean(.data[[vals_column]])) %>%
    distinct() %>%
    pivot_wider(names_from = context, values_from = vals) %>%
    column_to_rownames(var = cols_column) %>%
    as.matrix()

  label_mat = round(mat_df, 2)
  label_mat[is.na(label_mat)] = ""

  colors = get_d_hm_colors_context(kingdom = kingdom, add_mz_class = T)

  pheatmap(mat = mat_df,
           color = colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(100),
           display_numbers = label_mat,annotation_col = anno_df,
           show_colnames = F,cluster_rows = F, cluster_cols = T,
           annotation_colors = colors,number_color = "black",
           legend = T,annotation_legend = T,border_color = "white",
           fontsize_number = 12,fontsize = 12,
           main = paste0("Context size vs ", vals_column))
}

plot_context_hm_adjacency = function(df,kingdom = "Animal",
                                     vals_column = "AUC",
                                     row_column = "group1",
                                     col_column = "group2",
                                     label_column = "p.adju.signif"){

  row_anno_df = df %>%
    dplyr::select({{row_column}}) %>%
    distinct() %>%
    as.data.frame(check.names = F) %>%
    tidyr::separate(col = .data[[row_column]],
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    column_to_rownames(var = row_column)

  col_anno_df = df %>%
    dplyr::select({{col_column}}) %>%
    distinct() %>%
    as.data.frame(check.names = F) %>%
    tidyr::separate(col = .data[[col_column]],
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    column_to_rownames(var = col_column)


  mat_df = df %>% dplyr::select({{row_column}}, {{col_column}},{{vals_column}}) %>%
    pivot_wider(names_from = .data[[col_column]], values_from = .data[[vals_column]]) %>%
    column_to_rownames(var = row_column) %>%
    as.matrix()

  label_mat = df %>% dplyr::select({{row_column}}, {{col_column}},{{label_column}}) %>%
    pivot_wider(names_from = .data[[col_column]], values_from = .data[[label_column]]) %>%
    column_to_rownames(var = row_column) %>%
    as.matrix()

  label_mat[is.na(label_mat)] = ""
  mat_df[is.na(mat_df)] = 0

  colors = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = T)

  mat_colors = colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(100)
  mat_colors[1] = "#FFFFFF"

  p = ComplexHeatmap::pheatmap(mat = mat_df,
           color = mat_colors,
           display_numbers = label_mat,annotation_col = col_anno_df,
           annotation_row = row_anno_df,show_rownames = F,
           show_colnames = F,cluster_rows = F, cluster_cols = F,
           annotation_colors = colors,number_color = "black",
           legend = T,annotation_legend = T,border_color = "black",
           fontsize_number = 12,fontsize = 12, na_col = "white",heatmap_legend_param = list(title = "-Log10 p-value"))

  p@top_annotation@anno_list = lapply(p@top_annotation@anno_list, function(x) {
    x@show_legend = FALSE
    x
  })
  return(p)
}

plot_enrich_hm_context = function(enrich_res, filter_by_adjpval = T,
                                  min_TP = 3, use_FE = T,
                                  kingdom = c("Animal", "Plant"),
                                  min_dss_prop = 0.1){
  enrich_res = enrich_res %>%
    dplyr::filter(TP > min_TP, term != "") %>%
    dplyr::group_by(ds_id) %>%
    dplyr::mutate(adj_p_val = p.adjust(pval, method = "BH"))

  if (filter_by_adjpval){
    enrich_res = enrich_res %>%
      dplyr::filter(adj_p_val < 0.05)
  }
  else{
    enrich_res = enrich_res %>%
      dplyr::filter(pval < 0.05)
  }
  if (use_FE){
    enrich_res$OR = log2(enrich_res$FE)
  }
  else{
    enrich_res$OR = log2(enrich_res$OR)
  }

  context_df = testing_dss[[match.arg(kingdom)]] %>%
    as.data.frame(check.names = F) %>%
    gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F)

  plot_df = enrich_res %>% dplyr::left_join(context_df, by = "ds_id") %>%
    dplyr::select(ds_id, context, term, OR) %>%
    distinct() %>%
    dplyr::group_by(term, context) %>%
    dplyr::summarise(OR = mean(OR),
              n_sig_ds = n_distinct(ds_id),.groups = "drop")

  wide_mat = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>%
    dplyr::select(term, context, OR) %>%
    pivot_wider(names_from = context, values_from = OR) %>%
    column_to_rownames(var = "term") %>%
    as.matrix()

  anno_df = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>%
    dplyr::select(context) %>%
    distinct() %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    column_to_rownames(var = "context")

  label_mat = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>%
    dplyr::select(term, context, n_sig_ds) %>%
    pivot_wider(names_from = context, values_from = n_sig_ds) %>%
    column_to_rownames(var = "term") %>%
    as.matrix()

  label_mat[is.na(label_mat)] = ""

  colors = get_d_hm_colors_context(kingdom = match.arg(kingdom),
                                   add_mz_class = T)

  wide_mat[is.na(wide_mat)] = 0

  ComplexHeatmap::pheatmap(mat = wide_mat,
           color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
           display_numbers = label_mat,annotation_col = anno_df,
           show_colnames = F,cluster_rows = F, cluster_cols = F,
           annotation_colors = colors,number_color = "black",
           legend = T,annotation_legend = T,
           fontsize_number = 12,fontsize = 12,
           main = paste0("Enrichment per context\n", match.arg(kingdom)),
           na_col = "white", border_color = "white",
           heatmap_legend_param = list(title = "Log2 fold enrichment"))

}

plot_intens_hm_context = function(intens_res,
                                  kingdom = c("Animal", "Plant"),
                                  scale_intens = T){

    context_df = testing_dss[[match.arg(kingdom)]] %>%
    as.data.frame(check.names = F) %>%
    gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F)

  plot_df = intens_res %>%
    dplyr::left_join(context_df, by = "ds_id") %>%
    dplyr::select(ds_id, context, Group, median_intens) %>%
    distinct() %>%
    dplyr::group_by(Group, context) %>%
    dplyr::summarise(median_intens = median(median_intens),
                     n_sig_ds = n_distinct(ds_id),.groups = "drop")

  if (scale_intens){
    wide_mat = plot_df %>%
      dplyr::select(Group, context, median_intens) %>%
      dplyr::group_by(context) %>%
      dplyr::mutate(median_intens = minmax_intens(median_intens)) %>%
      pivot_wider(names_from = context, values_from = median_intens) %>%
      column_to_rownames(var = "Group") %>%
      as.matrix()

    col_scale = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(10)
    hm_title = "Scaled\nMedian\n(Log10 Intensity)"
  }
  else{
    wide_mat = plot_df %>%
      dplyr::select(Group, context, median_intens) %>%
      pivot_wider(names_from = context, values_from = median_intens) %>%
      column_to_rownames(var = "Group") %>%
      as.matrix()
    col_scale = colorRampPalette(brewer.pal(n = 6, name = "Blues"))(10)
    hm_title = "Median\n(Log10 Intensity)"
  }

  anno_df = plot_df %>%
    dplyr::select(context) %>%
    distinct() %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    column_to_rownames(var = "context")

  label_mat = plot_df %>%
    dplyr::select(Group, context, median_intens) %>%
    pivot_wider(names_from = context, values_from = median_intens) %>%
    column_to_rownames(var = "Group") %>%
    as.matrix()

  label_mat = round(label_mat, 1)

  label_mat[is.na(label_mat)] = ""

  colors = get_d_hm_colors_context(kingdom = match.arg(kingdom),
                                   add_mz_class = T)

  ComplexHeatmap::pheatmap(mat = wide_mat,
           color = col_scale,
           display_numbers = label_mat,annotation_col = anno_df,
           show_colnames = F,cluster_rows = F, cluster_cols = F,
           annotation_colors = colors,number_color = "black",
           legend = T,annotation_legend = T,
           fontsize_number = 6,fontsize = 12,
           main = paste0("Intensity per context\n", match.arg(kingdom)),
           na_col = "white", border_color = "white",cellheight = 18,
           heatmap_legend_param = list("title" = hm_title))

}

ggstats_wrapper = function(df, x_column, y_column,
                           paired = F, centrality.plotting = T,
                           results.subtitle = F,parametric = F,
                           sig_only = F){

  y_label = switch(y_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change ")
  x_label = switch(x_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change", "score_type" = "Score type")

  ggstats_func = ifelse(paired, ggstatsplot::ggwithinstats, ggstatsplot::ggbetweenstats)

  param_type = ifelse(parametric, "parametric", "nonparametric")

  p = ggstats_func(data = df,x = {{x_column}},y = {{y_column}},type = param_type,
                     p.adjust.method = "BH", pairwise.display = "none",
                     results.subtitle = results.subtitle ,pairwise.comparisons = T,
                   centrality.plotting = centrality.plotting)
  if(length(unique(df[[x_column]])) == 2){
    return(p)
  }

  p_comparison = extract_stats(p)[["pairwise_comparisons_data"]] %>%
    dplyr::mutate(groups = purrr::map2(group1, group2, function(x,y){c(x,y)})) %>%
    dplyr::mutate(p_astersik = sigFunc(p.value))

  p_comparison$p_format = ifelse(p_comparison$p.value < 0.05,
                                 format(p_comparison$p.value,
                                        scientific = T,digits = 2),
                                 round(p_comparison$p.value,2))
  p_comparison$p_format = paste0("p=",  p_comparison$p_format)

  if (sig_only){
    p_comparison = p_comparison %>%
      dplyr::filter(p_astersik != "ns")
  }

  g = ggstats_func(data = df,x = {{x_column}},y = {{y_column}},type = param_type,
                     p.adjust.method = "BH", pairwise.display = "none",
                     results.subtitle = F,pairwise.comparisons = F,
                     centrality.plotting = T, plot.type = "box")

  final_p = g + ggsignif::geom_signif(
    comparisons = p_comparison$groups,
    map_signif_level = T,
    annotations = p_comparison$p_format,
    test = NULL,
    y_position = ggstatsplot:::.ggsignif_xy(pull(df,.data[[x_column]])
                                            , pull(df,.data[[y_column]])),
    na.rm = T,
    textsize = 5
  ) +
    scale_y_continuous(n.breaks = 10) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 12)) +
    xlab(x_label) +
    ylab(y_label)

  return(final_p)

}



# Nat comms Utility functions -------------------------------------------------------
get_d_hm_colors_context = function(kingdom = c("Animal", "Plant"), add_mz_class = F){
  if (kingdom == "Animal"){
    cols = list("polarity" = c("Positive" = "#A6CEE3", "Negative" = "#1F78B4"),
                "analyzer" = c("Orbitrap" = "#B2DF8A", "FTICR" = "#33A02C"),
                "source" = c("MALDI" = "#FB9A99", "DESI" = "#E31A1C"),
                "Species" = c("Homo sapiens" = "#8DD3C7", "Mus musculus" = "#FFFFB3", "OTHER" = "#BC80BD"),
                "Sample_type" = c("Tissue" = "#FDBF6F", "Uncurated" = "#CAB2D6", "Whole Organism" = "#DFC497", "Cells" = "#664098"))
  }
  else{
    cols = list("polarity" = c("Positive" = "#A6CEE3", "Negative" = "#1F78B4"),
                "analyzer" = c("Orbitrap" = "#B2DF8A", "FTICR" = "#33A02C"),
                "source" = c("MALDI" = "#FB9A99", "DESI" = "#E31A1C"),
                "Species" = c("Populus" = "#8DD3C7", "Sorghum" = "#FFFFB3", "OTHER" = "#BC80BD"),
                "Sample_type" = c("Tissue" = "#FDBF6F", "Uncurated" = "#CAB2D6"))
  }
  if (add_mz_class){
    cols[["mz_class"]] = c("Lipids" = "#D1EAC6", "Lipids_and_small_moleclues" = "#62B5C2")
  }

  return(cols)
}

fix_species_context = function(context){
  spec = c("Homo sapiens", "Mus musculus", "OTHER", "Whole Organism", "DESI")
  names(spec) = c("Homo,sapiens", "Mus,musculus", "OTHER", "Whole,Organism", "ESI")

  str_replace_all(context, spec)
}


sigFunc <- function(x) {
  ifelse(x < 0.001, "***",
         ifelse(x < 0.01, "**",
                ifelse(x < 0.05, "*", "ns")))
}

mismatch_df = function (x, y, on = NULL)
{
  if (is.null(on)) {
    on <- intersect(names(x), names(y))
    message("Matching on: ", paste(on, collapse = ", "))
  }
  keys <- plyr::join.keys(x, y, on)
  x[keys$x %nin% keys$y, , drop = FALSE]
}

minmax_intens = function(x,target_min = -1, target_max = 1){
  target_min + (((x-min(x)) * (target_max - target_min))/(max(x) - min(x)))
}


# Nat Comms wrapper functions ---------------------------------------------

Performance_plots_pipeline = function(eval_df, annot_df,
                                      data_type = c("Training", "Testing"),
                                      context_size = 50,
                                      kingdom = c("Animal", "Plant"),
                                      mean_eval_context = F,
                                      metrics_per_grp = F,
                                      FDR.pct = NULL,
                                      plot_type = c("box", "bar", "scatter",
                                        "intens_box", "paired_box"),
                                      x_column = "score_type",
                                      y_column = "metric_value",
                                      color_column = "score_type",
                                      hide_x_axis = F,hide_x_axis_label = T,
                                      context_specific = T,
                                      in_plot_text_size = 2,
                                      plot_title = NULL){

  y_label = switch(y_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                               "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                               "LFC" = "Log2 Fold Change ")
  x_label = switch(x_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change", "score_type" = "Score type")

  if(!context_specific){
    plot_df = prepare_plot_df(eval_df = eval_df,
                              annot_df = annot_df,
                              data_type = data_type, context_size = context_size,
                              kingdom = match.arg(kingdom),metrics_per_grp = metrics_per_grp,
                              mean_eval_context = F,FDR.pct = FDR.pct)
    plot_df = plot_df[,!colnames(plot_df) %in% c("context","polarity","mz_class",
                                                 "analyzer", "source",
                                                 "Species", "Sample_type", "source_analyzer")] %>%
      distinct() %>%
      filter(!is.na(.data[[x_column]]), !is.na(.data[[y_column]]))
    if (plot_type == "bar"){
      plot_df = plot_df %>%
        dplyr::group_by(.data[[x_column]], metric_type) %>%
        dplyr::summarise(metric_sd = sd(metric_value,na.rm = T),
                         metric_value = mean(metric_value, na.rm = T)) %>%
        dplyr::ungroup()
    }

    if (plot_type == "paired_box"){
      PR_Data = plot_df %>%
        dplyr::select(ds_id,{{x_column}},{{y_column}}) %>%
        distinct() %>%
        spread(key = .data[[x_column]], value = .data[[y_column]])

      PR_Data$Delta = ifelse(PR_Data$`METASPACE-ML` > PR_Data$MSM, "Higher", "Lower")
      PR_Data = PR_Data[!is.na(PR_Data$Delta),]

      Delta_count = PR_Data %>%
        group_by(Delta) %>%
        summarise(n_ds = n_distinct(ds_id))

      PR_Data = PR_Data %>%
        left_join(Delta_count, by = "Delta") %>%
        mutate(Delta = paste(Delta, "\n(", n_ds, ")", sep = ""))

      p = PR_Data %>% ggpaired(cond1 = "METASPACE-ML", cond2 = "MSM",
                               color = "condition", linetype = "dotted",
                               line.size = 0.5, point.size = 2,
                               palette = c("#A6CEE3", "#FB9A99"),
                               ggtheme = theme_pubr(legend = "bottom"),
                               facet.by = "Delta",line.color = "#F4C08B",
                               xlab = FALSE, ylab = "PR-AUC", title = "Paired PR-AUC per dataset",
                               short.panel.labs = T) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 20),
              legend.title = element_blank(),
              strip.text = element_text(size = 20),
              strip.background = element_rect(fill = "white")) +
        stat_compare_means(aes(group = condition), label = "p.format" ,paired = F,
                           size = 8, label.x.npc = "center", label.y.npc = 0.95)
    }
    else{
      gg_base = get_ggmat_genbase(df = plot_df, plot_type = plot_type, x_column = x_column,
                                  y_column = y_column, color_column = color_column,
                                  in_plot_text_size = in_plot_text_size)

      p = gg_base +
        theme_pubr() +
        ylab(y_label) +
        xlab(x_label) +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14))
    }
  }
  else{
    plot_df = prepare_plot_df(eval_df = eval_df,
                              annot_df = annot_df,
                              data_type = data_type, context_size = 50,
                              kingdom = match.arg(kingdom),metrics_per_grp = metrics_per_grp,
                              mean_eval_context = mean_eval_context,FDR.pct = FDR.pct) %>%
      filter(!is.na(.data[[x_column]]), !is.na(.data[[y_column]]))

    p = ggmat_wrapper(plot_df = plot_df, plot_type = plot_type,
                  x_column = x_column, y_column = y_column,
                  color_column = color_column ,hide_x_axis = hide_x_axis,
                  hide_x_axis_label = hide_x_axis_label,
                  in_plot_text_size = in_plot_text_size)
  }
  if(hide_x_axis){
    p = p + theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
  }
  if(hide_x_axis_label){
    p = p + theme(axis.title.x = element_blank())
  }
  if(!is.null(plot_title)){
    p = p + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}

# Older plotting functions ------------------------------------------------------------
make_data_fig_D = function(eval_df_testing,model_param_comb = NULL,
                           metric_type = c("MAP", "nDCG"), add_err_bar = T){
  if (!is.null(model_param_comb)){
    eval_df_testing = eval_df_testing %>% dplyr::filter(model_params == model_param_comb)
  }
  if (any(str_detect(colnames(eval_df_testing), "AP"))){
    colnames(eval_df_testing)[str_which(colnames(eval_df_testing), "AP")] = "metric_value"
  }
  plot_data = eval_df_testing
  colnames(plot_data)[which(colnames(plot_data) == "metric_value")] = "AP"
  plot_data = plot_data[!is.nan(plot_data$AP),]
  all_db = plot_data %>%
    group_by(score_type) %>%
    dplyr::summarise(metric_value = mean(AP),
                     metric_val_sd = sd(AP)) %>%
    dplyr::mutate(db = "All\ndatabases")
  per_db = plot_data %>%
    group_by(score_type, db) %>%
    dplyr::summarise(metric_value = mean(AP),
                     metric_val_sd = sd(AP))
  plot_data = rbind.data.frame(all_db, per_db)
  plot_data$db = sub("-.*", "", plot_data$db)
  p = plot_eval_metrics_testing(df = plot_data, add_err_bar = add_err_bar)
  return(p)
} #12x10 portrait

make_data_fig_E = function(annot_df, FDR_pct = 10){
  plot_data = annot_df
  plot_data = plot_data %>%
    dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots, db) %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",") %>%
    dplyr::group_by(ds_id, FDR_pct,db) %>%
    dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                     pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
    as.data.frame()
  plot_data$Difference = plot_data$pred_fdr_per_ds - plot_data$msm_fdr_per_ds
  plot_data$diff_sign = ifelse(plot_data$Difference < 0,-1,1)
  plot_data$Difference = plot_data$diff_sign * log10(abs(plot_data$Difference) + 1)
  plot_data = plot_data[which(plot_data$FDR_pct == as.character(FDR_pct)),]
  x = plot_data %>%
    group_by(ds_id) %>%
    dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_per_ds),
                     pred_fdr_per_ds = sum(pred_fdr_per_ds)) %>%
    mutate(db = "All\ndatabases", Difference = (pred_fdr_per_ds - msm_fdr_per_ds),
           FDR_pct = as.character(FDR_pct))
  x$diff_sign = ifelse(x$Difference < 0,-1,1)
  x$Difference = x$diff_sign * log10(abs(x$Difference) + 1)

  x = x[,colnames(plot_data)] %>% as.data.frame()
  x = x[!duplicated(x),]
  plot_data = rbind.data.frame(x, plot_data)
  plot_data$db = sub("-.*", "", plot_data$db)

  y_break_range = c(floor(min(plot_data$Difference)),
                    ceiling(max(plot_data$Difference)))

  counts = plot_data %>%
    dplyr::group_by(db) %>%
    dplyr::summarise(n=n())
  n_labels <- paste0(counts[["db"]], "\n(n=", counts$n, ")")
  names(n_labels) <- counts[["db"]]

  p = ggboxplot(plot_data, x = "db", y = "Difference", add = "jitter",
                fill = "db", ggtheme = theme_pubr(),
                add.params = list(alpha = 0.4, color = 'black', size = 1),
                width = 0.9, palette = RColorBrewer::brewer.pal(length(unique(plot_data$db)),
                                                                "Set2")) +
    scale_y_continuous(breaks = seq(y_break_range[1], y_break_range[2], 1)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_blank())+
    stat_summary(fun=mean, geom="errorbar", aes(ymax = ..y..,
                                                ymin = ..y..), size = 1,
                 linetype = "dashed") +
    ylab("(+/-) Log 10 (Absolute Difference)") +
    xlab("") +
    scale_x_discrete(labels = n_labels)
  return(p)

} #12x10 portrait
make_data_fig_F = function(annot_df, y_axis_score = c("LFC", "LogDiff"), per_db = T){
  plot_data = annot_df
  plot_data = plot_data %>%
    dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots, db) %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",")
  if (per_db){
    plot_data = plot_data %>%
      dplyr::group_by(ds_id, FDR_pct,db) %>%
      dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                       pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
      as.data.frame()
    plot_data$db = sub("-.*", "", plot_data$db)
  }
  else{
    plot_data = plot_data %>%
      dplyr::group_by(ds_id, FDR_pct) %>%
      dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                       pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
      as.data.frame()
  }
  plot_data$Difference = plot_data$pred_fdr_per_ds - plot_data$msm_fdr_per_ds
  plot_data$diff_sign = ifelse(plot_data$Difference < 0,-1,1)
  plot_data$LogDiff = plot_data$diff_sign * log10(abs(plot_data$Difference) + 1)

  plot_data$FC = (plot_data$pred_fdr_per_ds + 1) / (plot_data$msm_fdr_per_ds + 1)
  plot_data$LFC = log2(plot_data$FC)

  y_break_range = c(floor(min(plot_data[,y_axis_score])),
                    ceiling(max(plot_data[,y_axis_score])))

  counts = plot_data %>%
    dplyr::group_by(db) %>%
    dplyr::summarise(n=n())
  n_labels <- paste0(counts[["db"]], "\n(n=", counts$n, ")")
  names(n_labels) <- counts[["db"]]

  plot_data[["db"]] = n_labels

  ggboxplot(plot_data, x = "FDR_pct", y = y_axis_score, add = ifelse(per_db, "none", "jitter"),
            color = ifelse(per_db, "db", "FDR_pct"), ggtheme = theme_pubr(),size = 1,
            add.params = ifelse(per_db, list(),
                                list(alpha = 0.5, color = 'black', size = 1.25)),
            width = 0.8, palette = RColorBrewer::brewer.pal(4, "Set2")) +
    scale_y_continuous(breaks = seq(y_break_range[1], y_break_range[2], 1)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.title = element_blank()) +
    ylab(ifelse(y_axis_score == "LFC",
                "Log2 Fold Change",
                "(+/-) Log 10 (Absolute Difference)"))

} #10x8 portrait

make_data_fig_G = function(raw_res, ds_id = "2021-03-22_15h30m15s", color_by_each_feat = F){
  umap_data = prepare_data_umap(raw_res = raw_res, ds_id =ds_id )

  umap_config = umap::umap.defaults
  umap_config$n_neighbors = 20

  x = modified_M3C_umap(mydata = umap_data$feature_mat, labels = as.factor(umap_data$labels), dotsize = 2, scale = 3,
                        controlscale = T, colvec = c("#FBB4AE", "#33A02C"), legendtitle = "Ion type",
                        umap_config = umap_config)

  y = modified_M3C_umap(mydata = umap_data$feature_mat, labels=as.numeric(umap_data$score_mat[row.names(umap_data$score_mat)=='MSM',]),
                        controlscale = TRUE,scale=1,legendtitle = 'MSM', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config) +
    xlab("UMAP 1") +
    ylab("UMAP 2")

  z = modified_M3C_umap(mydata = umap_data$feature_mat, labels=as.numeric(umap_data$score_mat[row.names(umap_data$score_mat)=='CatBoost',]),
                        controlscale = TRUE,scale=1,legendtitle = 'METASPACE-ML Score', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
    xlab("UMAP 1") +
    ylab("UMAP 2")

  if (!color_by_each_feat){
    return(list("Global_UMAPS" = list(x,y,z)))
  }
  else{
    a = modified_M3C_umap(mydata = umap_data$feature_mat, labels=as.numeric(umap_data$feature_mat[row.names(umap_data$feature_mat)=='spectral',]),
                          controlscale = TRUE,scale=1,legendtitle = 'spectral', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
      xlab("UMAP 1") +
      ylab("UMAP 2")
    b = modified_M3C_umap(mydata = umap_data$feature_mat, labels=as.numeric(umap_data$feature_mat[row.names(umap_data$feature_mat)=='chaos',]),
                          controlscale = TRUE,scale=1,legendtitle = 'chaos', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
      xlab("UMAP 1") +
      ylab("UMAP 2")
    c = modified_M3C_umap(mydata = umap_data$feature_mat, labels=as.numeric(umap_data$feature_mat[row.names(umap_data$feature_mat)=='spatial',]),
                          controlscale = TRUE,scale=1,legendtitle = 'spatial', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
      xlab("UMAP 1") +
      ylab("UMAP 2")
    d = modified_M3C_umap(mydata = umap_data$feature_mat, labels= -log10(abs(as.numeric(umap_data$feature_mat[row.names(umap_data$feature_mat)=='mz_err_abs',]))),
                          controlscale = TRUE,scale=1,legendtitle = '-Log10\nmz error rel', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
      xlab("UMAP 1") +
      ylab("UMAP 2")
    e = modified_M3C_umap(mydata = umap_data$feature_mat, labels= -log10(abs(as.numeric(umap_data$feature_mat[row.names(umap_data$feature_mat)=='mz_err_rel',]))),
                          controlscale = TRUE,scale=1,legendtitle = '-Log10\nmz error rel', low = "#450958", high = "#FAE726", dotsize = 2,umap_config = umap_config)+
      xlab("UMAP 1") +
      ylab("UMAP 2")
    return(list("Global_UMAPS" = list(x,y,z),
                "Feature_UMAPS" = list(a,b,c,d,e)))
  }

  #umap_compare_scores = cowplot::plot_grid(x,y,z)
  return(umap_compare_scores)
} #10x12 Landscape

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

  counts = plot_data %>%
    dplyr::group_by(db) %>%
    dplyr::summarise(n=n())
  n_labels <- paste0(counts[["db"]], "\n(n=", counts$n, ")")
  names(n_labels) <- counts[["db"]]

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
    scale_x_discrete(labels = n_labels) +
    xlab("") +
    ylab("MAP")
  return(p)
}

plot_enrich_res_global = function(enrich_result,use_FE = T, min_dss_prop = 0.1){
  plot_data = enrich_result %>% dplyr::filter(pval < 0.05, TP >= 3)
  n_sig_dss = length(unique(plot_data$ds_id))
  if (use_FE){
    plot_data$OR = log2(plot_data$FE)
  }
  else{
    plot_data$OR = log2(plot_data$OR)
  }

  filtered = plot_data %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(n_ds = n(), avg_LOR = mean(OR))
  filtered = filtered[which(filtered$n_ds >= ceiling(min_dss_prop * n_sig_dss)),]
  filtered = filtered[which(filtered$term != ""),]
  #filtered = filtered[which(abs(filtered$avg_LOR) >= 1),]

  x = HMDB_taxo_info[which(HMDB_taxo_info$sub_class %in% filtered$term),]
  x = x %>% dplyr::select(class, sub_class)
  x = x[!duplicated(x),]
  filtered = filtered %>% left_join(x, by = c("term" = "sub_class"))

  plot_data = plot_data[which(plot_data$term %in% filtered$term),]
  plot_data = plot_data %>%
    left_join(filtered) %>%
    select(term, class, ds_id, n_ds, OR)

  plot_data = plot_data[order(plot_data$class),]

  cols = RColorBrewer::brewer.pal(length(unique(plot_data$class)),"Paired")

  plot_data$myaxis = paste0(plot_data$term, " (", plot_data$n_ds, "/",
                            n_sig_dss,")")
  group_ordered = with(plot_data, reorder(myaxis, OR, median))

  plot_data$myaxis = factor(plot_data$myaxis, levels = levels(group_ordered))

  p = plot_data %>%
    ggplot(aes(x=myaxis, y=OR, fill=class)) +
    # geom_violin(width=1.4) +
    geom_boxplot(width= 1, alpha=0.6) +
    # geom_jitter(size= 0.5, alpha=0.5) +
    scale_fill_manual(values= cols) +
    scale_y_continuous(n.breaks = 16) +
    theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size=18, hjust = 0.5),
          plot.subtitle = element_text(size=16, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=18),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
    coord_flip() +
    geom_hline(yintercept =  0, linetype="dashed", color = "red") +
    ggtitle("Enrichment of Metabolite classes",
            subtitle = "Annotations only captured by METASPACE-ML") +
    ylab("Log2 Fold Enrichment")

  return(p)
}
