library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
library(ggsankey)
library(ggrepel)
library(ggstatsplot)
library(RColorBrewer)
library(DataExplorer)
library(cowplot)
library(cutpointr)
library(GGally)
library(rstatix)
library(ggbiplot)

source("../Scripts/Final_figures_functions.R")
"%nin%" <- Negate("%in%")


# Load data ---------------------------------------------------------------
all_public_ds = readRDS("../Data/Datasets_metadata.rds")
training_dss = readRDS("../Data/Training_datasets.rds")
testing_dss = readRDS("../Data/Testing_datasets.rds")
all_eval_res = readRDS("../Data/Evaluation_results.rds")
db_comp_new_res = readRDS("../Data/DB_comparison_eval.rds")
rel_scores = readRDS("../Data/Reliability_scores.rds")
ROC_PR_data = readRDS("../Data/ROC_PR_data.rds")
intens_data = readRDS("../Data/Intensity_plots_data.rds")
enrich_res = readRDS("../Data/context_enrich_results.rds")
HMDB_taxo_info = readRDS("../Data/HMDB_taxo_info.rds")
all_iter_res = readRDS("../Data/Learning_curves_data.rds")
decoy_diff_means = readRDS("../Data/decoy_mass_comparison_data.rds")

bulk_data = read.delim("../Data/Bulk_validation/RawPeakIntensity_curated.txt",
                       skip = 4, header = T) %>%
  dplyr::filter(Adduct.type %in% c("[M+H]+", "[M+Na]+"))

eval_ions = readRDS("../Data/Bulk_validation/all_FDR_annots_new_model_testing_w_ions.rds")

eval_df_bulk = readRDS("../Data/Bulk_validation/all_eval_new_model_testing.rds")
annot_df_bulk = readRDS("../Data/Bulk_validation/all_FDR_annots_new_model_testing.rds")


raw_res_bulk = read.csv("../Data/Bulk_validation/CoreMetabolome-v3_metrics_df_w_pred_score.csv.gz")


tmp = tempfile()
download.file(url = "https://www.ebi.ac.uk/biostudies/files/S-BIAD1283/Biostudies_files/Plant_Testing_Eval_feature_importance_per_context.csv", destfile = tmp)
plant_feat_imp = read.csv(tmp)

chosen_c_size = 30


# Prepare context data ----------------------------------------------------
MAP_bar_animal = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                 annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                 data_type = "Testing",context_size = chosen_c_size,
                                 kingdom = "Animal",metrics_per_grp = F,
                                 mean_eval_context = T)
MAP_bar_plant = prepare_plot_df(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",context_size = chosen_c_size,
                                kingdom = "Plant",metrics_per_grp = F,
                                mean_eval_context = T)

AP_box_animal = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",context_size = chosen_c_size,
                                kingdom = "Animal",metrics_per_grp = F,
                                mean_eval_context = F)

AP_box_plant = prepare_plot_df(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                               annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                               data_type = "Testing",context_size = chosen_c_size,
                               kingdom = "Plant",metrics_per_grp = F,
                               mean_eval_context = F)

Annot_MAP_animal = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                   annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                   data_type = "Testing",context_size = chosen_c_size,
                                   kingdom = "Animal",metrics_per_grp = F,
                                   mean_eval_context = F,FDR.pct = 10)

Annot_MAP_plant = prepare_plot_df(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                  annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                  data_type = "Testing",context_size = chosen_c_size,
                                  kingdom = "Plant",metrics_per_grp = F,
                                  mean_eval_context = F,FDR.pct = 10)

Annot_MAP_animal_nofdr = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                         annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                         data_type = "Testing",context_size = chosen_c_size,
                                         kingdom = "Animal",metrics_per_grp = F,
                                         mean_eval_context = F)

Annot_MAP_plant_nofdr = prepare_plot_df(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                        annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                        data_type = "Testing",context_size = chosen_c_size,
                                        kingdom = "Plant",metrics_per_grp = F,
                                        mean_eval_context = F)


# Prepare bulk validation data --------------------------------------------
perf_bulk_data = prepare_bulk_plot_df(eval_df = eval_df_bulk,
                                      annot_df = annot_df_bulk,
                                      mean_eval_context = F,
                                      metrics_per_grp = F,
                                      FDR.pct = NULL)

x_column = "score_type"
y_column = "metric_value"
color_column = "score_type"

PR_Data = perf_bulk_data %>%
  dplyr::select(ds_id,{{x_column}},{{y_column}}) %>%
  distinct() %>%
  spread(key = .data[[x_column]], value = .data[[y_column]])

PR_Data$Delta = ifelse(PR_Data$`METASPACE-ML` > PR_Data$MSM, "Higher", "Lower")
PR_Data = PR_Data[!is.na(PR_Data$Delta),]

Delta_count = PR_Data %>%
  dplyr::group_by(Delta) %>%
  dplyr::summarise(n_ds = n_distinct(ds_id))

PR_Data = PR_Data %>%
  left_join(Delta_count, by = "Delta") %>%
  mutate(Delta = paste(Delta, "\n(", n_ds, ")", sep = ""))

bulk_dss = unique(eval_ions$ds_id)
bulk_conting = list()
for (i in c(T,F)){
  for (s in c("METASPACE-ML", "MSM")){
    for (ds in bulk_dss){
      res = preprocess_bulk_conting(eval_ions = eval_ions,
                                    bulk_data = bulk_data,
                                    fdr = 0.1, exclusive = i,
                                    approach = s,
                                    ds_id = ds)
      res$ds_id = ds
      bulk_conting[[paste(i,s,ds,sep = "_")]] = res
    }
  }
}
bulk_conting = bulk_conting %>% dplyr::bind_rows()
bulk_conting$TPR = cutpointr::tpr(bulk_conting$TP, bulk_conting$FN) %>% as.numeric()
bulk_conting$FPR = cutpointr::fpr(bulk_conting$FP, bulk_conting$TN) %>% as.numeric()
bulk_conting$FNR = cutpointr::fnr(bulk_conting$TP, bulk_conting$FN) %>% as.numeric()


# Supplementary Figures ---------------------------------------------------

S1 = plot_context_dataexplorer(meta_df = all_public_ds, kingdom = "Plant")

S2_data = list()
counter = 0
for (k in c("Animal", "Plant")){
  for (c in seq(10,50,10)){
    counter = counter + 1
    S2_data[[counter]] = data.frame(kingdom = k,
                                                context_size = c,
                                                n_contexts = length(training_dss[[k]][[as.character(c)]]))
  }
}
S2_data = dplyr::bind_rows(S2_data)
S2_data$context_size = as.character(S2_data$context_size)

S2 = S2_data %>%
  ggplot(aes(x = context_size, y = n_contexts, group = kingdom,fill = kingdom)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = .data[["n_contexts"]]),
            position = position_dodge2(width = 0.9),
            vjust = 1, size = 4) +
  scale_fill_manual(values = c("#EC8377", "#BEBBD9")) +
  theme_pubr() +
  ylab("Number of contexts") +
  xlab("Number of datasets per context")

S3 = ggplot(data = all_iter_res, aes(x = iter, group = error_type)) +
  geom_line(aes(y = mean_error, color = error_type), size = 1) +
  geom_ribbon(aes(y = mean_error, ymin = mean_error - sd_error,
                  ymax = mean_error + sd_error, fill = error_type), alpha = .2) +
  xlab("Iteration") +
  ylab("PairLogit") +
  theme_pubr() +
  theme(legend.key = element_blank()) +
  theme(legend.title = element_blank()) +
  scale_x_continuous(n.breaks = 5) +
  facet_grid(kingdom ~ context_size)

S4 = plot_new_sankey(
  meta_df = all_public_ds,
  kingdom = "Plant",
  ds_context_list = training_dss$Plant[[as.character(chosen_c_size)]],
  sel_dss = unique(unlist(training_dss$Plant[[as.character(chosen_c_size)]])),
  plot_title = paste0("Training", "\n", "Plant")
)

S5 = plot_new_sankey(
  meta_df = all_public_ds,
  kingdom = "Plant",
  ds_context_list = testing_dss$Plant,
  sel_dss = unique(unlist(testing_dss$Plant)),
  plot_title = paste0("Testing", "\n", "Plant")
)

S6A = Performance_plots_pipeline(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Training$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Training$annot_df,
                                data_type = "Training",
                                context_size = chosen_c_size,
                                kingdom = "Animal",
                                mean_eval_context = T,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "bar",
                                x_column = "score_type",
                                y_column = "metric_value",
                                color_column = "score_type",
                                hide_x_axis = T,hide_x_axis_label = T,
                                context_specific = F,
                                in_plot_text_size = 5,
                                plot_title = paste0("Cross validation", "\n", "Animal"))
S6B = Performance_plots_pipeline(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Training$eval_df,
                                annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Training$annot_df,
                                data_type = "Training",
                                context_size = chosen_c_size,
                                kingdom = "Plant",
                                mean_eval_context = T,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "bar",
                                x_column = "score_type",
                                y_column = "metric_value",
                                color_column = "score_type",
                                hide_x_axis = T,hide_x_axis_label = T,
                                context_specific = F,
                                in_plot_text_size = 5,
                                plot_title = paste0("Cross validation", "\n", "Plant"))
S6 = cowplot::plot_grid(S6A, S6B, nrow = 1, ncol = 2, labels = c("A", "B"))

S7 = Performance_plots_pipeline(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Plant",
                                mean_eval_context = F,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "paired_box",
                                x_column = "score_type",
                                y_column = "metric_value",
                                color_column = "score_type",
                                hide_x_axis = T,hide_x_axis_label = T,
                                context_specific = F,
                                in_plot_text_size = 3,
                                plot_title = paste0("Testing", "\n", "Plant"))

S8 = ggmat_wrapper(plot_df = MAP_bar_animal, plot_type = "bar",
                  x_column = "score_type", y_column = "metric_value",
                  color_column = "score_type",hide_x_axis_label = T,
                  hide_x_axis = T,in_plot_text_size = 2)

S9 = ggmat_wrapper(plot_df = MAP_bar_plant, plot_type = "bar",
                  x_column = "score_type", y_column = "metric_value",
                  color_column = "score_type",hide_x_axis_label = T,
                  hide_x_axis = T,in_plot_text_size = 2)

S10 = make_density_hm_context(df = AP_box_plant, FDR = 10,
                            vals = "MAP", kingdom = "Plant",
                            plot_title = paste0("Plant", "\n", "Testing"))

S11_data_1 = AP_box_animal %>% dplyr::select(metric_value, score_type, context, ds_id) %>%
  filter(score_type == "METASPACE-ML") %>%
  distinct()

S11_data_2 = rstatix::wilcox_test(metric_value ~ context,data = S11_data_1,
                                  p.adjust.method = "fdr", paired = F)

S11_data_2$log10_p = -log10(S11_data_2$p.adj)
S11 = plot_context_hm_adjacency(S11_data_2,kingdom = "Animal",
                              col_column = "group1", row_column = "group2",
                              vals_column = "log10_p",
                              label_column = "p.adj.signif")

S12_data_1 = AP_box_plant %>% dplyr::select(metric_value, score_type, context, ds_id) %>%
  filter(score_type == "METASPACE-ML") %>%
  distinct()

S12_data_2 = rstatix::wilcox_test(metric_value ~ context,data = S12_data_1,
                                  p.adjust.method = "fdr", paired = F)

S12_data_2$log10_p = -log10(S12_data_2$p.adj)
S12 = plot_context_hm_adjacency(S12_data_2,kingdom = "Plant",
                              col_column = "group1", row_column = "group2",
                              vals_column = "log10_p",
                              label_column = "p.adj.signif")

S13 = make_figure_S8(db_comp_new_res$eval)

S14_data = Performance_plots_pipeline(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Animal",
                                mean_eval_context = F,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "box",
                                x_column = "FDR_pct",
                                y_column = "LogDiff",
                                color_column = "FDR_pct",
                                hide_x_axis = F,hide_x_axis_label = F,
                                context_specific = F,
                                in_plot_text_size = 5,
                                plot_title = paste0("Testing", "\n", "Animal"))
S14 = ggstats_wrapper(df = S14_data$data, x_column = "FDR_pct", y_column = "LFC", paired = F,
                     centrality.plotting = T, results.subtitle = F, parametric = F) +
  ggtitle("Animal") + theme(plot.title = element_text(hjust = 0.5))


S15_data = Performance_plots_pipeline(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Plant",
                                mean_eval_context = F,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "box",
                                x_column = "FDR_pct",
                                y_column = "LogDiff",
                                color_column = "FDR_pct",
                                hide_x_axis = F,hide_x_axis_label = F,
                                context_specific = F,
                                in_plot_text_size = 5,
                                plot_title = paste0("Testing", "\n", "Plant"))
S15 = ggstats_wrapper(df = S15_data$data, x_column = "FDR_pct", y_column = "LFC", paired = F,
                      centrality.plotting = T, results.subtitle = F, parametric = F) +
  ggtitle("Plant") + theme(plot.title = element_text(hjust = 0.5))

S16 = make_density_hm_context(df = Annot_MAP_plant,FDR = 10,vals = "LogDiff",
                            kingdom = "Plant", plot_title = NULL)

S17 = ggmat_wrapper(plot_df = Annot_MAP_animal_nofdr, plot_type = "box",
                  x_column = "FDR_pct", y_column = "LFC",
                  color_column = "FDR_pct",hide_x_axis = F,
                  hide_x_axis_label = F)

S18 = ggmat_wrapper(plot_df = Annot_MAP_plant_nofdr, plot_type = "box",
                  x_column = "FDR_pct", y_column = "LFC",
                  color_column = "FDR_pct",hide_x_axis = F,
                  hide_x_axis_label = F)

S19 = make_data_fig_F(annot_df = db_comp_new_res$annot, y_axis_score = "LFC",
                    per_db = T)

S20 = Performance_plots_pipeline(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Plant",
                                mean_eval_context = F,
                                metrics_per_grp = F,
                                FDR.pct = 10,
                                plot_type = "2d_bin",
                                x_column = "LogDiff",
                                y_column = "metric_value",
                                color_column = "score_type",
                                hide_x_axis = F,hide_x_axis_label = F,
                                context_specific = F,
                                in_plot_text_size = 5)

S21_data = S20$data
S21_data$Annot_diff = ifelse(S21_data$LogDiff > 0, "Higher", "Lower")
S21_data$Annot_diff[S21_data$LogDiff == 0] = "Equal"

S21 = ggstats_wrapper(df = S21_data, x_column = "Annot_diff", y_column = "metric_value",paired = F,
                     centrality.plotting = T, results.subtitle = F, parametric = F) +
  ylab("MAP") +
  xlab("") +
  ggtitle("Plant") +
  theme(plot.title = element_text(hjust = 0.5))

S22_data = merge_rel_annot_FDR(max_rel_df = rel_scores,
                            clean_annot_df = Annot_MAP_plant_nofdr,
                            kingdom = "Plant")

S22_data$FDR_pct = factor(S22_data$FDR_pct, levels = c("Optim FDR","5", "10", "20", "50"))

S22 = ggstats_wrapper(df = S22_data, x_column = "FDR_pct",
                     y_column = "LFC",paired = F,
                     centrality.plotting = T,
                     results.subtitle = F, parametric = F,
                     sig_only = T) +
  ggtitle("Plant") +
  theme(plot.title = element_text(hjust = 0.5))

S23_data = merge_rel_annot_FDR(max_rel_df = rel_scores,
                            clean_annot_df = Annot_MAP_animal_nofdr,
                            kingdom = "Animal")

S23_data$FDR_pct = factor(S23_data$FDR_pct, levels = c("Optim FDR","5", "10", "20", "50"))

S23 = get_ggmat_genbase(df = S23_data[S23_data$FDR_pct == "Optim FDR",],
                            plot_type = "2d_bin",
                            x_column = "LogDiff",
                            y_column = "metric_value",
                            color_column = "score_type",
                            in_plot_text_size = 5) +
  theme_pubr() +
  ylab("MAP") +
  xlab("(+-)Log10(Absolute Difference)")

S24_data = merge_rel_annot_FDR(max_rel_df = rel_scores,
                            clean_annot_df = Annot_MAP_plant_nofdr,
                            kingdom = "Plant")

S24_data$FDR_pct = factor(S24_data$FDR_pct, levels = c("Optim FDR","5", "10", "20", "50"))


S24 = get_ggmat_genbase(df = S24_data[S24_data$FDR_pct == "Optim FDR",],
                            plot_type = "2d_bin",
                            x_column = "LogDiff",
                            y_column = "metric_value",
                            color_column = "score_type",
                            in_plot_text_size = 5) +
  theme_pubr() +
  ylab("MAP") +
  xlab("(+-)Log10(Absolute Difference)")

S25 = make_density_hm_feat_imp(df = plant_feat_imp, kingdom = "Plant",
                              shap_abs_prop = T)

S26_data_1 = PCA_per_ds(raw_res = raw_res_bulk,ds_id = "2016-09-22_11h16m09s",
                     fix_mz_err = F,scale_pca = T, center_pca = T)


S26_data_2 = S26_data_1$orig_data %>%
  as.data.frame() %>%
  mutate(groups = S26_data_1$decoy_target$data$groups)

S26 = ggpairs(S26_data_2, mapping = aes(color = groups,alpha = 0.5),
             columns = c("chaos", "spectral", "spatial", "mz_err_abs", "mz_err_rel")) +
  theme_pubr()


all_PR_ROC = ROC_PR_data$AUC_results
PR_label_data = ROC_PR_data$PR_label_data
PR_label_data$score_type = gsub("_abserr", "", PR_label_data$score_type)
PR_label_data$label = gsub("_abserr", "", PR_label_data$label)
score_colors = c("#7AC1A5", "#EE8F68", "#92A0C9",
                 "#ADD65B", "#4C7EB6", "#9C572F",
                 "#D3282C")

S27 = ggplot(all_PR_ROC$PR, aes(x = Recall, y = Precision, color = score_type)) +
  geom_line(lwd = 0.8,stat = "identity") +
  theme_pubr() +
  xlab("Recall") +
  ylab("Precision") +
  geom_text(data = PR_label_data, aes(x = pos_x, y = pos_y, label = label,
                                      color = score_type),size = 4,
            show.legend = F) +
  scale_color_manual("Score type",
                     values = score_colors) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

S28 = ggbetweenstats(data = decoy_diff_means, x = mean_diff_sign, y = mass,
                   type = "nonparametric", p.adjust.method = "BH",
                   pairwise.comparisons = F) +
  xlab("") +
  ylab("Mass (Da)")

S29_data = readRDS("../Data/Fig_S29_data.rds")
S29 = ggbetweenstats(data = S29_data, x =Group, y = total_intensity,
                   type = "nonparametric",p.adjust.method = "BH",
                   title = "2019-08-19_11h28m42s", xlab = "", ylab = "Log10(Intensity)",
                   results.subtitle = F) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
          plot.subtitle = element_text(size=12, hjust = 0.5),
          axis.title.y = element_text(size=18, hjust = 0.5),
          axis.title.x = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=16),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
  scale_color_manual(values = c("#A6CEE3","#A6CEE3","#FB9A99","#FB9A99"))

intens_med_fdr10_plant = intens_data$Plant
S30 = plot_intens_hm_context(intens_med_fdr10_plant, "Plant")

S31 = plot_enrich_res_global(enrich_result = enrich_res$Plant)
S32 = plot_enrich_hm_context(enrich_res = enrich_res$Plant, filter_by_adjpval = T,
                            min_TP = 3, use_FE = T, kingdom = "Plant")

S33_data = bulk_conting %>%
  dplyr::select(-TP, -FP, -FN, -TN) %>%
  tidyr::gather(key = "metric", value = "metric_value", -ds_id, -incl_label, -approach)

S33 = ggplot(S33_data, aes(x=approach,
                      y=ds_id,
                      fill=metric_value)) +
  geom_tile(color="white")+
  geom_text(aes(label = round(metric_value,2)),size = 5) +
  scale_fill_gradient2(low = "lightblue", high = "red",midpoint = 0.5)+
  # coord_equal()+
  facet_grid(incl_label~metric, scales = "free") +
  theme_pubr(legend = "right") +
  rotate_x_text() +
  xlab("") +
  ylab("")
