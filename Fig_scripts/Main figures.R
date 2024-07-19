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

bulk_data = read.delim("../Data/Bulk_validation/RawPeakIntensity_curated.txt",
                       skip = 4, header = T) %>%
  dplyr::filter(Adduct.type %in% c("[M+H]+", "[M+Na]+"))

eval_ions = readRDS("../Data/Bulk_validation/all_FDR_annots_new_model_testing_w_ions.rds")

eval_df_bulk = readRDS("../Data/Bulk_validation/all_eval_new_model_testing.rds")
annot_df_bulk = readRDS("../Data/Bulk_validation/all_FDR_annots_new_model_testing.rds")

#TODO Add link from Biostudies to metrics file
tmp = tempfile()
download.file(url = "link_to_metrics_file", destfile = tmp)
raw_res_bulk = read.csv(tmp)

#TODO Add link from Biostudies to feat_imp file
tmp = tempfile()
download.file(url = "link_to_feat_imp", destfile = tmp)
animal_feat_imp_dss = read.csv(tmp)

chosen_c_size = 30


# Prepare context data ----------------------------------------------------
AP_box_animal = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",context_size = chosen_c_size,
                                kingdom = "Animal",metrics_per_grp = F,
                                mean_eval_context = F)
Annot_MAP_animal = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                   annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                   data_type = "Testing",context_size = chosen_c_size,
                                   kingdom = "Animal",metrics_per_grp = F,
                                   mean_eval_context = F,FDR.pct = 10)

Annot_MAP_animal_nofdr = prepare_plot_df(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                         annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                         data_type = "Testing",context_size = chosen_c_size,
                                         kingdom = "Animal",metrics_per_grp = F,
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

# Figure 2 ----------------------------------------------------------------

Fig_2A = plot_context_dataexplorer(meta_df = all_public_ds,
                                   kingdom = "Animal")

Fig_2B_data = all_public_ds %>%
  dplyr::filter(Kingdom %nin% c("N/A","Protista")) %>%
  dplyr::select(ds_id, Kingdom) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Kingdom) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(fraction = n / sum(n),
                ymax = cumsum(fraction),
                ymin = c(0, head(ymax, n=-1)),
                labelPosition = (ymax + ymin) / 2,
                label =  paste0(Kingdom, "\n", round(fraction * 100, 1), "%"))

Fig_2B = ggplot(Fig_2B_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Kingdom)) +
  geom_rect() +
  geom_label_repel(x=3.5, aes(y=labelPosition, label=label), size=4,
                   max.overlaps = 15,box.padding = 0.1,show.legend = F) +
  scale_fill_brewer(palette="Set2") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("Breakdown of METASPACE datasets by kingdom", subtitle = "Total = 7713") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

Fig_2C = plot_new_sankey(
  meta_df = all_public_ds,
  kingdom = "Animal",
  ds_context_list = training_dss$Animal[[as.character(chosen_c_size)]],
  sel_dss = unique(unlist(training_dss$Animal[[as.character(chosen_c_size)]])),
  plot_title = paste0("Training", "\n", "Animal")
)

Fig_2D = plot_new_sankey(
  meta_df = all_public_ds,
  kingdom = "Animal",
  ds_context_list = testing_dss$Animal,
  sel_dss = unique(unlist(testing_dss$Animal)),
  plot_title = paste0("Testing", "\n", "Animal")
)


# Figure 3 ----------------------------------------------------------------
Fig_3A = Performance_plots_pipeline(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Animal",
                                mean_eval_context = F,
                                metrics_per_grp = F,
                                FDR.pct = NULL,
                                plot_type = "paired_box",
                                x_column = "score_type",
                                y_column = "metric_value",
                                color_column = "score_type",
                                hide_x_axis = T,
                                hide_x_axis_label = T,
                                context_specific = F,
                                in_plot_text_size = 3,
                                plot_title = paste0("Testing", "\n", "Animal"))

Fig_3B = make_density_hm_context(df = AP_box_animal, FDR = 10,
                                 vals = "MAP", kingdom = "Animal",
                                 plot_title = paste0("Animal", "\n", "Testing"))

Fig3C = make_data_fig_D(eval_df_testing = db_comp_new_res$eval,
                        metric_type = "MAP",add_err_bar = F)

Fig3D_data = db_comp_new_res$eval %>%
  dplyr::filter(score_type == "METASPACE-ML")

x = Fig3D_data %>%
  dplyr::group_by(group_name) %>%
  dplyr::summarise(n = n())

Fig3D_data = Fig3D_data[Fig3D_data$group_name %in% x$group_name[x$n  == 3],]

Fig3D = ggwithinstats(data = Fig3D_data, x = db, y = metric_value,
                  type = "nonparametric", p.adjust.method = "fdr",
                  pairwise.display = "s", results.subtitle = F,
                  ylab = "MAP",xlab = "") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14))




# Figure 4 ----------------------------------------------------------------

Fig_4A = Performance_plots_pipeline(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
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
                                hide_x_axis = F,
                                hide_x_axis_label = F,
                                context_specific = F,
                                in_plot_text_size = 5,
                                plot_title = paste0("Testing", "\n", "Animal"))

Fig_4B = Performance_plots_pipeline(eval_df = all_eval_res$Plant[[as.character(chosen_c_size)]]$Testing$eval_df,
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

Fig4C = make_density_hm_context(df = Annot_MAP_animal,FDR = 10,vals = "LogDiff",
                            kingdom = "Animal", plot_title = NULL)

Fig_4D_data = Annot_MAP_animal[Annot_MAP_animal$score_type == "METASPACE-ML",]

Fig_4D = ggstats_wrapper(df = Fig_4D_data, x_column = "analyzer", paired = F,
                    y_column = "LogDiff", centrality.plotting = T,
                    results.subtitle = T, parametric = F) +
  ylab("(+-)Log10(Absolute Difference)") +
  xlab("Analyzer") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

Fig_4E = make_data_fig_E(annot_df = db_comp_new_res$annot, FDR_pct = 10)


# Figure 5 ----------------------------------------------------------------

Fig_5A = Performance_plots_pipeline(eval_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$eval_df,
                                annot_df = all_eval_res$Animal[[as.character(chosen_c_size)]]$Testing$annot_df,
                                data_type = "Testing",
                                context_size = chosen_c_size,
                                kingdom = "Animal",
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

Fig5B_data = Fig_5A$data
Fig5B_data$Annot_diff = ifelse(Fig5B_data$LogDiff > 0, "Higher", "Lower")
Fig5B_data$Annot_diff[Fig5B_data$LogDiff == 0] = "Equal"

Fig_5B = ggstats_wrapper(df = Fig5B_data, x_column = "Annot_diff",
                         y_column = "metric_value",paired = F,
                         centrality.plotting = T,
                         results.subtitle = F, parametric = F) +
  ylab("MAP") +
  xlab("") +
  ggtitle("Animal") +
  theme(plot.title = element_text(hjust = 0.5))


Fig_5C_data = merge_rel_annot_FDR(max_rel_df = rel_scores,
                            clean_annot_df = Annot_MAP_animal_nofdr,
                            kingdom = "Animal")

Fig_5C_data$FDR_pct = factor(Fig_5C_data$FDR_pct,
                             levels = c("Optim FDR","5", "10", "20", "50"))

Fig_5C = ggstats_wrapper(df = Fig_5C_data, x_column = "FDR_pct",
                     y_column = "LFC",paired = F,
                     centrality.plotting = T,
                     results.subtitle = F, parametric = F,
                     sig_only = T) +
  ggtitle("Animal") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

Fig_5D = ggstats_wrapper(df = rel_scores$Animal, x_column = "FDR_pct",
                    y_column = "rel_score",paired = F,
                    centrality.plotting = F,
                    results.subtitle = F, parametric = F) +
  ylab("Reliability Score") +
  xlab("FDR (%)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14)) +
  ylim(c(0,1))

# Figure 6 ----------------------------------------------------------------

animal_feat_imp_dss = animal_feat_imp
animal_feat_imp_dss$abs_pred = rowSums(abs(animal_feat_imp_dss[,c("chaos", "spatial", "spectral",
                                                                  "mz_err_abs_abserr", "mz_err_rel_abserr")]))
animal_feat_imp_dss = animal_feat_imp_dss %>%
  mutate(chaos = abs(chaos)/abs_pred,
         spatial = abs(spatial)/abs_pred,
         spectral = abs(spectral)/abs_pred,
         mz_err_abs_abserr = abs(mz_err_abs_abserr)/abs_pred,
         mz_err_rel_abserr = abs(mz_err_rel_abserr)/abs_pred)

animal_feat_imp_dss = animal_feat_imp_dss %>%
  dplyr::group_by(context, ds_id) %>%
  dplyr::summarise(chaos = mean(chaos, na.rm = T),
                   spatial = mean(spatial, na.rm = T),
                   spectral = mean(spectral, na.rm = T),
                   mz_err_abs_abserr = mean(mz_err_abs_abserr, na.rm = T),
                   mz_err_rel_abserr = mean(mz_err_rel_abserr, na.rm = T)) %>%
  gather(key = "feature", value = "Shap_contrib", -ds_id, -context)


animal_feat_imp_dss$feature = gsub("_abserr", "", animal_feat_imp_dss$feature)
animal_feat_imp_dss$feature = gsub("_", " ", animal_feat_imp_dss$feature)


Fig_6A = ggplot(animal_feat_imp_dss, aes(x = Shap_contrib, y = feature,
                                    fill = stat(quantile))) +
  stat_density_ridges(scale = 1, rel_min_height = 0.01, calc_ecdf = T,
                      geom = "density_ridges_gradient", quantile_lines = T) +
  scale_x_continuous(n.breaks = 10) +
  scale_fill_viridis_d(name = "Quartiles") +
  theme_pubr() +
  ylab("") +
  xlab("SHAP impact contribution") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

Fig_6B = make_density_hm_feat_imp(df = animal_feat_imp, kingdom = "Animal",
                              shap_abs_prop = T,col_split = 2)

Fig_6C_data = make_data_fig_G(raw_res = raw_res_bulk, ds_id = "2016-09-22_11h16m09s",
                           color_by_each_feat = F)

for (i in 1:length(Fig_6C_data$Global_UMAPS)) {
  Fig_6C_data$Global_UMAPS[[i]] = Fig_6C_data$Global_UMAPS[[i]] +
    ggpubr::theme_pubr()
}

Fig_6C = cowplot::plot_grid(plotlist = Fig_6C_data$Global_UMAPS)

score_colors = c("#7AC1A5", "#EE8F68", "#92A0C9",
                 "#ADD65B", "#4C7EB6", "#9C572F",
                 "#D3282C")
all_PR_ROC = ROC_PR_data$AUC_results
ROC_label_data = ROC_PR_data$ROC_label_data

Fig_6D = ggplot(all_PR_ROC$ROC, aes(x = FPR, y = Sensitivity,
                                color = score_type)) +
  geom_line(lwd = 0.8,stat = "identity") +
  theme_pubr() +
  xlab("FPR") +
  ylab("Sensitivity") +
  geom_text(data = ROC_label_data, aes(x = pos_x, y = pos_y, label = label,
                                       color = score_type),size = 4,
            show.legend = F) +
  scale_color_manual("Score type",
                     values = score_colors) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))



# Figure 7 ----------------------------------------------------------------
intens_med_fdr10_animal = intens_data$Animal

Fig_7A = plot_intensity_boxplot(intens_med_fdr10_animal)
Fig_7B = plot_intens_hm_context(intens_med_fdr10_animal,
                                "Animal",scale_intens = T)

Fig_7C = plot_enrich_res_global(enrich_result = enrich_res$Animal)

Fig_7D = plot_enrich_hm_context(enrich_res = enrich_res$Animal, filter_by_adjpval = T,
                            min_TP = 3, use_FE = T, kingdom = "Animal")



# Figure 8 ----------------------------------------------------------------

Fig_8A = PR_Data %>% ggpaired(cond1 = "METASPACE-ML", cond2 = "MSM",
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
  stat_compare_means(aes(group = condition), label = "p.signif" ,paired = F,
                     size = 8, label.x.npc = "center", label.y.npc = 0.95) +
  theme(legend.position = "none")

Fig_8B = get_ggmat_genbase(df = perf_bulk_data, plot_type = "box", x_column = "FDR_pct",
                      y_column = "LogDiff", color_column = "FDR_pct",
                      in_plot_text_size = 5) +
  theme_pubr() +
  geom_point(size = 3, alpha = 0.5) +
  ylab("(+-)Log10(Absolute Difference)") +
  xlab("FDR (%)") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))


zoom_ds_data = bulk_conting %>%
  dplyr::filter(ds_id == "2016−09−22_11h16m17s")
zoom_ds_data$approach = ifelse(zoom_ds_data$incl_label == "Exclusive",
                               paste(zoom_ds_data$approach, "only", sep = "_"),
                               zoom_ds_data$approach)

Fig_8C_data = zoom_ds_data[zoom_ds_data$approach == "METASPACE-ML",c(3:6)] %>%
  as.numeric() %>%
  matrix(nrow = 2, ncol = 2, byrow = F)

TPR = round((Fig_8C_data[1,1] / sum(Fig_8C_data[,1])),2) * 100

dimnames( Fig_8C_data ) <- list(`METASPACE-ML` = c("Predicted", "Not predicted"),
                          `LC-MS/MS` = c("Found", "Not found"))
fourfoldplot(Fig_8C_data,
             color = c("#B22222", "#2E8B57"),
             main = paste0("TPR : ",TPR ,"%"),
             conf.level = 0, std = "all.max",space = 0.2) +
  text(-0.2,0.2, "TP", cex=1.5) +
  text(0.2, -0.2, "TN", cex=1.5) +
  text(0.2,0.2, "FP", cex=1.5) +
  text(-0.2, -0.2, "FN", cex=1.5)
Fig_8C = recordPlot()


Fig_8D_data = zoom_ds_data[zoom_ds_data$approach == "METASPACE-ML_only",c(3:6)] %>%
  as.numeric() %>%
  matrix(nrow = 2, ncol = 2, byrow = F)

TPR = round((Fig_8D_data[1,1] / sum(Fig_8D_data[,1])),2) * 100

dimnames( Fig_8D_data ) <- list(`METASPACE-ML_only` = c("Predicted", "Not predicted"),
                          `LC-MS/MS` = c("Found", "Not found"))
fourfoldplot(Fig_8D_data,
             color = c("#B22222", "#2E8B57"),
             main = paste0("TPR : ",TPR ,"%"),
             conf.level = 0, std = "all.max",space = 0.2) +
  text(-0.2,0.2, "TP", cex=1.5) +
  text(0.2, -0.2, "TN", cex=1.5) +
  text(0.2,0.2, "FP", cex=1.5) +
  text(-0.2, -0.2, "FN", cex=1.5)
Fig_8D = recordPlot()

Fig_8E_data = zoom_ds_data[zoom_ds_data$approach == "MSM_only",c(3:6)] %>%
  as.numeric() %>%
  matrix(nrow = 2, ncol = 2, byrow = F)

TPR = round((Fig_8E_data[1,1] / sum(Fig_8E_data[,1])),2) * 100

dimnames( Fig_8E_data ) <- list(`MSM_only` = c("Predicted", "Not predicted"),
                          `LC-MS/MS` = c("Found", "Not found"))
fourfoldplot(Fig_8E_data,
             color = c("#B22222", "#2E8B57"),
             main = paste0("TPR : ",TPR ,"%"),
             conf.level = 0, std = "all.max",space = 0.2) +
  text(-0.2,0.2, "TP", cex=1.5) +
  text(0.2, -0.2, "TN", cex=1.5) +
  text(0.2,0.2, "FP", cex=1.5) +
  text(-0.2, -0.2, "FN", cex=1.5)
Fig_8E = recordPlot()
