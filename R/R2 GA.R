# ==========================================
# LIBRARY PREPARATION
# ==========================================
required_packages <- c("dplyr", "tidyr", "car", "lme4", "lmerTest", "ggplot2", "effectsize", "sjPlot", "minqa", "mediation")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(tidyr)
library(car)      
library(lme4)     
library(lmerTest) 
library(ggplot2)  
library(effectsize) 
library(sjPlot)   
library(mediation) 

# ==========================================
# 1. READ & PREPARE ORIGINAL DATA
# ==========================================
data <- read.csv("OGAT.csv", sep = ";")

# [UPDATE]: Sesuai dataset terbaru (0 = Control, 1 = Intervention)
# Mengubahnya menjadi factor dengan label yang jelas, 
# 0 (Control) otomatis akan menjadi level pertama (Baseline mutlak)
data$GROUP <- factor(data$GROUP, levels = c(0, 1), labels = c("Control", "Intervention"))

data$TIME <- as.numeric(data$TIME) 
data$RESPONDENT_ID <- as.factor(data$RESPONDENT_ID)
data$Gender <- as.factor(data$Gender)

score_columns <- c(
  "SDQ_Emotion", "SDQ_Conduct", "SDQ_Hyperactivity", "SDQ_Peer", "SDQ_Protective", "Total_SDQ",
  "GO_Bijak", "Assertiveness",
  "GAS_Salience", "GAS_Mood", "GAS_Tolerance", "GAS_Conflict", "GAS_Relapse", "GAS_Withdrawal", "GAS_Problems", "Total_GAS",
  "BPAQ_Physical", "BPAQ_Verbal", "BPAQ_Anger", "BPAQ_Hostility", "Total_BPAQ"
)

# Calculate Overall Change Score (Post - Pre) for Assumptions and MWU
data_change <- data %>%
  arrange(RESPONDENT_ID, TIME) %>%
  group_by(RESPONDENT_ID, GROUP, Gender) %>%
  summarise(across(all_of(score_columns), ~ last(.) - first(.), .names = "change_{.col}"), .groups = 'drop')

change_score_cols <- paste0("change_", score_columns)

# ==========================================
# STAGE 1 & 2: NORMALITY AND HOMOGENEITY
# ==========================================
cat("\n[Info] Performing Assumption Checks...\n")
assumption_results <- list()

for (col in change_score_cols) {
  subset_data <- data_change %>% drop_na(all_of(c(col, "GROUP", "Gender")))
  
  data_interv <- subset_data[[col]][subset_data$GROUP == "Intervention"] 
  data_ctrl <- subset_data[[col]][subset_data$GROUP == "Control"]
  
  p_norm_interv <- ifelse(length(data_interv) >= 3, shapiro.test(data_interv)$p.value, NA)
  p_norm_ctrl <- ifelse(length(data_ctrl) >= 3, shapiro.test(data_ctrl)$p.value, NA)
  normality_status <- ifelse(!is.na(p_norm_interv) & !is.na(p_norm_ctrl) & p_norm_interv > 0.05 & p_norm_ctrl > 0.05, "Normal", "Non-Normal")
  
  levene_group <- suppressWarnings(car::leveneTest(y = subset_data[[col]], group = subset_data$GROUP))
  levene_gender <- suppressWarnings(car::leveneTest(y = subset_data[[col]], group = subset_data$Gender))
  
  assumption_results[[col]] <- data.frame(
    `Variable` = col, `Normality Interv (p)` = round(p_norm_interv, 4), `Normality Control (p)` = round(p_norm_ctrl, 4), `Distribution` = normality_status,
    `Homogeneity Group (p)` = round(levene_group$`Pr(>F)`[1], 4), `Homogeneity Gender (p)` = round(levene_gender$`Pr(>F)`[1], 4), check.names = FALSE)
}
table_assumptions <- bind_rows(assumption_results)

# ==========================================
# STAGE 3: MANN-WHITNEY U TEST (OVERALL)
# ==========================================
cat("[Info] Performing Mann-Whitney U Test (Overall Change)...\n")
mwu_results <- list()

for (col in change_score_cols) {
  subset_data <- data_change %>% drop_na(all_of(c(col, "GROUP")))
  if(nrow(subset_data) < 2) next 
  
  mwu_test <- wilcox.test(subset_data[[col]] ~ subset_data$GROUP, exact = FALSE, correct = FALSE)
  z_score <- qnorm(mwu_test$p.value / 2)
  effect_size_r <- abs(z_score) / sqrt(nrow(subset_data))
  
  mwu_results[[col]] <- data.frame(
    `Change Score` = col, `U Statistic` = round(mwu_test$statistic, 2), `P-value` = round(mwu_test$p.value, 4),
    `Effect Size (r)` = round(effect_size_r, 4), `Significance` = ifelse(mwu_test$p.value < 0.05, "*", ""), check.names = FALSE)
}
table_mwu <- bind_rows(mwu_results)

# ==========================================
# STAGE 4 & 5: LMM (PRIMARY AND SUBSCALES)
# ==========================================
cat("[Info] Running Linear Mixed-Effects Models (LMM) with robust BOBYQA optimizer...\n")

run_lmm_models <- function(vars) {
  res_list <- list()
  mod_list <- list()
  for (var_y in vars) {
    lmm_formula <- as.formula(paste(var_y, "~ GROUP * TIME + Gender + (1 | RESPONDENT_ID)"))
    
    model_lmm <- lmer(lmm_formula, data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    mod_list[[var_y]] <- model_lmm 
    
    anova_lmm <- anova(model_lmm)
    coefs <- summary(model_lmm)$coefficients
    
    # Karena Control = Reference, nama koefisiennya otomatis "GROUPIntervention"
    beta_interaction <- coefs["GROUPIntervention:TIME", "Estimate"]
    p_val_interaction <- anova_lmm["GROUP:TIME", "Pr(>F)"]
    f_val_interaction <- anova_lmm["GROUP:TIME", "F value"]
    
    eta_sq <- suppressWarnings(effectsize::eta_squared(model_lmm, partial = TRUE))
    interaction_eta <- eta_sq[eta_sq$Parameter == "GROUP:TIME", "Eta2_partial"]
    cohens_d <- 2 * sqrt(interaction_eta / (1 - interaction_eta))
    correction_factor <- 1 - (3 / (4 * (nrow(data) - 2) - 1))
    
    res_list[[var_y]] <- data.frame(
      Variable = var_y, `Main Effect GROUP (p)` = round(anova_lmm["GROUP", "Pr(>F)"], 4), `Main Effect TIME (p)` = round(anova_lmm["TIME", "Pr(>F)"], 4),
      `Interaction F-val` = round(f_val_interaction, 2), `Interaction p-val` = round(p_val_interaction, 4), `Beta` = round(beta_interaction, 2),
      `Hedge's g` = round(cohens_d * correction_factor, 4), `Significance` = ifelse(p_val_interaction < 0.05, "*", ""), check.names = FALSE)
  }
  return(list(results = bind_rows(res_list), models = mod_list))
}

primary_variables <- c("Total_GAS", "Total_BPAQ", "Total_SDQ", "Assertiveness", "GO_Bijak")
subscale_variables <- setdiff(score_columns, primary_variables)

lmm_primary <- run_lmm_models(primary_variables)
lmm_subscales <- run_lmm_models(subscale_variables)

table_lmm <- lmm_primary$results
table_lmm_sub <- lmm_subscales$results

# ==========================================
# STAGE 6: EXPORT WIDE REGRESSION TABLE TO MS WORD (.doc)
# ==========================================
cat("\n[Info] Exporting WIDE Regression Table directly to Microsoft Word...\n")
tryCatch({
  # Fungsi sjPlot otomatis membuat tabel format Wide dan langsung jadi .doc
  tab_model(
    lmm_primary$models[[1]], lmm_primary$models[[2]], lmm_primary$models[[3]], lmm_primary$models[[4]], lmm_primary$models[[5]],
    file = "Table_3_LMM_Regression_Primary.doc", 
    dv.labels = names(lmm_primary$models),
    title = "Table 1. Linear Mixed-Effects Models Predicting Primary Outcomes Over Time",
    pred.labels = c("(Intercept)", "Intervention (vs Control)", "Time", "Female (vs Male)", "Intervention × Time"),
    show.se = TRUE, show.stat = FALSE, show.ci = 0.95, collapse.ci = TRUE, p.style = "numeric_stars"
  )
  
  tab_model(
    lmm_subscales$models[[1]], lmm_subscales$models[[2]], lmm_subscales$models[[3]], lmm_subscales$models[[4]], lmm_subscales$models[[5]],
    file = "Table_4_LMM_Regression_Subscales_Part1.doc",
    dv.labels = names(lmm_subscales$models)[1:5],
    title = "Table 2. Linear Mixed-Effects Models Predicting Subscale Outcomes (Part 1)",
    show.se = TRUE, show.stat = FALSE, show.ci = 0.95, collapse.ci = TRUE, p.style = "numeric_stars"
  )
  cat("[Success] Tabel Word format WIDE berhasil dibuat (Silakan buka file .doc yang baru)\n")
}, error = function(e) { cat("[Error] Could not generate Word tables.\n") })

# ==========================================
# STAGE 7: LONGITUDINAL LINE GRAPHS
# ==========================================
cat("\n[Info] Generating Longitudinal Line Graphs...\n")
significant_lmm_vars <- table_lmm$Variable[table_lmm$`Significance` == "*"]

if(length(significant_lmm_vars) > 0) {
  data_long_plot <- data %>% dplyr::select(RESPONDENT_ID, GROUP, TIME, all_of(significant_lmm_vars)) %>%
    pivot_longer(cols = all_of(significant_lmm_vars), names_to = "Variable", values_to = "Score")
  data_long_plot$Score <- as.numeric(data_long_plot$Score)
  
  p_trend <- ggplot(data_long_plot, aes(x = TIME, y = Score, color = GROUP, fill = GROUP)) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 2.5, shape = 21, color = "white", stroke = 1) +
    facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Intervention" = "#00BFC4", "Control" = "#F8766D")) +
    scale_fill_manual(values = c("Intervention" = "#00BFC4", "Control" = "#F8766D")) +
    labs(title = "Linear Trend of Score Trajectories Over Time", x = "Measurement Time (TIME 0, 1, 2, 3)", y = "Score") +
    theme_minimal() + theme(legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA))
  
  ggsave("Figure_1_LMM_Line_Graphs.png", plot = p_trend, width = 12, height = 7, dpi = 300, bg = "white")
}

# ==========================================
# STAGE 8 & 9: TRUE FOREST PLOTS (PRIMARY & SUBSCALES)
# ==========================================
cat("\n[Info] Calculating Specific Time Point Differences for Forest Plots...\n")
forest_results <- list()

for (var in primary_variables) {
  base_d <- data %>% filter(TIME == 0) %>% dplyr::select(RESPONDENT_ID, Base_Score = all_of(var))
  for (t in c(1, 2, 3)) {
    post_d <- data %>% filter(TIME == t) %>% dplyr::select(RESPONDENT_ID, GROUP, Post_Score = all_of(var))
    merged_d <- post_d %>% inner_join(base_d, by = "RESPONDENT_ID") %>% mutate(Change = Post_Score - Base_Score) %>% drop_na(Change, GROUP)
    if(nrow(merged_d) < 2) next
    
    # [KRUSIAL UNTUK FOREST PLOT]: 
    # Urutan ini menentukan rumusan "Intervention minus Control"
    merged_d$GROUP <- factor(merged_d$GROUP, levels = c("Intervention", "Control")) 
    
    eff <- suppressWarnings(effectsize::hedges_g(Change ~ GROUP, data = merged_d))
    tt <- t.test(Change ~ GROUP, data = merged_d)
    
    forest_results[[paste(var, t)]] <- data.frame(
      Variable = var, Timepoint = paste("Post-test", t), Hedges_g = eff$Hedges_g,
      CI_low = eff$CI_low, CI_high = eff$CI_high, p_value = tt$p.value, stringsAsFactors = FALSE)
  }
}

forest_data <- bind_rows(forest_results)
forest_data$Timepoint <- factor(forest_data$Timepoint, levels = rev(c("Post-test 1", "Post-test 2", "Post-test 3")))
forest_data$Variable <- factor(forest_data$Variable, levels = primary_variables)

forest_data <- forest_data %>%
  mutate(
    p_label = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
    stat_label = sprintf("g = %.2f [%.2f, %.2f]  |  p = %s", Hedges_g, CI_low, CI_high, p_label),
    Direction = case_when(
      Hedges_g < 0 & grepl("GAS|BPAQ|SDQ", Variable) & !grepl("Protective", Variable) ~ "Favorable",
      Hedges_g > 0 & grepl("Assertiveness|GO_Bijak|Protective", Variable) ~ "Favorable",
      TRUE ~ "Not Favorable"
    )
  )

max_ci_vals <- forest_data %>% group_by(Variable) %>% summarise(max_x = max(CI_high, na.rm = TRUE))
forest_data <- forest_data %>% left_join(max_ci_vals, by = "Variable")

p_forest <- ggplot(forest_data, aes(y = Timepoint, x = Hedges_g)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = Direction), height = 0.2, linewidth = 1) +
  geom_point(aes(color = Direction), size = 3, shape = 15) +
  geom_text(aes(x = max_x + 0.3, label = stat_label), hjust = 0, size = 3.5, fontface = "italic") +
  facet_wrap(~ Variable, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Favorable" = "#00BFC4", "Not Favorable" = "#F8766D")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.7))) +
  labs(
    title = "Forest Plot of Intervention Effects at Each Follow-up (Primary)",
    subtitle = "Standardized Mean Differences (Hedges' g) and 95% CIs comparing Intervention vs Control change scores",
    x = "Hedges' g (Negative = Greater Reduction in Intervention)", y = "Measurement Timepoint"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none", strip.text = element_text(face = "bold", size = 11, hjust = 0, color = "#333333"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "#666666"), plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("Figure_2_True_Forest_Plot_Primary.png", plot = p_forest, width = 11, height = 12, dpi = 300, bg = "white")
cat("\n[Success] Forest Plot successfully generated and saved as 'Figure_2_True_Forest_Plot_Primary.png'.\n")

# (Visualisasi Forest plot Subskala dihapus dari tampilan CLI agar fokus, file tetap sama logikanya)

# ==========================================
# STAGE 10: CAUSAL MEDIATION ANALYSIS
# ==========================================
cat("\n[Info] Performing Causal Mediation Analysis...\n")
mediation_data <- data_change %>% dplyr::select(RESPONDENT_ID, GROUP, Gender, change_Total_GAS, change_Assertiveness, change_GO_Bijak) %>% drop_na()
# 1 = Intervensi, 0 = Kontrol
mediation_data$GROUP_bin <- ifelse(mediation_data$GROUP == "Intervention", 1, 0)

model_m_assert <- lm(change_Assertiveness ~ GROUP_bin + Gender, data = mediation_data)
model_y_assert <- lm(change_Total_GAS ~ GROUP_bin + change_Assertiveness + Gender, data = mediation_data)

set.seed(123)
med_out_assert <- mediate(model_m_assert, model_y_assert, treat = "GROUP_bin", mediator = "change_Assertiveness", boot = TRUE, sims = 1000)

med_assert <- data.frame(
  `Mediator` = "Assertiveness", `ACME (Indirect Effect)` = round(med_out_assert$d1, 3), `ACME p-value` = round(med_out_assert$d1.p, 4),
  `ADE (Direct Effect)` = round(med_out_assert$z0, 3), `ADE p-value` = round(med_out_assert$z0.p, 4),
  `Total Effect` = round(med_out_assert$tau.coef, 3), `Total Effect p-value` = round(med_out_assert$tau.p, 4),
  `Prop. Mediated` = round(med_out_assert$n0, 3), `Prop. Mediated p-value` = round(med_out_assert$n0.p, 4), check.names = FALSE)

model_m_bijak <- lm(change_GO_Bijak ~ GROUP_bin + Gender, data = mediation_data)
model_y_bijak <- lm(change_Total_GAS ~ GROUP_bin + change_GO_Bijak + Gender, data = mediation_data)

set.seed(123)
med_out_bijak <- mediate(model_m_bijak, model_y_bijak, treat = "GROUP_bin", mediator = "change_GO_Bijak", boot = TRUE, sims = 1000)

med_bijak <- data.frame(
  `Mediator` = "Wise Gaming (GO_Bijak)", `ACME (Indirect Effect)` = round(med_out_bijak$d1, 3), `ACME p-value` = round(med_out_bijak$d1.p, 4),
  `ADE (Direct Effect)` = round(med_out_bijak$z0, 3), `ADE p-value` = round(med_out_bijak$z0.p, 4),
  `Total Effect` = round(med_out_bijak$tau.coef, 3), `Total Effect p-value` = round(med_out_bijak$tau.p, 4),
  `Prop. Mediated` = round(med_out_bijak$n0, 3), `Prop. Mediated p-value` = round(med_out_bijak$n0.p, 4), check.names = FALSE)

table_mediation <- bind_rows(med_assert, med_bijak)


# ==========================================
# STAGE 11: EXPORT TEXT SUMMARY
# ==========================================
sink("Publication_Results_Summary.txt")
cat("--- Linear Mixed-Effects Models (LMM) Summary - Primary ---\n")
print(table_lmm, row.names = FALSE)
cat("\n\n--- Causal Mediation Analysis Results ---\n")
print(table_mediation, row.names = FALSE)
cat("\n\n--- Specific Timepoint Differences (Primary Outcomes) ---\n")
forest_out <- forest_data %>% dplyr::select(Variable, Timepoint, Hedges_g, CI_low, CI_high, p_value, Direction) %>% mutate(across(c(Hedges_g, CI_low, CI_high, p_value), ~round(., 4)))
print(forest_out, row.names = FALSE)
sink()


# ==========================================
# STAGE 13: DISSECTING AGGRESSION SUB-VARIABLES (FIXED DIRECTION)
# ==========================================
cat("\n[Info] Membedah Sub-variabel Agresivitas (Sinkronisasi Arah)...\n")

bpaq_subscales <- c("BPAQ_Physical", "BPAQ_Verbal", "BPAQ_Anger", "BPAQ_Hostility")
forest_bpaq_results <- list()

for (var in bpaq_subscales) {
  base_d <- data %>% filter(TIME == 0) %>% dplyr::select(RESPONDENT_ID, Base_Score = all_of(var))
  
  for (t in c(1, 2, 3)) {
    post_d <- data %>% filter(TIME == t) %>% dplyr::select(RESPONDENT_ID, GROUP, Post_Score = all_of(var))
    merged_d <- post_d %>% inner_join(base_d, by = "RESPONDENT_ID") %>%
      mutate(Change = Post_Score - Base_Score) %>% drop_na(Change, GROUP)
    
    if(nrow(merged_d) < 2) next
    
    # [SINKRONISASI]: Memaksa urutan (Intervention, Control) 
    # agar hasil Negatif (-) berarti Intervention lebih rendah (BAGUS)
    merged_d$GROUP <- factor(merged_d$GROUP, levels = c("Intervention", "Control"))
    
    eff <- suppressWarnings(effectsize::hedges_g(Change ~ GROUP, data = merged_d))
    tt <- t.test(Change ~ GROUP, data = merged_d)
    
    # Tentukan Status Berdasarkan Angka
    # Negatif (-) = Favorable (Turun), Positif (+) = Spike (Naik)
    res_g <- eff$Hedges_g
    stat_label <- ifelse(res_g < 0, "REDUCTION (Favorable)", "SPIKE (Not Favorable)")
    
    forest_bpaq_results[[paste(var, t)]] <- data.frame(
      Sub_Variable = var,
      Timepoint = paste("Post-test", t),
      Hedges_g = round(res_g, 4),
      p_value = round(tt$p.value, 4),
      Status = stat_label
    )
  }
}

table_bpaq_sub_final <- bind_rows(forest_bpaq_results)

# Cetak hasil akhir yang sudah sinkron
cat("\n--- Rincian Perubahan Sub-Variabel Agresivitas (SINKRON) ---\n")
print(table_bpaq_sub_final, row.names = FALSE)


# ==========================================
# STAGE 14: ITEM-LEVEL CHANGE ANALYSIS FOR WISE GAMING
# ==========================================
cat("\n[Info] Menganalisis Perubahan Skor per Item pada Wise Gaming (GO_Bijak)...\n")

# Menyesuaikan dengan nama kolom terbaru (GO_Bijak1 hingga GO_Bijak15)
wisegaming_items <- paste0("GO_Bijak", 1:15)

missing_items <- setdiff(wisegaming_items, colnames(data))
if(length(missing_items) > 0) {
  cat("[Peringatan] Kolom item Wise Gaming berikut tidak ditemukan di dataset Anda:\n")
  print(missing_items)
} else {
  
  item_change_results <- list()
  
  for (item in wisegaming_items) {
    base_d <- data %>% filter(TIME == 0) %>% dplyr::select(RESPONDENT_ID, Base_Score = all_of(item))
    post_d <- data %>% filter(TIME == 3) %>% dplyr::select(RESPONDENT_ID, GROUP, Post_Score = all_of(item))
    
    merged_d <- post_d %>% inner_join(base_d, by = "RESPONDENT_ID") %>%
      mutate(Change = Post_Score - Base_Score) %>% drop_na(Change, GROUP)
    
    if(nrow(merged_d) < 2) next
    
    interv_data <- merged_d %>% filter(GROUP == "Intervention")
    
    mean_change <- mean(interv_data$Change, na.rm = TRUE)
    sd_change <- sd(interv_data$Change, na.rm = TRUE)
    
    item_change_results[[item]] <- data.frame(
      Item = item,
      `Mean Change (Post3 - Pre)` = round(mean_change, 3),
      `SD Change` = round(sd_change, 3),
      check.names = FALSE
    )
  }
  
  table_item_change <- bind_rows(item_change_results)
  table_item_change <- table_item_change %>% arrange(desc(`Mean Change (Post3 - Pre)`))
  
  cat("\n--- Rata-Rata Kenaikan Skor per Item Wise Gaming (Kelompok Intervensi) ---\n")
  print(table_item_change, row.names = FALSE)
  write.csv(table_item_change, "WiseGaming_Item_Change_Analysis.csv", row.names = FALSE)
  
  # ==========================================
  # STAGE 15: SUB-DIMENSION MEDIATION ANALYSIS
  # ==========================================
  cat("\n[Info] Melakukan Analisis Mediasi berdasarkan Sub-Dimensi Wise Gaming...\n")
  
  # Agregasi skor untuk masing-masing Sub-dimensi berdasarkan kolom GO_Bijak1 - 15
  data <- data %>%
    mutate(
      WG_SelfRegulation = rowMeans(dplyr::select(., all_of(paste0("GO_Bijak", c(1:6, 9)))), na.rm = TRUE),
      WG_HealthyContent = rowMeans(dplyr::select(., all_of(paste0("GO_Bijak", 7:8))), na.rm = TRUE),
      WG_LifeBalance    = rowMeans(dplyr::select(., all_of(paste0("GO_Bijak", 10:15))), na.rm = TRUE)
    )
  
  data_change_wg <- data %>%
    arrange(RESPONDENT_ID, TIME) %>%
    group_by(RESPONDENT_ID, GROUP, Gender) %>%
    summarise(
      change_Total_GAS = last(Total_GAS) - first(Total_GAS),
      change_SelfReg   = last(WG_SelfRegulation) - first(WG_SelfRegulation),
      change_Healthy   = last(WG_HealthyContent) - first(WG_HealthyContent),
      change_LifeBal   = last(WG_LifeBalance) - first(WG_LifeBalance),
      .groups = 'drop'
    ) %>% drop_na()
  
  data_change_wg$GROUP_bin <- ifelse(data_change_wg$GROUP == "Intervention", 1, 0)
  
  # --- Explicit Models to avoid scoping issues during bootstrapping ---
  
  # 1. Mediation for Self-Regulation
  cat("   -> Running mediation for Self-Regulation (Items 1-6, 9)...\n")
  model_m_selfreg <- lm(change_SelfReg ~ GROUP_bin + Gender, data = data_change_wg)
  model_y_selfreg <- lm(change_Total_GAS ~ GROUP_bin + change_SelfReg + Gender, data = data_change_wg)
  set.seed(123)
  med_out_selfreg <- mediate(model_m_selfreg, model_y_selfreg, treat = "GROUP_bin", mediator = "change_SelfReg", boot = TRUE, sims = 1000)
  
  med_selfreg <- data.frame(
    `Sub-Dimension` = "Self-Regulation (Items 1-6, 9)",
    `ACME (Indirect Effect)` = round(med_out_selfreg$d1, 3),
    `ACME p-value` = round(med_out_selfreg$d1.p, 4),
    `Prop. Mediated (%)` = round(med_out_selfreg$n0 * 100, 1),
    check.names = FALSE
  )
  
  # 2. Mediation for Healthy Content
  cat("   -> Running mediation for Healthy Content (Items 7-8)...\n")
  model_m_healthy <- lm(change_Healthy ~ GROUP_bin + Gender, data = data_change_wg)
  model_y_healthy <- lm(change_Total_GAS ~ GROUP_bin + change_Healthy + Gender, data = data_change_wg)
  set.seed(123)
  med_out_healthy <- mediate(model_m_healthy, model_y_healthy, treat = "GROUP_bin", mediator = "change_Healthy", boot = TRUE, sims = 1000)
  
  med_healthy <- data.frame(
    `Sub-Dimension` = "Healthy Content (Items 7-8)",
    `ACME (Indirect Effect)` = round(med_out_healthy$d1, 3),
    `ACME p-value` = round(med_out_healthy$d1.p, 4),
    `Prop. Mediated (%)` = round(med_out_healthy$n0 * 100, 1),
    check.names = FALSE
  )
  
  # 3. Mediation for Life Balance
  cat("   -> Running mediation for Life Balance (Items 10-15)...\n")
  model_m_lifebal <- lm(change_LifeBal ~ GROUP_bin + Gender, data = data_change_wg)
  model_y_lifebal <- lm(change_Total_GAS ~ GROUP_bin + change_LifeBal + Gender, data = data_change_wg)
  set.seed(123)
  med_out_lifebal <- mediate(model_m_lifebal, model_y_lifebal, treat = "GROUP_bin", mediator = "change_LifeBal", boot = TRUE, sims = 1000)
  
  med_lifebal <- data.frame(
    `Sub-Dimension` = "Life Balance (Items 10-15)",
    `ACME (Indirect Effect)` = round(med_out_lifebal$d1, 3),
    `ACME p-value` = round(med_out_lifebal$d1.p, 4),
    `Prop. Mediated (%)` = round(med_out_lifebal$n0 * 100, 1),
    check.names = FALSE
  )
  
  table_sub_mediation <- bind_rows(med_selfreg, med_healthy, med_lifebal)
  
  cat("\n--- Hasil Analisis Mediasi per Sub-Dimensi Wise Gaming ---\n")
  print(table_sub_mediation, row.names = FALSE)
  write.csv(table_sub_mediation, "WiseGaming_SubDimension_Mediation.csv", row.names = FALSE)
  
  cat("\n[Success] Analisis Wise Gaming selesai. Output tersimpan di CSV.\n")
}

# ==========================================
# STAGE 16: INITIAL RECOVERY FRICTION ANALYSIS (POST-TEST 1)
# ==========================================
cat("\n[Info] Menganalisis Hubungan antara Perbaikan Wise Gaming dan Lonjakan Agresi di Bulan ke-1...\n")

# 1. Menghitung Perubahan (Change Score) Khusus untuk Post-test 1 (Bulan ke-1)
data_t1 <- data %>%
  filter(TIME %in% c(0, 1)) %>%
  arrange(RESPONDENT_ID, TIME) %>%
  group_by(RESPONDENT_ID, GROUP, Gender) %>%
  summarise(
    diff_GO_Bijak = last(GO_Bijak) - first(GO_Bijak),
    diff_Assertiveness = last(Assertiveness) - first(Assertiveness),
    diff_Total_BPAQ = last(Total_BPAQ) - first(Total_BPAQ),
    diff_Hostility = last(BPAQ_Hostility) - first(BPAQ_Hostility),
    .groups = 'drop'
  ) %>%
  filter(GROUP == "Intervention") # Kita hanya melihat dinamika di kelompok Intervensi

# 2. Regresi: Apakah Perbaikan Wise Gaming memprediksi kenaikan Agresi di Bulan ke-1?
# Model: Perubahan Agresi ~ Perubahan Wise Gaming + Perubahan Assertiveness
model_friction <- lm(diff_Total_BPAQ ~ diff_GO_Bijak + diff_Assertiveness, data = data_t1)
summary_friction <- summary(model_friction)

# 3. Regresi Khusus untuk Hostility (karena ini yang tadi paling signifikan lonjakannya)
model_hostility_friction <- lm(diff_Hostility ~ diff_GO_Bijak + diff_Assertiveness, data = data_t1)
summary_hostility <- summary(model_hostility_friction)

cat("\n--- HASIL REGRESI: PREDIKTOR AGRESI PADA BULAN KE-1 (POST-TEST 1) ---\n")
cat("Variabel Dependen: Perubahan Total Agresi (BPAQ)\n")
print(coef(summary_friction))

cat("\nVariabel Dependen: Perubahan Hostility (Faset Kognitif Agresi)\n")
print(coef(summary_hostility))

# 4. Korelasi Sederhana untuk memperkuat narasi
cor_val <- cor(data_t1$diff_GO_Bijak, data_t1$diff_Total_BPAQ)
cat(sprintf("\nKorelasi antara Perbaikan Wise Gaming & Kenaikan Agresi di Bulan 1: r = %.3f\n", cor_val))

# ==========================================
# STAGE END: EXPORT AND PRINT TO CONSOLE
# ==========================================
cat("\n[Info] Finalizing output... Printing to console and saving to file.\n")

# Menggunakan split = TRUE agar output tampil di RStudio Console DAN tersimpan ke file .txt
sink("Publication_Results_Summary.txt", split = TRUE)

cat("\n=======================================================\n")
cat("          STATISTICAL SUMMARY FOR PUBLICATION          \n")
cat("=======================================================\n\n")

cat("--- Linear Mixed-Effects Models (LMM) Summary - Primary ---\n")
print(table_lmm, row.names = FALSE)
cat("\nNote for manuscript:\nSignificant 'Intervention p-val' indicates the differential effect of the intervention over time compared to control.\n")

cat("\n--- Causal Mediation Analysis Results ---\n")
print(table_mediation, row.names = FALSE)
cat("\nNote for manuscript:\nACME represents the indirect effect. A significant ACME (p < 0.05) indicates that the variable successfully mediates the intervention's effect.\n")

cat("\n--- Specific Timepoint Differences (Forest Plot Data) ---\n")
forest_out <- forest_data %>% dplyr::select(Variable, Timepoint, Hedges_g, CI_low, CI_high, p_value, Direction) %>% mutate(across(c(Hedges_g, CI_low, CI_high, p_value), ~round(., 4)))
print(forest_out, row.names = FALSE)
cat("\nNote for manuscript:\nHedges' g effect size comparing Intervention vs Control at specific timepoints.\n")

cat("\n=======================================================\n")
cat("          ALL ANALYSES COMPLETED SUCCESSFULLY          \n")
cat("=======================================================\n")

sink()
