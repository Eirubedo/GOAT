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

# [UPDATE]: Based on the latest dataset (0 = Control, 1 = Intervention)
# Convert to factor with clear labels, 
# 0 (Control) automatically becomes the first level (Absolute baseline)
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
# STAGE 1.1: BASELINE DEMOGRAPHIC PROFILING (TABLE 1)
# ==========================================
cat("\n[Info] Generating Baseline Demographic Profiling (Intervention vs Control)...\n")

baseline_data <- data %>% filter(TIME == 0)

# Define variables to profile
demo_vars_cat <- c("Gender", "Live_With", "Parent_Occupation", "SocEco_Status", "Communication_Parent_Style", "Gaming_Duration")
baseline_summary <- list()

# Age (Continuous variable -> Independent T-test)
age_t <- t.test(Age ~ GROUP, data = baseline_data)
baseline_summary[["Age"]] <- data.frame(
  Variable = "Age (Mean ± SD)",
  Category = "-",
  Control = sprintf("%.2f ± %.2f", mean(baseline_data$Age[baseline_data$GROUP=="Control"], na.rm=TRUE), sd(baseline_data$Age[baseline_data$GROUP=="Control"], na.rm=TRUE)),
  Intervention = sprintf("%.2f ± %.2f", mean(baseline_data$Age[baseline_data$GROUP=="Intervention"], na.rm=TRUE), sd(baseline_data$Age[baseline_data$GROUP=="Intervention"], na.rm=TRUE)),
  P_Value = as.character(round(age_t$p.value, 4)),
  stringsAsFactors = FALSE
)

# Categorical variables -> Chi-Square Test
for (v in demo_vars_cat) {
  if(v %in% colnames(baseline_data)) {
    tab <- table(baseline_data[[v]], baseline_data$GROUP)
    
    # Suppress warnings for approximations with small sample sizes
    test_res <- suppressWarnings(chisq.test(tab))
    prop_tab <- prop.table(tab, margin = 2) * 100
    
    for(cat_level in rownames(tab)) {
      ctrl_val <- sprintf("%d (%.1f%%)", tab[cat_level, "Control"], prop_tab[cat_level, "Control"])
      intv_val <- sprintf("%d (%.1f%%)", tab[cat_level, "Intervention"], prop_tab[cat_level, "Intervention"])
      
      # Show p-value only on the first category row for a clean table look
      pval_display <- ifelse(cat_level == rownames(tab)[1], as.character(round(test_res$p.value, 4)), "")
      
      baseline_summary[[paste(v, cat_level)]] <- data.frame(
        Variable = ifelse(cat_level == rownames(tab)[1], v, ""),
        Category = cat_level,
        Control = ctrl_val,
        Intervention = intv_val,
        P_Value = pval_display,
        stringsAsFactors = FALSE
      )
    }
  }
}
table_baseline <- bind_rows(baseline_summary)
cat("--- Baseline Characteristics (Table 1) ---\n")
print(table_baseline, row.names = FALSE)


# ==========================================
# STAGE 1.5: BASELINE DEMOGRAPHIC CORRELATION HEATMAP
# ==========================================
cat("\n[Info] Generating Demographic vs Baseline Outcomes Correlation Heatmap...\n")

# Filter for baseline data (TIME == 0) regardless of group
baseline_data <- data %>% filter(TIME == 0)

# 1. Dummy/Ordinal Coding for Categorical Demographics
# Gender: Male = 0, Female = 1
baseline_data$Gender_Num <- ifelse(tolower(as.character(baseline_data$Gender)) %in% c("perempuan", "female"), 1, 0)

# Communication Style: Closed = 0, Open = 1
if("Communication_Parent_Style" %in% colnames(baseline_data)) {
  baseline_data$Comm_Style_Num <- ifelse(tolower(as.character(baseline_data$Communication_Parent_Style)) == "terbuka", 1, 0)
}

# Socioeconomic Status: Lower/Same = 0, Higher = 1
if("SocEco_Status" %in% colnames(baseline_data)) {
  baseline_data$SocEco_Num <- ifelse(grepl("lebih besar", tolower(as.character(baseline_data$SocEco_Status))), 1, 0)
}

# Gaming Duration (Ordinal): <1h = 1, 1-2h = 2, 3-4h = 3, >4h = 4
if("Gaming_Duration" %in% colnames(baseline_data)) {
  baseline_data$Gaming_Dur_Num <- case_when(
    grepl("<1", as.character(baseline_data$Gaming_Duration)) ~ 1,
    grepl("1-2", as.character(baseline_data$Gaming_Duration)) ~ 2,
    grepl("3-4", as.character(baseline_data$Gaming_Duration)) ~ 3,
    grepl(">4", as.character(baseline_data$Gaming_Duration)) ~ 4,
    TRUE ~ NA_real_
  )
}

# Select relevant numeric columns for correlation (ADDED BPAQ_Anger)
cor_vars <- c("Age", "Gender_Num", "Comm_Style_Num", "SocEco_Num", "Gaming_Dur_Num", "Total_SDQ", "GO_Bijak", "Assertiveness", "Total_GAS", "Total_BPAQ", "BPAQ_Anger")

# Check if columns exist
valid_cor_vars <- cor_vars[cor_vars %in% colnames(baseline_data)]

if(length(valid_cor_vars) > 1) {
  cor_data <- baseline_data %>% dplyr::select(all_of(valid_cor_vars))
  
  # Calculate correlation matrix using Spearman (better for ordinal data like Gaming_Dur_Num)
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Format for ggplot2
  cor_df <- as.data.frame(as.table(cor_matrix))
  colnames(cor_df) <- c("Var1", "Var2", "Correlation")
  
  # Plot Heatmap
  p_heatmap <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
    scale_fill_gradient2(low = "#F8766D", high = "#00BFC4", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = "Correlation Heatmap: Demographics vs Baseline Outcomes (All Participants)",
         subtitle = "Positive value indicates variables move in the same direction.") +
    coord_fixed()
  
  ggsave("Figure_0_Baseline_Correlation_Heatmap.png", plot = p_heatmap, width = 9, height = 8, dpi = 300, bg = "white")
  cat("[Success] Heatmap successfully generated and saved as 'Figure_0_Baseline_Correlation_Heatmap.png'.\n")
} else {
  cat("[Warning] Could not generate heatmap: insufficient variables.\n")
}


# ==========================================
# STAGE 1.6: GAMING DURATION VS BASELINE OUTCOMES
# ==========================================
cat("\n[Info] Computing Relationship between Gaming Duration and Baseline Outcomes...\n")

if("Gaming_Dur_Num" %in% colnames(baseline_data)) {
  baseline_data_valid <- baseline_data %>% filter(!is.na(Gaming_Dur_Num))
  
  cor_gas <- cor.test(baseline_data_valid$Gaming_Dur_Num, baseline_data_valid$Total_GAS, method = "spearman", exact = FALSE)
  cor_bpaq <- cor.test(baseline_data_valid$Gaming_Dur_Num, baseline_data_valid$Total_BPAQ, method = "spearman", exact = FALSE)
  cor_sdq <- cor.test(baseline_data_valid$Gaming_Dur_Num, baseline_data_valid$Total_SDQ, method = "spearman", exact = FALSE)
  cor_gobijak <- cor.test(baseline_data_valid$Gaming_Dur_Num, baseline_data_valid$GO_Bijak, method = "spearman", exact = FALSE)
  
  cat("--- Spearman Rank Correlation (Gaming Duration vs Baseline Outcomes) ---\n")
  cat(sprintf("1. Gaming Duration vs Total GAS  : rho = %.3f, p-value = %.4f %s\n", cor_gas$estimate, cor_gas$p.value, ifelse(cor_gas$p.value < 0.05, "*", "")))
  cat(sprintf("2. Gaming Duration vs Total BPAQ : rho = %.3f, p-value = %.4f %s\n", cor_bpaq$estimate, cor_bpaq$p.value, ifelse(cor_bpaq$p.value < 0.05, "*", "")))
  cat(sprintf("3. Gaming Duration vs Total SDQ  : rho = %.3f, p-value = %.4f %s\n", cor_sdq$estimate, cor_sdq$p.value, ifelse(cor_sdq$p.value < 0.05, "*", "")))
  cat(sprintf("4. Gaming Duration vs Wise Gaming: rho = %.3f, p-value = %.4f %s\n", cor_gobijak$estimate, cor_gobijak$p.value, ifelse(cor_gobijak$p.value < 0.05, "*", "")))
  cat("Note: Positive 'rho' for GAS/BPAQ/SDQ means longer duration is associated with higher problem scores. Negative 'rho' for Wise Gaming means longer duration is associated with lower wise behavior.\n")
}

# ==========================================
# STAGE 1.7: SPECIFIC BASELINE CORRELATIONS FOR MANUSCRIPT CLAIMS
# ==========================================
cat("\n[Info] Computing exact rho and p-values for manuscript claims (Gender, SES, Gaming Duration vs Outcomes)...\n")

manuscript_cors <- list()
predictors <- c("Gender_Num", "SocEco_Num", "Gaming_Dur_Num")
pred_labels <- c("Gender (1=Female)", "SES (1=Higher)", "Gaming Duration")
outcomes <- c("Total_BPAQ", "BPAQ_Anger", "Total_SDQ", "Total_GAS", "GO_Bijak")

for (i in seq_along(predictors)) {
  pred <- predictors[i]
  pred_lbl <- pred_labels[i]
  if(pred %in% colnames(baseline_data)) {
    valid_data <- baseline_data %>% filter(!is.na(.[[pred]]))
    for (out in outcomes) {
      if(out %in% colnames(valid_data)) {
        ctest <- cor.test(valid_data[[pred]], valid_data[[out]], method = "spearman", exact = FALSE)
        manuscript_cors[[paste(pred, out)]] <- data.frame(
          Predictor = pred_lbl,
          Outcome = out,
          rho = round(ctest$estimate, 3),
          p_value = round(ctest$p.value, 4),
          Significance = ifelse(ctest$p.value < 0.05, "*", ""),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
table_manuscript_cors <- bind_rows(manuscript_cors)

cat("--- Specific Baseline Correlations for Manuscript ---\n")
print(table_manuscript_cors, row.names = FALSE)


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
cat("\n[Info] Performing Mann-Whitney U Test (Overall Change)...\n")
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
cat("\n[Info] Running Linear Mixed-Effects Models (LMM) with robust BOBYQA optimizer...\n")

run_lmm_models <- function(vars) {
  res_list <- list()
  mod_list <- list()
  # Correction fix for Hedges' g calculation ensuring accurate df for repeated measures
  n_subj <- length(unique(data$RESPONDENT_ID))
  correction_factor <- 1 - (3 / (4 * (n_subj - 2) - 1))
  
  for (var_y in vars) {
    lmm_formula <- as.formula(paste(var_y, "~ GROUP * TIME + Gender + (1 | RESPONDENT_ID)"))
    
    model_lmm <- lmer(lmm_formula, data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    mod_list[[var_y]] <- model_lmm 
    
    anova_lmm <- anova(model_lmm)
    coefs <- summary(model_lmm)$coefficients
    
    # Because Control = Reference, the coefficient name is automatically "GROUPIntervention:TIME"
    beta_interaction <- coefs["GROUPIntervention:TIME", "Estimate"]
    p_val_interaction <- anova_lmm["GROUP:TIME", "Pr(>F)"]
    f_val_interaction <- anova_lmm["GROUP:TIME", "F value"]
    
    eta_sq <- suppressWarnings(effectsize::eta_squared(model_lmm, partial = TRUE))
    interaction_eta <- eta_sq[eta_sq$Parameter == "GROUP:TIME", "Eta2_partial"]
    cohens_d <- 2 * sqrt(interaction_eta / (1 - interaction_eta))
    
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
  # sjPlot automatically creates a Wide format table and directly outputs to .doc
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
  cat("[Success] WIDE format Word tables successfully created (Please open the new .doc files)\n")
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
    
    # [CRUCIAL FOR FOREST PLOT]: 
    # This order determines the "Intervention minus Control" formulation
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


# ==========================================
# STAGE 10: CAUSAL MEDIATION ANALYSIS
# ==========================================
cat("\n[Info] Performing Causal Mediation Analysis...\n")
mediation_data <- data_change %>% dplyr::select(RESPONDENT_ID, GROUP, Gender, change_Total_GAS, change_Assertiveness, change_GO_Bijak) %>% drop_na()
# 1 = Intervention, 0 = Control
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

# (Visualisasi terpisah dihapus, dipindahkan menjadi satu panel di STAGE 15)

# ==========================================
# STAGE 13: DISSECTING AGGRESSION SUB-VARIABLES (FIXED DIRECTION)
# ==========================================
cat("\n[Info] Dissecting Aggression Sub-variables (Direction Synchronization)...\n")

bpaq_subscales <- c("BPAQ_Physical", "BPAQ_Verbal", "BPAQ_Anger", "BPAQ_Hostility")
forest_bpaq_results <- list()

for (var in bpaq_subscales) {
  base_d <- data %>% filter(TIME == 0) %>% dplyr::select(RESPONDENT_ID, Base_Score = all_of(var))
  
  for (t in c(1, 2, 3)) {
    post_d <- data %>% filter(TIME == t) %>% dplyr::select(RESPONDENT_ID, GROUP, Post_Score = all_of(var))
    merged_d <- post_d %>% inner_join(base_d, by = "RESPONDENT_ID") %>%
      mutate(Change = Post_Score - Base_Score) %>% drop_na(Change, GROUP)
    
    if(nrow(merged_d) < 2) next
    
    # [SYNCHRONIZATION]: Forcing the order (Intervention, Control) 
    # so Negative (-) means Intervention is lower (FAVORABLE)
    merged_d$GROUP <- factor(merged_d$GROUP, levels = c("Intervention", "Control"))
    
    eff <- suppressWarnings(effectsize::hedges_g(Change ~ GROUP, data = merged_d))
    tt <- t.test(Change ~ GROUP, data = merged_d)
    
    # Determine Status Based on Value
    # Negative (-) = Favorable (Reduction), Positive (+) = Spike (Increase)
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

# Print the final synchronized results
cat("\n--- Detailed Changes in Aggression Sub-variables (SYNCHRONIZED) ---\n")
print(table_bpaq_sub_final, row.names = FALSE)


# ==========================================
# STAGE 14: ITEM-LEVEL CHANGE ANALYSIS FOR WISE GAMING
# ==========================================
cat("\n[Info] Analyzing Item-Level Score Changes in Wise Gaming (GO_Bijak)...\n")

# Adjusting to the latest column names (GO_Bijak1 to GO_Bijak15)
wisegaming_items <- paste0("GO_Bijak", 1:15)

missing_items <- setdiff(wisegaming_items, colnames(data))
if(length(missing_items) > 0) {
  cat("[Warning] The following Wise Gaming item columns were not found in your dataset:\n")
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
  
  cat("\n--- Average Item Score Increase for Wise Gaming (Intervention Group) ---\n")
  print(table_item_change, row.names = FALSE)
  write.csv(table_item_change, "WiseGaming_Item_Change_Analysis.csv", row.names = FALSE)
  
  # ==========================================
  # STAGE 15: SUB-DIMENSION MEDIATION ANALYSIS
  # ==========================================
  cat("\n[Info] Performing Mediation Analysis based on Wise Gaming Sub-Dimensions...\n")
  
  # Aggregate scores for each Sub-dimension based on columns GO_Bijak1 - 15
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
  
  cat("\n[Info] Generating Combined Single-Panel Mediation Plot (Figure 3)...\n")
  
  # Ekstrak data model mediasi ke dalam data.frame untuk ggplot2
  extract_med <- function(med_obj, name) {
    data.frame(
      Mediator = name,
      Effect = c("Total Effect", "ADE (Direct Effect)", "ACME (Indirect Effect)"),
      Estimate = c(med_obj$tau.coef, med_obj$z0, med_obj$d1),
      CI_lower = c(med_obj$tau.ci[1], med_obj$z0.ci[1], med_obj$d1.ci[1]),
      CI_upper = c(med_obj$tau.ci[2], med_obj$z0.ci[2], med_obj$d1.ci[2]),
      stringsAsFactors = FALSE
    )
  }
  
  med_plot_data <- bind_rows(
    extract_med(med_out_bijak, "1. Overall Wise Gaming"),
    extract_med(med_out_selfreg, "2. Sub: Self-Regulation"),
    extract_med(med_out_healthy, "3. Sub: Healthy Content"),
    extract_med(med_out_lifebal, "4. Sub: Life Balance")
  )
  
  # Urutkan faktor agar tampil dari atas ke bawah
  med_plot_data$Effect <- factor(med_plot_data$Effect, levels = rev(c("Total Effect", "ADE (Direct Effect)", "ACME (Indirect Effect)")))
  med_plot_data$Mediator <- factor(med_plot_data$Mediator, levels = rev(c("1. Overall Wise Gaming", "2. Sub: Self-Regulation", "3. Sub: Healthy Content", "4. Sub: Life Balance")))
  
  p_med <- ggplot(med_plot_data, aes(x = Estimate, y = Effect, color = Mediator)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper), position = position_dodge(width = 0.6), size = 1, linewidth = 1.2) +
    scale_color_manual(values = c("1. Overall Wise Gaming" = "#333333", 
                                  "2. Sub: Self-Regulation" = "#F8766D", 
                                  "3. Sub: Healthy Content" = "#00BA38", 
                                  "4. Sub: Life Balance" = "#619CFF")) +
    theme_minimal(base_size = 12) +
    labs(title = "Causal Mediation Analysis (Primary & Sub-dimensions)",
         x = "Effect Size (with 95% Confidence Intervals)", y = "", color = "Mediator Model") +
    theme(legend.position = "right", 
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(face = "bold", size = 11))
  
  ggsave("Figure_3_Combined_Mediation_Plot.png", plot = p_med, width = 10, height = 5, dpi = 300, bg = "white")
  cat("[Success] Combined Mechanism Plot saved as 'Figure_3_Combined_Mediation_Plot.png'.\n")
  
  cat("\n[Success] Analisis Wise Gaming selesai. Output tersimpan di CSV.\n")
}

# ==========================================
# STAGE 16: INITIAL RECOVERY FRICTION ANALYSIS (POST-TEST 1)
# ==========================================
cat("\n[Info] Analyzing the Relationship Between Wise Gaming Improvement and Aggression Spike in Month 1...\n")

# 1. Calculate Change Score Specifically for Post-test 1 (Month 1)
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
  filter(GROUP == "Intervention") # We are only observing the dynamics within the Intervention group

# 2. Regression: Does Wise Gaming Improvement predict the Aggression spike in Month 1?
# Model: Aggression Change ~ Wise Gaming Change + Assertiveness Change
model_friction <- lm(diff_Total_BPAQ ~ diff_GO_Bijak + diff_Assertiveness, data = data_t1)
summary_friction <- summary(model_friction)

# 3. Specific Regression for Hostility (as this had the most significant spike)
model_hostility_friction <- lm(diff_Hostility ~ diff_GO_Bijak + diff_Assertiveness, data = data_t1)
summary_hostility <- summary(model_hostility_friction)

cat("\n--- REGRESSION RESULTS: PREDICTORS OF AGGRESSION IN MONTH 1 (POST-TEST 1) ---\n")
cat("Dependent Variable: Change in Total Aggression (BPAQ)\n")
print(coef(summary_friction))

cat("\nDependent Variable: Change in Hostility (Cognitive Facet of Aggression)\n")
print(coef(summary_hostility))

# 4. Simple Correlation to strengthen the narrative
cor_val <- cor(data_t1$diff_GO_Bijak, data_t1$diff_Total_BPAQ)
cat(sprintf("\nCorrelation between Wise Gaming Improvement & Aggression Spike in Month 1: r = %.3f\n", cor_val))

# ==========================================
# STAGE END: EXPORT AND PRINT TO CONSOLE
# ==========================================
cat("\n[Info] Finalizing output... Printing to console and saving to file.\n")

# Using split = TRUE so output appears in the RStudio Console AND is saved to the .txt file
sink("Publication_Results_Summary.txt", split = TRUE)

cat("\n=======================================================\n")
cat("          STATISTICAL SUMMARY FOR PUBLICATION          \n")
cat("=======================================================\n\n")

cat("--- 0. Baseline Characteristics (Table 1) ---\n")
print(table_baseline, row.names = FALSE)
cat("Note: P-values > 0.05 indicate successful randomization / baseline equivalence between groups.\n")

cat("\n--- 1. Assumption Checks (Normality & Homogeneity) ---\n")
print(table_assumptions, row.names = FALSE)

cat("\n--- 2. Mann-Whitney U Test (Overall Change) ---\n")
print(table_mwu, row.names = FALSE)

if(exists("cor_gas")) {
  cat("\n--- 3. Spearman Rank Correlation (Gaming Duration vs Baseline Outcomes) ---\n")
  cat(sprintf("   Gaming Duration vs Total GAS  : rho = %.3f, p-value = %.4f %s\n", cor_gas$estimate, cor_gas$p.value, ifelse(cor_gas$p.value < 0.05, "*", "")))
  cat(sprintf("   Gaming Duration vs Total BPAQ : rho = %.3f, p-value = %.4f %s\n", cor_bpaq$estimate, cor_bpaq$p.value, ifelse(cor_bpaq$p.value < 0.05, "*", "")))
  cat(sprintf("   Gaming Duration vs Total SDQ  : rho = %.3f, p-value = %.4f %s\n", cor_sdq$estimate, cor_sdq$p.value, ifelse(cor_sdq$p.value < 0.05, "*", "")))
  cat(sprintf("   Gaming Duration vs Wise Gaming: rho = %.3f, p-value = %.4f %s\n", cor_gobijak$estimate, cor_gobijak$p.value, ifelse(cor_gobijak$p.value < 0.05, "*", "")))
}

cat("\n--- 3.5 Specific Baseline Correlations for Manuscript ---\n")
if(exists("table_manuscript_cors")) {
  print(table_manuscript_cors, row.names = FALSE)
} else {
  cat("Data not available.\n")
}

cat("\n--- 4. Linear Mixed-Effects Models (LMM) Summary - Primary ---\n")
print(table_lmm, row.names = FALSE)
cat("Note: Significant 'Intervention p-val' indicates the differential effect of the intervention over time compared to control.\n")

cat("\n--- 5. Linear Mixed-Effects Models (LMM) Summary - Subscales ---\n")
print(table_lmm_sub, row.names = FALSE)

cat("\n--- 6. Specific Timepoint Differences (Primary Outcomes / Forest Plot Data) ---\n")
forest_out <- forest_data %>% dplyr::select(Variable, Timepoint, Hedges_g, CI_low, CI_high, p_value, Direction) %>% mutate(across(c(Hedges_g, CI_low, CI_high, p_value), ~round(., 4)))
print(forest_out, row.names = FALSE)
cat("Note: Hedges' g effect size comparing Intervention vs Control at specific timepoints.\n")

cat("\n--- 7. Causal Mediation Analysis Results (Primary) ---\n")
print(table_mediation, row.names = FALSE)
cat("Note: ACME represents the indirect effect. A significant ACME (p < 0.05) indicates that the variable successfully mediates the intervention's effect.\n")

cat("\n--- 8. Causal Mediation Analysis Results (Wise Gaming Sub-Dimensions) ---\n")
if(exists("table_sub_mediation")) {
  print(table_sub_mediation, row.names = FALSE)
} else {
  cat("Wise Gaming sub-dimension mediation data not available.\n")
}

cat("\n--- 9. Detailed Changes in Aggression Sub-variables (BPAQ - Synchronized) ---\n")
if(exists("table_bpaq_sub_final")) {
  print(table_bpaq_sub_final, row.names = FALSE)
} else {
  cat("BPAQ sub-variable data not available.\n")
}

cat("\n--- 10. Average Item Score Increase for Wise Gaming (Intervention Group) ---\n")
if(exists("table_item_change")) {
  print(table_item_change, row.names = FALSE)
} else {
  cat("Item-level change data not available.\n")
}

cat("\n--- 11. Initial Recovery Friction Analysis (Post-test 1) ---\n")
if(exists("summary_friction")) {
  cat("Dependent Variable: Change in Total Aggression (BPAQ)\n")
  print(coef(summary_friction))
  cat("\nDependent Variable: Change in Hostility (Cognitive Facet of Aggression)\n")
  print(coef(summary_hostility))
  cat(sprintf("\nCorrelation between Wise Gaming Improvement & Aggression Spike in Month 1: r = %.3f\n", cor_val))
} else {
  cat("Friction analysis data not available.\n")
}

cat("\n=======================================================\n")
cat("          ALL ANALYSES COMPLETED SUCCESSFULLY          \n")
cat("=======================================================\n")

sink()
