# =========================================================
# NHANES mini-RWE project:
# Treatment gaps and residual non-HDL risk in high-risk adults
# Cycles: 2011-2012, 2013-2014, 2015-2016, 2017-2018
# =========================================================

# -----------------------------
# 0) Packages
# -----------------------------
packages <- c("tidyverse", "haven", "janitor", "survey", "srvyr",
  "stringr", "forcats", "scales", "broom")

to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(haven)
library(janitor)
library(survey)
library(srvyr)
library(stringr)
library(forcats)
library(scales)
library(broom)

options(survey.lonely.psu = "adjust")

# -----------------------------
# 1) User inputs
# -----------------------------
DATA_DIR <- "C:/Users/Emerald/Documents/Career_Projects"
OUT_DIR  <- "C:/Users/Emerald/Documents/Career_Projects/output"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cycle_map <- c(G = "2011-2012",H = "2013-2014",I = "2015-2016",
  J = "2017-2018")

# -----------------------------
# 2) Helper functions
# -----------------------------
read_cycle_file <- function(stem, suffix, data_dir = DATA_DIR) {
  # Tries both upper and lower case extensions
  candidates <- c(
    file.path(data_dir, paste0(stem, "_", suffix, ".XPT")),
    file.path(data_dir, paste0(stem, "_", suffix, ".xpt"))
  )
  
  path <- candidates[file.exists(candidates)][1]
  
  if (is.na(path)) {
    stop("File not found for: ", stem, "_", suffix)
  }
  
  read_xpt(path) |>
    clean_names() |>
    mutate(cycle = cycle_map[[suffix]])
}

read_multi_cycle <- function(stem) {
  map_dfr(names(cycle_map), ~ read_cycle_file(stem = stem, suffix = .x))
}

pick_first_existing <- function(df, candidates) {
  nm <- names(df)
  found <- candidates[candidates %in% nm]
  if (length(found) == 0) return(NULL)
  found[1]
}

coalesce_existing <- function(df, candidates) {
  existing <- candidates[candidates %in% names(df)]
  if (length(existing) == 0) return(rep(NA, nrow(df)))
  out <- df[[existing[1]]]
  if (length(existing) > 1) {
    for (v in existing[-1]) {
      out <- dplyr::coalesce(out, df[[v]])
    }
  }
  out
}

weighted_prev <- function(design, var) {
  f <- as.formula(paste0("~", var))
  svymean(f, design, na.rm = TRUE) |>
    broom::tidy()
}

# -----------------------------
# 3) Read NHANES files
# -----------------------------
demo   <- read_multi_cycle("DEMO")
bpx    <- read_multi_cycle("BPX")
bmx    <- read_multi_cycle("BMX")
tchol  <- read_multi_cycle("TCHOL")
hdl    <- read_multi_cycle("HDL")
diq    <- read_multi_cycle("DIQ")
hiq    <- read_multi_cycle("HIQ")
huq    <- read_multi_cycle("HUQ")
smq    <- read_multi_cycle("SMQ")
rxq_rx <- read_multi_cycle("RXQ_RX")

# Optional: not required for main analysis, but kept if you want it later
# rxq_drug is not cycle-suffixed in 2017-2018 documentation; skip unless needed
# rxq_drug <- read_xpt(file.path(DATA_DIR, "RXQ_DRUG.XPT")) |> clean_names()

# -----------------------------
# 4) Medication classification
# -----------------------------
# We classify statins / lipid-lowering therapy from text in RXQ_RX.
# This avoids fragile joins and works well for a portfolio project.

# Identify likely text columns
char_cols <- names(rxq_rx)[vapply(rxq_rx, is.character, logical(1))]

rxq_rx <- rxq_rx |>
  mutate(
    med_text = pmap_chr(across(all_of(char_cols)), ~ {
      vals <- c(...)
      vals <- vals[!is.na(vals)]
      paste(vals, collapse = " | ")
    }),
    med_text_lower = str_to_lower(med_text),
    
    any_lipid_lowering = str_detect(
      med_text_lower,
      "statin|atorvastatin|simvastatin|rosuvastatin|pravastatin|lovastatin|fluvastatin|pitavastatin|ezetimibe|fenofibrate|gemfibrozil|niacin|pcsk9|alirocumab|evolocumab|bempedoic"
    ),
    
    statin_use = str_detect(
      med_text_lower,
      "atorvastatin|simvastatin|rosuvastatin|pravastatin|lovastatin|fluvastatin|pitavastatin|statin"
    ),
    
    nonstatin_llt = str_detect(
      med_text_lower,
      "ezetimibe|fenofibrate|gemfibrozil|niacin|pcsk9|alirocumab|evolocumab|bempedoic"
    )
  )

rx_person <- rxq_rx |>
  group_by(seqn, cycle) |>
  summarise(
    any_lipid_lowering = as.integer(any(any_lipid_lowering, na.rm = TRUE)),
    statin_use         = as.integer(any(statin_use, na.rm = TRUE)),
    nonstatin_llt      = as.integer(any(nonstatin_llt, na.rm = TRUE)),
    n_reported_meds    = n(),
    .groups = "drop")

# -----------------------------
# 5) Harmonise core person-level files
# -----------------------------
demo_h <- demo |>
  transmute(
    seqn,
    cycle,
    age = ridageyr,
    sex = factor(riagendr, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth = case_when(
      ridreth3 == 1 ~ "Mexican American",
      ridreth3 == 2 ~ "Other Hispanic",
      ridreth3 == 3 ~ "Non-Hispanic White",
      ridreth3 == 4 ~ "Non-Hispanic Black",
      ridreth3 == 6 ~ "Non-Hispanic Asian",
      ridreth3 == 7 ~ "Other / Multiracial",
      TRUE ~ NA_character_
    ),
    pir = indfmpir,
    education = case_when(
      dmdeduc2 == 1 ~ "<9th grade",
      dmdeduc2 == 2 ~ "9-11th grade",
      dmdeduc2 == 3 ~ "High school/GED",
      dmdeduc2 == 4 ~ "Some college/AA",
      dmdeduc2 == 5 ~ "College graduate+",
      TRUE ~ NA_character_
    ),
    sdmvpsu,
    sdmvstra,
    wtmec2yr
  ) |>
  mutate(
    race_eth = factor(race_eth),
    education = factor(education))

bpx_h <- bpx |>
  mutate(
    sbp = rowMeans(dplyr::pick(starts_with("bpxsy")), na.rm = TRUE),
    dbp = rowMeans(dplyr::pick(starts_with("bpxdi")), na.rm = TRUE)
  ) |>
  transmute(
    seqn,
    cycle,
    sbp = ifelse(is.nan(sbp), NA, sbp),
    dbp = ifelse(is.nan(dbp), NA, dbp))

bmx_h <- bmx |>
  transmute(
    seqn,
    cycle,
    bmi = bmxbmi)

tchol_h <- tchol |>
  transmute(
    seqn,
    cycle,
    total_chol = lbxtc)

hdl_h <- hdl |>
  transmute(
    seqn,
    cycle,
    hdl = lbdhdd)

insulin_var <- "diq050" %in% names(diq)
oral_var    <- "diq070" %in% names(diq)

diq_h <- diq

diq_h$diabetes_dx <- dplyr::case_when(
  diq_h$diq010 == 1 ~ 1,
  diq_h$diq010 %in% c(2, 3) ~ 0,
  TRUE ~ NA_real_)

diq_h$insulin_use <- if ("diq050" %in% names(diq_h)) {
  as.integer(diq_h$diq050 == 1)
} else {
  rep(NA_integer_, nrow(diq_h))
}

diq_h$oral_dm_rx <- if ("diq070" %in% names(diq_h)) {
  as.integer(diq_h$diq070 == 1)
} else {
  rep(NA_integer_, nrow(diq_h))
}

diq_h <- diq_h |>
  transmute(
    seqn,
    cycle,
    diabetes_dx,
    insulin_use,
    oral_dm_rx)

hiq_h <- hiq |>
  transmute(
    seqn,
    cycle,
    insured = case_when(
      hiq011 == 1 ~ 1,
      hiq011 == 2 ~ 0,
      TRUE ~ NA_real_))

huq_h <- huq

huq_h$routine_place <- if ("huq030" %in% names(huq_h)) {
  dplyr::case_when(
    huq_h$huq030 == 1 ~ 1,
    huq_h$huq030 == 2 ~ 0,
    TRUE ~ NA_real_)
} else {
  rep(NA_real_, nrow(huq_h))
}

visit_var <- c("huq090", "huq100", "huq071")
visit_var <- visit_var[visit_var %in% names(huq_h)]

huq_h$n_healthcare_visits <- if (length(visit_var) > 0) {
  huq_h[[visit_var[1]]]
} else {
  rep(NA_real_, nrow(huq_h))
}

huq_h <- huq_h |>
  transmute(
    seqn,
    cycle,
    routine_place,
    n_healthcare_visits)

smq_h <- smq |>
  transmute(
    seqn,
    cycle,
    ever_100_cigs = case_when(
      smq020 == 1 ~ 1,
      smq020 == 2 ~ 0,
      TRUE ~ NA_real_),
    current_smoker = case_when(
      smq040 %in% c(1, 2) ~ 1,
      smq040 == 3 ~ 0,
      TRUE ~ NA_real_))

# -----------------------------
# 6) Merge person-level analytic dataset
# -----------------------------
analytic <- demo_h |>
  left_join(bpx_h,   by = c("seqn", "cycle")) |>
  left_join(bmx_h,   by = c("seqn", "cycle")) |>
  left_join(tchol_h, by = c("seqn", "cycle")) |>
  left_join(hdl_h,   by = c("seqn", "cycle")) |>
  left_join(diq_h,   by = c("seqn", "cycle")) |>
  left_join(hiq_h,   by = c("seqn", "cycle")) |>
  left_join(huq_h,   by = c("seqn", "cycle")) |>
  left_join(smq_h,   by = c("seqn", "cycle")) |>
  left_join(rx_person, by = c("seqn", "cycle")) |>
  mutate(
    any_lipid_lowering = replace_na(any_lipid_lowering, 0L),
    statin_use         = replace_na(statin_use, 0L),
    nonstatin_llt      = replace_na(nonstatin_llt, 0L),
    n_reported_meds    = replace_na(n_reported_meds, 0L),
    
    non_hdl = total_chol - hdl,
    
    obese = if_else(!is.na(bmi) & bmi >= 30, 1, 0, missing = NA_real_),
    hypertensive = if_else(!is.na(sbp) & sbp >= 130 | !is.na(dbp) & dbp >= 80, 1, 0, missing = NA_real_),
    
    # High-risk cohort for this portfolio project:
    # age 40+ AND at least one major cardiometabolic risk marker
    high_risk = if_else(
      age >= 40 &
        (
          diabetes_dx == 1 |
            hypertensive == 1 |
            obese == 1 |
            current_smoker == 1),
      1, 0, missing = NA_real_
    ),
    
    # Residual-risk threshold using non-HDL cholesterol
    non_hdl_uncontrolled = if_else(!is.na(non_hdl) & non_hdl >= 130, 1, 0, missing = NA_real_),
    
    untreated_high_risk = if_else(high_risk == 1 & any_lipid_lowering == 0, 1, 0, missing = NA_real_),
    undertreated_high_risk = if_else(high_risk == 1 & statin_use == 0, 1, 0, missing = NA_real_),
    
    residual_risk_on_treatment = if_else(high_risk == 1 & any_lipid_lowering == 1 & non_hdl_uncontrolled == 1, 1, 0, missing = NA_real_),
    
    cycle = factor(cycle, levels = unname(cycle_map)))

# -----------------------------
# 7) Restrict to high-risk adults with core data
# -----------------------------
analytic_hr <- analytic |>
  filter(
    high_risk == 1,
    !is.na(wtmec2yr),
    !is.na(sdmvpsu),
    !is.na(sdmvstra),
    !is.na(non_hdl)
  ) |>
  mutate(
    wtmec8yr = wtmec2yr / 4,
    age_band = cut(
      age,
      breaks = c(40, 50, 60, 70, Inf),
      right = FALSE,
      labels = c("40-49", "50-59", "60-69", "70+")
    ),
    bmi_group = case_when(
      bmi < 25 ~ "Normal/Underweight",
      bmi < 30 ~ "Overweight",
      bmi >= 30 ~ "Obesity",
      TRUE ~ NA_character_
    ),
    insured_f = factor(if_else(insured == 1, "Insured", "Uninsured", missing = "Unknown")),
    routine_place_f = factor(if_else(routine_place == 1, "Usual source of care", "No usual source", missing = "Unknown"))
  )

# -----------------------------
# 8) Survey design
# -----------------------------
nhanes_design <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec8yr,
  nest = TRUE,
  data = analytic_hr)

nhanes_srvyr <- as_survey_design(
  analytic_hr,
  ids = sdmvpsu,
  strata = sdmvstra,
  weights = wtmec8yr,
  nest = TRUE)

# -----------------------------
# 9) Weighted descriptive results
# -----------------------------
weighted_prev <- function(design, var) {
  f <- as.formula(paste0("~", var))
  est <- svymean(f, design, na.rm = TRUE)
  ci  <- confint(est)
  
  tibble(
    term = var,
    estimate = as.numeric(est),
    std.error = as.numeric(SE(est)),
    conf.low = ci[1],
    conf.high = ci[2]
  )
}

overall_key <- bind_rows(
  weighted_prev(nhanes_design, "any_lipid_lowering") |> mutate(metric = "Any lipid-lowering therapy"),
  weighted_prev(nhanes_design, "statin_use") |> mutate(metric = "Statin use"),
  weighted_prev(nhanes_design, "untreated_high_risk") |> mutate(metric = "Untreated among high-risk"),
  weighted_prev(nhanes_design, "undertreated_high_risk") |> mutate(metric = "No statin among high-risk"),
  weighted_prev(nhanes_design, "residual_risk_on_treatment") |> mutate(metric = "Residual non-HDL risk on treatment"),
  weighted_prev(nhanes_design, "non_hdl_uncontrolled") |> mutate(metric = "Non-HDL >=130 mg/dL")
) |>
  select(metric, estimate, std.error, conf.low, conf.high) |>
  mutate(
    across(c(estimate, conf.low, conf.high), ~ scales::percent(.x, accuracy = 0.1))
  )

write.csv(overall_key, file.path(OUT_DIR, "overall_weighted_prevalence.csv"), row.names = FALSE)

# By subgroup
subgroup_table <- nhanes_srvyr |>
  group_by(age_band, sex, race_eth) |>
  summarise(
    untreated = survey_mean(untreated_high_risk, vartype = c("ci"), na.rm = TRUE),
    statin_use = survey_mean(statin_use, vartype = c("ci"), na.rm = TRUE),
    residual_risk = survey_mean(residual_risk_on_treatment, vartype = c("ci"), na.rm = TRUE),
    .groups = "drop")

write.csv(subgroup_table, file.path(OUT_DIR, "subgroup_weighted_results.csv"), row.names = FALSE)

# -----------------------------
# 10) Weighted regression models
# -----------------------------
# Model 1: undertreatment (no statin) among high-risk adults
m1 <- svyglm(
  undertreated_high_risk ~ age_band + sex + race_eth + insured_f + routine_place_f +
    diabetes_dx + hypertensive + obese + current_smoker,
  design = nhanes_design,
  family = quasibinomial())

m1_or <- broom::tidy(m1, conf.int = TRUE, exponentiate = TRUE)
write.csv(m1_or, file.path(OUT_DIR, "model1_undertreatment_or.csv"), row.names = FALSE)

# Model 2: uncontrolled non-HDL among those on any lipid-lowering therapy
treated_design <- subset(nhanes_design, any_lipid_lowering == 1)

m2 <- svyglm(
  non_hdl_uncontrolled ~ age_band + sex + race_eth + insured_f + routine_place_f +
    diabetes_dx + hypertensive + obese + current_smoker + statin_use,
  design = treated_design,
  family = quasibinomial())

m2_or <- broom::tidy(m2, conf.int = TRUE, exponentiate = TRUE)
write.csv(m2_or, file.path(OUT_DIR, "model2_nonhdl_uncontrolled_or.csv"), row.names = FALSE)

# -----------------------------
# 11) Simple figures
# -----------------------------
# Figure 1: overall weighted prevalence
fig1_dat <- overall_key |>
  mutate(
    estimate_num = readr::parse_number(estimate) / 100,
    metric = fct_reorder(metric, estimate_num))

p1 <- fig1_dat |>
  ggplot(aes(x = estimate_num, y = metric)) +
  geom_col() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Weighted prevalence of treatment gaps and residual risk",
    x = "Weighted prevalence",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "fig1_overall_prevalence.png"), p1, width = 9, height = 5, dpi = 300)

# Figure 2: undertreatment by subgroup
fig2_dat <- nhanes_srvyr |>
  group_by(race_eth) |>
  summarise(
    undertreated = survey_mean(undertreated_high_risk, vartype = c("ci"), na.rm = TRUE),
    .groups = "drop")

p2 <- fig2_dat |>
  ggplot(aes(x = fct_reorder(race_eth, undertreated), y = undertreated)) +
  geom_col() +
  geom_errorbar(aes(ymin = undertreated_low, ymax = undertreated_upp), width = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Undertreatment among high-risk adults by race/ethnicity",
    x = NULL,
    y = "Weighted prevalence"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "fig2_undertreatment_by_race.png"), p2, width = 8, height = 5, dpi = 300)

# Figure 3: residual risk among treated by insurance
fig3_dat <- nhanes_srvyr |>
  filter(any_lipid_lowering == 1) |>
  group_by(insured_f) |>
  summarise(
    residual = survey_mean(non_hdl_uncontrolled, vartype = c("ci"), na.rm = TRUE),
    .groups = "drop")

p3 <- fig3_dat |>
  ggplot(aes(x = insured_f, y = residual)) +
  geom_col() +
  geom_errorbar(aes(ymin = residual_low, ymax = residual_upp), width = 0.15) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Residual non-HDL risk among treated high-risk adults by insurance",
    x = NULL,
    y = "Weighted prevalence"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "fig3_residual_risk_by_insurance.png"), p3, width = 8, height = 5, dpi = 300)

# -----------------------------
# 12) Save analytic dataset
# -----------------------------
write.csv(analytic_hr, file.path(OUT_DIR, "analytic_high_risk_dataset.csv"), row.names = FALSE)

# -----------------------------
# 13) Quick console output
# -----------------------------
cat("\nAnalysis complete.\n")
cat("Rows in high-risk analytic sample:", nrow(analytic_hr), "\n")
cat("Outputs saved to:", OUT_DIR, "\n")

print(overall_key)
print(summary(m1))
print(summary(m2))