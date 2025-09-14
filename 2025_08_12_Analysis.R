## =================== HIV Lifetime Screening: Intervention-Based Analysis (CDC age bins) ===================

## ---- Libraries ----
library(dplyr)
library(tidyr)
library(lubridate)
library(binom)
library(broom)
library(forcats)
library(ggplot2)
library(scales)
library(purrr)
library(stringr)
library(DiagrammeR)
library(DiagrammeRsvg)


suppressPackageStartupMessages({
  if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
  if (!requireNamespace("performance", quietly = TRUE)) install.packages("performance")
  if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
  if (!requireNamespace("gtsummary", quietly = TRUE)) install.packages("gtsummary")
  if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
  if (!requireNamespace("corrplot", quietly = TRUE)) install.packages("corrplot")
})
library(car)
library(performance)
library(emmeans)
library(gtsummary)
library(gt)
library(corrplot)

#### Load cleaned data ####
load("R File/HIP_Cleaned.RData")

### Helper
table_na <- function(..., useNA = "ifany") table(..., useNA = useNA)

### Exclusion Criteria 
# Start with total
  n_total <- nrow(HIP)
  print(n_total)
    n_minors <- 13
    n_known_hiv <- 540
    n_wrong_team <- 1789
    n_missing_demo <- 8
    rem_after_minors   <- 9911
    rem_after_hiv      <- 9371
    rem_after_team     <- 7582
    final_n            <- 7574
    pre_n  <- 1228  # set to 1128 if that's your final pre-intervention N
    edu_n  <- 3900
    emr_n  <- 2446
    
    # ---- 3) Build a DiagrammeR flow (vertical, labeled arrows) ----
      label_box <- function(title, n) {
        paste0(title, "\\nN = ", comma(n))
      }
      
    
    g <- grViz(sprintf('
      digraph flow {
        graph [rankdir = TB, nodesep = 0.5, splines = ortho, fontname = "Helvetica"]
        node  [shape = box, style = "rounded", fontsize = 12, width = 4, penwidth = 1, fontname = "Helvetica"]
        edge  [fontsize = 12]
      
        // Main vertical flow
        a [label = "%s"];
        b [label = "%s"];
        c [label = "%s"];
        d [label = "%s"];
        e [label = "%s"];
      
        a -> b [label = "Excluded: Age <18 (n = %s)", fontname = "Helvetica"];
        b -> c [label = "Excluded: Known HIV prior to admission (n = %s)", fontname = "Helvetica"];
        c -> d [label = "Excluded: Non-eligible primary team (n = %s)", fontname = "Helvetica"];
        d -> e [label = "Excluded: Missing key demographics (n = %s)", fontname = "Helvetica"];
      
        // Bottom layer: split of final analytic cohort by intervention period
        subgraph cluster_bottom {
          rank = same;
          style = invis;
      
          f1 [label = "Pre-intervention\\n(Jul 1–Sep 30, 2020)\\nN = %s", width=4];
          f2 [label = "Education-only intervention\\n(Oct 1, 2020–Jun 30, 2021)\\nN = %s", width=4];
          f3 [label = "Education + EMR intervention\\n(Jul 1–Dec 31, 2021)\\nN = %s", width=4];
        }
      
        e -> f1;
        e -> f2;
        e -> f3;
      }
      ',
      # main boxes
      label_box("Admissions identified (Jul 1, 2020–Jun 30, 2021)", n_total),
      label_box("After minors excluded",                              rem_after_minors),
      label_box("After known HIV excluded",                           rem_after_hiv),
      label_box("After non-eligible team excluded",                   rem_after_team),
      label_box("Final analytic cohort",                              final_n),
      
      # exclusion counts
      comma(n_minors), comma(n_known_hiv), comma(n_wrong_team), comma(n_missing_demo),
      
      # bottom-layer counts
      comma(pre_n), comma(edu_n), comma(emr_n)
          ))
    
          g

          
          g2 <- grViz(sprintf('
            digraph flow {
              graph [rankdir = TB, nodesep = 0.5, splines = ortho, fontname = "Helvetica"]
              node  [shape = box, style = "rounded", fontsize = 12, width = 4, penwidth = 1, fontname = "Helvetica"]
            
              a [label = "%s"];
              b [label = "%s"];
              c [label = "%s"];
              d [label = "%s"];
              e [label = "%s"];
            
              edge [fontsize = 10]
              a -> b [label = "Excluded: Age <18 (n = %s)", fontname = "Helvetica", fontsize=12];
              b -> c [label = "Excluded: Known HIV prior to admission (n = %s)", fontname = "Helvetica", fontsize=12];
              c -> d [label = "Excluded: Non-eligible primary team (n = %s)", fontname = "Helvetica", fontsize=12];
              d -> e [label = "Excluded: Missing key demographics (n = %s)", fontname = "Helvetica", fontsize=12];
            }
            ',
            label_box("Admissions identified (Jul 1, 2020–Jun 30, 2021)", n_total),
            label_box("After minors excluded",                              rem_after_minors),
            label_box("After known HIV excluded",                           rem_after_hiv),
            label_box("After non-eligible team excluded",                   rem_after_team),
            label_box("Final analytic cohort",                              final_n),
            comma(n_minors), comma(n_known_hiv), comma(n_wrong_team), comma(n_missing_demo)
                      ))
          
          # Render in RStudio/Quarto viewer
          g2

## Exclude minors
n_before <- nrow(HIP)
  table_na(HIP$age)
    HIP <- HIP %>% filter(age >= 18) 
      table_na(HIP$age)
        n_after <- nrow(HIP)
          print(n_before - n_after)
  
#Exclude prior history of hiv
n_before <- nrow(HIP)
  table_na(HIP$HIV_Hx)
    HIP <- HIP %>% filter(HIV_Hx == "1")
      table_na(HIP$HIV_Hx)
        n_after <- nrow(HIP)
          print(n_before - n_after)

#Eclude wrong team 
n_before <- nrow(HIP)
  table_na(HIP$team_5)
    HIP <- HIP %>% filter(!!HIP$team_5 %in% c("Resident Medicine","Faculty Medicine","Cardiology"))
      table_na(HIP$team_5)
        n_after <- nrow(HIP)
          print(n_before - n_after)

#Missing Values 
n_before <- nrow(HIP)
  table_na(is.na(HIP$sex) | is.na(HIP$race) | is.na(HIP$language_major) | is.na(HIP$housing))
    HIP <- HIP %>% filter(! (is.na(HIP$sex) | is.na(HIP$race) | is.na(HIP$language_major) | is.na(HIP$housing)) )
      table_na(is.na(HIP$sex) | is.na(HIP$race) | is.na(HIP$language_major) | is.na(HIP$housing))
        n_after <- nrow(HIP)
          print(n_before - n_after)
          
### Creating diagram 
          
table(HIP$HIV_lifetime)

## Analysis dataset: teams 1–3 only
table_na(HIP$team_5)
dat <- HIP %>%
  filter(team_5 %in% levels(team_5)[1:3]) %>%
  mutate(
    team3    = fct_drop(team_5),
    month    = floor_date(admit_date, "month"),
    screened = as.integer(HIV_lifetime == "Yes"),
    # --- CDC-era age categories ---
    age_cat  = cut(
      age,
      breaks = c(-Inf, 25, 35, 45, 55, 65, Inf), right = FALSE,
      labels = c("≤24","25–34","35–44","45–54","55–64","≥65")
    ),
    age_cat  = fct_relevel(age_cat, "25–34") # reference
  ) 
  table_na(dat$team_5)

stopifnot(all(c("Intervention_Group","team3","month","screened","age_cat") %in% names(dat)))

# ----------------------------- PRIMARY: by Intervention -----------------------------
intervention_overall <- dat %>%
  group_by(Intervention_Group) %>%
  summarise(n = n(), y = sum(screened), .groups = "drop") %>%
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper),
    pct   = 100*mean,
    ci95  = sprintf("%.1f%% (%.1f–%.1f)", 100*mean, 100*lower, 100*upper)
  ) %>%
  select(Intervention_Group, n, y, pct, ci95)

intervention_overall

intervention_overall_nogroup <- dat %>%
  summarise(n = n(), y = sum(screened), .groups = "drop") %>%
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper),
    pct   = 100*mean,
    ci95  = sprintf("%.1f%% (%.1f–%.1f)", 100*mean, 100*lower, 100*upper)
  ) 

intervention_overall_nogroup


# Global chi-square
tab_int <- dat %>%
  count(Intervention_Group, screened = if_else(screened==1L, "Yes", "No")) %>%
  pivot_wider(names_from = screened, values_from = n, values_fill = 0)
chisq_int <- chisq.test(as.matrix(data.frame(row.names = tab_int$Intervention_Group, tab_int[, -1])))
chisq_int

# Pairwise (Holm) on unadjusted proportions
pairwise_int <- with(intervention_overall,
                     pairwise.prop.test(x = y, n = n, p.adjust.method = "holm"))
pairwise_int


# Recompute CIs to get numeric bounds for plotting
io <- intervention_overall %>%
  mutate(
    prop = y / n,
    ci = binom::binom.confint(y, n, methods = "wilson"),
    pct  = percent(prop, accuracy = 0.1),
    lcl  = ci$lower,
    ucl  = ci$upper,
    pct_ci = sprintf("%s (%s–%s)",
                     pct,
                     percent(lcl, accuracy = 0.1),
                     percent(ucl, accuracy = 0.1))
  )

## --------- 1) Publication-ready table with gt ---------
tab <- io %>%
  transmute(
    `Intervention period` = as.character(Intervention_Group),
    `Screened (n/N)`      = sprintf("%s / %s", y, n),
    `% (95% CI)`          = pct_ci
  ) %>%
  gt() %>%
  tab_header(title = md("**Lifetime HIV screening by intervention period**")) %>%
  cols_align(align = "right", columns = c(`Screened (n/N)`, `% (95% CI)`)) %>%
  tab_source_note(md("Global χ²(2) = 590.15, p < 0.001; all pairwise Holm-adjusted p < 0.001."))

tab

## --------- 2) Bar plot with 95% CIs ---------
# Order periods as they occurred
io$Intervention_Group <- factor(io$Intervention_Group,
                                levels = c("Pre-intervention","Education Intervention","Education and EMR Intervention")
)

ggplot(io, aes(x = Intervention_Group, y = prop)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.15, size = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Lifetime HIV screening (%)",
    title = "Screening increased stepwise across intervention phases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# ---------------------------- SECONDARY: by Team & Intervention ----------------------------
intervention_by_team <- dat %>%
  group_by(Intervention_Group, team3) %>%
  summarise(n = n(), y = sum(screened), .groups = "drop") %>%
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper),
    pct   = 100 * mean,
    ci95  = sprintf("%.1f%% (%.1f–%.1f)", 100*mean, 100*lower, 100*upper)
  ) %>%
  select(Intervention_Group, team3, n, y, pct, ci95) %>%
  arrange(Intervention_Group, team3)

intervention_by_team

library(dplyr)
library(ggplot2)
library(binom)
library(scales)

# Rebuild summary WITH CI columns (mean, lower, upper)
intervention_by_team_ci <- dat %>%
  group_by(Intervention_Group, team3) %>%
  summarise(n = n(), y = sum(screened), .groups = "drop") %>%
  bind_cols(binom::binom.confint(x = .$y, n = .$n, methods = "wilson")[, c("mean","lower","upper")]) %>%
  mutate(
    Intervention_Group = factor(
      Intervention_Group,
      levels = c("Pre-intervention", "Education Intervention", "Education and EMR Intervention")
    )
  )

nejm_cols <- c(
  "Pre-intervention"                 = "#6C8EBF",  # muted blue
  "Education Intervention"           = "#91C46C",  # muted green
  "Education and EMR Intervention"   = "#E39D3D"   # muted orange
)

ggplot(intervention_by_team_ci,
       aes(x = team3, y = mean, fill = Intervention_Group)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.68, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8),
                width = 0.18, linewidth = 0.4) +
  # Add N labels inside bars
  geom_text(aes(label = paste0("N=", n)),
            position = position_dodge(width = 0.8),
            vjust = 5, color = "black", size = 4, family = "Helvetica") +
  scale_fill_manual(values = nejm_cols, guide = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1.05), expand = c(0.02, 0)) +
  labs(
    x = NULL,
    y = "Lifetime HIV Screening (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_blank(),
    axis.ticks.length = unit(3, "pt"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.key.size  = unit(10, "pt"),
    axis.text.x      = element_text(size = 12, color="black", family = "Helvetica", margin = margin(t = 6)),
    axis.text.y      = element_text(size = 12, color="black", family = "Helvetica", margin = margin(r = 6)),
    axis.title.y     = element_text(size = 12, color="black", family = "Helvetica"),
    legend.text      = element_text(size = 12, color="black", family = "Helvetica")
  )

## ================= PRE-MODEL PREP: collapse levels + collinearity =================
dat <- dat %>%
  mutate(
    sex3 = fct_collapse(sex,
                        Female = "Female",
                        Male   = "Male",
                        `Other/Unknown` = c("Nonbinary","Other","Unknown")),
    Intervention_Group = fct_relevel(Intervention_Group, "Pre-intervention"),
    team3              = fct_relevel(team3, "Resident Medicine"),
    sex3               = fct_relevel(sex3, "Female"),
    race               = fct_relevel(race, "White, Non-Hispanic"),
    homeless           = fct_relevel(homeless, "No")     # Housing removed from analysis
  )

# ================= CONSISTENT LABELING (single source of truth) =================
make_labeled <- function(df) {
  df %>%
    mutate(
      Period = fct_recode(Intervention_Group,
                          "Pre-intervention" = "Pre-intervention",
                          "Education"        = "Education Intervention",
                          "Education + EMR"  = "Education and EMR Intervention") %>%
        fct_relevel("Pre-intervention"),
      Team = fct_recode(team3,
                        "Medicine (resident)" = "Resident Medicine",
                        "Medicine (faculty)"  = "Faculty Medicine",
                        "Cardiology"          = "Cardiology") %>%
        fct_relevel("Medicine (resident)"),
      Age  = fct_relevel(age_cat, "25–34"),
      Sex  = fct_recode(sex3, "Other" = "Other/Unknown") %>% fct_relevel("Female"),
      Race = fct_recode(race,
                        "White (non-Hispanic)" = "White, Non-Hispanic",
                        "White (Hispanic)"     = "White, Hispanic",
                        "Black" = "Black", "Asian" = "Asian", "Other" = "Other") %>%
        fct_relevel("White (non-Hispanic)"),
      Homelessness = fct_relevel(homeless, "No")
    )
}

dat_model <- make_labeled(dat)

dat_first_model <- dat %>%
  arrange(MRN, admit_date) %>%
  group_by(MRN) %>% slice(1L) %>% ungroup() %>%
  make_labeled()

# ========================= ADJUSTED ANALYSES =====================
fit_adj       <- glm(screened ~ Period + Team + Age + Sex + Race + Homelessness,
                     data = dat_model, family = binomial())
fit_adj_first <- glm(screened ~ Period + Team + Age + Sex + Race + Homelessness,
                     data = dat_first_model, family = binomial())

# Optional: Period x Team interaction
fit_intx <- glm(screened ~ Period + Team + Age + Sex + Race + Homelessness + Period:Team,
                data = dat_model, family = binomial())
anova(fit_adj, fit_intx, test = "LRT")

fit_intx_first <- glm(screened ~ Period + Team + Age + Sex + Race + Homelessness + Period:Team,
                      data = dat_first_model, family = binomial())
anova(fit_adj_first, fit_intx_first, test = "LRT")

# ================== Pairwise Period table (optional; shows all three) ==================
pairwise_period <- function(fit) {
  em  <- emmeans(fit, specs = ~ Period)
  con <- contrast(em, method = "pairwise", adjust = "holm")
  sm  <- summary(con, infer = TRUE, type = "response")
  df  <- as.data.frame(sm)
  if (!"odds.ratio" %in% names(df)) {
    if ("ratio" %in% names(df)) df$odds.ratio <- df$ratio
    else if ("estimate" %in% names(df)) df$odds.ratio <- df$estimate
  }
  lcl_candidates <- c("lower.CL","asymp.LCL","LCL")
  ucl_candidates <- c("upper.CL","asymp.UCL","UCL")
  lcl_col <- intersect(lcl_candidates, names(df))[1]
  ucl_col <- intersect(ucl_candidates, names(df))[1]
  if (is.na(lcl_col) || is.na(ucl_col)) {
    ci <- as.data.frame(confint(con, type = "response", adjust = "holm"))
    if (!"odds.ratio" %in% names(ci)) {
      if ("ratio" %in% names(ci)) ci$odds.ratio <- ci$ratio
      else if ("estimate" %in% names(ci)) ci$odds.ratio <- ci$estimate
    }
    lcl_col_ci <- intersect(lcl_candidates, names(ci))[1]
    ucl_col_ci <- intersect(ucl_candidates, names(ci))[1]
    df <- dplyr::left_join(df, ci[, c("contrast","odds.ratio", lcl_col_ci, ucl_col_ci)],
                           by = "contrast", suffix = c("", ".ci"))
    df$CI_low  <- df[[lcl_col_ci]]
    df$CI_high <- df[[ucl_col_ci]]
  } else {
    df$CI_low  <- df[[lcl_col]]
    df$CI_high <- df[[ucl_col]]
  }
  dplyr::transmute(dplyr::as_tibble(df), contrast, OR = odds.ratio, CI_low, CI_high, p_adj = p.value)
}

# ================== Coefficient-based contrasts (robust to labels) ==================
# EMR vs Education within Period
period_emr_vs_edu_from_coef <- function(fit, group_levels = c("Period","Team","Age","Sex","Race","Homelessness")) {
  b  <- coef(fit); V <- vcov(fit)
  levs <- levels(model.frame(fit)$Period); ref <- levs[1]; ed <- levs[2]; emr <- levs[3]
  name_ed  <- paste0("Period", ed)
  name_emr <- paste0("Period", emr)
  
  idx_ed  <- match(name_ed,  names(b), nomatch = 0)
  idx_emr <- match(name_emr, names(b), nomatch = 0)
  
  b_ed  <- if (idx_ed  > 0) b[idx_ed]  else 0
  b_emr <- if (idx_emr > 0) b[idx_emr] else 0
  
  v_ed  <- if (idx_ed  > 0) V[idx_ed,  idx_ed]  else 0
  v_emr <- if (idx_emr > 0) V[idx_emr, idx_emr] else 0
  cov12 <- if (idx_ed  > 0 && idx_emr > 0) V[idx_emr, idx_ed] else 0
  
  logOR <- b_emr - b_ed
  se    <- sqrt(v_emr + v_ed - 2 * cov12)
  z     <- 1.96
  
  tibble::tibble(
    pretty    = "Education + EMR vs Education",
    estimate  = exp(logOR),
    conf.low  = exp(logOR - z*se),
    conf.high = exp(logOR + z*se),
    group     = factor("Period", levels = group_levels)
  )
}

# General: contrast two levels within a factor (e.g., Team)
contrast_levels_from_coef <- function(fit, factor_name, level_A, level_B, label,
                                      group_levels = c("Period","Team","Age","Sex","Race","Homelessness")) {
  b  <- coef(fit); V <- vcov(fit)
  name_A <- paste0(factor_name, level_A)
  name_B <- paste0(factor_name, level_B)
  
  idx_A <- match(name_A, names(b), nomatch = 0)
  idx_B <- match(name_B, names(b), nomatch = 0)
  
  b_A  <- if (idx_A > 0) b[idx_A] else 0
  b_B  <- if (idx_B > 0) b[idx_B] else 0
  
  v_A  <- if (idx_A > 0) V[idx_A, idx_A] else 0
  v_B  <- if (idx_B > 0) V[idx_B, idx_B] else 0
  covAB<- if (idx_A > 0 && idx_B > 0) V[idx_A, idx_B] else 0
  
  logOR <- b_A - b_B
  se    <- sqrt(v_A + v_B - 2 * covAB)
  z     <- 1.96
  
  tibble::tibble(
    pretty    = label,
    estimate  = exp(logOR),
    conf.low  = exp(logOR - z*se),
    conf.high = exp(logOR + z*se),
    group     = factor(factor_name, levels = group_levels)
  )
}

# ================= UNADJUSTED OVERALL (Period only) =================
fit_unadj <- glm(screened ~ Period, data = dat_model, family = binomial())

# ------------------------------ PLOTTING HELPERS ------------------------------
# Build tidy frame for forests from a fitted model
prep_forest <- function(fit,
                        var_order  = c("Period","Team","Age","Sex","Race","Homelessness"),
                        var_labels = c(
                          Period       = "Period",
                          Team         = "Team",
                          Age          = "Age",
                          Sex          = "Sex",
                          Race         = "Race/Ethnicity",
                          Homelessness = "Homelessness"
                        ),
                        add_emr_edu = TRUE,
                        add_team_pair = FALSE, # NEW
                        period_only = FALSE) {
  
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)")
  
  mf <- model.frame(fit)
  
  if (period_only) {
    levs <- levels(mf$Period); ref  <- levs[1]
    td <- td %>%
      dplyr::filter(grepl("^Period", term)) %>%
      dplyr::mutate(level = sub("^Period", "", term),
                    pretty = paste0(level, " vs ", ref),
                    group  = factor("Period", levels = "Period")) %>%
      dplyr::select(estimate, conf.low, conf.high, pretty, group)
  } else {
    refs <- setNames(
      sapply(var_order, function(v)
        if (v %in% names(mf) && is.factor(mf[[v]])) levels(mf[[v]])[1] else NA_character_),
      var_order
    )
    td <- td %>%
      dplyr::mutate(
        var = purrr::map_chr(term, function(tt) {
          hit <- var_order[startsWith(tt, var_order)]
          if (length(hit)) hit[[1]] else NA_character_
        })
      ) %>%
      dplyr::filter(!is.na(var)) %>%
      dplyr::mutate(level_raw = trimws(sub(paste0("^", var), "", term)))
    
    drop_prefixes <- var_order
    prefix_re <- paste0("^(", paste(drop_prefixes, collapse = "|"), ")")
    
    td <- td %>%
      dplyr::mutate(
        level_clean = stringr::str_trim(stringr::str_remove(level_raw, prefix_re)),
        pretty      = paste0(level_clean, " vs ", refs[var])
      ) %>%
      dplyr::mutate(
        pretty = dplyr::if_else(var == "Homelessness" & grepl("^Yes\\s*vs\\s*No$", pretty), "Homeless", pretty)
      ) %>%
      dplyr::mutate(group = factor(var, levels = var_order, labels = unname(var_labels[var_order]))) %>%
      dplyr::arrange(group) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(pretty = forcats::fct_inorder(pretty)) %>%
      dplyr::ungroup() %>%
      dplyr::select(estimate, conf.low, conf.high, pretty, group)
  }
  
  # Add EMR+Education vs Education (Period)
  if (add_emr_edu) {
    td <- dplyr::bind_rows(td, period_emr_vs_edu_from_coef(fit))
  }
  
  # NEW: Add Cardiology vs Medicine (faculty) (Team)
  if (!period_only && add_team_pair) {
    td <- dplyr::bind_rows(
      td,
      contrast_levels_from_coef(fit,
                                factor_name = "Team",
                                level_A = "Cardiology",
                                level_B = "Medicine (faculty)",
                                label   = "Cardiology vs Medicine (faculty)")
    )
  }
  
  td
}

# Forest plot with right-side OR column
plot_forest <- function(fit, title = "",
                        add_emr_edu = TRUE, add_team_pair = FALSE, period_only = FALSE) {
  df <- prep_forest(fit, add_emr_edu = add_emr_edu, add_team_pair = add_team_pair, period_only = period_only) %>%
    dplyr::mutate(or_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high))
  
  max_hi      <- max(df$conf.high[is.finite(df$conf.high)], na.rm = TRUE)
  x_min_graph <- 0.06
  x_max_graph <- max(4, max_hi) * 1.55
  x_right_lab <- max(4, max_hi) * 1.35
  
  ggplot(df, aes(x = estimate, y = pretty)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4, color = "black") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 0.6, color = "black") +
    geom_point(shape = 15, size = 2.6, color = "black") +
    geom_text(aes(x = x_right_lab, y = pretty, label = or_ci),
              inherit.aes = FALSE, hjust = 0, size = 3.5) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8),
      labels = scales::label_number(accuracy = 0.01),
      limits = c(x_min_graph, x_max_graph),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      x = "Odds ratio (log scale)", y = NULL, title = title
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", margin = margin(b = 8)),
      strip.background = element_blank(),
      strip.text.y     = element_blank(),
      panel.border     = element_blank(),
      axis.text.y      = element_text(size = 11),
      plot.margin      = margin(10, 140, 10, 16)
    ) +
    coord_cartesian(clip = "off")
}

# ----------------------- Draw forest plots -----------------------
p_forest_main   <- plot_forest(fit_adj,
                               add_emr_edu = TRUE, add_team_pair = TRUE, period_only = FALSE)

p_forest_first  <- plot_forest(fit_adj_first,
                               add_emr_edu = TRUE, add_team_pair = TRUE, period_only = FALSE)

p_forest_unadj  <- plot_forest(fit_unadj,
                               add_emr_edu = TRUE, add_team_pair = FALSE, period_only = TRUE)

p_forest_main
p_forest_first
p_forest_unadj

### Removing age less than 66
table_na(dat_model$Age)
fit_adj_lt66 <- glm(
  screened ~ Period + Team + Age + Sex + Race + Homelessness,
  data   = subset(dat_model, Age!= "≥65"),
  family = binomial()
)

fit_adj_first_lt66 <- glm(
  screened ~ Period + Team + Age + Sex + Race + Homelessness,
  data   = subset(dat_first_model, Age!= "≥65"),
  family = binomial()
)

# Optional: Period x Team interaction (Age < 66)
fit_intx_lt66 <- glm(
  screened ~ Period + Team + Age + Sex + Race + Homelessness + Period:Team,
  data   = subset(dat_model, Age!= "≥65"),
  family = binomial()
)
anova(fit_adj_lt66, fit_intx_lt66, test = "LRT")

fit_intx_first_lt66 <- glm(
  screened ~ Period + Team + Age + Sex + Race + Homelessness + Period:Team,
  data   = subset(dat_first_model, Age!= "≥65"),
  family = binomial()
)
anova(fit_adj_first_lt66, fit_intx_first_lt66, test = "LRT")

# ================== Pairwise Period table (Age < 66) ==================
pairwise_period_lt66_all   <- pairwise_period(fit_adj_lt66)
pairwise_period_lt66_first <- pairwise_period(fit_adj_first_lt66)

# ================= UNADJUSTED (Period only), Age < 66 =================
fit_unadj_lt66 <- glm(screened ~ Period,
                      data = subset(dat_model, Age!= "≥65"),
                      family = binomial())

# ----------------------- Forest plots (Age < 66) -----------------------
p_forest_main_lt66  <- plot_forest(fit_adj_lt66,
                                   add_emr_edu = TRUE, add_team_pair = TRUE, period_only = FALSE)

p_forest_first_lt66 <- plot_forest(fit_adj_first_lt66,
                                   add_emr_edu = TRUE, add_team_pair = TRUE, period_only = FALSE)

p_forest_unadj_lt66 <- plot_forest(fit_unadj_lt66,
                                   add_emr_edu = TRUE, add_team_pair = FALSE, period_only = TRUE)

# Print or save as needed
p_forest_main_lt66
p_forest_first_lt66
p_forest_unadj_lt66

# ---------- Table 1: Overall and by Intervention Period (descriptive; Housing can remain) ----------
library(gtsummary)
library(gt)
library(dplyr)
library(stringr)

#Changing sex 
  #Current 
    table_na(dat$sex3)

# ---------- Table 1: Overall and by Intervention Period ----------
dat_tbl1 <- dat %>%
  mutate(
    HIV_history = case_when(
      "hiv_Hx" %in% names(.) ~ as.factor(hiv_Hx),
      "HIV_Hx" %in% names(.) ~ factor(ifelse(HIV_Hx == 2, "Yes", "No"),
                                      levels = c("No","Yes")),
      TRUE ~ NA
    ),
    Intervention_Group = forcats::fct_relevel(
      Intervention_Group,
      "Pre-intervention", "Education Intervention", "Education and EMR Intervention"
    )
  )

tbl1 <-
  dat_tbl1 %>%
  select(Intervention_Group, age, age_cat, sex3,race,
         homeless, team3) %>%
  tbl_summary(
    by       = Intervention_Group,
    missing  = "no",
    type     = list(age ~ "continuous2", age_cat ~ "categorical"),
    statistic = list(
      age ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 0, # age in whole years
    label = list(
      age ~ "Age (years)",
      age_cat ~ "Age category (CDC)",
      sex3 ~ "Sex",
      race ~ "Race/Ethnicity",
      homeless ~ "Homeless",
      team3 ~ "Primary Team"
    )
  ) %>%
  add_overall(last = FALSE) %>%
  add_n() %>%
  add_p(
    test = list(
      age ~ "kruskal.test",
      age_cat ~ "chisq.test",
      all_categorical() ~ "chisq.test"
    ),
    test.args   = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 5000)),
    pvalue_fun  = ~ style_pvalue(.x, digits = 3)
  ) %>%
  bold_labels() %>%
  modify_header(
    label ~ md("**Baseline Characteristic**"),
    n ~ md("**N**"),
    all_stat_cols() ~ "**n (%) or median (Q1, Q3)**",
    p.value ~ md("**p-value**")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ md("**Intervention Period**")) %>%
  modify_column_alignment(
    columns = c(n, all_stat_cols(), p.value),
    align = "right"
  ) %>%
  modify_footnote(
    update = list(
      all_stat_cols() ~ "Percentages calculated within column.",
      p.value ~ "Kruskal–Wallis for age; Pearson χ² for categorical (with Monte Carlo simulation)."
    ),
    abbreviation = TRUE
  ) %>%
  as_gt() %>%
  # tidy the p-values and numbers
  fmt_number(columns = c(p.value), decimals = 3) %>%
  text_transform(
    locations = cells_body(columns = c(p.value)),
    fn = function(x) ifelse(suppressWarnings(as.numeric(x) < 0.001), "<0.001", x)
  ) %>%
  # visual polish
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(4),
    heading.align = "left",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.background.color = "#F6F7F9",
    column_labels.font.weight = "bold"
  ) %>%
  cols_width(
    label ~ pct(28),
    everything() ~ pct(12),
    p.value ~ px(90)
  ) %>%
  tab_header(
    title = md("**Table 1. Baseline Characteristics Overall and by Intervention Period**")
  ) %>%
  tab_source_note(md("*Values are **n (%)** unless noted. Medians reported with **Q1, Q3**.*"))

tbl1

    ### Doing my own test 
      #Age categorical age_cat by Intervention Group using chisq.test
        table_na(dat$age_cat)
          #Chi-square test 
          chisq.test(table(dat$age_cat, dat$Intervention_Group))
      #Age continous using kruskal wall list test (interventiong group) 
        table_na(dat$age)
          kruskal.test(age ~ Intervention_Group, data = dat)
      #Race
        table_na(dat$race)
          chisq.test(table(dat$race, dat$Intervention_Group))
      #Homless
        table_na(dat$homeless)
          chisq.test(table(dat$homeless, dat$Intervention_Group))
      #Team
        table_na(dat$team3)
          chisq.test(table(dat$team3, dat$Intervention_Group))
      #Sex
        table_na(dat$sex3)
          #Fischer Exact test 
          fisher.test(table(dat$sex3, dat$Intervention_Group), workspace = 2e7)
          
library(gtsummary)
library(gt)
library(openxlsx)

# Convert the gtsummary object to a data frame
tbl1_df <- as.data.frame(
  dat_tbl1 %>%
    select(Intervention_Group, age, age_cat, sex, gender, language_major, race,
           homeless, housing, team_5, HIV_history) %>%
    tbl_summary(
      by       = Intervention_Group,
      missing  = "no",
      type     = list(age ~ "continuous2", age_cat ~ "categorical"),
      statistic = list(
        age ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 0,
      label = list(
        age ~ "Age (years)",
        age_cat ~ "Age category (CDC)",
        sex ~ "Sex",
        gender ~ "Gender",
        language_major ~ "Language (major)",
        race ~ "Race/Ethnicity",
        homeless ~ "Homeless",
        housing ~ "Housing",
        team_5 ~ "Primary Team",
        HIV_history ~ "History of HIV"
      )
    ) %>%
    add_overall(last = FALSE) %>%
    add_n() %>%
    add_p(
      test = list(
        age ~ "kruskal.test",
        age_cat ~ "chisq.test",
        all_categorical() ~ "chisq.test"
      ),
      test.args   = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 5000)),
      pvalue_fun  = ~ style_pvalue(.x, digits = 3)
    )
)

# Write to Excel
write.xlsx(tbl1_df, file = "HIP Analysis/Table1_Baseline_Characteristics.xlsx")

library(dplyr)
library(tidyr)
library(broom)



## 1) Inappropriate HIV Testing
inapp_data <- HIP %>%
  filter(hiv_Hx == "Yes") %>%
  group_by(Intervention_Group) %>%
  summarise(
    n_total = n(),
    n_yes = sum(Inappropriate_HIV_Testing == "Yes", na.rm = TRUE),
    pct = 100 * n_yes / n_total,
    .groups = "drop"
  )

load("R File/HIP_Cleaned.RData")
table_na(HIP$hiv_Hx)

# Global chi-square
tab_inapp <- HIP %>%
  filter(hiv_Hx == "Yes")

table_na(tab_inapp$Inappropriate_HIV_Testing)


# Create Table: Overall and by Intervention Group
tbl_inapp <- tab_inapp %>%
  select(Intervention_Group, Inappropriate_HIV_Testing) %>%
  tbl_summary(
    by = Intervention_Group,
    missing = "no",
    statistic = all_categorical() ~ "{n} ({p}%)",
    label = list(Inappropriate_HIV_Testing ~ "Inappropriate HIV Testing")
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(test = list(Inappropriate_HIV_Testing ~ "chisq.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  modify_header(
    label ~ md("**Characteristic**"),
    all_stat_cols() ~ "**n (%)**",
    p.value ~ md("**p-value**")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ md("**Intervention Group**")) %>%
  as_gt()

tbl_inapp
  
# Manual Chi-Square 
chisq_inapp <- fisher.test(table(tab_inapp$Intervention_Group, tab_inapp$Inappropriate_HIV_Testing))
 print(chisq_inapp)
 
 table_na(tab_inapp$Intervention_Group)

 pairs <- list(
   c("Pre-intervention", "Education Intervention"),
   c("Pre-intervention", "Education and EMR Intervention"),
   c("Education Intervention", "Education and EMR Intervention")
 )
 
 pairwise_results <- lapply(pairs, function(g) {
   subdat <- tab_inapp %>%
     filter(Intervention_Group %in% g) %>%
     droplevels() %>%                       # drop unused factor levels
     filter(!is.na(Inappropriate_HIV_Testing))
   
   # 2x2 contingency table for this pair only (after droplevels)
   xtab <- table(subdat$Intervention_Group, subdat$Inappropriate_HIV_Testing)
   
   # Exact Fisher for 2x2 (no simulation)
   ft <- fisher.test(xtab)
   
   # Optional: return risks/counts for clarity
   yes <- xtab[, "Yes", drop = TRUE]
   n   <- rowSums(xtab)
   
   data.frame(
     group1  = g[1],
     group2  = g[2],
     yes_g1  = as.integer(yes[1]),
     n_g1    = as.integer(n[1]),
     risk_g1 = yes[1] / n[1],
     yes_g2  = as.integer(yes[2]),
     n_g2    = as.integer(n[2]),
     risk_g2 = yes[2] / n[2],
     p.value = ft$p.value,
     stringsAsFactors = FALSE
   )
 }) %>%
   bind_rows() %>%
   mutate(
     p.adj_holm = p.adjust(p.value, method = "holm"),
     p.adj_bonf = p.adjust(p.value, method = "bonferroni")
   )
 
 pairwise_results

 library(epitools)   # for oddsratio()
 
 pairwise_results <- lapply(pairs, function(g) {
   subdat <- tab_inapp %>%
     filter(Intervention_Group %in% g) %>%
     droplevels() %>%
     filter(!is.na(Inappropriate_HIV_Testing))
   
   # 2x2 contingency
   xtab <- table(subdat$Intervention_Group, subdat$Inappropriate_HIV_Testing)
   
   # Exact Fisher
   ft <- fisher.test(xtab)
   
   # Odds ratio (Wald method for CI)
   or_res <- oddsratio(xtab, method = "wald")
   
   # Pull counts, risks, OR + CI
   yes <- xtab[, "Yes", drop = TRUE]
   n   <- rowSums(xtab)
   
   data.frame(
     group1   = g[1],
     group2   = g[2],
     yes_g1   = as.integer(yes[1]),
     n_g1     = as.integer(n[1]),
     risk_g1  = yes[1] / n[1],
     yes_g2   = as.integer(yes[2]),
     n_g2     = as.integer(n[2]),
     risk_g2  = yes[2] / n[2],
     odds_ratio = or_res$measure[2,1],
     or_lci     = or_res$measure[2,2],
     or_uci     = or_res$measure[2,3],
     p.value    = ft$p.value,
     stringsAsFactors = FALSE
   )
 }) %>%
   bind_rows() %>%
   mutate(
     p.adj_holm = p.adjust(p.value, method = "holm"),
     p.adj_bonf = p.adjust(p.value, method = "bonferroni")
   )
 
 pairwise_results

## 2) New HIV Diagnosis
 
table_na(dat$New_HIV_Diagnosis)

# Global chi-square
tab_newdx <- dat %>%
  filter(!is.na(New_HIV_Diagnosis)) %>% 
  count(Intervention_Group, New_HIV_Diagnosis) %>%
  pivot_wider(names_from = New_HIV_Diagnosis, values_from = n, values_fill = 0)
chisq_newdx <- chisq.test(as.matrix(tab_newdx[,-1]))

tab_newdx

new_dx <- fisher.test(table(dat$Intervention_Group, dat$New_HIV_Diagnosis))
print(new_dx)

pairs <- list(
  c("Pre-intervention", "Education Intervention"),
  c("Pre-intervention", "Education and EMR Intervention"),
  c("Education Intervention", "Education and EMR Intervention")
)

pairwise_results2 <- lapply(pairs, function(g) {
  subdat <- dat %>%
    filter(Intervention_Group %in% g) %>%
    droplevels()
  
  # 2x2 contingency for this pair only
  xtab <- table(subdat$Intervention_Group, subdat$New_HIV_Diagnosis)
  
  # Exact Fisher
  ft <- fisher.test(xtab)
  
  # Odds ratio + CI (Wald)
  or_res <- oddsratio(xtab, method = "wald")
  
  # Pull counts + risks
  yes <- xtab[, "Yes", drop = TRUE]
  n   <- rowSums(xtab)
  
  data.frame(
    group1     = g[1],
    group2     = g[2],
    yes_g1     = as.integer(yes[1]),
    n_g1       = as.integer(n[1]),
    risk_g1    = yes[1] / n[1],
    yes_g2     = as.integer(yes[2]),
    n_g2       = as.integer(n[2]),
    risk_g2    = yes[2] / n[2],
    odds_ratio = or_res$measure[2, 1],
    or_lci     = or_res$measure[2, 2],
    or_uci     = or_res$measure[2, 3],
    p.value    = ft$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    p.adj_holm = p.adjust(p.value, method = "holm"),
    p.adj_bonf = p.adjust(p.value, method = "bonferroni")
  )

pairwise_results2


# =========================
# FIGURE 1 (center title + lower labels, bigger axis text, 3-month ticks)
# =========================
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(purrr)
library(binom)

monthly <- dat %>%
  group_by(month, Intervention_Group) %>%
  summarise(n = n(),
            y = sum(screened),
            .groups = "drop") %>%
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper)
  )
monthly_all <- monthly %>% arrange(month)

cuts <- monthly %>%
  group_by(Intervention_Group) %>%
  summarise(start_month = min(month), .groups = "drop") %>%
  arrange(start_month) %>%
  filter(row_number() != 1) %>%
  mutate(label = c("Education Intervention", "Education and EMR Intervention"))

p_time_nejm <-
  ggplot() +
  # 95% CI ribbon
  geom_ribbon(
    data = monthly_all,
    aes(x = month, ymin = lower, ymax = upper),
    fill = "grey70", alpha = 0.35
  ) +
  # Main line
  geom_line(
    data = monthly_all,
    aes(x = month, y = mean, group = 1),
    linewidth = 0.9, color = "black"
  ) +
  geom_point(
    data = monthly_all,
    aes(x = month, y = mean),
    size = 2, color = "black", stroke = 0
  ) +
  # Phase cut lines
  geom_vline(
    data = cuts,
    aes(xintercept = as.numeric(start_month)),
    linetype = "dashed", linewidth = 0.5, color = "grey40"
  ) +
  # Phase labels placed lower
  geom_text(
    data = cuts,
    aes(x = start_month, y = 0.05, label = label),
    vjust = 1, hjust = 0, size = 4, family = "Helvetica"
  ) +
  # Axes
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    limits = c(0, 1.04), expand = c(0.01, 0)
  ) +
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b %Y",
    expand = c(0.01, 0)
  ) +
  labs(
    x = NULL,
    y = "Lifetime HIV Screening (%)"
  ) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(
      size = 14, family = "Helvetica",
      angle = 45, hjust = 1, margin = margin(t = 6)
    ),
    axis.text.y      = element_text(size = 14, family = "Helvetica", margin = margin(r = 6)),
    axis.title.y     = element_text(size = 14, family = "Helvetica"),
    axis.line        = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    legend.position  = "none",
    plot.margin      = margin(15, 16, 15, 15)
  ) +
  coord_cartesian(clip = "off")

p_time_nejm


# =========================
# FIGURE 1 (by Team: colored lines per team, bigger axis text, 3-month ticks)
# =========================
monthly_team <- dat %>%
  group_by(month, team3) %>%
  summarise(n = n(),
            y = sum(screened),
            .groups = "drop") %>%
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper)
  )

cuts <- dat %>%
  group_by(Intervention_Group) %>%
  summarise(start_month = min(month), .groups = "drop") %>%
  arrange(start_month) %>%
  filter(row_number() != 1) %>%
  mutate(label = c("Education Intervention", "Education and EMR Intervention"))

p_team_lines <-
  ggplot() +
  geom_ribbon(
    data = monthly_team,
    aes(x = month, ymin = lower, ymax = upper, fill = team3, group = team3),
    alpha = 0.3, colour = NA
  ) +
  geom_line(
    data = monthly_team,
    aes(x = month, y = mean, color = team3, group = team3),
    linewidth = 1
  ) +
  geom_point(
    data = monthly_team,
    aes(x = month, y = mean, color = team3),
    size = 2
  ) +
  geom_vline(
    data = cuts,
    aes(xintercept = as.numeric(start_month)),
    linetype = "dashed", linewidth = 0.6, color = "grey40"
  ) +
  geom_text(
    data = cuts,
    aes(x = start_month, y = 0.05, label = label),
    angle = 90, vjust = -0.5, hjust = .15, size = 4.5
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.0)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  scale_color_brewer(palette = "Dark2", name = "Team") +
  scale_fill_brewer(palette = "Pastel2", name = "Team") +
  labs(
    x = "Time (Months)", y = "Proportion screened (95% CI)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_team_lines

# =============================================================================
# Table 3. Period effects within Team (Education vs Pre; Education+EMR vs Pre)
# =============================================================================
suppressPackageStartupMessages({ library(emmeans); library(dplyr); library(stringr) })

em_team <- emmeans(fit_adj, ~ Period, by = "Team")
con_team <- contrast(em_team, method = "trt.vs.ctrl", ref = 1) # vs Pre-intervention
sm_team <- summary(con_team, type = "response", infer = TRUE, adjust = "holm")
team_df <- as.data.frame(sm_team)

# harmonize names
if (!"odds.ratio" %in% names(team_df)) {
  if ("ratio" %in% names(team_df)) team_df$odds.ratio <- team_df$ratio
  else if ("estimate" %in% names(team_df)) team_df$odds.ratio <- team_df$estimate
}
lcl <- intersect(c("lower.CL","asymp.LCL","LCL"), names(team_df))[1]
ucl <- intersect(c("upper.CL","asymp.UCL","UCL"), names(team_df))[1]

tab3_df <- team_df |>
  transmute(
    Team,
    Contrast = factor(contrast,
                      levels = c("Education - Pre-intervention",
                                 "Education and EMR - Pre-intervention"),
                      labels = c("Education vs Pre", "Education + EMR vs Pre")),
    OR = odds.ratio,
    LCL = .data[[lcl]],
    UCL = .data[[ucl]],
    `p (Holm)` = p.value
  ) |>
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", OR, LCL, UCL)) |>
  select(Team, Contrast, `OR (95% CI)`, `p (Holm)`)

tab3_gt <- gt(tab3_df) |>
  fmt_number(columns = "p (Holm)", decimals = 3) |>
  text_transform(
    locations = cells_body(columns = "p (Holm)"),
    fn = function(x) ifelse(suppressWarnings(as.numeric(x) < 0.001), "<0.001", x)
  ) |>
  cols_align(align = "right", columns = c(`OR (95% CI)`, `p (Holm)`)) |>
  tab_header(title = md("**Table 3. Period Effects Within Team**")) |>
  tab_source_note(md("Odds ratios from emmeans contrasts on the adjusted model; multiplicity adjusted (Holm)."))

tab3_gt


# =============================================================================
# Table 4. Period effects within Age and Sex strata (subgroups)
# =============================================================================
extract_period_vs_pre <- function(fit, by_var, nice_by) {
  em  <- emmeans(fit, ~ Period, by = by_var)
  con <- contrast(em, method = "trt.vs.ctrl", ref = 1)
  sm  <- summary(con, type = "response", infer = TRUE, adjust = "holm")
  df  <- as.data.frame(sm)
  if (!"odds.ratio" %in% names(df)) {
    if ("ratio" %in% names(df)) df$odds.ratio <- df$ratio
    else if ("estimate" %in% names(df)) df$odds.ratio <- df$estimate
  }
  lcl <- intersect(c("lower.CL","asymp.LCL","LCL"), names(df))[1]
  ucl <- intersect(c("upper.CL","asymp.UCL","UCL"), names(df))[1]
  df |>
    transmute(
      Stratum = .data[[by_var]],
      Facet = nice_by,
      Contrast = factor(contrast,
                        levels = c("Education - Pre-intervention",
                                   "Education and EMR - Pre-intervention"),
                        labels = c("Education vs Pre", "Education + EMR vs Pre")),
      OR = odds.ratio,
      LCL = .data[[lcl]],
      UCL = .data[[ucl]],
      `p (Holm)` = p.value,
      `OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", OR, LCL, UCL)
    )
}

age_df <- extract_period_vs_pre(fit_adj, "Age", "Age")
sex_df <- extract_period_vs_pre(fit_adj, "Sex", "Sex")

tab4_df <- dplyr::bind_rows(age_df, sex_df) |>
  select(Facet, Stratum, Contrast, `OR (95% CI)`, `p (Holm)`)

tab4_gt <- gt(tab4_df) |>
  fmt_number(columns = "p (Holm)", decimals = 3) |>
  text_transform(
    locations = cells_body(columns = "p (Holm)"),
    fn = function(x) ifelse(suppressWarnings(as.numeric(x) < 0.001), "<0.001", x)
  ) |>
  tab_spanner(label = "Age strata", columns = c()) |>
  cols_align(align = "right", columns = c(`OR (95% CI)`, `p (Holm)`)) |>
  tab_header(title = md("**Table 4. Period Effects Within Subgroups (Age, Sex)**")) |>
  tab_source_note(md("Odds ratios from emmeans contrasts on the adjusted model; multiplicity adjusted (Holm)."))

tab4_gt


# =============================================================================
# Table 5. Unadjusted screening rates by Period × Team (Wilson 95% CI)
# =============================================================================
suppressPackageStartupMessages({ library(purrr); library(binom) })

tab5_df <- dat |>
  group_by(Intervention_Group, team3) |>
  summarise(n = n(), y = sum(screened), .groups = "drop") |>
  mutate(
    mean  = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$mean),
    lower = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$lower),
    upper = map2_dbl(y, n, ~ binom::binom.wilson(.x, .y)$upper),
    `Proportion (95% CI)` = sprintf("%.1f%% (%.1f–%.1f)", 100*mean, 100*lower, 100*upper)
  ) |>
  arrange(Intervention_Group, team3) |>
  select(Intervention_Group, Team = team3, n, y, `Proportion (95% CI)`)

tab5_gt <- gt(tab5_df) |>
  cols_label(
    Intervention_Group = "Intervention Period",
    n = "N",
    y = "Yes"
  ) |>
  cols_align(align = "right", columns = c(n, y, `Proportion (95% CI)`)) |>
  tab_header(title = md("**Table 5. Unadjusted Screening Rates by Period and Team**")) |>
  tab_source_note(md("Wilson score intervals shown."))

tab5_gt


# =========================================================
# Table S1. Sensitivity: first admission only (adjusted model)
# =========================================================
tblS1 <-
  tbl_regression(
    fit_adj_first,
    exponentiate = TRUE,
    label = var_labels_tbl2
  ) |>
  add_glance_table(include = c(AIC, BIC, nobs)) |>
  bold_labels() |>
  modify_header(label ~ md("**Covariate**")) |>
  modify_spanning_header(starts_with("estimate_") ~ md("**Adjusted Association (OR, 95% CI)**")) |>
  modify_footnote(
    update = everything() ~ "Odds ratios from logistic regression using only each patient's first admission."
  )

as_gt(tblS1) |>
  tab_header(title = md("**Table S1. Adjusted Odds Ratios — First Admission Only**"))


