
#  Chennai Reservoir Levels & Rainfall Analysis

#Initial import from github
levels_url <- "https://raw.githubusercontent.com/marishaabell88/Shortage-Warning-/refs/heads/main/chennai_reservoir_levels.csv"
rainfall_url <- "https://raw.githubusercontent.com/marishaabell88/Shortage-Warning-/refs/heads/main/chennai_reservoir_rainfall.csv"

#Download both scripts automatically when running the script
download.file(
  url = levels_url,
  destfile = "chennai_reservoir_levels.csv",
  mode = "wb"
)

download.file(
  url = rainfall_url,
  destfile = "chennai_reservoir_rainfall.csv",
  mode = "wb"
)

#  1. Libraries
req <- c("readr","dplyr","ggplot2","tidyr","lubridate","broom","randomForest","pROC")
for (p in req) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(lubridate); library(broom); library(randomForest); library(pROC)
})

#  1b. SETTINGS
TOTAL_CAP     <- 11257          # in mcft
THRESH        <- 20             # shortage = storage <= THRESH % of capacity
SENS_THRESH   <- sort(unique(c(15, THRESH, 25)))   # sensitivity check
BASE_MARGIN   <- 5              # baseline rule: warn if storage <= THRESH + 5 (i.e. 25% when THRESH = 20)
COST_RATIOS   <- c(5, 10, 20, 50)
TIMELINE_RATIO <- 20            # cost ratio used for the warning timeline / months-of-notice
CUTS          <- seq(0.05, 0.95, by = 0.05)        # candidate probability thresholds
SPLIT_DATE    <- as.Date("2017-01-01")             # train = before, test = on/after. Adjust to your data
RES           <- c("POONDI","CHOLAVARAM","REDHILLS","CHEMBARAMBAKKAM")
set.seed(2280)

#2. READ + CLEAN: reservoir LEVELS and RAINFALL 
read_res <- function(file) {
  raw <- readr::read_csv(file, col_types = readr::cols(.default = readr::col_character())) %>%
    rename(Date = 1)                                   # first column is the date/label
  cat("\n", file, "- columns:", paste(names(raw), collapse = ", "), "\n")
  cat("First-column values:\n"); print(head(raw$Date, 5))
  out <- raw %>%
    filter(grepl("^\\d{2}-\\d{2}-\\d{4}$", Date)) %>%
    mutate(Date = dmy(Date), across(all_of(RES), as.numeric)) %>%
    select(Date, all_of(RES)) %>%
    arrange(Date)
  if (nrow(out) == 0) stop("No rows matched DD-MM-YYYY in ", file, ": check the date format above")
  cat("Rows kept:", nrow(out), "| range:", format(min(out$Date)), "to", format(max(out$Date)), "\n")
  out
}
levels_data <- read_res("chennai_reservoir_levels.csv")
rainfall    <- read_res("chennai_reservoir_rainfall.csv")
if (nrow(levels_data) < 365)
  warning("Levels file has fewer than 365 rows: a monthly/seasonal analysis of storage needs the daily series.")

# 3. DATA-QUALITY CHECKS
date_gaps <- function(d) d %>% arrange(Date) %>%
  mutate(gap = as.numeric(Date - lag(Date))) %>% filter(gap > 1) %>% select(Date, gap)
qc <- list(
  missing_levels = levels_data %>% summarise(across(all_of(RES), ~ sum(is.na(.)))),
  missing_rain   = rainfall    %>% summarise(across(all_of(RES), ~ sum(is.na(.)))),
  dup_levels     = sum(duplicated(levels_data$Date)),
  dup_rain       = sum(duplicated(rainfall$Date)),
  gaps_levels    = date_gaps(levels_data),
  gaps_rain      = date_gaps(rainfall)
)
cat("\n[QC] Missing LEVELS per reservoir:\n");   print(qc$missing_levels)
cat("[QC] Missing RAINFALL per reservoir (set to 0 below):\n"); print(qc$missing_rain)
cat("[QC] Duplicate dates  levels:", qc$dup_levels, " rainfall:", qc$dup_rain, "\n")
cat("[QC] Date gaps > 1 day  levels:", nrow(qc$gaps_levels), " rainfall:", nrow(qc$gaps_rain), "\n")

#  4. merge, capacity %, MONTHLY table 
# 4a. Merge by date ; missing rainfall -> 0
daily <- inner_join(
  levels_data %>% rename_with(~ paste0(.x, "_lvl"),  all_of(RES)),
  rainfall    %>% rename_with(~ paste0(.x, "_rain"), all_of(RES)),
  by = "Date") %>%
  mutate(
    across(ends_with("_rain"), ~ coalesce(.x, 0)),
    Total_Storage = POONDI_lvl + CHOLAVARAM_lvl + REDHILLS_lvl + CHEMBARAMBAKKAM_lvl,
    Pct_Capacity  = Total_Storage / TOTAL_CAP * 100,
    # average over the 4 reservoirs each day; summed over the month below
    Rain_Mean     = (POONDI_rain + CHOLAVARAM_rain + REDHILLS_rain + CHEMBARAMBAKKAM_rain) / 4
  )
cat("\nMerged daily rows:", nrow(daily), "\n")

# 4b. Monthly table. Months with missing days are set to NA (not analysed),
# and every calendar month is kept so lag()/lead() always mean "previous/next month".
monthly_raw <- daily %>%
  mutate(Month_Start = floor_date(Date, "month")) %>%
  group_by(Month_Start) %>%
  summarise(Storage_Pct = mean(Pct_Capacity, na.rm = TRUE),   # monthly MEAN storage (% of capacity)
            Rain        = sum(Rain_Mean),                     # monthly TOTAL rainfall (mm)
            Days        = n(), .groups = "drop") %>%
  mutate(Complete = Days == days_in_month(Month_Start),
         Storage_Pct = ifelse(Complete, Storage_Pct, NA_real_),
         Rain        = ifelse(Complete, Rain, NA_real_))
cat("Incomplete months set to NA:", sum(!monthly_raw$Complete), "of", nrow(monthly_raw), "\n")

monthly <- tibble(Month_Start = seq(min(monthly_raw$Month_Start), max(monthly_raw$Month_Start), by = "month")) %>%
  left_join(monthly_raw, by = "Month_Start") %>%
  mutate(
    Year      = year(Month_Start),
    MonthNum  = month(Month_Start),
    Cal_Month = factor(month.abb[MonthNum], levels = month.abb),
    Monsoon   = MonthNum %in% 10:12,                                   # report: Oct-Dec = NE monsoon
    Season    = ifelse(Monsoon, "Monsoon (Oct-Dec)", "Non-monsoon"),
    Change    = Storage_Pct - lag(Storage_Pct),                        # monthly change in storage
    Rain_L1   = lag(Rain, 1),                                          # rainfall 1 month ago
    Rain_L2   = lag(Rain, 2),                                          # rainfall 2 months ago
    Shortage  = Storage_Pct <= THRESH
  )
cat("Months analysed (non-NA storage):", sum(!is.na(monthly$Storage_Pct)),
    "| shortage months:", sum(monthly$Shortage, na.rm = TRUE), "\n")


#  RESEARCH QUESTIONS: tabulated against the data and methods
rq_table <- tibble(
  Q = c("RQ1","RQ2","RQ3"),
  Question = c(
    "Is a seasonal pattern evident in the water storage levels?",
    "Does historical reservoir data help forecast future water shortages?",
    "How might climate change affect long-term water availability?"),
  Data_Used = c(
    "Monthly mean storage (% capacity), monthly rainfall, monsoon indicator",
    "Storage, monsoon, rainfall and 2 lags -> next-month shortage (storage <= 20%)",
    "Annual rainfall (full years only): trend and year-to-year variability"),
  Method = c(
    "Boxplots by month/season; Welch t-test, Wilcoxon, Spearman",
    "Time-based split; logistic regression + random forest; cost-based warning rule vs baseline",
    "lm(annual rain ~ year); coefficient of variation (interpret cautiously)"),
  Output = c(
    "storage_by_month.png, storage_by_season.png, seasonal_means table",
    "metrics, odds ratios, cost tables, warning timeline",
    "annual_summary.csv, trend and CV")
)
cat("\n================ RESEARCH-QUESTION PLAN ================\n")
print(rq_table, n = Inf, width = Inf)


#  5. SUMMARY + STATISTICS
ok <- monthly %>% filter(!is.na(Storage_Pct))

# 5.1 Seasonal summary table -> RQ1
seasonal_means <- ok %>%
  group_by(Season) %>%
  summarise(n = n(), Mean_Storage = mean(Storage_Pct), SD_Storage = sd(Storage_Pct),
            Mean_Rain = mean(Rain, na.rm = TRUE), Mean_Change = mean(Change, na.rm = TRUE),
            Pct_Shortage_Months = mean(Shortage) * 100, .groups = "drop")
cat("\n[5.1] Monsoon vs non-monsoon summary:\n"); print(seasonal_means, width = Inf)

# 5.2 Calendar-month summary -> RQ1
monthly_means <- ok %>%
  group_by(MonthNum, Cal_Month) %>%
  summarise(Mean_Storage = mean(Storage_Pct), Mean_Rain = mean(Rain, na.rm = TRUE),
            Pct_Shortage = mean(Shortage) * 100, .groups = "drop") %>% arrange(MonthNum)
cat("\n[5.2] Calendar-month means:\n"); print(monthly_means, n = 12, width = Inf)

# 5.3 Storage summary
storage_summary <- ok %>%
  summarise(Mean_Pct = mean(Storage_Pct), Min_Pct = min(Storage_Pct), Max_Pct = max(Storage_Pct),
            Shortage_Months = sum(Shortage), Total_Months = n())
cat("\n[5.3] Storage summary:\n"); print(storage_summary, width = Inf)

# 5.4 Hypothesis tests -> RQ1 / Objective 1
#   H0: mean storage is the same in monsoon and non-monsoon months; H1: it differs; alpha = 0.05
tt <- t.test(Storage_Pct ~ Monsoon, data = ok)                       # Welch (unequal variances)
cat("\n[5.4a] Welch t-test: storage, monsoon vs non-monsoon\n")
cat(sprintf("  mean non-monsoon = %.2f, mean monsoon = %.2f\n", tt$estimate[1], tt$estimate[2]))
cat(sprintf("  t = %.3f, df = %.1f, p = %s, 95%% CI for difference = [%.2f, %.2f]\n",
            tt$statistic, tt$parameter, format.pval(tt$p.value), tt$conf.int[1], tt$conf.int[2]))
wt <- wilcox.test(Storage_Pct ~ Monsoon, data = ok, exact = FALSE)   # non-parametric check
cat("[5.4b] Wilcoxon rank-sum: W =", wt$statistic, " p =", format.pval(wt$p.value), "\n")

tt_chg <- t.test(Change ~ Monsoon, data = ok %>% filter(!is.na(Change)))   # change in storage by season
cat(sprintf("[5.4c] Welch t-test on CHANGE in storage: t = %.3f, p = %s\n",
            tt_chg$statistic, format.pval(tt_chg$p.value)))

# 5.5 Rainfall vs change in storage (Spearman) -> Objective 1
cc <- ok %>% filter(!is.na(Change), !is.na(Rain))
sp <- cor.test(cc$Rain, cc$Change, method = "spearman", exact = FALSE)
cat(sprintf("\n[5.5] Spearman rho (rain vs change in storage) = %.3f, p = %s\n", sp$estimate, format.pval(sp$p.value)))
cc1 <- ok %>% filter(!is.na(Change), !is.na(Rain_L1))
sp1 <- cor.test(cc1$Rain_L1, cc1$Change, method = "spearman", exact = FALSE)
cat(sprintf("      Spearman rho (LAST month's rain vs change) = %.3f, p = %s\n", sp1$estimate, format.pval(sp1$p.value)))
ac1 <- acf(monthly$Storage_Pct, plot = FALSE, na.action = na.pass)$acf[2]
cat(sprintf("      Lag-1 autocorrelation of monthly storage = %.2f (observations are NOT independent: p-values are optimistic)\n", ac1))

# 5.6 Warning-rule models (logistic + random forest) 
eval_rule <- function(warn, y) {
  TP <- sum(warn & y); FN <- sum(!warn & y); FP <- sum(warn & !y); TN <- sum(!warn & !y)
  tibble(Accuracy    = (TP + TN) / length(y),
         Sensitivity = ifelse(TP + FN > 0, TP / (TP + FN), NA_real_),
         Specificity = ifelse(TN + FP > 0, TN / (TN + FP), NA_real_),
         Precision   = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
         TP = TP, FN = FN, FP = FP, TN = TN)
}
cost_of <- function(warn, y, ratio) ratio * sum(!warn & y) + sum(warn & !y)   # cost = ratio*missed + 1*false alarms
auc_of  <- function(p, y) as.numeric(pROC::auc(pROC::roc(as.integer(y), p, levels = c(0, 1),
                                                          direction = "<", quiet = TRUE)))
pick_cut <- function(p, y, ratio) {
  costs <- vapply(CUTS, function(ct) cost_of(p >= ct, y, ratio), numeric(1))
  CUTS[which.min(costs)]                       # ties -> lowest cut (warns more)
}
# months of advance notice before the first shortage in the test period
lead_months <- function(warn, next_short) {
  onset <- which(next_short & !dplyr::lag(next_short, default = TRUE))[1]
  if (is.na(onset)) return(NA_integer_)
  if (!warn[onset]) return(0L)
  i <- onset
  while (i > 1 && warn[i - 1]) i <- i - 1
  as.integer(onset - i + 1)
}

run_pipeline <- function(thresh) {
  set.seed(2280)
  mdf <- monthly %>%
    mutate(Next_Shortage = lead(Storage_Pct) <= thresh,        # target: is NEXT month a shortage?
           Monsoon = as.integer(Monsoon)) %>%
    select(Month_Start, Storage_Pct, Monsoon, Rain, Rain_L1, Rain_L2, Next_Shortage) %>%
    drop_na()
  train <- mdf %>% filter(Month_Start <  SPLIT_DATE)           # time-based split, NOT random
  test  <- mdf %>% filter(Month_Start >= SPLIT_DATE)
  counts <- tibble(Thresh = thresh, Train_Months = nrow(train), Train_Shortages = sum(train$Next_Shortage),
                   Test_Months = nrow(test), Test_Shortages = sum(test$Next_Shortage))
  if (min(sum(train$Next_Shortage), sum(!train$Next_Shortage)) < 2 ||
      min(sum(test$Next_Shortage),  sum(!test$Next_Shortage))  < 1) {
    message("Threshold ", thresh, "%: not enough shortage / non-shortage months in train or test. Try another SPLIT_DATE.")
    return(list(ok = FALSE, counts = counts))
  }
  y <- test$Next_Shortage

  logit <- glm(Next_Shortage ~ Storage_Pct + Monsoon + Rain + Rain_L1 + Rain_L2,
               data = train, family = binomial)
  rf <- randomForest(factor(Next_Shortage, levels = c(FALSE, TRUE)) ~
                       Storage_Pct + Monsoon + Rain + Rain_L1 + Rain_L2,
                     data = train, ntree = 500)

  p_test_log  <- predict(logit, newdata = test, type = "response")
  p_test_rf   <- predict(rf,    newdata = test, type = "prob")[, "TRUE"]
  p_train_log <- fitted(logit)
  p_train_rf  <- rf$votes[, "TRUE"]                              # out-of-bag probabilities (honest)
  base_warn   <- test$Storage_Pct <= thresh + BASE_MARGIN        # baseline rule

  metrics <- bind_rows(
    eval_rule(p_test_log >= 0.5, y) %>% mutate(Model = "Logistic (cut 0.5)", AUC = auc_of(p_test_log, y)),
    eval_rule(p_test_rf  >= 0.5, y) %>% mutate(Model = "Random forest (cut 0.5)", AUC = auc_of(p_test_rf, y)),
    eval_rule(base_warn, y)         %>% mutate(Model = paste0("Baseline (storage <= ", thresh + BASE_MARGIN, "%)"), AUC = NA_real_)
  ) %>% mutate(Thresh = thresh) %>% select(Thresh, Model, everything())

  best <- bind_rows(lapply(COST_RATIOS, function(r) {
    cl <- pick_cut(p_train_log, train$Next_Shortage, r)
    cr <- pick_cut(p_train_rf,  train$Next_Shortage, r)
    wl <- p_test_log >= cl; wr <- p_test_rf >= cr
    tibble(Thresh = thresh, Ratio = r, Rule = c("Logistic", "Random forest", "Baseline"),
           Cut = c(cl, cr, NA_real_),
           Missed       = c(sum(!wl & y), sum(!wr & y), sum(!base_warn & y)),
           False_Alarms = c(sum(wl & !y), sum(wr & !y), sum(base_warn & !y))) %>%
      mutate(Cost = Ratio * Missed + False_Alarms)
  }))

  cost_grid <- expand_grid(Ratio = COST_RATIOS, Cut = CUTS) %>%
    mutate(Logistic        = mapply(function(r, ct) cost_of(p_test_log >= ct, y, r), Ratio, Cut),
           `Random forest` = mapply(function(r, ct) cost_of(p_test_rf  >= ct, y, r), Ratio, Cut))

  cl <- pick_cut(p_train_log, train$Next_Shortage, TIMELINE_RATIO)
  cr <- pick_cut(p_train_rf,  train$Next_Shortage, TIMELINE_RATIO)
  warn_table <- test %>%
    mutate(Prob_Logistic = p_test_log, Prob_Forest = p_test_rf,
           Warn_Logistic = p_test_log >= cl, Warn_Forest = p_test_rf >= cr, Warn_Baseline = base_warn)
  notice <- tibble(Thresh = thresh, Rule = c("Logistic", "Random forest", "Baseline"),
                   Months_Notice = c(lead_months(warn_table$Warn_Logistic, y),
                                     lead_months(warn_table$Warn_Forest,   y),
                                     lead_months(warn_table$Warn_Baseline, y)))

  list(ok = TRUE, counts = counts, metrics = metrics, best = best, cost_grid = cost_grid,
       warn_table = warn_table, notice = notice, logit = logit, rf = rf)
}

runs <- setNames(lapply(SENS_THRESH, run_pipeline), as.character(SENS_THRESH))
main <- runs[[as.character(THRESH)]]

cat("\n[5.6] Train/test shortage counts at each threshold (sensitivity check):\n")
print(bind_rows(lapply(runs, `[[`, "counts")), width = Inf)

if (isTRUE(main$ok)) {
  cat("\n[5.7] Logistic regression odds ratios (trained on months before", format(SPLIT_DATE), ")\n")
  print(tidy(main$logit, conf.int = TRUE, exponentiate = TRUE) %>%
          select(term, Odds_Ratio = estimate, conf.low, conf.high, p.value), width = Inf)
  cat("\n[5.8] Test-set performance (threshold", THRESH, "%): accuracy alone is misleading, use sensitivity/precision/AUC\n")
  print(main$metrics, width = Inf)
  cat("\n[5.9] Cost of the warning rule: threshold picked on TRAIN, evaluated on TEST\n")
  print(main$best, n = Inf, width = Inf)
  cat("\n[5.10] Months of advance notice before the first test-period shortage (cost ratio", TIMELINE_RATIO, "):\n")
  print(main$notice, width = Inf)
  cat("\n[5.11] Random forest variable importance:\n"); print(importance(main$rf))
}

# Sensitivity tables across 15 / 20 / 25 %
good <- Filter(function(s) isTRUE(s$ok), runs)
sens_metrics <- bind_rows(lapply(good, `[[`, "metrics"))
sens_best    <- bind_rows(lapply(good, `[[`, "best"))
sens_notice  <- bind_rows(lapply(good, `[[`, "notice"))
cat("\n[5.12] Sensitivity: metrics by threshold\n"); print(sens_metrics, n = Inf, width = Inf)
cat("\n[5.13] Sensitivity: cost by threshold, ratio and rule\n"); print(sens_best, n = Inf, width = Inf)

# 5.14 Climate-variability trend -> RQ3 
annual <- monthly %>%
  group_by(Year) %>%
  summarise(Months = sum(!is.na(Rain)), Annual_Rain = sum(Rain, na.rm = TRUE), .groups = "drop") %>%
  mutate(Full_Year = Months == 12)
full <- annual %>% filter(Full_Year)
if (nrow(full) >= 4) {
  cv <- sd(full$Annual_Rain) / mean(full$Annual_Rain)
  tr <- lm(Annual_Rain ~ Year, data = full)
  cat(sprintf("\n[5.14] Full years used: %d | annual rain CV = %.2f | linear trend = %.1f mm/year (p = %s)\n",
              nrow(full), cv, unname(coef(tr)[2]), format.pval(summary(tr)$coefficients[2, 4])))
  cat("       Only ~15 years of data: describe this as variability, not proof of a climate trend.\n")
}
print(annual)


#  6. VISUALISATIONS
ymax <- max(100, max(monthly$Storage_Pct, na.rm = TRUE))   # if this exceeds 100, check TOTAL_CAP

# Figure 1: storage over time with the shortage line
p_storage <- ggplot(monthly, aes(Month_Start, Storage_Pct)) +
  geom_line(colour = "steelblue", linewidth = 1, na.rm = TRUE) +
  geom_hline(yintercept = THRESH, linetype = "dashed", colour = "red") +
  annotate("text", x = min(monthly$Month_Start), y = THRESH, label = paste0(THRESH, "% shortage line"),
           hjust = -0.05, vjust = -0.6, colour = "red") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(title = "Total reservoir storage over time (monthly mean, % of capacity)",
       x = "Date", y = "% of total capacity") + theme_minimal(base_size = 13)
ggsave("storage_over_time.png", p_storage, width = 10, height = 5); print(p_storage)

# Figure 2: storage by calendar month
p_month <- ggplot(ok, aes(Cal_Month, Storage_Pct)) +
  geom_boxplot(fill = "lightblue", outlier.size = 1) +
  labs(title = "Storage by calendar month", x = "Month", y = "% of total capacity") +
  theme_minimal(base_size = 13)
ggsave("storage_by_month.png", p_month, width = 10, height = 5); print(p_month)

# Figure 3: monsoon vs non-monsoon
p_season <- ggplot(ok, aes(Season, Storage_Pct, fill = Season)) +
  geom_boxplot(outlier.size = 1) +
  labs(title = "Storage: monsoon (Oct-Dec) vs non-monsoon", x = NULL, y = "% of total capacity") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")
ggsave("storage_by_season.png", p_season, width = 7, height = 5); print(p_season)

# Figure 4: rainfall vs change in storage
p_rain <- ggplot(cc, aes(Rain, Change)) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
  labs(title = "Monthly rainfall vs change in storage", x = "Monthly total rainfall (mm)",
       y = "Change in storage (percentage points)") + theme_minimal(base_size = 13)
ggsave("rain_vs_change.png", p_rain, width = 8, height = 5); print(p_rain)

# Figure 5: monthly rainfall by year
p_trend <- monthly %>% filter(!is.na(Rain)) %>%
  ggplot(aes(MonthNum, Rain, colour = factor(Year))) +
  geom_line(alpha = 0.7) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Monthly total rainfall by year", x = "Month", y = "Total rainfall (mm)", colour = "Year") +
  theme_minimal(base_size = 13)
ggsave("rainfall_by_year.png", p_trend, width = 12, height = 6); print(p_trend)

if (isTRUE(main$ok)) {
  # Figure 6: random forest variable importance
  imp <- as.data.frame(importance(main$rf)); imp$Variable <- rownames(imp)
  p_imp <- ggplot(imp, aes(reorder(Variable, MeanDecreaseGini), MeanDecreaseGini)) +
    geom_col(fill = "steelblue") + coord_flip() +
    labs(title = "Random forest variable importance", x = NULL, y = "Mean decrease in Gini") +
    theme_minimal(base_size = 13)
  ggsave("rf_importance.png", p_imp, width = 8, height = 4); print(p_imp)

  # Figure 7: cost vs warning threshold
  p_cost <- main$cost_grid %>%
    pivot_longer(c(Logistic, `Random forest`), names_to = "Model", values_to = "Cost") %>%
    ggplot(aes(Cut, Cost, colour = Model)) + geom_line() +
    facet_wrap(~ Ratio, scales = "free_y", labeller = label_both) +
    labs(title = "Total cost vs warning threshold (test set)", x = "Warning threshold (probability)", y = "Total cost") +
    theme_minimal(base_size = 13)
  ggsave("cost_vs_threshold.png", p_cost, width = 10, height = 6); print(p_cost)

  # Figure 8: warning timeline in the test period
  wt_tab <- main$warn_table
  warn_long <- wt_tab %>%
    select(Month_Start, Warn_Logistic, Warn_Forest, Warn_Baseline) %>%
    pivot_longer(-Month_Start, names_to = "Rule", values_to = "Warn") %>%
    filter(Warn) %>% mutate(Rule = sub("Warn_", "", Rule),
                            Y = unname(c(Logistic = 2, Forest = 5, Baseline = 8)[Rule]))
  p_time <- ggplot(wt_tab, aes(Month_Start, Storage_Pct)) +
    geom_line(colour = "grey30") +
    geom_hline(yintercept = THRESH, linetype = "dashed", colour = "red") +
    geom_point(data = warn_long, aes(y = Y, colour = Rule), shape = 15, size = 2) +
    labs(title = "Warnings issued in the test period (squares = warning that month)",
         x = "Date", y = "% of capacity") + theme_minimal(base_size = 13)
  ggsave("warning_timeline.png", p_time, width = 10, height = 5); print(p_time)
}

if (isTRUE(main$ok)) {
  wt_tab <- main$warn_table
  
  # Figure 9: logistic curve (other predictors held at training medians)
  tr <- main$logit$model
  med <- lapply(tr[, c("Monsoon","Rain","Rain_L1","Rain_L2")], median)
  grid <- tibble(Storage_Pct = seq(0, 60, 0.5),
                 Monsoon = med$Monsoon, Rain = med$Rain,
                 Rain_L1 = med$Rain_L1, Rain_L2 = med$Rain_L2)
  grid$Prob <- predict(main$logit, newdata = grid, type = "response")
  
  p_logit <- ggplot() +
    geom_line(data = grid, aes(Storage_Pct, Prob), colour = "steelblue", linewidth = 1) +
    geom_point(data = wt_tab, aes(Storage_Pct, as.integer(Next_Shortage)),
               alpha = 0.6, position = position_jitter(height = 0.02, width = 0)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
    labs(title = "Logistic regression: probability of next-month shortage",
         x = "Current storage (% of capacity)",
         y = "Predicted probability (points = test months)") +
    theme_minimal(base_size = 13)
  ggsave("logistic_curve.png", p_logit, width = 8, height = 5); print(p_logit)
  
  # Figure 10: ROC curves, logistic vs random forest (test set)
  y <- as.integer(wt_tab$Next_Shortage)
  roc_log <- pROC::roc(y, wt_tab$Prob_Logistic, levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_rf  <- pROC::roc(y, wt_tab$Prob_Forest,   levels = c(0, 1), direction = "<", quiet = TRUE)
  p_roc <- pROC::ggroc(list(Logistic = roc_log, "Random forest" = roc_rf)) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "grey60") +
    labs(title = "ROC curves on the test set", colour = "Model") +
    theme_minimal(base_size = 13)
  ggsave("roc_curves.png", p_roc, width = 7, height = 5); print(p_roc)
}

#  7. EXPORT
write_csv(monthly, "chennai_monthly.csv")
write_csv(annual, "annual_summary.csv")
write_csv(seasonal_means, "seasonal_means.csv")
if (nrow(sens_metrics) > 0) {
  write_csv(sens_metrics, "model_metrics_by_threshold.csv")
  write_csv(sens_best, "cost_by_threshold_ratio.csv")
  write_csv(sens_notice, "months_of_notice.csv")
}
if (isTRUE(main$ok)) write_csv(main$warn_table, "warning_table_test_period.csv")
cat("\nSaved CSVs and PNGs in the working directory.\n")
