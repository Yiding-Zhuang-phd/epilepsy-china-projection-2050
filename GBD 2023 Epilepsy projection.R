setwd("E:/研究生文件/研0/2023癫痫GBD/Projection")
library(tidyverse)
library(forecast)
library(MASS)
library(tidyverse)
library(tidyr)
## ==================== 1. 导入三个维度 ====================
EDU <- read.csv("E:/研究生文件/研0/2023癫痫GBD/China - Subnational/China - Subnational/CHN_edu_years_above_15.csv")
LDI <- read.csv("E:/研究生文件/研0/2023癫痫GBD/China - Subnational/China - Subnational/CHN_LDI_10_yrs.csv")
fert <- read.csv("E:/研究生文件/研0/2023癫痫GBD/China - Subnational/China - Subnational/CHN_tot_fert_under_25.csv")
fert <- fert %>%
  filter(sex == "Female")
edu_hist <- EDU %>% dplyr::select(location_name, year_id, mean_value)
ldi_hist <- LDI %>% dplyr::select(location_name, year_id, mean_value, lower_value, upper_value)
fert_hist <- fert %>% dplyr::select(location_name, year_id, mean_value, lower_value, upper_value)
locations <- sort(unique(edu_hist$location_name))
## ==================== 2. 预测函数（保持原逻辑） ====================
forecast_component <- function(df, end_year = 2050, n_draws = 1000) {
  trend_model <- lm(mean_value ~ year_id, data = df)
  resid_ts <- ts(residuals(trend_model), start = min(df$year_id), frequency = 1)
  resid_fit <- Arima(resid_ts, order = c(0,1,0))
  
  future_years <- (max(df$year_id) + 1):end_year
  h <- length(future_years)
  resid_fc <- forecast(resid_fit, h = h)
  
  set.seed(123)
  draws <- matrix(NA, nrow = h, ncol = n_draws)
  
  for (i in 1:n_draws) {
    coef_i <- mvrnorm(1, coef(trend_model), vcov(trend_model))
    trend_pred <- coef_i[1] + coef_i[2] * future_years
    resid_i <- rnorm(h, mean = resid_fc$mean, sd = sd(resid_fc$residuals, na.rm = TRUE))
    draws[, i] <- trend_pred + resid_i
  }
  
  list(years = future_years, draws = draws)
}
predict_all <- function(df) split(df, df$location_name) %>% lapply(forecast_component)
edu_pred <- predict_all(edu_hist)
ldi_pred <- predict_all(ldi_hist)
fert_pred <- predict_all(fert_hist)
## ==================== 3. 历史draws函数 ====================
make_hist_draws <- function(mean, lower = NULL, upper = NULL, n_draws = 1000) {
  if (!is.null(lower) && !is.null(upper)) {
    sd <- (upper - lower) / (2 * 1.96)
    matrix(rnorm(length(mean) * n_draws, mean, sd), ncol = n_draws)
  } else {
    matrix(rep(mean, n_draws), ncol = n_draws)
  }
}
## ==================== 4. 合并历史+预测draws ====================
make_full_draws <- function(hist_df, pred_list, n_draws = 1000) {
  map(locations, ~ {
    hist_loc <- hist_df %>% filter(location_name == .x)
    hist_mat <- if ("lower_value" %in% names(hist_loc)) {
      make_hist_draws(hist_loc$mean_value, hist_loc$lower_value, hist_loc$upper_value, n_draws)
    } else {
      make_hist_draws(hist_loc$mean_value, n_draws = n_draws)
    }
    rbind(hist_mat, pred_list[[.x]]$draws)
  }) %>% set_names(locations)
}
edu_full <- make_full_draws(edu_hist, edu_pred)
ldi_full <- make_full_draws(ldi_hist, ldi_pred)
fert_full <- make_full_draws(fert_hist, fert_pred)
## ==================== 5. min/max（保留原逻辑） ====================
edu_min <- min(unlist(edu_full), na.rm = TRUE)
edu_max <- max(unlist(edu_full), na.rm = TRUE)
ldi_min <- min(unlist(ldi_full), na.rm = TRUE)
ldi_max <- max(unlist(ldi_full), na.rm = TRUE)
fert_min <- min(unlist(fert_full), na.rm = TRUE)
fert_max <- max(unlist(fert_full), na.rm = TRUE)
## ==================== 6. SDI合成 ====================
compute_sdi <- function(edu, ldi, fert) {
  edu_i <- (edu - edu_min) / (edu_max - edu_min)
  ldi_i <- (ldi - ldi_min) / (ldi_max - ldi_min)
  fert_i <- (fert_max - fert) / (fert_max - fert_min)
  pmin(pmax((edu_i * ldi_i * fert_i)^(1/3), 0), 1)
}
## ==================== 7. 生成各省SDI draws ====================
years_common <- c(sort(unique(edu_hist$year_id)), edu_pred[[1]]$years)
sdi_draws <- map(locations, ~ compute_sdi(edu_full[[.x]], ldi_full[[.x]], fert_full[[.x]])) %>%
  set_names(locations)
## ==================== 8. 汇总 ====================
summarise_draws <- function(mat) {
  tibble(
    sdi = rowMeans(mat),
    sdi_lower = apply(mat, 1, quantile, 0.025),
    sdi_upper = apply(mat, 1, quantile, 0.975)
  )
}
sdi_all <- map_dfr(locations, ~ {
  summarise_draws(sdi_draws[[.x]]) %>%
    mutate(location_name = .x, year_id = years_common)
})
## ==================== 9. 疾病数据读入（完全保留你原来的方式） ====================
all_cause1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 1/IHME-GBD_2023_DATA-522f4c75-1.csv")
all_cause2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 2/IHME-GBD_2023_DATA-9f89cdeb-1.csv")
all_cause3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 3/IHME-GBD_2023_DATA-5b64a7e5-1.csv")
Idiopathic1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 1/IHME-GBD_2023_DATA-6c53952c-1.csv")
Idiopathic2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 2/IHME-GBD_2023_DATA-65e64f67-1.csv")
Idiopathic3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 3/IHME-GBD_2023_DATA-5132aaa5-1.csv")
Secondary1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 1/IHME-GBD_2023_DATA-bd4ebaf5-1.csv")
Secondary2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 2/IHME-GBD_2023_DATA-fb6253af-1.csv")
Secondary3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 3/IHME-GBD_2023_DATA-0b56cca8-1.csv")
all_cause <- bind_rows(all_cause1, all_cause2, all_cause3)
Idiopathic <- bind_rows(Idiopathic1, Idiopathic2, Idiopathic3)
Secondary <- bind_rows(Secondary1, Secondary2, Secondary3)
all_cause$epilepsy_type <- "Total_Epilepsy"
Idiopathic$epilepsy_type <- "Idiopathic"
Secondary$epilepsy_type <- "Secondary"
epilepsy_all <- bind_rows(all_cause, Idiopathic, Secondary)
sdi_data_for_join <- sdi_all %>% rename(sdi_value = sdi)
epilepsy_all <- epilepsy_all %>%
  dplyr::select(measure_name, location_name, sex_name, age_name, metric_name, year, val, upper, lower, epilepsy_type) %>%
  filter(metric_name == "Rate") %>%
  mutate(location_name = gsub("People's Republic of China", "China", location_name)) %>%
  left_join(sdi_data_for_join, by = c("location_name", "year" = "year_id"))
rm(list = setdiff(ls(), c("sdi_all", "sdi_data_for_join", "epilepsy_all")))
library(tidyverse)
library(glmmTMB) # 混合效应
library(forecast) # ARIMA
library(foreach) # 循环
library(zoo) # na.approx
# 参数
future_years <- 2024:2050
n_draws <- 1 # 测试用 5，正式跑改成 1000
expo_val <- 100000 # 虽然不再用 count，但保留以备后用
# 准备 SDI
sdi_data <- sdi_all %>%
  dplyr::rename(location = location_name, year = year_id,
                sdi_mean = sdi, lower = sdi_lower, upper = sdi_upper) %>%
  tidyr::complete(location, year = 1990:2050) %>%
  dplyr::group_by(location) %>%
  dplyr::mutate(dplyr::across(c(sdi_mean, lower, upper),
                              ~ zoo::na.approx(., na.rm = FALSE, rule = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(sdi_mean))
# 主预测函数（已修改为按 location × sex 分别计算 BMA 权重）
run_forecast <- function(data, epilepsy_type_filter, measure_filter = "Prevalence") {
  bma_weight_store <- list() # 仍然保留，但现在会存更多分组的权重信息
  
  # 数据准备（保持不变）
  df <- data %>%
    dplyr::filter(epilepsy_type == epilepsy_type_filter,
                  measure_name == measure_filter,
                  metric_name == "Rate") %>%
    filter(location_name != "China") %>%
    dplyr::mutate(
      location = location_name,
      sex = sex_name,
      age_group = age_name,
      year = as.integer(year),
      rate = val,
      sdi = sdi_value
    ) %>%
    dplyr::select(location, sex, age_group, year, rate, lower, upper, sdi)
  
  if (nrow(df) == 0) {
    cat("No data for", epilepsy_type_filter, "-", measure_filter, "Rate\n")
    return(NULL)
  }
  
  train_df <- df %>% dplyr::filter(year <= 2010)
  withhold_df <- df %>% dplyr::filter(year >= 2011 & year <= 2023)
  
  cat("Training rows (≤2010):", nrow(train_df), "\n")
  cat("Withhold rows (2011-2023):", nrow(withhold_df), "\n")
  cat("Unique locations:", toString(unique(df$location)), "\n")
  
  # 未来网格（包括全国）
  future_grid <- tidyr::expand_grid(
    location = unique(df$location),
    sex = unique(df$sex),
    age_group = unique(df$age_group),
    year = future_years
  ) %>%
    dplyr::left_join(sdi_data, by = c("location", "year")) %>%
    dplyr::mutate(sdi = sdi_mean) %>%
    dplyr::filter(!is.na(sdi))
  
  cat("Future grid rows:", nrow(future_grid), "\n")
  
  # Monte Carlo 抽样（串行）
  mc_results <- foreach(draw = 1:n_draws, .combine = dplyr::bind_rows) %do% {
    cat("\n=== Starting Draw", draw, "/", n_draws, "===\n")
    
    # 1. 抽样 SDI（全网格）
    sdi_sd <- (sdi_data$upper - sdi_data$lower) / 3.92
    sdi_draw_all <- rnorm(nrow(sdi_data), sdi_data$sdi_mean, sdi_sd)
    sdi_draw_df <- sdi_data %>% dplyr::mutate(sdi_draw = sdi_draw_all)
    
    # 2. 抽样历史 rate（仅训练集）
    rate_sd_train <- (train_df$upper - train_df$lower) / 3.92
    rate_draw_train <- rnorm(nrow(train_df), train_df$rate, rate_sd_train)
    rate_draw_train <- pmax(rate_draw_train, 0)
    
    train_draw <- train_df %>%
      dplyr::mutate(rate = rate_draw_train)
    
    cat("train_draw rows after sampling:", nrow(train_draw), "\n")
    
    # 3. 为未来网格匹配本次抽样的 SDI
    future_grid_draw <- future_grid %>%
      dplyr::select(-sdi) %>%
      dplyr::left_join(sdi_draw_df %>% dplyr::select(location, year, sdi_draw),
                       by = c("location", "year")) %>%
      dplyr::rename(sdi = sdi_draw)
    
    if (any(is.na(future_grid_draw$sdi))) {
      future_grid_draw <- future_grid_draw %>%
        dplyr::mutate(sdi = tidyr::replace_na(sdi, mean(sdi, na.rm = TRUE)))
    }
    
    # ──────────────────────────────── 核心修改：按 location × sex 分组拟合 & BMA ───────
    result_list <- list()
    weight_records <- list()
    
    # 创建分组标识
    train_draw <- train_draw %>%
      mutate(group_id = paste(location, sex, sep = "___"))
    
    future_grid_draw <- future_grid_draw %>%
      mutate(group_id = paste(location, sex, sep = "___"))
    
    unique_groups <- unique(train_draw$group_id)
    
    for (grp in unique_groups) {
      cat("Processing group:", grp, "\n")
      
      train_sub <- train_draw %>% filter(group_id == grp)
      future_sub <- future_grid_draw %>% filter(group_id == grp)
      # 先准备两种 SDI：地区-年份特异 vs 全国年份平均
      base_models_sub <- list()
      tryCatch({
        # M1: 最简单模型
        base_models_sub$M1 <- glmmTMB::glmmTMB(
          rate ~ sdi,
          family = gaussian, data = train_sub
        )
        
        # M2: 和 M1 一样
        base_models_sub$M2 <- glmmTMB::glmmTMB(
          rate ~ sdi + age_group,
          family = gaussian, data = train_sub
        )
        
        # M3: 交互项
        base_models_sub$M3 <- glmmTMB::glmmTMB(
          rate ~ sdi * age_group,
          family = gaussian, data = train_sub
        )
      }, error = function(e) {
        cat("Model fit error in group", grp, ":", e$message, "\n")
      })
      
      valid_models <- names(base_models_sub)[!sapply(base_models_sub, is.null)]
      if (length(valid_models) < 2) {
        cat("Too few models in group", grp, "\n")
        next
      }
      
      # 计算 BIC 和 BMA 权重（不变）
      bics_sub <- sapply(base_models_sub[valid_models], BIC)
      delta <- bics_sub - min(bics_sub)
      weights_sub <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
      
      weight_records[[grp]] <- tibble(
        draw = draw,
        group_id = grp,
        model = valid_models,
        weight = weights_sub
      )
      
      # 基础预测（直接用 future_sub，不区分）
      pred_base_sub <- lapply(base_models_sub[valid_models], function(m) {
        predict(m, newdata = future_sub, type = "response", allow.new.levels = TRUE)
      })
      
      # ──────────────── 残差 ARIMA（按 group_id 做更精细） ────────────────
      resid_list_sub <- lapply(base_models_sub[valid_models], function(m) residuals(m, type = "response"))
      
      arima_resid_sub <- lapply(resid_list_sub, function(res) {
        train_sub %>%
          mutate(resid = res) %>%
          group_by(year) %>% # 同一个 group 内按年聚合残差（数据少时可不分组）
          summarise(resid_mean = mean(resid, na.rm = TRUE), .groups = "drop") %>%
          mutate(ts_res = ts(resid_mean, start = min(year), frequency = 1)) %>%
          { if (nrow(.) >= 3) {
            fit <- try(Arima(.$ts_res, order = c(0,1,0)), silent = TRUE)
            if (!inherits(fit, "try-error")) {
              fc <- forecast(fit, h = length(future_years))$mean
              tibble(year = future_years, resid_fc = as.numeric(fc))
            } else {
              tibble(year = future_years, resid_fc = 0)
            }
          } else {
            tibble(year = future_years, resid_fc = 0)
          }}
      })
      
      # 残差调整预测
      pred_resid_sub <- mapply(function(base_p, resid_df) {
        if (nrow(resid_df) == 0) return(base_p)
        future_sub %>%
          left_join(resid_df, by = "year") %>%
          mutate(resid_fc = replace_na(resid_fc, 0)) %>%
          pull(resid_fc) +
          base_p
      }, pred_base_sub, arima_resid_sub, SIMPLIFY = FALSE)
      
      # 合并基础 + 残差预测（这里假设基础和残差用相同权重）
      all_preds_sub <- c(pred_base_sub, pred_resid_sub)
      weights_all_sub <- rep(weights_sub, 2)
      weights_all_sub <- weights_all_sub / sum(weights_all_sub)
      
      # BMA 加权
      bma_rate_sub <- rowSums(do.call(cbind, Map(`*`, all_preds_sub, weights_all_sub)), na.rm = TRUE)
      
      result_list[[grp]] <- future_sub %>%
        mutate(
          pred_rate = as.numeric(bma_rate_sub),
          draw = draw
        ) %>%
        dplyr::select(location, sex, age_group, year, pred_rate, draw)
    }
    
    # 合并本 draw 所有分组的结果
    result_this_draw <- bind_rows(result_list)
    
    # 存权重记录（可选，用于事后分析不同组的权重差异）
    bma_weight_store[[draw]] <- bind_rows(weight_records)
    
    cat("Draw", draw, "completed - rows:", nrow(result_this_draw), "\n")
    return(result_this_draw)
  }
  
  if (is.null(mc_results) || nrow(mc_results) == 0) {
    cat("All draws failed\n")
    return(NULL)
  }
  
  cat("Monte Carlo completed:", nrow(mc_results), "rows\n")
  
  # 汇总权重（现在是按组的，可进一步 group_by group_id 分析）
  bma_weights_df <- bind_rows(bma_weight_store)
  bma_weights_summary <- bma_weights_df %>%
    group_by(group_id, model) %>%
    summarise(
      mean_weight = mean(weight, na.rm = TRUE),
      sd_weight = sd(weight, na.rm = TRUE),
      .groups = "drop"
    )%>%
    # 关键：拆分 group_id
    tidyr::separate(group_id, into = c("location", "sex"), sep = "___") %>%
    dplyr::select(location, sex, model, mean_weight, sd_weight)
  
  print("BMA weight summary by group:")
  print(bma_weights_summary)
  
  # 后续汇总、加全国平均、验证部分保持原样（或后续再同步修改验证为分组）
  predictions <- mc_results %>%
    group_by(location, sex, age_group, year) %>%
    summarise(
      pred_rate_mean = mean(pred_rate, na.rm = TRUE),
      lower = quantile(pred_rate, 0.025, na.rm = TRUE),
      upper = quantile(pred_rate, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      epilepsy_type = epilepsy_type_filter,
      measure_name = measure_filter,
      metric_name = "Rate"
    )
  
  # 全国平均（简单均值，注意现在已按组预测，全国需小心是否重复）
  
  national <- predictions %>%
    filter(location != "China") %>%
    group_by(sex, age_group, year, measure_name, metric_name, epilepsy_type) %>%
    summarise(
      pred_rate_mean = mean(pred_rate_mean, na.rm = TRUE),
      lower = mean(lower, na.rm = TRUE),
      upper = mean(upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(location = "China")
  
  predictions <- bind_rows(predictions, national)
  
  # ==================== 验证部分（关键修复）====================
  cat("\n=== 开始验证 ===\n")
  
  # 验证的核心问题：我们需要预测历史验证期（2011-2023），不仅仅是未来
  # 修改策略：重新运行模型来预测验证期
  
  # 准备验证期的 SDI 数据
  validation_years <- 2011:2023
  sdi_validation <- sdi_data %>%
    dplyr::filter(year %in% validation_years) %>%
    dplyr::select(location, year, sdi_mean, lower, upper)
  
  # 为验证期创建预测网格
  validation_grid <- tidyr::expand_grid(
    location = unique(df$location),
    sex = unique(df$sex),
    age_group = unique(df$age_group),
    year = validation_years
  ) %>%
    dplyr::left_join(sdi_validation, by = c("location", "year")) %>%
    dplyr::rename(sdi = sdi_mean) %>%
    dplyr::filter(!is.na(sdi))
  
  cat("Validation grid rows:", nrow(validation_grid), "\n")
  
  # 使用 Monte Carlo 样本来预测验证期
  validation_preds <- foreach(draw = 1:min(100, n_draws), .combine = dplyr::bind_rows) %do% {
    cat("Validation draw", draw, "/", min(100, n_draws), "\n")
    
    # 抽样 SDI 用于验证期
    sdi_sd_val <- (sdi_validation$upper - sdi_validation$lower) / 3.92
    sdi_draw_val <- rnorm(nrow(sdi_validation), sdi_validation$sdi_mean, sdi_sd_val)
    sdi_val_df <- sdi_validation %>% dplyr::mutate(sdi_draw = sdi_draw_val)
    
    # 为验证网格匹配抽样的 SDI
    validation_grid_draw <- validation_grid %>%
      dplyr::select(-sdi) %>%
      dplyr::left_join(sdi_val_df %>% dplyr::select(location, year, sdi_draw),
                       by = c("location", "year")) %>%
      dplyr::rename(sdi = sdi_draw)
    
    if (any(is.na(validation_grid_draw$sdi))) {
      validation_grid_draw <- validation_grid_draw %>%
        dplyr::mutate(sdi = tidyr::replace_na(sdi, mean(sdi, na.rm = TRUE)))
    }
    
    # 使用训练数据抽样
    rate_sd_train <- (train_df$upper - train_df$lower) / 3.92
    rate_draw_train <- rnorm(nrow(train_df), train_df$rate, rate_sd_train)
    rate_draw_train <- pmax(rate_draw_train, 0)
    
    train_draw_val <- train_df %>%
      dplyr::mutate(rate = rate_draw_train)
    
    # 拟合模型（使用相同的模型公式）
    tryCatch({
      m1_val <- glmmTMB::glmmTMB(rate ~ sdi + (1 | location),
                                 family = gaussian, data = train_draw_val)
      m2_val <- glmmTMB::glmmTMB(rate ~ sdi + age_group + (1 | location),
                                 family = gaussian, data = train_draw_val)
      m3_val <- glmmTMB::glmmTMB(rate ~ sdi * age_group + (1 | location),
                                 family = gaussian, data = train_draw_val)
      
      pred1 <- predict(m1_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      pred2 <- predict(m2_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      pred3 <- predict(m3_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      
      bics_val <- c(BIC(m1_val), BIC(m2_val), BIC(m3_val))
      bics_adj_val <- bics_val - min(bics_val, na.rm = TRUE)
      weights_val <- exp(-0.5 * bics_adj_val) / sum(exp(-0.5 * bics_adj_val))
      
      if(any(is.na(weights_val)) || sum(weights_val, na.rm = TRUE) == 0) {
        weights_val <- rep(1/3, 3)
      }
      
      bma_pred <- weights_val[1] * pred1 + weights_val[2] * pred2 + weights_val[3] * pred3
      
      validation_grid_draw %>%
        dplyr::mutate(pred_rate = as.numeric(bma_pred), draw = draw) %>%
        dplyr::select(location, sex, age_group, year, pred_rate, draw)
      
    }, error = function(e) {
      cat("Validation model error in draw", draw, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  if (!is.null(validation_preds) && nrow(validation_preds) > 0) {
    cat("Validation predictions generated:", nrow(validation_preds), "rows\n")
    
    # 汇总验证预测
    validation_summary <- validation_preds %>%
      dplyr::group_by(location, sex, age_group, year) %>%
      dplyr::summarise(
        pred_rate_mean = mean(pred_rate, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 与实际数据比较
    validation_comparison <- withhold_df %>%
      dplyr::mutate(
        location = as.character(location),
        sex = as.character(sex),
        age_group = as.character(age_group),
        year = as.integer(year)
      ) %>%
      dplyr::left_join(
        validation_summary %>%
          dplyr::mutate(
            location = as.character(location),
            sex = as.character(sex),
            age_group = as.character(age_group),
            year = as.integer(year)
          ),
        by = c("location", "sex", "age_group", "year")
      )
    
    cat("Matched rows for validation:", sum(!is.na(validation_comparison$pred_rate_mean)),
        "/", nrow(validation_comparison), "\n")
    
    # 计算验证指标
    if (sum(!is.na(validation_comparison$pred_rate_mean)) > 10) {
      validation_metrics <- validation_comparison %>%
        dplyr::filter(!is.na(pred_rate_mean)) %>%
        dplyr::summarise(
          rmse = sqrt(mean((rate - pred_rate_mean)^2, na.rm = TRUE)),
          mae = mean(abs(rate - pred_rate_mean), na.rm = TRUE),
          mape = mean(abs((rate - pred_rate_mean) / pmax(rate, 0.001)) * 100, na.rm = TRUE),
          n = n(),
          correlation = cor(rate, pred_rate_mean, use = "complete.obs"),
          .groups = "drop"
        ) %>%
        dplyr::mutate(epilepsy_type = epilepsy_type_filter,
                      measure_name = measure_filter)
      
      cat(sprintf("Validation results - RMSE: %.4f, MAE: %.4f, MAPE: %.2f%%, Corr: %.3f, n: %d\n",
                  validation_metrics$rmse, validation_metrics$mae,
                  validation_metrics$mape, validation_metrics$correlation,
                  validation_metrics$n))
      
      validation <- validation_metrics
    } else {
      cat("Insufficient matched data for validation\n")
      validation <- tibble(
        rmse = NA_real_, mae = NA_real_, mape = NA_real_,
        n = 0, correlation = NA_real_, epilepsy_type = epilepsy_type_filter,
        measure_name = measure_filter
      )
    }
  } else {
    cat("No validation predictions generated\n")
    validation <- tibble(
      rmse = NA_real_, mae = NA_real_, mape = NA_real_,
      n = 0, correlation = NA_real_, epilepsy_type = epilepsy_type_filter
    )
  }
  
  # ==================== 返回结果 ====================
  list(
    predictions = predictions,
    validation = validation, # 来自原代码
    bma_weights = bma_weights_summary # 现在包含 group_id
  )
}
# 执行：分别跑 Prevalence Rate 和 YLDs Rate
epilepsy_types <- c("Total_Epilepsy", "Idiopathic", "Secondary")
measures <- c("Prevalence", "YLDs (Years Lived with Disability)")
all_results <- list()
for (meas in measures) {
  cat("\n===== 正在预测:", meas, "Rate =====\n")
  res_list <- lapply(epilepsy_types, function(type) {
    run_forecast(epilepsy_all, type, measure_filter = meas)
  })
  all_results[[meas]] <- res_list
}
# 提取所有预测结果
final_predictions <- bind_rows(
  lapply(unlist(all_results, recursive = FALSE), `[[`, "predictions")
)
final_validations <- bind_rows(
  lapply(unlist(all_results, recursive = FALSE), `[[`, "validation")
)
# 合并权重（可选，按 measure 添加区分）
all_weights <- lapply(seq_along(all_results), function(i) {
  meas <- names(all_results)[i]
  bind_rows(lapply(all_results[[meas]], `[[`, "bma_weights")) %>%
    mutate(measure_name = meas)
}) %>% bind_rows()
# 保存
write_csv(final_predictions, "test epilepsy_2050_prevalence_ylds_rate_predictions.csv")
write_csv(final_validations, "test epilepsy_2050_prevalence_ylds_rate_validation.csv")
if (nrow(all_weights) > 0) {
  write_csv(all_weights, "test epilepsy_2050_prevalence_ylds_model_weights_by_group.csv")
  cat("Saved weights to: epilepsy_2050_prevalence_ylds_model_weights_by_group.csv\n")
}



########!!!Prevalence的number的预测############
all_cause1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 1/IHME-GBD_2023_DATA-522f4c75-1.csv")
all_cause2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 2/IHME-GBD_2023_DATA-9f89cdeb-1.csv")
all_cause3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/All cause 3/IHME-GBD_2023_DATA-5b64a7e5-1.csv")
Idiopathic1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 1/IHME-GBD_2023_DATA-6c53952c-1.csv")
Idiopathic2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 2/IHME-GBD_2023_DATA-65e64f67-1.csv")
Idiopathic3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/Non communicable 3/IHME-GBD_2023_DATA-5132aaa5-1.csv")
Secondary1 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 1/IHME-GBD_2023_DATA-bd4ebaf5-1.csv")
Secondary2 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 2/IHME-GBD_2023_DATA-fb6253af-1.csv")
Secondary3 <- read.csv("E:/研究生文件/研0/2023癫痫GBD/Materials/for projection/communicable 3/IHME-GBD_2023_DATA-0b56cca8-1.csv")
all_cause <- bind_rows(all_cause1, all_cause2, all_cause3)
Idiopathic <- bind_rows(Idiopathic1, Idiopathic2, Idiopathic3)
Secondary <- bind_rows(Secondary1, Secondary2, Secondary3)
all_cause$epilepsy_type <- "Total_Epilepsy"
Idiopathic$epilepsy_type <- "Idiopathic"
Secondary$epilepsy_type <- "Secondary"
epilepsy_all <- bind_rows(all_cause, Idiopathic, Secondary)
sdi_data_for_join <- sdi_all %>% rename(sdi_value = sdi)
epilepsy_all <- epilepsy_all %>%
  dplyr::select(measure_name, location_name, sex_name, age_name, metric_name, year, val, upper, lower, epilepsy_type) %>%
  filter(metric_name == "Number") %>%
  mutate(location_name = gsub("People's Republic of China", "China", location_name)) %>%
  left_join(sdi_data_for_join, by = c("location_name", "year" = "year_id")) %>%
  mutate(val = val/1000000,
         upper = upper/1000000,
         lower = lower/1000000)
library(tidyverse)
library(glmmTMB) # 混合效应
library(forecast) # ARIMA
library(foreach) # 循环
library(zoo) # na.approx
# 参数
future_years <- 2024:2050
n_draws <- 2 # 测试用 5，正式跑改成 1000
expo_val <- 100000 # 虽然不再用 count，但保留以备后用
# 准备 SDI
sdi_data <- sdi_all %>%
  dplyr::rename(location = location_name, year = year_id,
                sdi_mean = sdi, lower = sdi_lower, upper = sdi_upper) %>%
  tidyr::complete(location, year = 1990:2050) %>%
  dplyr::group_by(location) %>%
  dplyr::mutate(dplyr::across(c(sdi_mean, lower, upper),
                              ~ zoo::na.approx(., na.rm = FALSE, rule = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(sdi_mean))
# 主预测函数（已修改为按 location × sex 分别计算 BMA 权重）
run_forecast <- function(data, epilepsy_type_filter, measure_filter = "Prevalence") {
  bma_weight_store <- list() # 仍然保留，但现在会存更多分组的权重信息
  
  # 数据准备（保持不变）
  df <- data %>%
    dplyr::filter(epilepsy_type == epilepsy_type_filter,
                  measure_name == measure_filter,
                  metric_name == "Number") %>%
    filter(location_name != "China") %>%
    dplyr::mutate(
      location = location_name,
      sex = sex_name,
      age_group = age_name,
      year = as.integer(year),
      rate = val,
      sdi = sdi_value
    ) %>%
    dplyr::select(location, sex, age_group, year, rate, lower, upper, sdi)
  
  if (nrow(df) == 0) {
    cat("No data for", epilepsy_type_filter, "-", measure_filter, "Rate\n")
    return(NULL)
  }
  
  train_df <- df %>% dplyr::filter(year <= 2010)
  withhold_df <- df %>% dplyr::filter(year >= 2011 & year <= 2023)
  
  cat("Training rows (≤2010):", nrow(train_df), "\n")
  cat("Withhold rows (2011-2023):", nrow(withhold_df), "\n")
  cat("Unique locations:", toString(unique(df$location)), "\n")
  
  # 未来网格（包括全国）
  future_grid <- tidyr::expand_grid(
    location = unique(df$location),
    sex = unique(df$sex),
    age_group = unique(df$age_group),
    year = future_years
  ) %>%
    dplyr::left_join(sdi_data, by = c("location", "year")) %>%
    dplyr::mutate(sdi = sdi_mean) %>%
    dplyr::filter(!is.na(sdi))
  
  cat("Future grid rows:", nrow(future_grid), "\n")
  
  # Monte Carlo 抽样（串行）
  mc_results <- foreach(draw = 1:n_draws, .combine = dplyr::bind_rows) %do% {
    cat("\n=== Starting Draw", draw, "/", n_draws, "===\n")
    
    # 1. 抽样 SDI（全网格）
    sdi_sd <- (sdi_data$upper - sdi_data$lower) / 3.92
    sdi_draw_all <- rnorm(nrow(sdi_data), sdi_data$sdi_mean, sdi_sd)
    sdi_draw_df <- sdi_data %>% dplyr::mutate(sdi_draw = sdi_draw_all)
    
    # 2. 抽样历史 rate（仅训练集）
    rate_sd_train <- (train_df$upper - train_df$lower) / 3.92
    rate_draw_train <- rnorm(nrow(train_df), train_df$rate, rate_sd_train)
    rate_draw_train <- pmax(rate_draw_train, 0)
    
    train_draw <- train_df %>%
      dplyr::mutate(rate = rate_draw_train)
    
    cat("train_draw rows after sampling:", nrow(train_draw), "\n")
    
    # 3. 为未来网格匹配本次抽样的 SDI
    future_grid_draw <- future_grid %>%
      dplyr::select(-sdi) %>%
      dplyr::left_join(sdi_draw_df %>% dplyr::select(location, year, sdi_draw),
                       by = c("location", "year")) %>%
      dplyr::rename(sdi = sdi_draw)
    
    if (any(is.na(future_grid_draw$sdi))) {
      future_grid_draw <- future_grid_draw %>%
        dplyr::mutate(sdi = tidyr::replace_na(sdi, mean(sdi, na.rm = TRUE)))
    }
    
    # ──────────────────────────────── 核心修改：按 location × sex 分组拟合 & BMA ───────
    result_list <- list()
    weight_records <- list()
    
    # 创建分组标识
    train_draw <- train_draw %>%
      mutate(group_id = paste(location, sex, sep = "___"))
    
    future_grid_draw <- future_grid_draw %>%
      mutate(group_id = paste(location, sex, sep = "___"))
    
    unique_groups <- unique(train_draw$group_id)
    
    for (grp in unique_groups) {
      cat("Processing group:", grp, "\n")
      
      train_sub <- train_draw %>% filter(group_id == grp)
      future_sub <- future_grid_draw %>% filter(group_id == grp)
      # 先准备两种 SDI：地区-年份特异 vs 全国年份平均
      base_models_sub <- list()
      tryCatch({
        # M1: 最简单模型
        base_models_sub$M1 <- glmmTMB::glmmTMB(
          rate ~ sdi,
          family = gaussian, data = train_sub
        )
        
        # M2: 和 M1 一样
        base_models_sub$M2 <- glmmTMB::glmmTMB(
          rate ~ sdi + age_group,
          family = gaussian, data = train_sub
        )
        
        # M3: 交互项
        base_models_sub$M3 <- glmmTMB::glmmTMB(
          rate ~ sdi * age_group,
          family = gaussian, data = train_sub
        )
      }, error = function(e) {
        cat("Model fit error in group", grp, ":", e$message, "\n")
      })
      
      valid_models <- names(base_models_sub)[!sapply(base_models_sub, is.null)]
      if (length(valid_models) < 2) {
        cat("Too few models in group", grp, "\n")
        next
      }
      
      # 计算 BIC 和 BMA 权重（不变）
      bics_sub <- sapply(base_models_sub[valid_models], BIC)
      delta <- bics_sub - min(bics_sub)
      weights_sub <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
      
      weight_records[[grp]] <- tibble(
        draw = draw,
        group_id = grp,
        model = valid_models,
        weight = weights_sub
      )
      
      # 基础预测（直接用 future_sub，不区分）
      pred_base_sub <- lapply(base_models_sub[valid_models], function(m) {
        predict(m, newdata = future_sub, type = "response", allow.new.levels = TRUE)
      })
      
      # ──────────────── 残差 ARIMA（按 group_id 做更精细） ────────────────
      resid_list_sub <- lapply(base_models_sub[valid_models], function(m) residuals(m, type = "response"))
      
      arima_resid_sub <- lapply(resid_list_sub, function(res) {
        train_sub %>%
          mutate(resid = res) %>%
          group_by(year) %>% # 同一个 group 内按年聚合残差（数据少时可不分组）
          summarise(resid_mean = mean(resid, na.rm = TRUE), .groups = "drop") %>%
          mutate(ts_res = ts(resid_mean, start = min(year), frequency = 1)) %>%
          { if (nrow(.) >= 3) {
            fit <- try(Arima(.$ts_res, order = c(0,1,0)), silent = TRUE)
            if (!inherits(fit, "try-error")) {
              fc <- forecast(fit, h = length(future_years))$mean
              tibble(year = future_years, resid_fc = as.numeric(fc))
            } else {
              tibble(year = future_years, resid_fc = 0)
            }
          } else {
            tibble(year = future_years, resid_fc = 0)
          }}
      })
      
      # 残差调整预测
      pred_resid_sub <- mapply(function(base_p, resid_df) {
        if (nrow(resid_df) == 0) return(base_p)
        future_sub %>%
          left_join(resid_df, by = "year") %>%
          mutate(resid_fc = replace_na(resid_fc, 0)) %>%
          pull(resid_fc) +
          base_p
      }, pred_base_sub, arima_resid_sub, SIMPLIFY = FALSE)
      
      # 合并基础 + 残差预测（这里假设基础和残差用相同权重）
      all_preds_sub <- c(pred_base_sub, pred_resid_sub)
      weights_all_sub <- rep(weights_sub, 2)
      weights_all_sub <- weights_all_sub / sum(weights_all_sub)
      
      # BMA 加权
      bma_rate_sub <- rowSums(do.call(cbind, Map(`*`, all_preds_sub, weights_all_sub)), na.rm = TRUE)
      
      result_list[[grp]] <- future_sub %>%
        mutate(
          pred_rate = as.numeric(bma_rate_sub),
          draw = draw
        ) %>%
        dplyr::select(location, sex, age_group, year, pred_rate, draw)
    }
    
    # 合并本 draw 所有分组的结果
    result_this_draw <- bind_rows(result_list)
    
    # 存权重记录（可选，用于事后分析不同组的权重差异）
    bma_weight_store[[draw]] <- bind_rows(weight_records)
    
    cat("Draw", draw, "completed - rows:", nrow(result_this_draw), "\n")
    return(result_this_draw)
  }
  
  if (is.null(mc_results) || nrow(mc_results) == 0) {
    cat("All draws failed\n")
    return(NULL)
  }
  
  cat("Monte Carlo completed:", nrow(mc_results), "rows\n")
  
  # 汇总权重（现在是按组的，可进一步 group_by group_id 分析）
  bma_weights_df <- bind_rows(bma_weight_store)
  bma_weights_summary <- bma_weights_df %>%
    group_by(group_id, model) %>%
    summarise(
      mean_weight = mean(weight, na.rm = TRUE),
      sd_weight = sd(weight, na.rm = TRUE),
      .groups = "drop"
    )%>%
    # 关键：拆分 group_id
    tidyr::separate(group_id, into = c("location", "sex"), sep = "___") %>%
    dplyr::select(location, sex, model, mean_weight, sd_weight)
  
  print("BMA weight summary by group:")
  print(bma_weights_summary)
  
  # 后续汇总、加全国平均、验证部分保持原样（或后续再同步修改验证为分组）
  # 先算每个 draw 的全国总人数
  national_per_draw <- mc_results %>%
    filter(location != "China") %>%
    group_by(draw, sex, age_group, year) %>%
    summarise(
      national_number = sum(pred_rate, na.rm = TRUE),  # 每个 draw 的全国总人数
      .groups = "drop"
    ) %>%
    mutate(location = "China")
  
  # 汇总全国的预测统计量
  national_summary <- national_per_draw %>%
    group_by(location, sex, age_group, year) %>%
    summarise(
      pred_rate_mean = mean(national_number, na.rm = TRUE),
      lower = quantile(national_number, 0.025, na.rm = TRUE),
      upper = quantile(national_number, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      measure_name = measure_filter,
      metric_name = "Number",
      epilepsy_type = epilepsy_type_filter
    )
  
  # 最终 predictions = 省预测 + 全国预测
  predictions <- bind_rows(
    # 原省/地区预测（保持不变）
    mc_results %>%
      group_by(location, sex, age_group, year) %>%
      summarise(
        pred_rate_mean = mean(pred_rate, na.rm = TRUE),
        lower = quantile(pred_rate, 0.025, na.rm = TRUE),
        upper = quantile(pred_rate, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        epilepsy_type = epilepsy_type_filter,
        measure_name = measure_filter,
        metric_name = "Number"
      ),
    
    # 加上全国
    national_summary
  )
  
  # ==================== 验证部分（关键修复）====================
  cat("\n=== 开始验证 ===\n")
  
  # 验证的核心问题：我们需要预测历史验证期（2011-2023），不仅仅是未来
  # 修改策略：重新运行模型来预测验证期
  
  # 准备验证期的 SDI 数据
  validation_years <- 2011:2023
  sdi_validation <- sdi_data %>%
    dplyr::filter(year %in% validation_years) %>%
    dplyr::select(location, year, sdi_mean, lower, upper)
  
  # 为验证期创建预测网格
  validation_grid <- tidyr::expand_grid(
    location = unique(df$location),
    sex = unique(df$sex),
    age_group = unique(df$age_group),
    year = validation_years
  ) %>%
    dplyr::left_join(sdi_validation, by = c("location", "year")) %>%
    dplyr::rename(sdi = sdi_mean) %>%
    dplyr::filter(!is.na(sdi))
  
  cat("Validation grid rows:", nrow(validation_grid), "\n")
  
  # 使用 Monte Carlo 样本来预测验证期
  validation_preds <- foreach(draw = 1:min(100, n_draws), .combine = dplyr::bind_rows) %do% {
    cat("Validation draw", draw, "/", min(100, n_draws), "\n")
    
    # 抽样 SDI 用于验证期
    sdi_sd_val <- (sdi_validation$upper - sdi_validation$lower) / 3.92
    sdi_draw_val <- rnorm(nrow(sdi_validation), sdi_validation$sdi_mean, sdi_sd_val)
    sdi_val_df <- sdi_validation %>% dplyr::mutate(sdi_draw = sdi_draw_val)
    
    # 为验证网格匹配抽样的 SDI
    validation_grid_draw <- validation_grid %>%
      dplyr::select(-sdi) %>%
      dplyr::left_join(sdi_val_df %>% dplyr::select(location, year, sdi_draw),
                       by = c("location", "year")) %>%
      dplyr::rename(sdi = sdi_draw)
    
    if (any(is.na(validation_grid_draw$sdi))) {
      validation_grid_draw <- validation_grid_draw %>%
        dplyr::mutate(sdi = tidyr::replace_na(sdi, mean(sdi, na.rm = TRUE)))
    }
    
    # 使用训练数据抽样
    rate_sd_train <- (train_df$upper - train_df$lower) / 3.92
    rate_draw_train <- rnorm(nrow(train_df), train_df$rate, rate_sd_train)
    rate_draw_train <- pmax(rate_draw_train, 0)
    
    train_draw_val <- train_df %>%
      dplyr::mutate(rate = rate_draw_train)
    
    # 拟合模型（使用相同的模型公式）
    tryCatch({
      m1_val <- glmmTMB::glmmTMB(rate ~ sdi + (1 | location),
                                 family = gaussian, data = train_draw_val)
      m2_val <- glmmTMB::glmmTMB(rate ~ sdi + age_group + (1 | location),
                                 family = gaussian, data = train_draw_val)
      m3_val <- glmmTMB::glmmTMB(rate ~ sdi * age_group + (1 | location),
                                 family = gaussian, data = train_draw_val)
      
      pred1 <- predict(m1_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      pred2 <- predict(m2_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      pred3 <- predict(m3_val, newdata = validation_grid_draw, type = "response", allow.new.levels = TRUE)
      
      bics_val <- c(BIC(m1_val), BIC(m2_val), BIC(m3_val))
      bics_adj_val <- bics_val - min(bics_val, na.rm = TRUE)
      weights_val <- exp(-0.5 * bics_adj_val) / sum(exp(-0.5 * bics_adj_val))
      
      if(any(is.na(weights_val)) || sum(weights_val, na.rm = TRUE) == 0) {
        weights_val <- rep(1/3, 3)
      }
      
      bma_pred <- weights_val[1] * pred1 + weights_val[2] * pred2 + weights_val[3] * pred3
      
      validation_grid_draw %>%
        dplyr::mutate(pred_rate = as.numeric(bma_pred), draw = draw) %>%
        dplyr::select(location, sex, age_group, year, pred_rate, draw)
      
    }, error = function(e) {
      cat("Validation model error in draw", draw, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  if (!is.null(validation_preds) && nrow(validation_preds) > 0) {
    cat("Validation predictions generated:", nrow(validation_preds), "rows\n")
    
    # 汇总验证预测
    validation_summary <- validation_preds %>%
      dplyr::group_by(location, sex, age_group, year) %>%
      dplyr::summarise(
        pred_rate_mean = mean(pred_rate, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 与实际数据比较
    validation_comparison <- withhold_df %>%
      dplyr::mutate(
        location = as.character(location),
        sex = as.character(sex),
        age_group = as.character(age_group),
        year = as.integer(year)
      ) %>%
      dplyr::left_join(
        validation_summary %>%
          dplyr::mutate(
            location = as.character(location),
            sex = as.character(sex),
            age_group = as.character(age_group),
            year = as.integer(year)
          ),
        by = c("location", "sex", "age_group", "year")
      )
    
    cat("Matched rows for validation:", sum(!is.na(validation_comparison$pred_rate_mean)),
        "/", nrow(validation_comparison), "\n")
    
    # 计算验证指标
    if (sum(!is.na(validation_comparison$pred_rate_mean)) > 10) {
      validation_metrics <- validation_comparison %>%
        dplyr::filter(!is.na(pred_rate_mean)) %>%
        dplyr::summarise(
          rmse = sqrt(mean((rate - pred_rate_mean)^2, na.rm = TRUE)),
          mae = mean(abs(rate - pred_rate_mean), na.rm = TRUE),
          mape = mean(abs((rate - pred_rate_mean) / pmax(rate, 0.001)) * 100, na.rm = TRUE),
          n = n(),
          correlation = cor(rate, pred_rate_mean, use = "complete.obs"),
          .groups = "drop"
        ) %>%
        dplyr::mutate(epilepsy_type = epilepsy_type_filter,
                      measure_name = measure_filter)
      
      cat(sprintf("Validation results - RMSE: %.4f, MAE: %.4f, MAPE: %.2f%%, Corr: %.3f, n: %d\n",
                  validation_metrics$rmse, validation_metrics$mae,
                  validation_metrics$mape, validation_metrics$correlation,
                  validation_metrics$n))
      
      validation <- validation_metrics
    } else {
      cat("Insufficient matched data for validation\n")
      validation <- tibble(
        rmse = NA_real_, mae = NA_real_, mape = NA_real_,
        n = 0, correlation = NA_real_, epilepsy_type = epilepsy_type_filter,
        measure_name = measure_filter
      )
    }
  } else {
    cat("No validation predictions generated\n")
    validation <- tibble(
      rmse = NA_real_, mae = NA_real_, mape = NA_real_,
      n = 0, correlation = NA_real_, epilepsy_type = epilepsy_type_filter
    )
  }
  
  # ==================== 返回结果 ====================
  list(
    predictions = predictions,
    validation = validation, # 来自原代码
    bma_weights = bma_weights_summary # 现在包含 group_id
  )
}
# 执行：分别跑 Prevalence Rate 和 YLDs Rate
epilepsy_types <- c("Total_Epilepsy", "Idiopathic", "Secondary")
measures <- c("Prevalence", "YLDs (Years Lived with Disability)")
all_results <- list()
for (meas in measures) {
  cat("\n===== 正在预测:", meas, "Rate =====\n")
  res_list <- lapply(epilepsy_types, function(type) {
    run_forecast(epilepsy_all, type, measure_filter = meas)
  })
  all_results[[meas]] <- res_list
}
# 提取所有预测结果
final_predictions <- bind_rows(
  lapply(unlist(all_results, recursive = FALSE), `[[`, "predictions")
)
final_validations <- bind_rows(
  lapply(unlist(all_results, recursive = FALSE), `[[`, "validation")
)
# 合并权重（可选，按 measure 添加区分）
all_weights <- lapply(seq_along(all_results), function(i) {
  meas <- names(all_results)[i]
  bind_rows(lapply(all_results[[meas]], `[[`, "bma_weights")) %>%
    mutate(measure_name = meas)
}) %>% bind_rows()
# 保存
write_csv(final_predictions, "num epilepsy_2050_prevalence_ylds_rate_predictions.csv")
write_csv(final_validations, "num epilepsy_2050_prevalence_ylds_rate_validation.csv")
if (nrow(all_weights) > 0) {
  write_csv(all_weights, "num epilepsy_2050_prevalence_ylds_model_weights_by_group.csv")
  cat("Saved weights to: epilepsy_2050_prevalence_ylds_model_weights_by_group.csv\n")
}
