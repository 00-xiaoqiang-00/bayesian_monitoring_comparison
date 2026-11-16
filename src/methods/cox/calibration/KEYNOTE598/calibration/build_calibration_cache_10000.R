if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Cox / KEYNOTE-598 calibration-cache build only.
# Uses explicit interim enrollment counts 228/340/454 and writes only inside
# this version directory, so a cache built with a different sample-size layout cannot be reused.
# build_calibration_cache_10000.R (Cox version - Parallel MCMC Calibration)
# 目标：对 H0 (Type I error) 和 H1 (Power) 进行两阶段联合校准
# 扫描逻辑：在共同 KEYNOTE598 标准下校准 T1E 和 Power

# 1. 加载依赖并配置编译环境
# Optional compiler settings may be supplied through R_MAKEVARS_USER.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- Sys.getenv("R_ENV_BIN", unset = "")
if (nzchar(r_env_bin) && dir.exists(r_env_bin)) {
  Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))
}

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(MASS)
library(survival)
library(brms)
library(parallel)
library(posterior)

cal_version <- "CAL10000_V2_INTERIM228_20260805"
cal_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(cal_file_arg) != 1L) stop("Run this calibration-cache file with Rscript.")
cal_script_file <- normalizePath(sub("^--file=", "", cal_file_arg), mustWork = TRUE)
cal_version_dir <- dirname(cal_script_file)
setwd(cal_version_dir)

# 载入数据生成函数
source("data_simulation_final.R", local = TRUE)

# --- 校准与模拟参数设置 ---
n_sim_env    <- Sys.getenv("N_SIM")
n_sim_h0_env <- Sys.getenv("N_SIM_H0")
n_sim_h1_env <- Sys.getenv("N_SIM_H1")

N.sim_h0 <- if (n_sim_h0_env != "") as.integer(n_sim_h0_env) else (if (n_sim_env != "") as.integer(n_sim_env) else 10000)
N.sim_h1 <- if (n_sim_h1_env != "") as.integer(n_sim_h1_env) else (if (n_sim_env != "") as.integer(n_sim_env) else 10000)
cal_smoke_test <- identical(Sys.getenv("CAL_SMOKE_TEST", unset = "0"), "1")
if (cal_smoke_test) {
  if (!is.finite(N.sim_h0) || !is.finite(N.sim_h1) ||
      N.sim_h0 < 1L || N.sim_h1 < 1L || N.sim_h0 > 20L || N.sim_h1 > 20L) {
    stop("Smoke mode requires N_SIM_H0 and N_SIM_H1 between 1 and 20.")
  }
} else if (N.sim_h0 != 10000L || N.sim_h1 != 10000L) {
  stop("Formal calibration requires exactly 10,000 H0 and 10,000 H1-standard replicates.")
}

cores_env <- Sys.getenv("NUM_CORES")
num_cores <- if (cores_env != "") as.integer(cores_env) else 30     # 并行核数

trial_settings <- keynote598_trial_settings()
N.max        <- trial_settings$N.max
cens_upper   <- trial_settings$cens_upper
interim_fracs<- trial_settings$interim_fracs
interim_counts <- c(228L, 340L, 454L)
stopifnot(
  length(interim_counts) == length(interim_fracs),
  all(interim_counts > 0L),
  all(interim_counts < N.max)
)
recruit.int  <- trial_settings$recruit.int
alpha_final  <- 0.05      # 频率派初始显著性水平
bayes_iter   <- 2000      # 每链总迭代数
bayes_warmup <- 1000      # 每链保留 1000；4 链共 4000 个 post-warmup draws
target_t1e <- 0.045
t1e_tolerance <- 0.010
cache_tag <- "KEYNOTE598_v2_interim228_H0_H1standard_cox_normal_0_1p5_iter2000_warmup1000"
cache_h0_path <- sprintf("calib_stats_H0_%s.rds", cache_tag)
cache_h1_path <- sprintf("calib_stats_H1_%s.rds", cache_tag)

writeLines(
  c(
    paste0("version=", cal_version),
    "project=cox_KEYNOTE598",
    "method=Cox",
    "case=KEYNOTE-598",
    "purpose=parameter_calibration_cache_generation",
    "scenarios=H0,H1-standard",
    "N_max=568",
    "interim_counts=228,340,454",
    "futility_rule=PP <= b_k",
    paste0("N_sim_H0=", N.sim_h0),
    paste0("N_sim_H1_standard=", N.sim_h1),
    paste0("cores=", num_cores),
    "mcmc_chains=4",
    paste0("mcmc_iter=", bayes_iter),
    paste0("mcmc_warmup=", bayes_warmup),
    "data_generator=data_simulation_final.R",
    "incompatible_227_cache_reuse=false",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  "RUN_METADATA.txt"
)

# H0 and H1_calibration arm parameters are selected inside the case generator.

# --- MCMC 诊断与后验概率 ---
extract_valid_pp <- function(fit) {
  if (inherits(fit, "try-error")) return(NA_real_)

  draws <- try(posterior::as_draws_array(fit, variable = "b_group"), silent = TRUE)
  if (inherits(draws, "try-error")) return(NA_real_)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NA_real_)
  }

  draw_matrix <- posterior::as_draws_matrix(draws)
  mean(draw_matrix[, "b_group"] < 0)
}

# --- 模拟与提取统计量函数 ---
run_simulation_and_extract_stats <- function(sim_array, base_model, N.sim, label) {
  cat(sprintf("\n【MCMC】开始对 [%s] 模拟数据进行采样并提取指标 (总次数: %d, 核数: %d)...\n", label, N.sim, num_cores))
  entry_times <- recruit.int * (0:(N.max - 1))
  
  stats_list <- mclapply(seq_len(N.sim), function(j) {
    df_j <- as.data.frame(sim_array[, , j])
    
    # 复合终点
    time_to_first <- pmin(df_j$death_t, df_j$progression_t)
    final_censor_t <- pmin(
      df_j$censor_t,
      keynote598_admin_followup(
        entry_times, trial_settings$final_analysis_time
      )
    )
    df_j$time   <- pmin(time_to_first, final_censor_t)
    df_j$status <- as.integer(time_to_first <= final_censor_t)

    interim_CP <- numeric(length(interim_fracs))
    interim_PP <- numeric(length(interim_fracs))
    simulation_ok <- TRUE
    
    original_censor_t <- df_j$censor_t
    
    for (k in seq_along(interim_fracs)) {
    n_k <- interim_counts[k]

      # 行政随访截断
      Time.current <- recruit.int * (n_k - 1)
      obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
      censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
      
      time_interim   <- pmin(time_to_first[1:n_k], censor_t_interim)
      status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
      
      valid_rows <- censor_t_interim > 0
      df_interim <- data.frame(
        time   = time_interim[valid_rows],
        status = status_interim[valid_rows],
        group  = df_j$group[1:n_k][valid_rows]
      )
      
      # 1. 频率派 CP 计算
      fit_interim <- try(coxph(Surv(time, status) ~ group, data = df_interim), silent = TRUE)
      if (!inherits(fit_interim, "try-error") && fit_interim$nevent > 0) {
        s <- summary(fit_interim)
        beta_hat <- s$coefficients[1, "coef"]
        se       <- s$coefficients[1, "se(coef)"]
        
        t_k <- n_k / N.max
        Z_k <- -beta_hat / se
        CP_k <- 1 - pnorm((qnorm(1 - alpha_final) - sqrt(t_k) * Z_k) / sqrt(1 - t_k))
        interim_CP[k] <- CP_k
      } else {
        interim_CP[k] <- 0
      }
      
      # 2. 计算贝叶斯派真实的 MCMC 后验概率 PP_k
      fit_bayes_interim <- try(
        update(
          base_model,
          newdata = df_interim,
          chains  = 4,
          cores   = 1,
          iter    = bayes_iter,
          warmup  = bayes_warmup,
          seed    = 12345 + j,
          refresh = 0,
          silent  = 2
        ),
        silent = TRUE
      )
      
      interim_PP[k] <- extract_valid_pp(fit_bayes_interim)
      if (is.na(interim_PP[k])) {
        simulation_ok <- FALSE
        break
      }
    }

    if (!simulation_ok) {
      return(list(ok = FALSE, CP = interim_CP, PP = interim_PP,
                  pval = NA_real_, PP_final = NA_real_))
    }
    
    # 最终分析 - 频率派
    fit_final <- try(coxph(Surv(time, status) ~ group, data = df_j), silent = TRUE)
    if (!inherits(fit_final, "try-error")) {
      s <- summary(fit_final)
      pval_final <- s$coefficients[1, "Pr(>|z|)"]
    } else {
      pval_final <- 1.0
    }
    
    # 最终分析 - 贝叶斯派真实 MCMC
    df_final_clean <- df_j[, c("time", "status", "group")]
    fit_bayes_final <- try(
      update(
        base_model,
        newdata = df_final_clean,
        chains  = 4,
        cores   = 1,
        iter    = bayes_iter,
        warmup  = bayes_warmup,
        seed    = 12345 + j,
        refresh = 0,
        silent  = 2
      ),
      silent = TRUE
    )
    PP_final <- extract_valid_pp(fit_bayes_final)
    if (is.na(PP_final)) {
      return(list(ok = FALSE, CP = interim_CP, PP = interim_PP,
                  pval = pval_final, PP_final = NA_real_))
    }
    
    list(ok = TRUE, CP = interim_CP, PP = interim_PP,
         pval = pval_final, PP_final = PP_final)
  }, mc.cores = num_cores)
  
  return(stats_list)
}

# --- 缓存加载或重新计算 ---
use_cache <- FALSE
force_recalculate <- Sys.getenv("FORCE_RECALC") == "TRUE"

if (!force_recalculate && file.exists(cache_h0_path) && file.exists(cache_h1_path)) {
  stats_list_H0_temp <- readRDS(cache_h0_path)
  stats_list_H1_temp <- readRDS(cache_h1_path)
  if (length(stats_list_H0_temp) == N.sim_h0 && length(stats_list_H1_temp) == N.sim_h1) {
    cat("\nCompatible MCMC cache found for the current N_SIM; loading it...\n")
    stats_list_H0 <- stats_list_H0_temp
    stats_list_H1 <- stats_list_H1_temp
    use_cache <- TRUE
  }
}

if (!use_cache) {
  # 生成 H0 模拟数据
  cat(sprintf("\n【1】在 H0（原假设）下生成 %d 次模拟数据集...\n", N.sim_h0))
  set.seed(999)
  sim_array_H0 <- data.simulation(
    N.sim = N.sim_h0, N.max = N.max, 
    cens_upper = cens_upper,
    scenario = "H0"
  )

  # 生成 H1 模拟数据
  cat(sprintf("\n【2】在 H1（标准疗效）下生成 %d 次模拟数据集...\n", N.sim_h1))
  set.seed(888)
  sim_array_H1 <- data.simulation(
    N.sim = N.sim_h1, N.max = N.max, 
    cens_upper = cens_upper,
    scenario = "H1_calibration"
  )

  # 预编译基础贝叶斯模型
  cat("\n正在为所有分析预编译通用的基础贝叶斯 Cox 模型...\n")
  df_dummy <- data.frame(
    time   = runif(N.max, 0.5, cens_upper),
    status = sample(c(0, 1), N.max, replace = TRUE),
    group  = sample(c(0, 1), N.max, replace = TRUE)
  )
  base_model <- brm(
    formula = time | cens(1 - status) ~ group,
    data    = df_dummy,
    family  = brmsfamily("cox"),
    prior   = prior(normal(0, 1.5), class = "b", coef = "group"),
    chains  = 4,
    cores   = 1,
    iter    = bayes_iter,
    warmup  = bayes_warmup,
    seed    = 12345,
    refresh = 0,
    silent  = 2
  )
  cat("Base model precompilation completed.即将开始并行采样分析。\n\n")

  # 运行 MCMC 并提取结果
  stats_list_H0 <- run_simulation_and_extract_stats(sim_array_H0, base_model, N.sim_h0, "H0")
  saveRDS(stats_list_H0, cache_h0_path)

  stats_list_H1 <- run_simulation_and_extract_stats(sim_array_H1, base_model, N.sim_h1, "H1")
  saveRDS(stats_list_H1, cache_h1_path)
}

# 过滤并提取结果
successful_runs_H0 <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list_H0)
N.valid_H0 <- length(successful_runs_H0)

successful_runs_H1 <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list_H1)
N.valid_H1 <- length(successful_runs_H1)

if (N.valid_H0 < N.sim_h0) {
  cat(sprintf("警告：H0 仿真有 %d 次报错，实际成功 %d 次。\n", N.sim_h0 - N.valid_H0, N.valid_H0))
}
if (N.valid_H1 < N.sim_h1) {
  cat(sprintf("警告：H1 仿真有 %d 次报错，实际成功 %d 次。\n", N.sim_h1 - N.valid_H1, N.valid_H1))
}

# --- 辅助评价函数 ---
get_h0_metrics <- function(lam, pp_thresh, stats_list, N.valid) {
    c_primes <- pp_thresh * (interim_counts / N.max)^lam
  h0_stopped_count <- 0
  h0_type_I_count  <- 0
  ess_sum <- 0
  for (j in seq_len(N.valid)) {
    stats <- stats_list[[j]]
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      h0_stopped_count <- h0_stopped_count + 1
        ess_sum <- ess_sum + interim_counts[stop_idx[1]]
    } else {
      ess_sum <- ess_sum + N.max
      if (!is.na(stats$PP_final) && stats$PP_final > pp_thresh) {
        h0_type_I_count <- h0_type_I_count + 1
      }
    }
  }
  return(list(
    t1e = h0_type_I_count / N.valid, 
    stop_rate = h0_stopped_count / N.valid,
    ess = ess_sum / N.valid
  ))
}

evaluate_h1 <- function(lam, pp_thresh, stats_list_H1, N.valid_H1) {
    c_primes <- pp_thresh * (interim_counts / N.max)^lam
  h1_stopped_count <- 0
  h1_success_count <- 0
  ess_sum <- 0
  for (j in seq_len(N.valid_H1)) {
    stats <- stats_list_H1[[j]]
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      h1_stopped_count <- h1_stopped_count + 1
        ess_sum <- ess_sum + interim_counts[stop_idx[1]]
    } else {
      ess_sum <- ess_sum + N.max
      if (!is.na(stats$PP_final) && stats$PP_final > pp_thresh) {
        h1_success_count <- h1_success_count + 1
      }
    }
  }
  return(list(
    Power = h1_success_count / N.valid_H1,
    H1_stop_rate = h1_stopped_count / N.valid_H1,
    ess = ess_sum / N.valid_H1
  ))
}

# --- 网格扫描评估与输出 ---
cat("\n【3】开始执行决策参数网格扫描并保存校准结果...\n")

PP_thresh_grid <- seq(0.85, 0.99, by = 0.01)
lambda_grid <- seq(0.20, 2.00, by = 0.05)

candidates_df <- data.frame(
  PP_Thresh = numeric(),
  Lambda = numeric(),
  H0_T1E = numeric(),
  H0_StopRate = numeric(),
  H0_ESS = numeric(),
  H1_Power = numeric(),
  H1_StopRate = numeric(),
  H1_ESS = numeric(),
  stringsAsFactors = FALSE
)

for (pp_thresh in PP_thresh_grid) {
  for (lam in lambda_grid) {
    h0_res <- get_h0_metrics(lam, pp_thresh, successful_runs_H0, N.valid_H0)
    h1_res <- evaluate_h1(lam, pp_thresh, successful_runs_H1, N.valid_H1)
    
    candidates_df <- rbind(candidates_df, data.frame(
      PP_Thresh = pp_thresh,
      Lambda = lam,
      H0_T1E = h0_res$t1e,
      H0_StopRate = h0_res$stop_rate,
      H0_ESS = h0_res$ess,
      H1_Power = h1_res$Power,
      H1_StopRate = h1_res$H1_stop_rate,
      H1_ESS = h1_res$ess
    ))
  }
}

# 方法间公平性优先：先对齐共同的实际 T1E，再比较 H0 ESS。
matching_candidates <- subset(
  candidates_df,
  abs(H0_T1E - target_t1e) <= t1e_tolerance & H1_Power >= 0.80
)

# 进行自动参数推荐：
# 推荐原则：
# 1. 优先使 T1E 落在共同目标 target_t1e 的容差范围内，且 Power >= 80%
# 2. 先最小化与共同 T1E 目标的距离，再比较 H0_ESS 和 Power。
if (nrow(matching_candidates) > 0) {
  recommended_idx <- order(abs(matching_candidates$H0_T1E - target_t1e),
                           matching_candidates$H0_ESS,
                           -matching_candidates$H1_Power)[1]
  rec_param <- matching_candidates[recommended_idx, ]
  
  rec_text <- sprintf(
    "【系统推荐决策参数】\nPP_threshold = %.2f\nlambda'      = %.2f\n\n推荐理由：\n- 一类错误率（T1E）：%.1f%%，与方法比较的共同目标对齐。\n- 检验功效（Power）：%.1f%%。\n- H0 期中无效停试率（PET）：%.1f%%，H0 期望样本量（ESS）：%.1f 人，较最大样本量节省 %.1f%%。",
    rec_param$PP_Thresh, rec_param$Lambda, rec_param$H0_T1E * 100, rec_param$H1_Power * 100, rec_param$H0_StopRate * 100, rec_param$H0_ESS, (1 - rec_param$H0_ESS / N.max) * 100
  )
} else {
  rec_text <- sprintf("警告：未找到 T1E 在 %.3f +/- %.3f 且 Power >= 80%% 的候选参数组合。",
                      target_t1e, t1e_tolerance)
}

# 写入 txt 结果文件
output_path <- "results/top_candidates_KEYNOTE598_Cox.txt"
dir.create("results", showWarnings = FALSE)

# 构造表格文本
write_output <- function() {
  cat("=========================================================================\n")
  cat("                   自适应设计决策参数校准网格报告\n")
  cat("=========================================================================\n")
  cat(sprintf("有效模拟次数：H0 = %d, H1 = %d\n", N.valid_H0, N.valid_H1))
  cat(sprintf("筛选条件：H0 T1E = %.3f +/- %.3f 且 H1 Power >= 80%%\n\n",
              target_t1e, t1e_tolerance))
  
  cat("【所有候选参数表现列表】\n")
  cat(sprintf("%-10s %-8s %-10s %-12s %-10s %-10s %-12s %-10s\n", 
              "PP_Thresh", "Lambda", "H0_T1E", "H0_StopRate", "H0_ESS", "H1_Power", "H1_StopRate", "H1_ESS"))
  for (i in 1:nrow(candidates_df)) {
    row <- candidates_df[i, ]
    cat(sprintf("%-10.2f %-8.2f %-10s %-12s %-10.1f %-10s %-12s %-10.1f\n",
                row$PP_Thresh, row$Lambda, 
                sprintf("%.1f%%", row$H0_T1E * 100), 
                sprintf("%.1f%%", row$H0_StopRate * 100), 
                row$H0_ESS, 
                sprintf("%.1f%%", row$H1_Power * 100), 
                sprintf("%.1f%%", row$H1_StopRate * 100), 
                row$H1_ESS))
  }
  
  cat("\n=========================================================================\n")
  cat("                        最优候选参数及系统推荐\n")
  cat("=========================================================================\n")
  cat(rec_text)
  cat("\n=========================================================================\n")
}

# 同时输出到控制台和文件
sink(output_path)
write_output()
sink()

# 仅输出推荐结果到终端
write_output()
cat(sprintf("\n完整的校准评估列表已保存至：%s\n", output_path))
