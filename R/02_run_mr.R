# ==============================================================================
#  02_run_mr.R
#  用 TwoSampleMR 标准流程跑 血浆 CPN2 → 三个肝病结局 的孟德尔随机化
#
#  主分析:IVW(单 SNP 时退化为 Wald ratio)
#  敏感性:MR-Egger、加权中位数、加权众数、Cochran Q、Egger 截距、leave-one-out
#  另跑一版:剔除带蛋白改变变异(PAV)的 rs3732477
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR)
})

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
stopifnot(dir.exists("data"))
dir.create("results", showWarnings = FALSE)

exp_dat <- fread("data/exposure_dat.tsv")
out_dat <- fread("data/outcome_dat.tsv")

exp_cols <- c("SNP","beta.exposure","se.exposure","pval.exposure","eaf.exposure",
              "effect_allele.exposure","other_allele.exposure","samplesize.exposure",
              "exposure","id.exposure")
out_cols <- c("SNP","beta.outcome","se.outcome","pval.outcome","eaf.outcome",
              "effect_allele.outcome","other_allele.outcome","samplesize.outcome",
              "outcome","id.outcome")

pav_snps <- exp_dat[pav == TRUE, SNP]

run_one <- function(exposure, outcome, tag) {
  dat <- harmonise_data(as.data.frame(exposure[, ..exp_cols]),
                        as.data.frame(outcome[,  ..out_cols]),
                        action = 2)          # 2 = 用频率推断回文 SNP,推不出就剔除
  dat <- dat[dat$mr_keep, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  res <- mr(dat, method_list = if (nrow(dat) >= 3)
              c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")
            else if (nrow(dat) == 2) c("mr_ivw") else c("mr_wald_ratio"))
  res <- generate_odds_ratios(res)

  het  <- tryCatch(mr_heterogeneity(dat),  error = function(e) NULL)
  plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  loo  <- tryCatch(mr_leaveoneout(dat),     error = function(e) NULL)

  list(tag = tag, n_snp = nrow(dat), dat = dat,
       res = as.data.table(res), het = as.data.table(het),
       plei = as.data.table(plei), loo = as.data.table(loo))
}

all_res <- list(); all_het <- list(); all_plei <- list(); all_loo <- list()

for (ph in unique(out_dat$id.outcome)) {
  o <- out_dat[id.outcome == ph]
  lab <- o$outcome[1]

  for (tag in c("全部工具变量", "剔除 PAV 变异")) {
    e <- if (tag == "全部工具变量") exp_dat else exp_dat[!SNP %in% pav_snps]
    oo <- if (tag == "全部工具变量") o else o[!SNP %in% pav_snps]
    r <- run_one(e, oo, tag)
    if (is.null(r)) next

    cat(sprintf("\n%s\n### %s (%s)  —  %s  [%d 个工具变量]\n",
                strrep("=", 92), lab, ph, tag, r$n_snp))
    print(r$res[, .(method, nsnp, b = round(b, 4), se = round(se, 4),
                    OR = round(or, 3), lo = round(or_lci95, 3),
                    hi = round(or_uci95, 3), pval = signif(pval, 3))])
    if (nrow(r$het)) for (i in seq_len(nrow(r$het)))
      cat(sprintf("  Cochran Q [%s] = %.2f (df=%d), P = %.3f\n",
                  r$het$method[i], r$het$Q[i], r$het$Q_df[i], r$het$Q_pval[i]))
    if (nrow(r$plei)) cat(sprintf("  MR-Egger 截距 = %.4f (se %.4f), P = %.3f\n",
                                  r$plei$egger_intercept[1], r$plei$se[1], r$plei$pval[1]))

    add <- function(x) if (!is.null(x) && nrow(x)) x[, `:=`(outcome_code = ph, subset = tag)] else NULL
    all_res[[length(all_res)+1]]  <- add(r$res)
    all_het[[length(all_het)+1]]  <- add(r$het)
    all_plei[[length(all_plei)+1]]<- add(r$plei)
    all_loo[[length(all_loo)+1]]  <- add(r$loo)
  }
}

fwrite(rbindlist(all_res,  fill = TRUE), "results/mr_results.csv")
fwrite(rbindlist(all_het,  fill = TRUE), "results/mr_heterogeneity.csv")
fwrite(rbindlist(all_plei, fill = TRUE), "results/mr_pleiotropy.csv")
fwrite(rbindlist(all_loo,  fill = TRUE), "results/mr_leaveoneout.csv")
cat("\n\n结果已写入 results/\n")
