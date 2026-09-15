# ==============================================================================
#  04_reverse_mr.R
#  反向孟德尔随机化:非酒精性脂肪肝(暴露) → 血浆 CPN2 水平(结局)
#
#  动机:正向 MR + 共定位未发现 CPN2 因果影响肝病,但湿实验在脂肪肝病人的
#        EV 组分中观察到 CPN2 升高。若 CPN2 是疾病的【结果】而非原因,
#        反向 MR 应当出现信号。本脚本检验这一点。
#
#  暴露:FinnGen R12 NAFLD(3,504 例 / 496,844 对照),工具变量为全基因组显著位点
#  结局:deCODE SomaScan 血浆 CPN2,SeqId 6415_90,N ≈ 35,300
#
#  ⚠ 反向设计里 CPN2 是【结局】,结局端不需要工具变量,
#    因此正向分析中「只有 3 个 cis-pQTL 可用」的限制在此不适用。
#    NAFLD 的工具变量相对 CPN2 全部是 trans,且无需达到显著。
#
#  ⚠ 当前版本的工具变量来自文献已知的 NAFLD 主效应位点(经 FinnGen 验证显著),
#    不是全基因组无偏扫描。要做无偏版需下载 FinnGen 完整汇总统计,
#    见文件末尾 TODO。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR)
})

set.seed(42)   # 加权中位数/众数内部用 bootstrap,固定种子保证可重复

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
stopifnot(dir.exists("data"))
dir.create("results", showWarnings = FALSE)

FG_BASE  <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
PHENO    <- "NAFLD"
N_CASE   <- 3504L
N_CTRL   <- 496844L
P_THRESH <- 5e-8

# ------------------------------------------------------------------------------
# 1. 候选工具变量
#
#    位置均为 GRCh38。等位基因写法两库不一致,与 01_prepare_data.R 同样显式映射,
#    不按位置盲配:
#      rs58542926  FinnGen chr19:19268740  C / T
#                  deCODE  同位置 CGGAGCTGTATTTGCCTTCCATGGTGCA /
#                                TGGAGCTGTATTTGCCTTCCATGGTGCA  (右扩 27 碱基)
#                  仅首位不同,故 deCODE 的 T... 对应 FinnGen 的 T。
#      rs72613567  两库一致,TA / T(插入)。
# ------------------------------------------------------------------------------
LOCI <- data.table(
  SNP       = c("rs738409", "rs58542926", "rs72613567", "rs1260326"),
  gene      = c("PNPLA3",   "TM6SF2",     "HSD17B13",   "GCKR"),
  chr       = c(22L,        19L,          4L,           2L),
  pos       = c(43928847L,  19268740L,    87310240L,    27508073L),
  decode_ea = c("G", "TGGAGCTGTATTTGCCTTCCATGGTGCA", "TA", "C"),
  decode_oa = c("C", "CGGAGCTGTATTTGCCTTCCATGGTGCA", "T",  "T")
)

# ------------------------------------------------------------------------------
# 2. 暴露端(FinnGen NAFLD)—— tabix 远程取点,无需下载完整文件
# ------------------------------------------------------------------------------
fg_file <- "data/reverse_exposure_NAFLD.tsv"
if (!file.exists(fg_file)) {
  tbi <- file.path(tempdir(), sprintf("%s.tbi", PHENO))
  if (!file.exists(tbi))
    download.file(sprintf("%s/finngen_R12_%s.gz.tbi", FG_BASE, PHENO), tbi, quiet = TRUE)
  remote <- sprintf("%s/finngen_R12_%s.gz##idx##%s", FG_BASE, PHENO, tbi)

  fg <- rbindlist(lapply(seq_len(nrow(LOCI)), function(i) {
    reg <- sprintf("%d:%d-%d", LOCI$chr[i], LOCI$pos[i], LOCI$pos[i])
    txt <- system2("tabix", c(shQuote(remote), reg), stdout = TRUE)
    if (!length(txt)) return(NULL)
    d <- fread(text = paste(txt, collapse = "\n"), header = FALSE,
               col.names = c("chrom","pos","ref","alt","rsids","nearest_genes",
                             "pval","mlogp","beta","sebeta","af_alt",
                             "af_alt_cases","af_alt_controls"))
    d[rsids == LOCI$SNP[i]][1][, SNP := LOCI$SNP[i]][]
  }))
  fwrite(fg, fg_file, sep = "\t")
} 
fg <- fread(fg_file)

# ------------------------------------------------------------------------------
# 3. 结局端(deCODE CPN2)
#    909 MB 全基因组文件扫一次代价高,结果已缓存。缓存缺失时重扫。
# ------------------------------------------------------------------------------
DECODE_RAW <- "data/raw/deCODE_6415_90_CPN2_CPN2.txt.gz"

# 在 deCODE 全基因组文件里按 chr:pos 取指定变异。909 MB 扫一次约 1-2 分钟,
# 故结果落盘缓存;删掉缓存文件即可重扫。
scan_decode <- function(chrpos, out) {
  if (file.exists(out)) return(invisible(NULL))
  if (!file.exists(DECODE_RAW))
    stop("缺少 ", DECODE_RAW, "\n获取方式见 README「数据获取 · 暴露端 deCODE」。")
  message("扫描 deCODE 全基因组文件(约 1-2 分钟)…")
  pf <- tempfile(); writeLines(chrpos, pf)
  # $5 != $6 跳过 deCODE 中效应/非效应等位基因写成同一串的畸形行
  # $4 != "NA" 要求有 rsID,便于与暴露端按 rsID 合并
  prog <- sprintf('BEGIN{ OFS="\t"; while((getline l < "%s")>0) want[l]=1
      print "decode_chrpos","rsids","decode_ea","decode_oa","Beta","SE","Pval","N","ImpMAF" }
    NR==1{next}
    { k=$1":"$2; if (k in want && $5!=$6 && $4!="NA") print k,$4,$5,$6,$7,$10,$8,$11,$12 }', pf)
  cmd <- sprintf("gzcat %s | awk -F'\t' %s > %s",
                 shQuote(DECODE_RAW), shQuote(prog), shQuote(out))
  if (system(cmd) != 0L) { unlink(out); stop("deCODE 扫描失败") }
  invisible(NULL)
}

dc_file <- "data/reverse_decode_CPN2_lookup.tsv"
scan_decode(sprintf("chr%d:%d", LOCI$chr, LOCI$pos), dc_file)
dc <- fread(dc_file)
setnames(dc, c("decode_ea","decode_oa","rsids"), c("effectAllele","otherAllele","SNP"),
         skip_absent = TRUE)

# ------------------------------------------------------------------------------
# 4. 合并 + 工具变量筛选
# ------------------------------------------------------------------------------
dat <- Reduce(function(a, b) merge(a, b, by = "SNP"), list(
  LOCI[, .(SNP, gene, chr, pos)],
  fg[, .(SNP, fg_ref = ref, fg_alt = alt, beta.exposure = beta,
         se.exposure = sebeta, pval.exposure = pval, eaf.exposure = af_alt)],
  dc[, .(SNP, dc_ea = effectAllele, dc_oa = otherAllele,
         beta.outcome = Beta, se.outcome = SE, pval.outcome = Pval)]
))

# deCODE 效应等位基因对齐到 FinnGen alt。首位碱基一致即视为同一等位基因
# (TM6SF2 在 deCODE 中右扩了 27 个碱基,见上方说明)。
dat[, dc_matches_alt := substr(dc_ea, 1, nchar(fg_alt)) == fg_alt |
                        dc_ea == fg_alt]
stopifnot(all(dat$dc_matches_alt))   # 若为 FALSE 说明映射表需要修订

dat[, `:=`(
  effect_allele.exposure = fg_alt, other_allele.exposure = fg_ref,
  effect_allele.outcome  = fg_alt, other_allele.outcome  = fg_ref,
  samplesize.exposure    = N_CASE + N_CTRL,
  samplesize.outcome     = 35300L,
  exposure = "NAFLD", id.exposure = "NAFLD_FinnGenR12",
  outcome  = "plasma CPN2", id.outcome = "CPN2_deCODE"
)]

dat[, gw_sig := pval.exposure < P_THRESH]
cat("\n=== 候选工具变量 ===\n")
print(dat[, .(SNP, gene,
              NAFLD_beta = round(beta.exposure, 4), NAFLD_P = signif(pval.exposure, 3),
              CPN2_beta  = round(beta.outcome, 4),  CPN2_P  = signif(pval.outcome, 3),
              selected = fifelse(gw_sig, "yes", sprintf("no (P>%g)", P_THRESH)))])

MR_COLS <- c("SNP","beta.exposure","se.exposure","pval.exposure","eaf.exposure",
             "effect_allele.exposure","other_allele.exposure","samplesize.exposure",
             "exposure","id.exposure",
             "beta.outcome","se.outcome","pval.outcome",
             "effect_allele.outcome","other_allele.outcome","samplesize.outcome",
             "outcome","id.outcome")

# ------------------------------------------------------------------------------
# 5. MR
# ------------------------------------------------------------------------------
run_rev <- function(d, tag) {
  d <- as.data.frame(d[, ..MR_COLS])
  d$mr_keep <- TRUE
  k <- nrow(d)
  res <- mr(d, method_list = if (k >= 3)
              c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")
            else if (k == 2) "mr_ivw" else "mr_wald_ratio")
  res$analysis <- tag
  het <- if (k >= 3) mr_heterogeneity(d, method_list = "mr_ivw") else NULL
  plt <- if (k >= 3) mr_pleiotropy_test(d) else NULL
  loo <- if (k >= 3) mr_leaveoneout(d) else NULL
  list(res = res, het = het, plt = plt, loo = loo, n = k)
}

main <- run_rev(dat[gw_sig == TRUE], "全基因组显著位点")
cat(sprintf("\n=== 主分析:NAFLD → 血浆 CPN2 (n = %d) ===\n", main$n))
print(main$res[, c("method","nsnp","b","se","pval")], digits = 3)
if (!is.null(main$het))
  cat(sprintf("\nCochran Q = %.2f, df = %d, P = %.3f\n",
              main$het$Q, main$het$Q_df, main$het$Q_pval))
if (!is.null(main$plt))
  cat(sprintf("Egger 截距 = %.4f (SE %.4f), P = %.3f\n",
              main$plt$egger_intercept, main$plt$se, main$plt$pval))
if (!is.null(main$loo)) { cat("\n--- leave-one-out ---\n"); print(main$loo[, c("SNP","b","se","p")], digits = 3) }

# 敏感性:纳入 GCKR。GCKR 是著名的代谢多效性位点,对血浆蛋白组有广泛直接效应,
# 纳入它主要用于展示多效性会如何扭曲估计,不作为主结果。
sens <- run_rev(dat, "另加 GCKR(仅作多效性示例)")
cat(sprintf("\n=== 敏感性:纳入 GCKR (n = %d) ===\n", sens$n))
print(sens$res[, c("method","nsnp","b","se","pval")], digits = 3)
if (!is.null(sens$het))
  cat(sprintf("Cochran Q = %.2f, df = %d, P = %.4f\n",
              sens$het$Q, sens$het$Q_df, sens$het$Q_pval))

# 单 SNP Wald ratio,便于看每个位点各自说了什么
wald <- dat[, .(SNP, gene, gw_sig,
                ratio = beta.outcome / beta.exposure,
                ratio_se = se.outcome / abs(beta.exposure))]
wald[, `:=`(lo = ratio - 1.96 * ratio_se, hi = ratio + 1.96 * ratio_se)]
cat("\n=== 各位点单独的 Wald ratio(CPN2 SD / NAFLD log-odds)===\n")
print(wald[, .(SNP, gene, gw_sig, ratio = round(ratio, 4),
               CI = sprintf("%.3f, %.3f", lo, hi))])

fwrite(rbind(main$res, sens$res, fill = TRUE), "results/reverse_mr_results.csv")
fwrite(wald, "results/reverse_mr_wald.csv")
if (!is.null(main$loo)) fwrite(main$loo, "results/reverse_mr_leaveoneout.csv")
cat("\n已写出 results/reverse_mr_*.csv\n")

# ------------------------------------------------------------------------------
# TODO:全基因组无偏版
#   当前工具变量来自文献已知位点,可能遗漏 FinnGen 中其他显著的 NAFLD 位点。
#   无偏做法见 R/05_reverse_mr_liverfat.R —— 那一版改用 UKB 肝脏 MRI 脂肪含量
#   作暴露,做了全基因组扫描 + 距离剪枝,得到 8 个独立位点。
#   若要对 FinnGen NAFLD 本身做无偏扫描:
#     curl -O https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_NAFLD.gz   # 约 800 MB
# ------------------------------------------------------------------------------
