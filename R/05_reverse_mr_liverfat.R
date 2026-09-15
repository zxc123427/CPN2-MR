# ==============================================================================
#  05_reverse_mr_liverfat.R
#  反向 MR 之二:肝脏脂肪含量(暴露) → 血浆 CPN2 水平(结局)
#
#  相对 04_reverse_mr.R 的改进:
#    04 用 FinnGen NAFLD(3,504 例),工具变量取自文献已知位点,只有 3 个,
#       且全是肝细胞脂质处理基因。
#    05 改用 UK Biobank 肝脏 MRI 定量脂肪含量,连续表型、全基因组无偏扫描,
#       得到 8 个独立位点。
#
#  暴露:Haas et al. 2021, Cell Genomics(PMID 34957434)
#        GWAS Catalog GCST90029073,N = 32,974 欧洲人
#        harmonised 文件为 GRCh38 / 1-based(已实测 rs738408 = 22:43928850 确认)
#  结局:deCODE SomaScan 血浆 CPN2,SeqId 6415_90,N ≈ 35,300
#
#  ⚠ 暴露单位:GWAS Catalog 元数据未声明。从效应量推断为【标准化后的 SD】
#    (PNPLA3 每等位基因 0.195;若为绝对百分比则与已知效应量差一个数量级)。
#    正式引用前需回原文核实。
#
#  ⚠ 本机未装 plink,无法做 LD clumping。改用 ±1 Mb 距离剪枝取各位点 lead SNP。
#    距离剪枝在长 LD 区块可能把同一信号切成两个,chr19 的 rs56252442
#    (18.12 Mb)与 rs58542926(19.27 Mb)相距 1.15 Mb,需留意。
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })

set.seed(42)   # 加权中位数/众数内部用 bootstrap,固定种子保证可重复

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
stopifnot(dir.exists("data"))
dir.create("results", showWarnings = FALSE)

N_OUT <- 35300L
N_EXP <- 32974L

HAAS_RAW   <- "data/raw/GCST90029073_liverfat_Haas2021.h.tsv.gz"
DECODE_RAW <- "data/raw/deCODE_6415_90_CPN2_CPN2.txt.gz"
P_THRESH   <- 5e-8
PRUNE_BP   <- 1e6

# ------------------------------------------------------------------------------
# 0. 工具变量:全基因组显著位点 → ±1 Mb 距离剪枝
#
#    ⚠ 只取 hm_* 列(harmonise 后)定义等位基因、坐标与 beta。
#      同文件里不带 hm_ 前缀的是原始投稿值,两套的效应等位基因可能相反
#      (本文件首行:hm_effect_allele=AC/hm_beta=+0.00336 vs effect_allele=A/beta=-0.00336)。
#      SE 与 P 只有非前缀列,但二者与等位基因方向无关,可安全取用。
# ------------------------------------------------------------------------------
lead_file <- "data/reverse_exposure_liverfat_lead.tsv"
if (!file.exists(lead_file)) {
  if (!file.exists(HAAS_RAW))
    stop("缺少 ", HAAS_RAW, "\n获取方式见 README「反向暴露端 B · UKB 肝脏脂肪」。")
  message("扫描 Haas 肝脂汇总统计,取 P < ", P_THRESH, " …")
  prog <- sprintf('BEGIN{OFS="\t"; print "rsid","chr","pos","other_allele","effect_allele","beta","se","eaf","pval"}
    NR==1{next}
    $14!="NA" && $14+0 < %g && $7!="NA" && $21!="NA" && $3!="NA" && $4!="NA" { print $2,$3,$4,$5,$6,$7,$21,$11,$14 }', P_THRESH)
  gws <- fread(cmd = sprintf("gzcat %s | awk -F'\t' %s", shQuote(HAAS_RAW), shQuote(prog)))
  message("  全基因组显著变异: ", nrow(gws), " 个")

  # 贪心距离剪枝:按 P 排序,依次取 lead,排除同染色体 PRUNE_BP 内的其余变异。
  # 这是本机无 plink 时对 LD clumping 的替代。局限:长 LD 区块可能被切成多个
  # 「独立」位点而重复计权。核验办法是看相邻位点的单位点 Wald ratio 是否一致 ——
  # 同号同量级则可疑,反向则倾向确为独立信号。
  setorder(gws, pval)
  keep <- gws[0]
  for (i in seq_len(nrow(gws))) {
    r <- gws[i]
    if (nrow(keep) && keep[chr == r$chr & abs(pos - r$pos) < PRUNE_BP, .N]) next
    keep <- rbind(keep, r)
  }
  setorder(keep, chr, pos)
  fwrite(keep, lead_file, sep = "\t")
  message("  距离剪枝后独立位点: ", nrow(keep), " 个")
}
exp_raw <- fread(lead_file)

# ------------------------------------------------------------------------------
# 0b. 结局端:在 deCODE 全基因组文件里取这些位点
#
#     ⚠ 不按位置精确配。deCODE 对部分变异用左锚多碱基写法,坐标会偏移
#       (rs2250802:GWAS Catalog 112161596,deCODE 112161593)。
#       故按 ±25 bp 窗口粗取,再按 rsID 精确筛。
# ------------------------------------------------------------------------------
dc_file <- "data/reverse_decode_CPN2_liverfat.tsv"
if (!file.exists(dc_file)) {
  if (!file.exists(DECODE_RAW))
    stop("缺少 ", DECODE_RAW, "\n获取方式见 README「数据获取 · 暴露端 deCODE」。")
  message("扫描 deCODE 全基因组文件(约 1-2 分钟)…")
  wf <- tempfile()
  writeLines(sprintf("chr%s\t%d\t%d", exp_raw$chr, exp_raw$pos - 25L, exp_raw$pos + 25L), wf)
  prog <- sprintf('BEGIN{ OFS="\t"; while((getline l < "%s")>0){ split(l,a,"\t"); n++; C[n]=a[1]; L[n]=a[2]; H[n]=a[3] }
      print "decode_chrpos","rsids","decode_ea","decode_oa","Beta","SE","Pval","N","ImpMAF" }
    NR==1{next}
    $5==$6 || $4=="NA" { next }
    { for(i=1;i<=n;i++) if ($1==C[i] && $2>=L[i] && $2<=H[i]) { print $1":"$2,$4,$5,$6,$7,$10,$8,$11,$12; break } }', wf)
  cmd <- sprintf("gzcat %s | awk -F'\t' %s > %s",
                 shQuote(DECODE_RAW), shQuote(prog), shQuote(dc_file))
  if (system(cmd) != 0L) { unlink(dc_file); stop("deCODE 扫描失败") }
}
dc_all <- fread(dc_file)
# deCODE 的 rsids 字段可能是逗号分隔的多个 ID,逐个拆开后再与工具变量表匹配
dc_raw <- dc_all[vapply(strsplit(rsids, ","), function(v) any(v %in% exp_raw$rsid), TRUE)]
dc_raw[, rsids := vapply(strsplit(rsids, ","), function(v) v[v %in% exp_raw$rsid][1], "")]
stopifnot(!anyDuplicated(dc_raw$rsids))

GENE <- c(rs2642438 = "MTARC1", rs1229984 = "ADH1B",   rs112875651 = "TRIB1",
          rs2250802 = "GPAM",   rs56252442 = "MAST3",  rs58542926  = "TM6SF2",
          rs429358  = "APOE",   rs738408  = "PNPLA3")

# ------------------------------------------------------------------------------
# 1. 等位基因表示法核验
#
#    deCODE 对部分变异使用左锚/右扩的多碱基写法,与 GWAS Catalog 的单碱基写法
#    不一致(01_prepare_data.R 里 rs34225900 是同一类问题)。本处不按位置盲配,
#    而是定位两条串真正相异的那一位,再与暴露端等位基因比对。
#      rs2250802   deCODE TGAA / TGAG  → 第 4 位 A / G
#      rs58542926  deCODE T…… / C……(28 碱基)→ 第 1 位 T / C
# ------------------------------------------------------------------------------
discriminating <- function(ea, oa) {
  if (nchar(ea) != nchar(oa)) return(c(ea, oa))          # 真 indel,原样返回
  e <- strsplit(ea, "")[[1]]; o <- strsplit(oa, "")[[1]]
  i <- which(e != o)
  if (length(i) != 1L) return(c(ea, oa))                 # 非单点差异,原样返回
  c(e[i], o[i])
}
dc <- copy(dc_raw)
dc[, c("dc_ea1", "dc_oa1") := as.list(discriminating(decode_ea, decode_oa)), by = rsids]

dat <- merge(
  exp_raw[, .(SNP = rsid, chr, pos,
              effect_allele.exposure = effect_allele,
              other_allele.exposure  = other_allele,
              beta.exposure = beta, se.exposure = se,
              pval.exposure = pval, eaf.exposure = eaf)],
  dc[, .(SNP = rsids, dc_ea1, dc_oa1, decode_ea, decode_oa,
         beta.outcome = Beta, se.outcome = SE, pval.outcome = Pval)],
  by = "SNP")
stopifnot(nrow(dat) == nrow(exp_raw))

dat[, allele_state := fifelse(
  dc_ea1 == effect_allele.exposure & dc_oa1 == other_allele.exposure, "match",
  fifelse(dc_ea1 == other_allele.exposure & dc_oa1 == effect_allele.exposure,
          "flip", "nomatch"))]
cat("\n=== 等位基因核验 ===\n")
print(dat[, .(SNP, gene = GENE[SNP], exposure_EA_OA = paste0(effect_allele.exposure, "/", other_allele.exposure),
              decode_EA_OA = paste0(substr(decode_ea, 1, 8), "/", substr(decode_oa, 1, 8)),
              discrim = paste0(dc_ea1, "/", dc_oa1), allele_state)])
stopifnot(!any(dat$allele_state == "nomatch"))
dat[allele_state == "flip", beta.outcome := -beta.outcome]   # 对齐到暴露端效应等位基因

dat[, `:=`(
  effect_allele.outcome = effect_allele.exposure,
  other_allele.outcome  = other_allele.exposure,
  samplesize.exposure = N_EXP, samplesize.outcome = N_OUT,
  exposure = "liver fat", id.exposure = "LIVERFAT_Haas2021",
  outcome  = "plasma CPN2", id.outcome = "CPN2_deCODE",
  gene = GENE[SNP],
  F_stat = (beta.exposure / se.exposure)^2
)]
setorder(dat, chr, pos)

cat("\n=== 工具变量 ===\n")
print(dat[, .(SNP, gene, chr, F = round(F_stat), eaf = round(eaf.exposure, 3),
              liverfat_beta = round(beta.exposure, 4), liverfat_P = signif(pval.exposure, 2),
              CPN2_beta = round(beta.outcome, 4), CPN2_P = signif(pval.outcome, 2))])

MR_COLS <- c("SNP","beta.exposure","se.exposure","pval.exposure","eaf.exposure",
             "effect_allele.exposure","other_allele.exposure","samplesize.exposure",
             "exposure","id.exposure","beta.outcome","se.outcome","pval.outcome",
             "effect_allele.outcome","other_allele.outcome","samplesize.outcome",
             "outcome","id.outcome")

run_rev <- function(d, tag) {
  df <- as.data.frame(d[, ..MR_COLS]); df$mr_keep <- TRUE; k <- nrow(df)
  res <- mr(df, method_list = if (k >= 3)
              c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")
            else if (k == 2) "mr_ivw" else "mr_wald_ratio")
  res$analysis <- tag
  cat(sprintf("\n=== %s (n = %d) ===\n", tag, k))
  print(res[, c("method","nsnp","b","se","pval")], digits = 3)
  if (k >= 3) {
    h <- mr_heterogeneity(df, method_list = "mr_ivw")
    p <- mr_pleiotropy_test(df)
    cat(sprintf("Cochran Q = %.2f, df = %d, P = %.3g   |   Egger 截距 = %.4f (SE %.4f), P = %.3f\n",
                h$Q, h$Q_df, h$Q_pval, p$egger_intercept, p$se, p$pval))
    ci <- res$b[1] + c(-1.96, 1.96) * res$se[1]
    cat(sprintf("IVW 95%% CI: %.4f 到 %.4f (CPN2 SD / 肝脂 SD)\n", ci[1], ci[2]))
  }
  res
}

# ------------------------------------------------------------------------------
# 2. 主分析与敏感性
#
#    剔除理由:
#      ADH1B rs1229984 —— 乙醇脱氢酶。该位点通过【饮酒】影响肝脂,
#        用它做「非酒精性」脂肪肝的工具变量在概念上就不成立;
#        且冰岛人群 MAF 仅 0.004,deCODE 端 SE = 0.065,几乎不提供信息。
#      APOE rs429358 / TRIB1 rs112875651 —— 经典脂质多效性位点,
#        对血浆蛋白组有广泛直接效应(TRIB1 对 CPN2 P = 6.3e-10),
#        与 04 脚本里 GCKR 属同一类问题。
# ------------------------------------------------------------------------------
PLEIO <- c("rs429358", "rs112875651")
all_res <- rbind(
  run_rev(dat,                                          "主分析:全部 8 个位点"),
  run_rev(dat[SNP != "rs1229984"],                      "敏感性 A:剔除 ADH1B(酒精通路)"),
  run_rev(dat[SNP != "rs1229984" & !SNP %in% PLEIO],    "敏感性 B:再剔除 APOE + TRIB1"),
  fill = TRUE)

cat("\n=== leave-one-out(主分析)===\n")
loo <- mr_leaveoneout(transform(as.data.frame(dat[, ..MR_COLS]), mr_keep = TRUE))
loo$gene <- ifelse(loo$SNP == "All", "—", GENE[loo$SNP])
print(loo[, c("SNP","gene","b","se","p")], digits = 3)

wald <- dat[, .(SNP, gene, ratio = beta.outcome / beta.exposure,
                se = se.outcome / abs(beta.exposure))]
wald[, CI := sprintf("%.3f, %.3f", ratio - 1.96 * se, ratio + 1.96 * se)]
cat("\n=== 各位点单独 Wald ratio(CPN2 SD / 肝脂 SD)===\n")
print(wald[order(ratio), .(SNP, gene, ratio = round(ratio, 4), CI)])

fwrite(all_res, "results/reverse_mr_liverfat_results.csv")
fwrite(loo,     "results/reverse_mr_liverfat_leaveoneout.csv")
fwrite(wald,    "results/reverse_mr_liverfat_wald.csv")
fwrite(dat,     "results/reverse_mr_liverfat_instruments.csv")
cat("\n已写出 results/reverse_mr_liverfat_*.csv\n")
