# ==============================================================================
#  01_prepare_data.R
#  把 deCODE(暴露)和 FinnGen(结局)的 cis 区域数据整理成 TwoSampleMR 格式
#
#  暴露:血浆 CPN2 蛋白水平,deCODE SomaScan,N = 35,299,SeqId 6415_90
#  结局:FinnGen R12 三个肝病表型
#  区域:chr3 CPN2 TSS(194,351,387)± 1 Mb
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# 从项目根目录运行(cpn2/)。若当前在 R/ 下,自动上跳一级。
if (!dir.exists("data") && dir.exists("../data")) setwd("..")
stopifnot(dir.exists("data"))
dir.create("results", showWarnings = FALSE)

CPN2_TSS <- 194351387L
N_EXP    <- 35299L

OUTCOMES <- list(
  NAFLD                      = list(label = "非酒精性脂肪肝", ncase = 3504L, ncontrol = 496844L),
  CIRRHOSIS_BROAD            = list(label = "肝硬化",          ncase = 5545L, ncontrol = 494803L),
  C3_HEPATOCELLU_CARC_EXALLC = list(label = "肝细胞癌",        ncase =  947L, ncontrol = 378749L)
)

# ------------------------------------------------------------------------------
# 1. 读入
# ------------------------------------------------------------------------------
exp_raw <- fread("data/exposure_CPN2_cis.tsv")
message(sprintf("暴露端 cis 区域: %s 个变异", format(nrow(exp_raw), big.mark = ",")))

out_raw <- lapply(names(OUTCOMES), function(ph) {
  d <- fread(sprintf("data/outcome_%s_cis.tsv", ph))
  setnames(d, "#chrom", "chrom")
  d[, phenocode := ph]
  d
})
names(out_raw) <- names(OUTCOMES)

# ------------------------------------------------------------------------------
# 2. 工具变量
#
#    来源:deCODE 附表 ST08,位点 chr3_326,6 个条件独立的 cis-pQTL。
#    过 QC 只保留 MAF >= 0.01 的 3 个(其余 3 个 MAF 0.0003~0.007,
#    携带者过少,估计不稳,见 results/instrument_qc.csv)。
#
#    ⚠ 等位基因表示法在两个数据库里不一致,必须显式映射,不能按位置盲配:
#      rs34225900  deCODE  chr3:194350491  TAATG / TAATT   (左补齐 5 碱基)
#                  FinnGen chr3:194350495      G / T        (规范化后)
#      差异位就是 deCODE 串的第 5 位(194350491 + 4 = 194350495),
#      因此 deCODE 的 TAATG 对应 FinnGen 的 G。
# ------------------------------------------------------------------------------
INSTRUMENTS <- data.table(
  SNP        = c("rs3732477", "rs34225900", "rs1466733"),
  decode_pos = c(194341790L,  194350491L,   194400269L),
  decode_ea  = c("T",         "TAATG",      "G"),
  decode_oa  = c("C",         "TAATT",      "A"),
  fg_pos     = c(194341790L,  194350495L,   194400269L),
  fg_alt     = c("T",         "G",          "G"),   # = deCODE 效应等位基因
  fg_ref     = c("C",         "T",          "A"),
  rank       = c(1L, 2L, 3L),
  pav        = c(TRUE, FALSE, FALSE)                # 是否有相关的蛋白改变变异
)

# ------------------------------------------------------------------------------
# 3. 暴露端
# ------------------------------------------------------------------------------
exp_dat <- merge(
  INSTRUMENTS[, .(SNP, decode_pos, decode_ea, decode_oa, rank, pav)],
  exp_raw[, .(decode_pos = Pos, decode_ea = effectAllele, decode_oa = otherAllele,
              beta.exposure = Beta, se.exposure = SE, pval.exposure = Pval,
              eaf.exposure = ImpMAF)],
  by = c("decode_pos", "decode_ea", "decode_oa")
)
stopifnot(nrow(exp_dat) == nrow(INSTRUMENTS))

# 统一到 FinnGen 的等位基因写法,便于 harmonise
exp_dat <- merge(exp_dat, INSTRUMENTS[, .(SNP, fg_alt, fg_ref)], by = "SNP")
exp_dat[, `:=`(
  effect_allele.exposure = fg_alt,
  other_allele.exposure  = fg_ref,
  samplesize.exposure    = N_EXP,
  exposure               = "plasma CPN2",
  id.exposure            = "CPN2_deCODE",
  F_stat                 = (beta.exposure / se.exposure)^2
)]

# ------------------------------------------------------------------------------
# 4. 结局端
# ------------------------------------------------------------------------------
out_dat <- rbindlist(lapply(names(OUTCOMES), function(ph) {
  info <- OUTCOMES[[ph]]
  d <- merge(INSTRUMENTS[, .(SNP, pos = fg_pos, ref = fg_ref, alt = fg_alt)],
             out_raw[[ph]][, .(pos, ref, alt, beta, sebeta, pval, af_alt)],
             by = c("pos", "ref", "alt"))
  d[, .(SNP,
        beta.outcome           = beta,
        se.outcome             = sebeta,
        pval.outcome           = pval,
        eaf.outcome            = af_alt,
        effect_allele.outcome  = alt,
        other_allele.outcome   = ref,
        ncase.outcome          = info$ncase,
        ncontrol.outcome       = info$ncontrol,
        samplesize.outcome     = info$ncase + info$ncontrol,
        outcome                = info$label,
        id.outcome             = ph)]
}))
stopifnot(nrow(out_dat) == nrow(INSTRUMENTS) * length(OUTCOMES))

# ------------------------------------------------------------------------------
# 5. QC:等位基因频率交叉核对(两端频率应接近;相加≈1 说明基准反了)
# ------------------------------------------------------------------------------
qc <- merge(exp_dat[, .(SNP, rank, pav, F_stat, eaf.exposure)],
            out_dat[id.outcome == "NAFLD", .(SNP, eaf.outcome)], by = "SNP")
qc[, `:=`(freq_diff = abs(eaf.exposure - eaf.outcome),
          freq_sum  = eaf.exposure + eaf.outcome)]
qc[, flag := fifelse(abs(freq_sum - 1) < 0.1, "⚠ 可能基准相反",
              fifelse(freq_diff > 0.15, "⚠ 频率差异较大(人群差异?)", "OK"))]
setorder(qc, rank)
print(qc[, .(SNP, rank, PAV = pav, F = round(F_stat, 1),
             eaf_deCODE = round(eaf.exposure, 3),
             eaf_FinnGen = round(eaf.outcome, 3), flag)])
fwrite(qc, "results/instrument_qc.csv")

fwrite(exp_dat, "data/exposure_dat.tsv", sep = "\t")
fwrite(out_dat, "data/outcome_dat.tsv",  sep = "\t")
message("\n已写出 data/exposure_dat.tsv 与 data/outcome_dat.tsv")
