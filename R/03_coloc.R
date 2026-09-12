# ==============================================================================
#  03_coloc.R
#  贝叶斯共定位:血浆 CPN2 与三个肝病结局是否共享同一个因果变异
#
#  目的:排除「遗传连锁」造成的假阳性 —— MR 显著可能只是因为 CPN2 的 cis-pQTL
#        与真正致病的变异处于连锁不平衡,而非 CPN2 本身致病(见手册 §03、§10)。
#
#  判据:PPH4 > 0.8         共享同一因果变异
#        PPH3 高而 PPH4 低  两个信号但因果变异不同 → 就是要排除的情况
#        PPH3 + PPH4 都低   区域内证据不足,功效不够,不能下结论
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(coloc)
})

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
stopifnot(dir.exists("data"))
dir.create("results", showWarnings = FALSE)

N_EXP    <- 35299L
CPN2_TSS <- 194351387L
WINDOW   <- 5e5L          # coloc 用 ±500 kb,比 MR 的 ±1 Mb 窄,减少噪声

OUTCOMES <- list(
  NAFLD                      = list(label = "非酒精性脂肪肝", ncase = 3504L, ncontrol = 496844L),
  CIRRHOSIS_BROAD            = list(label = "肝硬化",          ncase = 5545L, ncontrol = 494803L),
  C3_HEPATOCELLU_CARC_EXALLC = list(label = "肝细胞癌",        ncase =  947L, ncontrol = 378749L)
)

# ------------------------------------------------------------------------------
# 暴露端:去重、去缺失、限定窗口
# 用 rsID 作为两端的连接键 —— 两个数据库的 indel 写法不同,按位置+等位基因会漏配
# ------------------------------------------------------------------------------
exp <- fread("data/exposure_CPN2_cis.tsv")
exp <- exp[!is.na(rsids) & rsids != "NA" & rsids != "" &
           !is.na(Beta) & !is.na(SE) & SE > 0 &
           ImpMAF > 0.01 & ImpMAF < 0.99 &
           abs(Pos - CPN2_TSS) <= WINDOW]
exp <- exp[!duplicated(rsids)]
message(sprintf("暴露端可用变异: %s", format(nrow(exp), big.mark = ",")))

res_all <- list()

for (ph in names(OUTCOMES)) {
  info <- OUTCOMES[[ph]]
  out <- fread(sprintf("data/outcome_%s_cis.tsv", ph))
  setnames(out, "#chrom", "chrom")
  out <- out[!is.na(rsids) & rsids != "NA" & rsids != "" &
             !is.na(beta) & !is.na(sebeta) & sebeta > 0 &
             af_alt > 0.01 & af_alt < 0.99 &
             abs(pos - CPN2_TSS) <= WINDOW]
  out <- out[!duplicated(rsids)]

  m <- merge(exp[, .(snp = rsids, b1 = Beta, v1 = SE^2, maf = ImpMAF)],
             out[, .(snp = rsids, b2 = beta, v2 = sebeta^2)],
             by = "snp")
  message(sprintf("%-28s 两端共有变异: %s", ph, format(nrow(m), big.mark = ",")))
  if (nrow(m) < 100) { warning("共有变异过少,跳过"); next }

  # deCODE 的 beta 单位就是标准差(秩逆正态变换后),故 sdY = 1,
  # 不需要 coloc 从 MAF 和 varbeta 反推(反推会有额外误差)
  d1 <- list(beta = m$b1, varbeta = m$v1, snp = m$snp, MAF = m$maf,
             type = "quant", N = N_EXP, sdY = 1)
  d2 <- list(beta = m$b2, varbeta = m$v2, snp = m$snp, MAF = m$maf,
             type = "cc", N = info$ncase + info$ncontrol,
             s = info$ncase / (info$ncase + info$ncontrol))

  r <- coloc.abf(d1, d2)
  s <- as.list(r$summary)

  cat(sprintf("\n%s\n### %s (%s)   共有变异 %s\n", strrep("=", 78),
              info$label, ph, format(nrow(m), big.mark = ",")))
  cat(sprintf("  PPH0 %.3f  无关联\n  PPH1 %.3f  仅蛋白\n  PPH2 %.3f  仅疾病\n",
              s$PP.H0.abf, s$PP.H1.abf, s$PP.H2.abf))
  cat(sprintf("  PPH3 %.3f  两个信号,因果变异【不同】 ← 要排除的\n", s$PP.H3.abf))
  cat(sprintf("  PPH4 %.3f  两个信号,【共享】同一因果变异 ← 要的\n", s$PP.H4.abf))
  min_p_out <- min(2 * pnorm(-abs(m$b2 / sqrt(m$v2))))
  verdict <- if (s$PP.H4.abf > 0.8) "✅ 支持共定位"
             else if (s$PP.H3.abf > 0.8) "❌ 不同因果变异(遗传连锁)"
             else if (s$PP.H1.abf > 0.5) "❌ 仅暴露端有信号,结局端在该区域无关联"
             else if (s$PP.H2.abf > 0.5) "❌ 仅结局端有信号"
             else "⚠ 两端信号均不足,无法判断"
  cat(sprintf("  区域内结局端最强关联 P = %.2e(全基因组显著线 5e-08)\n", min_p_out))
  cat(sprintf("  判定:%s\n", verdict))

  res_all[[ph]] <- data.table(outcome_code = ph, outcome = info$label,
                              nsnps = nrow(m),
                              PPH0 = s$PP.H0.abf, PPH1 = s$PP.H1.abf, PPH2 = s$PP.H2.abf,
                              PPH3 = s$PP.H3.abf, PPH4 = s$PP.H4.abf, verdict = verdict)
}

fwrite(rbindlist(res_all), "results/coloc_results.csv")
cat("\n结果已写入 results/coloc_results.csv\n")
