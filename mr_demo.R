# ============================================================
#  两样本 MR 手算演示  ——  CPN2 表达量  →  肝癌
#  纯 base R,不需要任何包。目的是看清每一步在算什么。
# ============================================================

# ---------- 表1:暴露(CPN2 的 cis-eQTL,来自 eQTLGen 之类) ----------
# beta = 每多带一个 effect_allele,CPN2 表达量升高几个标准差
exp_dat <- data.frame(
  SNP           = c("rs12345", "rs67890", "rs11111"),
  effect_allele = c("G",       "T",       "A"),
  other_allele  = c("A",       "C",       "G"),
  beta          = c( 0.35,      0.22,     -0.28),
  se            = c( 0.04,      0.05,      0.045),
  eaf           = c( 0.28,      0.41,      0.15),
  stringsAsFactors = FALSE
)

# ---------- 表2:结局(肝癌 GWAS,来自 FinnGen 之类) ----------
# 注意 rs67890:这里报的是 C/T,和暴露表的 T/C 反过来了 —— 需要翻转
out_dat <- data.frame(
  SNP           = c("rs12345", "rs67890", "rs11111"),
  effect_allele = c("G",       "C",       "A"),
  other_allele  = c("A",       "T",       "G"),
  beta          = c( 0.070,    -0.048,    -0.052),
  se            = c( 0.030,     0.032,     0.031),
  stringsAsFactors = FALSE
)

cat("\n===== 第0步:两张原始表 =====\n")
cat("\n[暴露 CPN2 表达量]\n"); print(exp_dat)
cat("\n[结局 肝癌]\n");        print(out_dat)


# ---------- 第1步:检查工具变量强度 F ----------
# F = (beta/se)^2 。F < 10 叫"弱工具变量",会让结果偏向观察性关联,必须剔除
cat("\n===== 第1步:工具变量强度 F 统计量 =====\n")
exp_dat$F <- (exp_dat$beta / exp_dat$se)^2
for (i in 1:nrow(exp_dat)) {
  cat(sprintf("  %-9s F = %6.1f   %s\n", exp_dat$SNP[i], exp_dat$F[i],
              ifelse(exp_dat$F[i] > 10, "OK 合格", "!! 弱工具,剔除")))
}


# ---------- 第2步:harmonize 对齐效应等位基因 ----------
cat("\n===== 第2步:harmonize(把两张表对齐到同一个字母) =====\n")
dat <- merge(exp_dat, out_dat, by = "SNP", suffixes = c(".exp", ".out"))
dat <- dat[match(exp_dat$SNP, dat$SNP), ]   # 保持原顺序

flip <- function(a) c(A="T", T="A", G="C", C="G")[a]

for (i in 1:nrow(dat)) {
  ea_e <- dat$effect_allele.exp[i]; oa_e <- dat$other_allele.exp[i]
  ea_o <- dat$effect_allele.out[i]; oa_o <- dat$other_allele.out[i]

  palindromic <- (ea_e == flip(oa_e))   # A/T 或 C/G:两条链读出来一样,无法判断

  if (palindromic) {
    cat(sprintf("  %-9s %s/%s  -> !! 回文 SNP,无法判断链方向,标准做法是剔除\n",
                dat$SNP[i], ea_e, oa_e))
  } else if (ea_e == ea_o && oa_e == oa_o) {
    cat(sprintf("  %-9s %s/%s vs %s/%s  -> 情况1:完全一致,不动\n",
                dat$SNP[i], ea_e, oa_e, ea_o, oa_o))
  } else if (ea_e == oa_o && oa_e == ea_o) {
    cat(sprintf("  %-9s %s/%s vs %s/%s  -> 情况2:等位基因互换(同一条链,选的效应字母不同)\n",
                dat$SNP[i], ea_e, oa_e, ea_o, oa_o))
  } else if (ea_e == flip(ea_o) && oa_e == flip(oa_o)) {
    # 互补链:字母换成互补碱基即可,beta 符号不变
    cat(sprintf("  %-9s %s/%s vs %s/%s  -> 情况3:链翻转,换成互补碱基,beta 不变号\n",
                dat$SNP[i], ea_e, oa_e, ea_o, oa_o))
    dat$effect_allele.out[i] <- flip(ea_o); dat$other_allele.out[i] <- flip(oa_o)
    ea_o <- dat$effect_allele.out[i]; oa_o <- dat$other_allele.out[i]
  }
  # 对齐后如果效应等位基因仍相反,则 beta 变号
  if (ea_e == dat$other_allele.out[i] && oa_e == dat$effect_allele.out[i]) {
    cat(sprintf("             效应等位基因相反 -> beta.out %+.3f 翻转为 %+.3f  ★\n",
                dat$beta.out[i], -dat$beta.out[i]))
    dat$beta.out[i] <- -dat$beta.out[i]
    dat$effect_allele.out[i] <- ea_e; dat$other_allele.out[i] <- oa_e
  }
}


# ---------- 第3步:每个 SNP 各算一个 Wald ratio ----------
cat("\n===== 第3步:Wald ratio(结局beta ÷ 暴露beta) =====\n")
dat$wald    <- dat$beta.out / dat$beta.exp
dat$wald_se <- dat$se.out   / abs(dat$beta.exp)   # 一阶 delta 法

for (i in 1:nrow(dat)) {
  cat(sprintf("  %-9s  %+.3f / %+.3f = %+.4f   (se %.4f, OR %.3f)\n",
              dat$SNP[i], dat$beta.out[i], dat$beta.exp[i],
              dat$wald[i], dat$wald_se[i], exp(dat$wald[i])))
}


# ---------- 第4步:IVW 把多个 SNP 合并 ----------
# 本质就是一个固定效应 meta 分析:权重 = 1/se^2 ,越准的 SNP 说话越算数
cat("\n===== 第4步:IVW 合并 =====\n")
w        <- 1 / dat$wald_se^2
beta_ivw <- sum(w * dat$wald) / sum(w)
se_ivw   <- sqrt(1 / sum(w))
z        <- beta_ivw / se_ivw
p        <- 2 * pnorm(-abs(z))
lo <- beta_ivw - 1.96 * se_ivw
hi <- beta_ivw + 1.96 * se_ivw

for (i in 1:nrow(dat))
  cat(sprintf("  %-9s 权重 %6.1f  (占 %4.1f%%)\n",
              dat$SNP[i], w[i], 100 * w[i] / sum(w)))

cat("\n----------------------------------------------------------\n")
cat(sprintf("  IVW beta = %+.4f   se = %.4f   Z = %.2f   P = %.4f\n",
            beta_ivw, se_ivw, z, p))
cat(sprintf("  OR = %.3f  (95%% CI %.3f - %.3f)\n", exp(beta_ivw), exp(lo), exp(hi)))
cat("----------------------------------------------------------\n")

cat("\n【结论怎么读】\n")
cat(sprintf("  CPN2 表达量每升高 1 个标准差,肝癌风险变为 %.2f 倍(P = %.4f)。\n",
            exp(beta_ivw), p))
cat("  95%CI 不跨过 1,统计学上显著。\n\n")
