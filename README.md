# 血浆 CPN2 与肝病的孟德尔随机化

检验血浆 CPN2 蛋白水平是否因果影响 非酒精性脂肪肝 → 肝硬化 → 肝细胞癌。

## 结论(2026-09-12)

**没有找到 CPN2 因果影响这三个肝病结局的证据。**

共定位进一步显示:**CPN2 位点上根本不存在肝病的遗传关联信号** —— ±500 kb
区域内结局端最强关联 P ≈ 1×10⁻³,离全基因组显著线(5×10⁻⁸)差五个数量级。
因此 MR 中 NAFLD 那个边缘显著(P = 0.040)不应解读为阳性。

## 数据来源

| 角色 | 来源 | 规模 |
|---|---|---|
| 暴露 | deCODE SomaScan 血浆蛋白组,SeqId `6415_90` | 35,299 人 |
| 结局 | FinnGen R12(tabix 远程取 chr3 区域) | NAFLD 3,504 / 肝硬化 5,545 / 肝癌 947 例 |

区域:chr3 CPN2 TSS(194,351,387,GRCh38)± 1 Mb。

**为什么用 deCODE**:UKB-PPP、Sun_2018、eQTLGen 均无 CPN1/CPN2。eQTLGen 是下载
完整显著 cis-eQTL 表(1,050 万条、16,923 个基因)逐行核对的,以 ALB / IL6R /
CLEC12A 作阳性对照确认流程有效。缺席的原因是生物学的:两基因在肝细胞转录、
蛋白分泌入血,全血本身不表达。

> ⚠ **CPN1 在 SomaScan 平台上命名为 `CBPN`**,按基因名检索会漏掉。

## 工具变量

deCODE 附表 ST08 列出 6 个条件独立 cis-pQTL,过 MAF ≥ 0.01 后保留 3 个:

| rsID | Beta (SD) | SE | F | MAF | 备注 |
|---|---|---|---|---|---|
| rs3732477 | −0.4994 | 0.0071 | 4942 | 0.230 | **有相关蛋白改变变异(PAV)** |
| rs34225900 | +0.0744 | 0.0132 | 32 | 0.118 | indel,两库写法不同 |
| rs1466733 | +0.0495 | 0.0100 | 24 | 0.232 | |

两个关键问题,分析中已处理:

1. **rs3732477 占 IVW 约 98% 的权重**(单独解释血浆 CPN2 的 8.83% 变异,另两个各约 0.1%),
   且带 PAV,存在表位效应风险 —— SomaScan 靠适配体结合定量,氨基酸改变可能使
   适配体结合变弱、读数偏低,而蛋白实际含量未变。故另跑一版剔除它的敏感性分析。
2. **rs34225900 的等位基因表示法两库不一致**:deCODE 为 `chr3:194350491 TAATG/TAATT`
   (左补齐),FinnGen 为 `chr3:194350495 T/G`。差异位即 deCODE 串第 5 位,
   因此 deCODE 的 `TAATG` 对应 FinnGen 的 `G`。脚本中显式映射,未按位置盲配。

## 数据获取

本仓库**不包含数据**。`data/` 已在 `.gitignore` 中排除,原因与获取方式如下。

### 暴露端 · deCODE

deCODE 的数据使用协议规定 **"You will not redistribute the downloaded files."**,
因此原始文件与其区域切片均不入库。请自行获取:

1. 访问 <https://download.decode.is>,阅读并接受数据使用协议
2. 搜索 `CPN` ,下载 `6415_90_CPN2_CPN2.txt.gz`(909 MB)
3. 放到 `data/raw/`,并用同目录的 md5 校验(校验的是**解压后**的文件):
   ```bash
   gunzip -c data/raw/6415_90_CPN2_CPN2.txt.gz | md5
   ```
4. 切出 cis 区域(CPN2 TSS 194,351,387 ± 1 Mb):
   ```bash
   gunzip -c data/raw/6415_90_CPN2_CPN2.txt.gz \
     | awk -F'\t' 'NR==1 || ($1=="chr3" && $2>=193351387 && $2<=195351387)' \
     > data/exposure_CPN2_cis.tsv
   ```

> ⚠ CPN1 在 SomaScan 平台命名为 `CBPN`,对应 `7142_5_CPN1_CBPN.txt.gz`。

### 结局端 · FinnGen R12

公开可取,**无需下载完整的 800 MB 文件** —— 用 tabix 配远程索引只取所需区域:

```bash
brew install htslib          # 提供 tabix

B=https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release
HDR=$'#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls'

for PH in NAFLD CIRRHOSIS_BROAD C3_HEPATOCELLU_CARC_EXALLC; do
  curl -sL -o "/tmp/$PH.tbi" "$B/finngen_R12_$PH.gz.tbi"
  { echo "$HDR"
    tabix "$B/finngen_R12_$PH.gz##idx##/tmp/$PH.tbi" 3:193351387-195351387
  } > "data/outcome_${PH}_cis.tsv"
done
```

`##idx##` 语法指定本地索引 + 远程数据,htslib 直接用 HTTP range 请求取对应字节区间。

## 运行

```bash
Rscript R/01_prepare_data.R   # 整理两端数据 + 等位基因频率交叉核对
Rscript R/02_run_mr.R         # TwoSampleMR 主分析与敏感性分析
Rscript R/03_coloc.R          # 贝叶斯共定位
```

依赖:`TwoSampleMR` `coloc` `data.table`(r-universe:`https://mrcieu.r-universe.dev`)。
区域数据取用了 `tabix`(htslib)+ FinnGen 的远程 `.tbi` 索引,无需下载完整的 800 MB 文件。

## 目录

```
data/raw/    deCODE 原始全基因组汇总统计(909 MB,已 md5 校验)
data/        两端 cis 区域切片、整理后的暴露/结局数据
R/           分析脚本
results/     MR、异质性、多效性、leave-one-out、共定位结果
```

## 已知限制

- 工具变量仅 3 个,MR-PRESSO 无法进行;剔除 PAV 后仅剩 2 个,MR-Egger 与加权中位数失效
- 剔除 PAV 后置信区间极宽(如 NAFLD:0.50–2.04),属功效不足,**不构成证伪**
- FinnGen 的肝硬化与肝癌 ICD 编码不区分病因,「非病毒非酒精」在公共数据层面无法切分
- NASH 在 FinnGen R12 仅 254 例,功效不足,已排除出分析
- 未纳入邻近的 **LRRC15**(距 CPN2 仅 3.8 kb,该位点信号强得多,是竞争解释),待补
