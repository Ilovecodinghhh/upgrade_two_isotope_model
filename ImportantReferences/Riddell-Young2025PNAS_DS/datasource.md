# 数据来源与处理流程总结

本文档总结 `Riddell-Young_2025_MassBalancePackage` 中，各排放源排放强度以及两种同位素组成时间序列是如何得到的。结论基于项目中的主反演脚本和其上游预处理脚本，而不是仅看最终 `Output/` 结果文件。

## 1. 总体流程

这个项目的数据链可以分成三层：

1. 原始/外部数据层
   - 大气观测：全球 CH4、全球 `d13C-CH4`、全球 `dD-CH4`
   - 底层排放清单或后验通量：CarbonTracker CH4、CTCH4 prior、EDGAR、GFED、EIA 等
   - 各源同位素端元数据库：煤、油气、湿地、反刍动物、植被 C3/C4 分布、水体氢同位素等

2. 项目内部预处理层
   - 先把不同来源的数据转成年度的源同位素时间序列
   - 再把这些时间序列写到 `Output/*.csv`

3. 主质量平衡反演层
   - 用全球大气 CH4 年均浓度推总源强
   - 用全球大气同位素年均值推总源同位素组成
   - 再用三端元质量平衡把总源拆成化石源、微生物源和生物质燃烧源

## 2. 主反演真正读取的输入

### 2.1 `d13C` 反演主脚本

主脚本：`Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/d13C_MassBalance_MC.py`

直接读取的数据包括：

| 类型 | 文件 | 用途 |
|---|---|---|
| 全球大气 `d13C-CH4` 年均值 | `data/ch4c13_nh_sh_mean.xlsx` | 先按年平均，得到全球 `d13C-CH4` 序列 |
| 全球大气 `d13C-CH4` Monte Carlo 迭代 | `data/d13C_dei_compiled.txt` | 传播大气观测不确定度 |
| 全球 CH4 年均浓度 | `data/GML_CH4_AnnualMean.xlsx` | 通过一箱模型反推出总源强 |
| CarbonTracker CH4 年度汇总 | `data/CarbonTracker_CH4.xlsx` | 提供对照用的微生物/化石/火排放，以及 BB 固定量 |
| 微生物源 `d13C` 年度序列 | `Output/Mic_d13C_annual.csv`、`Output/Mic_d13C_MC.csv` | 微生物端元及其 MC 扰动 |
| 生物质燃烧 `d13C` 年度序列 | `Output/BB_d13C_annual.csv` | BB 端元均值和不确定度 |
| 化石源 `d13C` 年度序列 | `Output/FF_d13C_GlobUnc.csv`、`Output/FF_d13C_GlobMC_EDGAR.csv`、`Output/FF_d13C_GlobMC_CTCH4.csv` | 化石端元均值和 MC 扰动 |

### 2.2 `dD` 反演主脚本

主脚本：`Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/dD_MassBalance_MC.py`

直接读取的数据包括：

| 类型 | 文件 | 用途 |
|---|---|---|
| 全球大气 `dD-CH4` 年均值 | `../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_dei_UmezawaCal_noBUDS.csv` | 全球 `dD-CH4` 年度序列 |
| 全球大气 `dD-CH4` Monte Carlo 迭代 | `../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx` | 传播大气观测不确定度 |
| 全球 CH4 年均浓度 | `data/GML_CH4_AnnualMean.xlsx` | 推总源强 |
| CarbonTracker CH4 年度汇总 | `data/CarbonTracker_CH4.xlsx` | 提供 BB 固定量和对照 |
| 微生物源 `dD` 年度序列 | `Output/Mic_dD_AnnGlob.csv`、`Output/Mic_dD_MC.csv` | 微生物端元及其 MC 扰动 |
| 生物质燃烧 `dD` 年度序列 | `Output/BB_dD_annual.csv` | BB 端元均值和不确定度 |
| 化石源 `dD` 年度序列 | `Output/FF_dD_GlobUnc.csv`、`Output/FF_dD_GlobMC_EDGAR.csv`、`Output/FF_dD_GlobMC_CTCH4.csv` | 化石端元均值和 MC 扰动 |

## 3. 排放强度时间序列是如何得到的

这里的“排放强度”在主脚本里实际是各大类源的年排放量时间序列，而不是先验清单直接给定的结果。

### 3.1 先由全球 CH4 年均浓度反推总源强

`d13C_MassBalance_MC.py` 和 `dD_MassBalance_MC.py` 都先读取 `data/GML_CH4_AnnualMean.xlsx`，然后用一箱模型计算：

- 年际大气 CH4 储量变化
- 给定寿命 `Lifetime = 9` 年时的总源强 `SumSource`
- 不同汇情景下的总源强，如 `SumSource_OH_inc`、`SumSource_Cl_inc`

因此：

- 总源强时间序列来自大气 CH4 浓度变化
- 不是直接从 EDGAR 或 CarbonTracker 拿来的

### 3.2 再由同位素质量平衡拆分为 FF / Mic / BB

主脚本先从全球大气同位素变化推总源同位素组成：

- `d13C_source`
- `dD_source`

再做三端元质量平衡：

- 生物质燃烧排放量 `BB` 默认取 `CarbonTracker_CH4.xlsx` 中 BB 年度序列的均值，因此在主反演中通常被当作常数
- 化石源端元随时间变化，来自 `FF_*` 预处理结果
- 微生物源端元随时间变化，来自 `Mic_*` 预处理结果

得到的就是：

- `FFS_ffvary` / `MicS_ffvary`：绝对排放量时间序列
- `FFS_ffvaryR` / `MicS_ffvaryR`：相对基准期的变化量

其中相对序列的基准期：

- `d13C` 脚本是 2005–2007 平均
- `dD` 脚本也是 2005–2007 平均

### 3.3 平滑与输出

主反演后续对每次 Monte Carlo 结果做 5 年平滑，再计算均值和标准差，写出：

- `Output/Results_d13C-MassBalance_*.csv`
- `Output/Results_Rd13C-MassBalance_*.csv`
- `Output/Results_dD-MassBalance_*.csv`
- `Output/Results_RdD-MassBalance_*.csv`

这些就是论文图中真正使用的各源排放时间序列。

## 4. 生物质燃烧源同位素时间序列

### 4.1 `d13C`：`BB_d13C.py`

脚本：`Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/BB_d13C.py`

原始数据：

- `data/C4_distribution_NUS_v2.2.nc`
  - Luo et al. 2023 的 C4 植被分布，2001–2019
- `data/Still2003_C4.xlsx`
  - Still et al. 2003 的旧版 C4 地图，用于比较
- `data/CTCH4_2023_flux3x2.nc`
  - CarbonTracker CH4 的 pyrogenic flux
- `data/GFED5_Beta/GFED5_Beta_monthly_*.nc`
  - GFED5 monthly CH4 fire emissions
- `data/prior_monthly_emission_kg_lei.nc`
  - CTCH4 prior 中的 biomass burning 通量

处理方式：

1. 把 C4 分布从更高分辨率压缩到 1x1 网格，并补齐到 1998–2021。
2. 将火源通量转成年度网格。
3. 设定 C3、C4 火源端元：
   - `C3_d13C = -26.8`
   - `C4_d13C = -12.7`
4. 用每年火排放在 C3/C4 区域中的权重，计算全球年度平均 `d13C-BB`。
5. 用 Monte Carlo 扰动 C3/C4 端元不确定度，得到年均值和标准差。

输出：

- `Output/BB_d13C_annual.csv`
- `Output/BB_d13C_1x1_YYYY.txt`

结论上，`BB d13C` 的时间变化主要来自：

- 年际 C3/C4 空间权重变化
- 所采用的火排放分布

### 4.2 `dD`

主脚本直接读取 `Output/BB_dD_annual.csv`。仓库里可见其最终产物，但当前目录下没有和 `BB_d13C.py` 对应的 `BB_dD.py` 源脚本。也就是说：

- `BB dD` 年度序列已经预先生成
- 但本仓库当前可直接回溯的处理链不如 `BB d13C` 完整

从输出文件命名看，项目至少保留了：

- `Output/BB_dD_annual.csv`
- `Output/BB_dD_1x1.txt`
- `Output/BB_dD_1x1_Unc.csv`

因此在主反演中，`BB dD` 被作为已完成的上游输入使用。

## 5. 微生物源同位素时间序列

### 5.1 微生物源排放分配的共同框架

`Mic_d13C.py` 和 `Mic_dD.py` 都用了同一类思路：

原始数据：

- `data/prior_monthly_emission_kg_lei.nc`
  - 读取 rice、ruminant、wetland、termite、waste_landfill、wild_animals 等 prior 通量
- `data/CTCH4_2023_flux3x2.nc`
  - 读取 `microbial_flux` 作为后验微生物总通量

处理方式：

1. 从 prior 中提取各微生物子类年总量和空间分布。
2. 从 CarbonTracker 后验中取微生物总通量。
3. 用 `Posterior_to_prior = Mic_ann / prior_mic` 形成“后验/先验”缩放因子。
4. 用这个缩放因子把各子类先验通量重权重化。
5. 把后验无法归到 rice/ruminant/termite/landfill/wild_animals 的剩余量并入 wetland。
6. 得到每年各微生物子类在总微生物中的比例 `Rice_Post_Frac`、`Wetland_Post_Frac` 等。

这一步很关键：  
微生物源同位素时间序列不是简单文献平均，而是“子类端元 × 后验约束后的子类占比”加权得到的。

### 5.2 `d13C`：`Mic_d13C.py`

额外原始数据：

- `data/C4_distribution_NUS_v2.2.nc`
  - 用于反刍动物和野生动物的 C3/C4 饮食分配
- `data/Still2003_C4.xlsx`
  - 旧版 C4 地图，对比用
- `data/isotem_wetland_d13C-CH4.nc4`
  - Oh et al. 2022 的湿地 `d13C-CH4` 时空场
- `data/Chang_2019_ruminants.xlsx`
  - 反刍动物 `d13C-CH4` 年度序列
- `data/Oh_2022_Wetlands.xlsx`
  - 湿地全球平均 `d13C-CH4` 年度序列

子类端元处理：

- Wetland：直接用 `Oh_2022_Wetlands.xlsx` 年度序列
- Ruminant：直接用 `Chang_2019_ruminants.xlsx` 年度序列
- Rice / Waste / Termite：使用固定端元，并叠加 Suess effect 线性趋势
- Wild animals：先按 C3/C4 饮食结构估算，再纳入总微生物混合

Monte Carlo 处理：

- 对各子类端元加高斯扰动
- 对各子类比例加 10% 扰动后归一化
- 生成 1000 组微生物 `d13C` 年度序列

输出：

- `Output/Mic_d13C_annual.csv`
- `Output/Mic_d13C_MC.csv`
- 以及湿地 1x1 地图 `Output/Wetland_d13C_1x1_YYYY.txt`

### 5.3 `dD`：`Mic_dD.py`

额外原始数据：

- `data/d2h_MA.tif`
  - 全球水体 `d2H` 地图

子类端元处理：

- Wetland、Rice、Ruminant、Termite、Wild animals 大都不是直接从单独观测时间序列读取
- 而是用全球水体 `d2H` 地图，通过经验线性关系推导源 `dD-CH4`
- 脚本中对 wetland、rice、landfill 分别设置了斜率和截距，并在 MC 中扰动
- 缺测海洋格点用纬向平均值回填，避免乘以通量时出现空值

随后：

1. 用各子类后验权重图计算每年的全球平均子类 `dD`
2. 再用各子类年比例混合成总微生物 `dD`
3. 用 MC 扰动经验关系和子类比例

输出：

- `Output/Mic_dD_AnnGlob.csv`
- `Output/Mic_dD_MC.csv`

## 6. 化石源同位素时间序列

### 6.1 共同思路

`FF_d13C_GlobMean.py` 和 `FF_dD_GlobMean.py` 的框架几乎一致：

原始数据：

- `data/coal_d13C.csv` / `data/coal_dD.csv`
  - 国家尺度煤源同位素统计
- `data/ONG_d13C.csv` / `data/ONG_dD.csv`
  - 国家尺度油气源同位素统计
- `data/EDGAR8_Coal.csv`
  - 各国煤源 CH4 排放
- `data/EDGAR8_ONG.csv`
  - 各国油气源 CH4 排放
- `data/US_ONG_trends.csv`
  - 美国油气端元时间趋势
- `data/China_Canada_ONG_Trends.csv`
  - 中国、加拿大油气端元时间趋势

处理方式：

1. 先把煤和油气的国家端元库标准化到统一国家名。
2. 用各国历史平均 EDGAR 排放，把国家端元混成国家级 fossil mean。
3. 对没有同位素数据的国家，用全球均值回填。
4. 对中国、美国、加拿大的油气源，覆盖为逐年变化的时间序列；其他国家默认常数。
5. 用 EDGAR 各国逐年煤/油气排放作权重，计算全球逐年：
   - coal mean
   - ONG mean
   - total fossil mean

### 6.2 不确定度处理

两套脚本都会：

1. 从国家端元库中读取 `std` 和 `n`
2. 计算国家级标准误 `StdErr = std / sqrt(n)`
3. 缺测国家用全球平均标准误回填
4. 进行 1000 次 Monte Carlo
   - 对煤端元按国家加扰动
   - 对油气端元按国家加扰动
   - 再按 EDGAR 年度排放权重求全球均值

输出：

- `Output/Coal_*_GlobUnc.csv`
- `Output/ONG_*_GlobUnc.csv`
- `Output/FF_*_GlobUnc.csv`
- `Output/R_FF_*_GlobUnc.csv`
- `Output/FF_*_GlobMC_EDGAR.csv`

其中 `*` 代表 `d13C` 或 `dD`。

### 6.3 CTCH4 版本的化石源时间序列

这两套脚本还额外生成一套 “CTCH4 posterior flux-weighted” 的化石源时间序列：

额外原始数据：

- `data/CTCH4_2023_flux3x2.nc`
  - fossil_flux
- `data/prior_monthly_emission_kg_lei.nc`
  - `flux_coal`、`flux_oil_gas`、`flux_geologic_seep`、`flux_other_industry`

处理方式：

1. 取 CarbonTracker posterior fossil flux 作为总化石后验通量。
2. 用 prior 中 coal / oil_gas / geologic seep / other industry 的比例构建后验/先验缩放。
3. 将未分配部分并入 ONG。
4. 用生成好的 1x1 化石同位素地图和 CTCH4 后验排放权重，计算全球年度平均化石端元。
5. 对 1x1 端元地图加上格点标准误的 MC 扰动，形成 `FF_*_GlobMC_CTCH4.csv`。

这套 CTCH4 结果在主反演里主要作为敏感性比较或替代方案，而默认主反演使用的是 EDGAR 加权的 `FF_*_GlobMC_EDGAR.csv`。

## 7. 大气同位素时间序列

### 7.1 全球 `d13C-CH4`

`d13C` 主脚本直接读取：

- `data/ch4c13_nh_sh_mean.xlsx`
- `data/d13C_dei_compiled.txt`

处理方式：

1. 先对 `ch4c13_nh_sh_mean.xlsx` 做年平均，形成全球年均 `d13C-CH4`
2. 再用 `d13C_dei_compiled.txt` 中的多组迭代传播不确定度

当前仓库里没有看到生成 `d13C_dei_compiled.txt` 的完整上游脚本，因此这一步对本项目而言属于“已预处理输入”。

### 7.2 全球 `dD-CH4`

`dD` 主脚本不是直接从观测原始文件读，而是依赖另一个子项目：

- `Riddell-Young_2025_dD_GlobMean`

关键脚本：

- `Riddell-Young_2025_dD_GlobMean/Riddell-Young_2025_dD_GlobMean/dD_globmean.py`

原始数据：

- `Riddell-Young_2025_dD_GlobMean/.../data/*.txt`
  - 多个站点的 `01D0` 观测文件
- `siteinfo_all_ch4h2.txt`
  - 台站纬度与站点信息
- 各站点的 `output/*_smoothedMC.txt`
  - 站点平滑曲线与 MC 结果

处理方式：

1. 对不同实验室尺度做统一校正（MPI / IMAU / NIPR / INSTAAR）。
2. 每次 MC 随机剔除 2 个站点，传播网络抽样不确定度。
3. 按纬带聚合，构造半球与全球平均。
4. 输出全球 `dD-CH4` 年度序列和迭代矩阵。

主反演实际使用的是该子项目已生成的：

- `GlobMean_dD_dei_UmezawaCal_noBUDS.csv`
- `GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx`

## 8. 项目里“底层清单”的真实角色

容易混淆的一点是：EDGAR、CarbonTracker、GFED、CT prior 这些数据在项目中承担的角色并不完全一样。

### 8.1 用来约束端元时间变化和空间权重

- EDGAR8：主要用于化石源国家权重
- CTCH4 posterior：主要用于微生物源和化石源的后验空间权重
- CT prior：主要用于把后验总通量拆回各子类的比例
- GFED5：主要用于火源空间/时间权重

### 8.2 不直接作为最终排放强度结论

最终发表的 FF / Mic 排放时间序列并不是这些 bottom-up 数据原样输出，而是：

- 大气 CH4 给总量
- 大气同位素给总源同位素组成
- 源端元时间序列给拆分依据
- 最后通过质量平衡反演得到

## 9. 可直接追踪的关键输入输出对应关系

### 9.1 `d13C`

| 中间结果 | 来源脚本 | 主反演读取文件 |
|---|---|---|
| 全球微生物 `d13C` | `Mic_d13C.py` | `Output/Mic_d13C_annual.csv` / `Output/Mic_d13C_MC.csv` |
| 全球生物质燃烧 `d13C` | `BB_d13C.py` | `Output/BB_d13C_annual.csv` |
| 全球化石源 `d13C` | `FF_d13C_GlobMean.py` | `Output/FF_d13C_GlobUnc.csv` / `Output/FF_d13C_GlobMC_EDGAR.csv` / `Output/FF_d13C_GlobMC_CTCH4.csv` |
| 全球大气 `d13C` | 上游已预处理 | `data/ch4c13_nh_sh_mean.xlsx` / `data/d13C_dei_compiled.txt` |

### 9.2 `dD`

| 中间结果 | 来源脚本 | 主反演读取文件 |
|---|---|---|
| 全球微生物 `dD` | `Mic_dD.py` | `Output/Mic_dD_AnnGlob.csv` / `Output/Mic_dD_MC.csv` |
| 全球生物质燃烧 `dD` | 上游已预处理 | `Output/BB_dD_annual.csv` |
| 全球化石源 `dD` | `FF_dD_GlobMean.py` | `Output/FF_dD_GlobUnc.csv` / `Output/FF_dD_GlobMC_EDGAR.csv` / `Output/FF_dD_GlobMC_CTCH4.csv` |
| 全球大气 `dD` | `Riddell-Young_2025_dD_GlobMean` 子项目 | `../Riddell-Young_2025_dD_GlobMean/output/GlobMean_dD_dei_UmezawaCal_noBUDS.csv` / `.xlsx` |

## 10. 一句话总结

这个项目不是“从某个排放清单直接读出各源排放时间序列”，而是：

- 先用底层清单和通量产品构建各源同位素端元的年度时间序列，
- 再用全球大气 CH4 和全球大气同位素观测做一箱质量平衡，
- 最后反推出化石源和微生物源的排放强度时间序列。

如果后续还需要，我可以继续把这份文档再补成“按脚本执行顺序”的版本，或者再加一张“数据流关系图”。

已删除项目中全部大于 100MB 的文件；删除后再次扫描，当前仓库里已经没有超过 100MB 的文件了。

删除的 3 个文件是：

Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/data/prior_monthly_emission_kg_lei.nc，约 2930.73 MB
Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/data/isotem_wetland_d13C-CH4.nc4，约 783.13 MB
Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/data/C4_distribution_NUS_v2.2.nc，约 237.38 MB
删除依据是：

这 3 个文件都在代码里被上游脚本处理过
且项目内已经保留了对应的派生输出，例如 Wetland_d13C_1x1_*.txt、BB_d13C_1x1_*.txt、Mic_d13C_annual.csv、Mic_dD_AnnGlob.csv 等
主反演脚本本身读取的是这些 Output/*.csv / Output/*.txt 结果，而不是这几个大原始文件
需要注意的是，如果你之后想重新运行这些上游预处理脚本，比如 Mic_d13C.py、Mic_dD.py、BB_d13C.py、FF_*_GlobMean.py，就需要把对应原始大文件重新放回 data/。