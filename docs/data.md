# 数据表

数据与引擎**分离**：`chemkit/data/*.json` 是唯一事实源，加载器在
`chemkit.data`，派生（Hess 焓、晶格阴离子电对、弱配形态池）在加载时完成。
所有数值都能追溯到文献或明确的计算约定，**没有"为了让某条用例变绿"选的值**
（纪律见文末）。

## 1. 文件与规模

| 文件 | 内容 | 条数 |
|---|---|---|
| `couples.json` | 氧化还原电对（半反应 `ox/red`、电子数 `n`、`E0`） | 163 → 加载后 **180**（含派生电对） |
| `pka.json` | 酸碱解离常数（`acid/base/pka/n`） | 89 |
| `ksp.json` | 溶度积（`solid/pair/pKsp`） | 230 |
| `beta.json` | 配合物累积稳定常数（`complex/center/ligand/nu/logb`） | 162 |
| `substance_ex.json` | 物质附加属性（`form`/`ions`/`conc_forms`/`conc_M`/`spec`） | 201 |
| `thermo.json` | 标准生成焓 ΔHf°（kJ/mol，扁平 `名称 → 值`） | 320 |
| `overrides.json` | OVERRIDE 逃生舱（`id/match/result`） | 27 |
| `tests.json` | 回归用例库 | 1176 |

元素覆盖：数据表共 **710 个物种名、79 种元素**（含稀土全族、锕系、铂族、
后过渡与类金属）；未收录元素由元素周期表级拆盐兜底电离入账，是否参与反应
仍由数据表决定。

## 2. 字段

### `couples.json`（每对一个半反应）

| 字段 | 含义 |
|---|---|
| `ox` / `red` | 氧化态 / 还原态物种（`n` 为转移电子数） |
| `E0` | 标准电极电势（V，298.15 K） |
| `oh_ref` | 碱性侧参考写法备注（15 条） |
| `kinetics` | 动力学闸门（见 §3；39 条） |
| `gate` | 浓度/形态闸门（见 §4；10 条） |
| `rate` | 速率/时间尺度维度（见 §5；1 条） |
| `consistency` | 与相邻电对的自洽性核算过程（14 条，供复核） |
| `note` | 出处与选值理由（134 条） |

### `pka.json` / `ksp.json` / `beta.json`

| 字段 | 含义 |
|---|---|
| `acid` / `base` / `pka` / `n` | 共轭酸碱对、pKa、质子数（多元酸 `n>1` 时 pKa 为该级**累计**值的 1/n 折算源） |
| `solid` / `pair` / `pKsp` | 固相、阳阴离子对、溶度积指数 |
| `phase` | **相声明**：该 Ksp 对应哪个相（15 条，见 §6） |
| `slight` / `film` / `unlock_T` | 微溶标记 / 钝化膜 / 膜解锁温度（6 条 `film`） |
| `complex` / `center` / `ligand` / `nu` / `logb` | 配离子、中心、配体、配体数、**累积** logβ |
| `calibrated` | 与 Ksp 联动校准的说明（3 条；如 `[Al(OH)₄]⁻` 与 pKsp(Al(OH)₃)=33 配对） |
| `note` / `consistency` | 出处与复核过程 |

### `substance_ex.json` / `thermo.json` / `overrides.json`

| 字段 | 含义 |
|---|---|
| `form: molecule` | 该物质在账本里是分子态（弱酸/弱碱，如 HF、HCN） |
| `ions` | 该物质电离出的离子（34 条） |
| `conc_forms` / `conc_M` | 浓酸的分子态形态与浓度标定点（4 条） |
| `spec` | 特殊形态说明 |
| ΔHf° | `thermo.json` 扁平表；**裸名 = 账本态**（水溶/凝聚），`"X(g)"` 键 = 逸出气相态 |
| `override.id/match/result` | 逃生舱：匹配条件（`species`/`T_min`）与返回结果（`degree`/`annotations`/`note`） |

## 3. 动力学字段（`couples.json` 的 `kinetics`，唯一书写形式）

`chemkit/data.py::_expand_kinetics` 是**唯一接线点**：只转换下列键，
未列出的键不会生效（只进 `kinetics_all` 留档，由 `tools/data_audit.py`
审计"注了却没被消费"）。

| 字段 | 含义 |
|---|---|
| `ox_closed` / `red_closed` | 该电对作氧化剂 / 还原剂方向封闭 |
| `below_T` | 低于该温度（K）双向封闭 |
| `below_T_only_red` | 温度闸门只对指定还原剂生效 |
| `closed_with_red` / `closed_with_ox` | 只对指定搭档封闭 |
| `closed_except_red` | 除指定还原剂外封闭 |
| `red_pH_min` / `ox_pH_max` | 还原/氧化通道的 pH 形态闸门（方向感知） |
| `h2_passivation` | 金属-水析氢的钝化窗 `[lo, hi]`（pH 落窗内则封闭析氢） |
| `h2o_red_oh_min` | 氧化水的浓碱解锁阈值（如浓碱制锰酸钾） |
| `rev_gate` | 逆向驱动闸门（产物选择性） |

语义：这些是"**无限时间也不发生**"的化学硬事实（钝化、过电位、催化依赖），
与热力学严格分离；`kinetics=False` 可回到纯热力学基线做对照。

## 4. 浓度/形态闸门（`gate`）

| 字段 | 含义 |
|---|---|
| `pH_max` / `T_min` | 形态窗口（超出即该电对不驱动） |
| `c_min` / `c_max` | 氧化剂浓度窗口 |
| `acid` + `c_basis: shadow` | 以**游离酸影子库存**（不含盐的共同离子）为基准 |
| `ligand_max` | 配位掩蔽：指定配体超过阈值时该电对不再作为氧化剂 |
| `vs_red_block_E_max` | 只对电势低于阈值的还原剂开放 |

## 5. 速率 / 时间尺度（`rate`，v0.5.2 新增维度）

```json
"rate": {"Eu^{2+}": {"tier": "slow", "pH_min": 3.0}}
```

`tier` 取 `slow`（η=0.1）或 `very_slow`（η=1e-3）；走步时对**该搭档**的
反应程度按 η 限幅（不动点不变，只是到不了）。用途：把"热力学允许但动力学
慢到看不见"的通道如实呈现（如近中性介质里 Eu²⁺ 不被水氧化），而不用
`closed_*` 一刀切封死。

## 6. 选值纪律

1. **`phase` 只做声明，不改数值**。文献 Ksp 的测定相与账本写法可能不一致
   （无水盐 vs 水合物）。换相是化学决策（温度、过饱和度、陈化时间），
   必须逐条定并写进 `note`，**不许为了让某条用例变绿而选值**。
2. **宁可如实记录反例，不许删真数据躲测试**。已知与实验不符或无条件
   活度模型外推失真的汇编值（如 `[AgI₃]²⁻`）在 `agents/log.md` 留档拒录理由，
   而不是悄悄不要。
3. **宁缺毋假**：没有 ΔHf 数据的条目就不做 van't Hoff 温度修正
   （`dH` 由 `thermo.json` 经 Hess 派生：电对 100/180、pKa 53/89、
   Ksp 91/230、logβ 9/160 有条目覆盖），而不是编一个斜率。
4. **`calibrated` 必须写清与谁联动**（如 `[Al(OH)₄]⁻` 的 logβ 与
   pKsp(Al(OH)₃)=33 联合标定，使 `Al(OH)₃ + OH⁻ → [Al(OH)₄]⁻` 的整体
   logK ≈ +1.5 符合两性溶解的实测）。
5. **准入即审计**：新增/改动的数据要过
   `tools/data_audit.py`（机械审计：死注解、值域、重复、单侧缺失）、
   `tools/hess_audit.py`（派生候选的 logK 必须等于基候选的线性组合，
   含变温 273.15/363.15 K）、`tools/eqcheck.py`（用例期望方程式的精确守恒）。

## 跨表完备性审计（`tools/db_matrix.py`）

`tools/data_audit.py` 检查**表内**自洽（守恒、电子数、重复、冲突），
`tools/db_matrix.py` 检查**跨表**完备：某物种在一张表里有、在另一张表里没有，
不报错也不违反守恒，只是化学少了一块。当前缺口家族（按测试用例触达数排序）：

- **有 pKa 无 thermo**（45 条）：硅酸/锗酸、草酸、过氧族 —— 缺温度依赖。
- **电对成员缺 thermo**（82 条）：CO₂/H₂C₂O₄、碱土金属、Al/Al(OH)₃ 等 —— 无 dH。
- **Ksp 阳离子缺一级水解 β**（55 条）：Fe²⁺、Cu²⁺、Zn²⁺、Mn²⁺、Mg²⁺ 等仍在
  Ksp 派生的复合近似上（Fe³⁺/Al³⁺ 已由显式 β 接管）。
- **化合物固相缺溶解度数据**（6 条）：MnO₂、Fe₂O₃、PbO₂、MnOOH。
- **化合物固相缺 thermo**（147 条）：Cr(OH)₃、SiO₂、Ga(OH)₃、镧系 M(OH)₃ 成片缺失。

```bash
python tools/db_matrix.py                 # 全部规则 + 规模摘要
python tools/db_matrix.py --rule solid_no_ksp --top 20
python tools/db_matrix.py --json .tmp_db.json
```
