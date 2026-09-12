# chemkit 0.4.4

水溶液反应判定与产物计算引擎（Python，**零第三方依赖**——纯标准库实现，
含精确有理数零空间配平）。
输入投料与条件，输出：是否发生变化 / 是否发生化学反应 / 反应程度 / 初态与终态组成 /
pH / 总净离子方程式 / 分步离子方程式 / 逸出气体 / 绝热耦合热效应。

设计哲学：**从热力学数据推出反应行为**（电对 E0、pKa、Ksp、logβ、Henry），而不是
维护一张庞大的反应事实表；热力学解释不了的硬事实才以动力学闸门标注。
完整公理、模块职责与不变量见下文「模块文档」，未完成项与失败实验记录见
`architecture.md`。

## 安装

Python ≥ 3.12。

```bash
git clone https://github.com/arlenWKX/chemkit.git
cd chemkit
python -m chemkit.testsuit   # 1173 例 / 1294 断言 + 酸碱锚点 + 环闭合（~40 s）
```

> **基线如实声明**：当前 HEAD 实测 **1290/1294** 断言通过（4 例长期红例
> 及其病根见下文「测试」与「已知局限」）。任何改动都必须跑全量复核——
> "全绿"也不等于"没变"，判定语义的位移要靠 `chemkit.converg` 的全量
> 差分与 `tools/xver.py` 的逐例对照才能看见。

## 快速上手（5 分钟）

### 统一入口：`Engine` 对象（v0.3.6 起）

进程内建一次、长期复用——judge / react / system 三个层级是它的**平级方法**，
没有 react()→System()→judge() 的层层包装：

```python
import chemkit

eng = chemkit.Engine()                              # 默认包内数据表（持表一次、缓存随行）
r = eng.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)   # 一步式（直接落求解管线，零中间对象）
sys = eng.system(V=1.0, T_C=25)                     # 连续投料体系（共享引擎的全部缓存）
r = sys.add("NaOH", 0.1)
raw = eng.judge([{"name": "HCl", "mol": 0.1}],      # 引擎 raw dict 直通（底层）
                {"V_L": 1.0})
```

函数式便捷入口 `chemkit.react / chemkit.System` 与 Engine 方法完全同语义
（隐式使用全局默认表）；需要换表或隔离缓存时用 Engine 对象。

### 一步式反应：`react()`

```python
r = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)

r.changed        # True  —— 体系发生显著净变化（宽口径：含溶解/电离/水解形态变化）
r.reacted        # True  —— 发生狭义化学反应（氧化还原/中和/跨投料沉淀配位）
r.degree         # 2     —— 2=完全反应 / 1=可逆（部分）反应 / 0=难反应或未反应
r.net_equation   # '2H^+ + Zn -> H_2 + Zn^{2+}'（总净离子反应方程式，两侧内部顺序无关）
r.equations      # 分步离子方程式列表（按贡献降序）
r.consumption    # {'Zn': 1.0, ...} 净消耗；r.production 净生成 {化学式: mol}
r.initial        # 初态组成（强电解质已电离、SO3 等已与水反应、酸碱已中和）
r.final          # 终态组成 {化学式: mol}（H2O 为溶剂不计入）
r.pH             # 终态 pH
```

参数：

| 参数 | 含义 | 默认 |
|---|---|---|
| `substances` | 投料 `{化学式: mol}`（必填） | — |
| `V` | 溶液体积（L） | 1.0 |
| `T` / `T_C` | 温度（K）/（°C），液态水域 273.15–373.15 K | 298.15 |
| `p` | 外界气压（kPa），影响气体逸出阈值 | 101.3 |
| `tables` | 自定义数据表（默认用包内 `TABLES` 单例） | None |
| `isothermal` | 恒温假设（默认 True）；False = 绝热耦合 | True |
| `kinetics` | 动力学层（默认 True）；False = 纯热力学基线 | True |
| `gas_escape` | 自产气体逸出（默认 True）；False = 闭口体系 | True |

### 连续投料：`System`

```python
sys = chemkit.System(V=1.0, T_C=25)
sys.add("NaOH", 0.1)            # 纯水 + NaOH
r = sys.add("HCl", 0.15)        # 再投 HCl —— 按累计投料整体重新平衡
r.pH, r.net_equation
sys.feeds                       # {'NaOH': 0.1, 'HCl': 0.15}
sys.history                     # 历次 Reaction 列表
sys.result                      # 最近一次 Reaction
```

### 底层一步：`judge()`

```python
from chemkit import load_tables, judge
T = load_tables()               # 加载包内 data/ 数据表（进程内单例）
r = judge([{"name": "HCl", "mol": 1.0}, {"name": "NaOH", "mol": 1.0}],
          {"V_L": 1.0, "T_K": 298.15, "p_kpa": 101.3}, T)
# 返回 raw dict；chemkit.Reaction(r) 包装为人类可读对象
# Engine 对象上：eng.judge(subs, cond)（自动带引擎的表）
```

`judge` 的条件键：`V_L` / `T_K` / `T_C` / `c_H`（初始强酸浓度 mol/L）/ `c_OH` / `pH` /
`p_kpa` / `isothermal` / `kinetics` / `gas_escape`；未知键一律报错不静默。温度超出
液态水温度域、气压非正均报错。

## Reaction 结果对象

数据分两层：**人类可读层**（化学习惯，含 H⁺/OH⁻/H₂O，系数为最简整数比）与
**引擎记账层**（`*_raw`：H₂O 为溶剂不入账，H⁺/OH⁻ 合记为带符号质子账本
`H_excess_raw`，正 = 残余游离强酸 mol，负 = 残余游离强碱 mol）。

### 判定层

| 属性 | 类型 | 含义 |
|---|---|---|
| `changed` | `bool` | 体系是否发生显著净变化（净账本判定；含纯溶解、弱酸弱碱电离、水解等单投料形态变化——宽口径） |
| `reacted` | `bool` | 是否发生狭义化学反应：氧化还原 / 酸碱中和 / 跨投料的沉淀与配位等；纯溶解、电离、水解为 False |
| `degree` | `int` | 2 = 完全反应；1 = 可逆（部分）反应；0 = 难反应或未反应 |

### 组成层

| 属性 | 类型 | 含义 |
|---|---|---|
| `consumption` | `dict[str, float]` | 净消耗 {化学式: mol} |
| `production` | `dict[str, float]` | 净生成 {化学式: mol} |
| `initial` | `dict[str, float]` | 初态组成（post-normalize：强电解质已完成电离、SO₃ 等气体已完成水合、酸碱中和已记账） |
| `final` | `dict[str, float]` | 终态组成（H₂O 为溶剂不计入） |
| `pH` | `float \| None` | 终态 pH（OVERRIDE 路径为 None） |
| `escaped` | `dict[str, float]` | 逸出气相 {化学式: mol}（自产气体超过 H(T)·p_ext 溶解上限即逸出；投料气体视为持续供给不逸出） |

### 方程式层

| 属性 | 类型 | 含义 |
|---|---|---|
| `net_equation` | `str \| None` | 总净离子反应方程式（全部显著步骤的净和：中间体自然抵消，H₂O 显式配平，OH⁻ 从 H⁺ 正则形还原，系数最简整数比）。无显著反应为 None |
| `equations` | `list[str]` | 分步离子方程式（按贡献降序）。许多反应用多步概括更贴近书写习惯（如 Ca(OH)₂+CO₂ 是 `CO_2 + 2OH^- -> CO_3^{2-} + H_2O` 与 `Ca^{2+} + CO_3^{2-} -> CaCO_3` 两步） |
| `steps` | `list[dict]` | 引擎逐步过程（kind / equation / logK / S / extent / conversion） |

### 标注层

| 属性 | 含义 |
|---|---|
| `annotations` | `["slow"]`（该体系存在显著但动力学缓慢的通道）、`["blocked"]`（金属表面膜封锁）等 |
| `override` | 命中的 OVERRIDE 通道 id（如 `"ko2_water"`），未命中为 None |
| `raw` | `judge()` 原始 dict（备用；含 `consumption_raw` / `production_raw` / `final_raw` / `initial_raw` / `H_excess_raw` 等记账层字段） |

### 热效应层（isothermal=False 时）

| 属性 | 类型 | 含义 |
|---|---|---|
| `heat_kJ` | `float \| None` | 放热为正（= −ΔH；None = ΔHf 数据不足，宁缺毋假） |
| `dT_K` | `float \| None` | 溶液温升（水比热容、V×1000 g 计；负 = 降温） |
| `T_final_K` | `float` | 终温（绝热耦合收敛值；相变时钳在 273.15/373.15 K） |
| `thermal` | `dict` | 完整分析：`water_mol` / `missing` / `flags` / `reason`，以及耦合模式特有的 `trace`（逐轮温度/热量记录）与 `converged`（外层不动点收敛标志） |

默认 `isothermal=True`：单遍求解不计热效应（判定引擎主用途，
`thermal={"mode": "isothermal"}`）。`isothermal=False` 为**绝热耦合**：

```python
r = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0, isothermal=False)
print(r.heat_kJ, r.dT_K, r.T_final_K)
# 5.584 1.335 299.48 —— 中和热 55.84 kJ/mol，温升与自洽终温
r.thermal["trace"]   # 逐轮 [{T_K, heat_kJ, dT_K, clamped}]（温度不动点轨迹）
```

温度反馈真实进入求解：独立温度模块在化学平衡与能量平衡之间迭代温度
不动点（锚定初温 T₀：T_f = T₀ + q(T_f)/(m·c_p)，放热反应 q 随 T 单调
减弱，收缩映射 2–4 轮收敛），化学结果在自洽终温下给出（K(T) 经
van't Hoff 随温度移动）。这是状态函数意义上的严格绝热终态，而非
"初温解完再补一个 ΔT 数字"的事后报告。

`bool(reaction)` 等价于 `changed`（有无净变化）。

## 化学式写法

类 LaTeX 语法，下标 `_n`、上标电荷 `^{n±}`、括号分组、水合点：

- `H_2SO_4`、`Ca(OH)_2`、`Ca_3(PO_4)_2`、`Cu_2(OH)_2CO_3`
- 离子：`Fe^{3+}`、`SO_4^{2-}`、`Cl^-`（一价电荷可写 `^+` / `^-`）
- 配合物：`[Fe(CN)_6]^{4-}`、`[Cu(NH_3)_4]^{2+}`、`K_3[Fe(CN)_6]`
- 水合：`CuSO_4·5H_2O`（按元素总量解析）

**离子无需注册即可投料**：任何元素合法的带电物种直接入账（如 `Eu^{2+}`、`TiO^{2+}`）。
未注册的离子为旁观离子（无反应数据即不参与反应）。

**元素周期表级拆盐兜底**：未收录进数据表的盐（如 `NdCl_3`、`Ra(NO_3)_2`、
`Tb(NO_3)_4`）由公式结构 + 电荷平衡 + 常见氧化态表自动电离（`NdCl_3 → Nd³⁺ + 3Cl⁻`），
保证电荷与质量守恒；离子是否参与反应仍由数据表决定。

## 模块文档

| 模块 | 职责 | 关键对象 |
|---|---|---|
| `chemkit.interfaces` | **用户侧唯一门面**（v0.3.3 起）：Engine 对象与三层 API（react/System/judge）与辅助对象的总入口，内部实现细节不从此导出 | `Engine`、全部公开对象 |
| `chemkit.core` | 化学式解析（元素/电荷/水合）、方程式精确配平（**手写免分数整数消元零空间**，v0.3.5 起取代 sympy）、热力学温度函数 | `elements_of`、`charge_of`、`balance`、`FormulaError` |
| `chemkit.data` | 数据表加载与派生：统一 kinetics 字段展开、晶格阴离子电对派生、van't Hoff 反应焓 Hess 派生、Henry 温度修正 | `load_tables()`、`Tables` |
| `chemkit.candidates` | 引擎地基（v0.3.4 拆出）：全局常数唯一事实源、Cand 候选对象与 logK(T) 三通道、半反应配平、元素周期表级拆盐兜底、Hess 派生候选 | `Cand`、`logK_T()`、`build_derived()` |
| `chemkit.normalize` | 投料规范化（v0.3.4 拆出）：投料 → 引擎账本（电离/中和/浓酸分子态曲线/电离映射） | `normalize()` |
| `chemkit.acidbase` | **精确酸碱求解地基（0.4.3 新增，未接入主路径）**：质子化族（pKa 连通块）+ 逐级 pKa 严格多级分布 + 电荷平衡闭式解 pH；作诊断/校验用，替代启发式分支机的施工入口 | `build_families()`、`charge_pH()`、`dist_charge()` |
| `chemkit.speciation` | pH 机器（v0.3.4 拆出）：缓冲滴定记账、四分支 pH 估计、多级形态分布、分子态强酸再电离 | `estimate_state()`、`estimate_pH()` |
| `chemkit.templates` | 候选宇宙（v0.3.4 拆出）：红氧模板四级缓存（温度无关静态层/温度过滤层/静态候选）与三层枚举 memo + 静态段 present 二级缓存（v0.3.6） | `enumerate_candidates()` |
| `chemkit.engine` | 判定引擎主循环：S = logK−logQ 评估、二分平衡程度求解（微步快通道 v0.3.6）、驱动力排序 walk（震荡冻结/膜封锁/气体扫气/爬行外推加速/微步排空 v0.3.6） | `judge()` |
| `chemkit.equations` | 方程式装配（v0.3.4 拆出）：步骤聚合、瞬态中间体折叠、H/O 配平、系数有理化、净离子方程式组 | `_build_equations()` |
| `chemkit.system` | 高层实现（v0.3.6 重构）：**Engine 对象（judge/react/system 平级方法）**、Reaction 结果包装、System 连续投料、共用求解管线 `_solve()`（三层 API 直落此处，无包装链）、三模式开关 | `Engine`、`Reaction`、`System`、`react()` |
| `chemkit.thermo` | 独立温度模块：反应热事后分析（analyze）与绝热耦合求解（coupled，平衡⇌能量外层温度不动点）；ΔHf 单一数据源 thermo.json + 逸出气相拆分；与平衡求解路径完全解耦 | `analyze()`、`coupled()` |
| `chemkit.testsuit` | 测试套件（独立模块，默认不加载）：1218 化学用例 + 温度域校验 + 高层 API 自检（含 Engine 三层级/绝热耦合/热覆盖）+ 热力学环闭合检查 | `python -m chemkit.testsuit` |

单向依赖链：core ← data ← candidates ← normalize ← speciation/templates
← engine ← system ← interfaces（thermo 延迟引 engine）。

### 数据表（`chemkit/data/`）

| 文件 | 内容 | 规模 |
|---|---|---|
| `couples.json` | 氧化还原电对（ox/red/E0/n/dH/kinetics 闸门；含配合物形态电对如 [Cu(NH₃)₂]⁺/Cu） | 164 条 |
| `pka.json` | 酸碱解离常数（强酸 pKa ≤ 0 在规范化即拆解） | 89 条 |
| `ksp.json` | 溶度积（pair 阳阴离子对；`slight: true` 微溶盐；含 AgSCN/SrCrO₄/KClO₄ 等 v0.3.6 新增） | 230 条 |
| `beta.json` | 配合物累积稳定常数 logβ（v0.3.6 新增 Cu(I) 氨/氰/碘/硫氰、Co(SCN)₄、Cd(OH)₄；v0.3.7 恢复 [PbCl₃]⁻/[ZnI₄]²⁻ 并以事件口径解除弱配合判定的幻影约束） | 159 条 |
| `substance_ex.json` | 物质附加属性（solid/gas/ions 形态、浓酸 conc_forms） | 201 条 |
| `thermo.json` | 物质 ΔHf°（单一事实源：Hess 派生 + 热分析；裸名 = 账本态，"X(g)" 键 = 逸出气相态；v0.3.6 补齐 9 种标准态元素） | 319 条 |
| `overrides.json` | OVERRIDE 逃生舱（原理覆盖不到的显式注册反应） | 27 条 |
| `tests.json` | 回归测试库 | 1218 例 |

元素覆盖：**78 种元素**有反应数据（H、碱金属/碱土全族、稀土全 17 种、锕系 9 种、
过渡金属全族、铂族全族、后过渡、类金属、卤素/氧族/氮族/碳族常见非金属）；
未覆盖元素（稀有气体、At/Fr/Pa/Po 等超稀有元素）经拆盐兜底仍可正确电离入账。

### 动力学字段（`couples.json` 唯一书写形式）

| 字段 | 含义 |
|---|---|
| `ox_closed` / `red_closed` | 该电对作氧化剂/还原剂方向封闭 |
| `below_T` / `below_T_only_red` | 低于该温度（K）冻结 / 仅对指定还原剂冻结 |
| `closed_with_red` / `closed_with_ox` / `closed_except_red` | 只对指定搭档封闭 / 除指定还原剂外封闭 |
| `red_pH_min` / `ox_pH_max` | pH 形态闸门（方向感知） |
| `h2o_red_oh_min` | 氧化水的浓碱解锁阈值（如浓碱制锰酸钾） |
| `rev_gate` | 逆向驱动闸门 |

`judge()` 的 conditions 传 `"kinetics": False`（或 react/System 的同名参数）可回到纯热力学基线（无限时间）做对照。

## 测试

```bash
python -m chemkit.testsuit              # 全部 1247 条断言通道
python -m chemkit.testsuit my.json      # 自定义用例库
```

用例断言维度：`changed` / `reacted` / `degree`（int 或 str）/ `has` / `has_not` /
`has_range` / `has_any` / `has_initial` / `has_in_initial_not` / `ph` / `ann` /
`override` / `eq`（净方程）/ `eq_has`（分步方程）。

## 注意点

- **热力学为体、动力学为例外**：反应性由平衡推出；热力学无法解释的硬事实
  （钝化、光/催化依赖、动力学冻结）以 kinetics 闸门标注，与热力学显式分离。
- **气体**：外界恒压惰性环境；自产气体超过 H(T)·p_ext 溶解上限即逸出（计入
  `production` 与 `escaped`，`final` 只剩溶解态）；投料气体按持续供给处理。
- **纯溶解/电离**：`changed=True` 但 `reacted=False`（形态变化不算化学反应）。
- **OVERRIDE**：极少数原理覆盖不到的反应（白磷热碱歧化、钝化保护等 27 例）
  显式注册，命中时 `r.override` 给出通道 id。

## 性能

引擎面向长期运行设计：进程内建一次 Engine（或首次 judge）触发静态层构建
（~3 s，电对配对/派生候选，之后全部驻内存跨判定复用）；同温度判定命中
热路径 ~8 ms。首个判定的构建成本一次性摊销，稳态单判定吞吐是优化主目标。

1218 条化学用例的实测分布（每例 3 遍取均值，预热后；0.4.2；共享沙箱 ±10% 抖动）：

| 分界线 | 用例数 | 占比 |
|--------|--------|------|
| > 10 ms | 397 | 32.6% |
| > 20 ms | 195 | 16.0% |
| > 50 ms | 76 | 6.2% |
| > 100 ms | 21 | 1.7% |
| > 500 ms | 0 | 0.0% |

- 均值 14.4 ms；中位数 6.9 ms；P90 31.2 ms；最值 0.27 s。
- 全套 1247 条断言约 19 s；跨 `PYTHONHASHSEED` 结果 digest 恒定
（v0.3.9 修复两处 set 迭代序非确定性：残留字典插入序影响金属守衡
浮点累加、`T.gases` set 迭代序影响逸出/产物呈现序——v0.3.8 曾在
N08/H25/UO02 三例跨进程掷骰子）。
- v0.3.9 性能：派生候选过滤的元素位预筛（pset miss 时 ~13k 候选
  逐个子集判定 → 逐侧元素掩码 AND，典型体系幸存者 ~1%）；
  均值 19.7→13.0 ms（−34%）、>50ms 92→66 例、>10ms 682→357 例。
- v0.3.10 性能：`_buffer_titration` 的 solve_extent 级堆条目缓存
  （二分探针间仅本平衡 changing 物种量变，其余物种堆条目整 tuple
  复用；平局 cnt 恒为首次构建序保证 bit 级等价，D32 实测命中率 91%）：
  单步走步热点 H46 184→104 ms（−43%）、X10 190→105 ms（−45%）、
  D32 210→174 ms（−17%）；全量均值 13.0→12.9 ms、>50ms 66→65。
- v0.4.0 成本口径：均值 12.9→14.2 ms（+10%）来自六新物种（Zn-Cl
  β1–β4 / Ni-Cl β1–β2）的枚举空间与平衡网络扩张（候选族 112→
  118），非引擎退化；深水区反向受益：O02 376→80 ms（pH 悬崖冻结
  兼带 −79%）、RX13 227→127 ms（−44%）。
- v0.4.1 成本口径：呈现层自洽闸的一次性 _buffer_titration 调用
  （finalize + 收敛探针各一次），均值 14.2 ms 持平（沙箱抖动内）；
  判定层零位移（31 例 pH 呈现诚实化，全部 deg/reacted 不变）。
- v0.4.2 成本口径：均值 14.2→14.4 ms（+1.4%，张力族联立失败尝试
  的 Newton 成本，黑名单限 2 次即止）换取鲁棒性实质提升——resid
  债净化（V25 1.20→0.0、N34 0.084→0.009、Ag32 0.023→0.009、
  Hg2_1 0.018→0.0）与慢例族大赢（N34 162.7→70.2 ms −57%、
  Ag32 44.1→32.7 −26%、Hg2_1 333.9→287.0 −14%）；5 例语义
  位移全部定向向好、判定层零翻案。
- 恒温（默认）单遍求解；`isothermal=False` 绝热耦合另跑外层温度不动点
  （近零热效应/数据不足时单遍直出，与恒温耗时相同）。
- 深水区（多平衡耦合慢收敛）v0.3.8 已由联立 Newton 求解器收编：
  AgNO₃+氨水 1:2（628→31 ms，走步 1786→101）、Cu+AgNO₃（291→10 ms）、
  AgNO₃+氨水 1:1.5（240→39 ms）、Hg₂²⁺ 体系（533→187 ms）、
  最值 541→348 ms、>500 ms 归零。数据张力型爬行（联立不动点在物理域外）
  以边界冻结提前仲裁（物种级周转冻结判据），仍 >50 ms 的残余为非 298K
  模板首建摊销与真实迭代开销。

0.3.8 的结构性变更（收敛质量基准制度化 + 联立求解 + 净方程美化）：

- **联立 Newton 求解器**（`chemkit/joint.py`）：耦合平衡集作为联立非
  线性方程组（变量=各平衡净程度 x_j，状态线性外推保元素/电荷守恒恒等，
  方程 F_j = logK_j − logQ_j = 0），阻尼 Newton + 数值 Jacobian + 回溯
  行搜索 + 非负投影，奇异走 Tikhonov 最小二乘。三态返回：内点不动点
  跳步 / boundary（不动点在物理域外——数据张力型爬行）→ 提前冻结循环键
  / fail 严格无操作。触发：爬行检测（it≥64 且近 24 步全微步且 ≥2 个
  循环键），只联立解实际在循环的平衡（Gauss-Seidel 循环坐标的 Newton
  加速）——全收会因候选集跨数据源 Hess 互斥无解。
- **收敛质量基准**（`chemkit/converg.py`）：1218 例残差画像
  dump/diff/top（live 残差=非冻结两侧平衡的 max|S|）——动 walk 语义前
  的差分基线制度（触发点② idle 精修即被此差分否决：Co32/Co33/T34/Fe33
  四例语义翻案，idle 点走步仲裁已完成、子集不动点≠仲裁点）。
- **净方程大系数美化**（equations.py）：混合通道净差（浓度依赖物种
  分配使真实比为无理数）→ 移除次要物种 + 元素/电荷平衡子空间投影
  （约束最小二乘）+ 小分母有理重构，配平硬保证；仅 ≤3 物种的教科书
  单通道形式替换（N10 296:133:105:28 → 2Br⁻+Pb²⁺→PbBr₂），≥4 物种
  混合物保留诚实定量比。
- 步执行记账重构为 `_exec` 闭包（walk 步与联立跳步共用，bit 级等价
  剪切）；judge 增 `_probe` 只读探针（退出画像：残差/冻结/迭代数）。
- **Web 工作台**（仓库外，沙箱前端）：mini-service（bun + 持久 Python
  桥，行协议）+ 响应式 UI——投料/条件/预设/净方程/物种账本/走步过程。

0.3.7 的结构性变更（全量差分 sha256 比对 bit 级一致，1247/1247）：

- **判定语义：事件量口径**（`changed` 三通道）：物种账本差 ÷ 计量系数
  ν̄ 为溶液相显著度（ν=4 弱配位 0.86% 转化不再被 4× 账本差误判显著）、
  固相产物量与新相/逸出/中和单列——数据与判定语义正交（详见
  architecture.md §7）；代码结构：`_apply_override`/`_finalize_result`
  从 654 行 judge() 单体中提出（判定语义有了唯一居所）。
- 数据：恢复 v0.3.6 为躲 NR 幻影撤回的 [PbCl₃]⁻（Luo 2007 实测）与
  [ZnI₄]²⁻；[AgI₃]²⁻/[ZnCl₄]²⁻/[NiCl₄]²⁻ 以反例如实记录拒录理由
  （条件常数外推失真 / β4-only 逐级失真）。

0.3.6 的性能手段（全部经 1205 例全量差分验证 bit 级一致，sha256 快照比对）：

- **微步二分早停**（solve_extent 微步快通道）：多平衡体系 walk 中六成以上
  pick 求解是微步（根 < 1e-4·x_max，结果只进废弃分支）——二分区间上界
  坍缩到微步阈值即止损（探针 ~47 → ~14），三重守卫保证与跑满迭代在
  分支决策上等价（非单调 f 的 lo < hi ≤ 阈值推理不依赖 f 形态）；
- **微步排空（drain）**：pick 被判微步禁用后，若账本冻结（无强酸分子影子
  → respeciate 结构性空转），本轮评估结果仍完全有效——直接剔除已禁用者
  重挑次优，免去每个空转迭代一整轮 pH/形态/枚举/评估重算；
- **枚举静态段 present 二级缓存**：静态+派生候选的过滤只依赖 present
  物种集（与 pH/浓度闸门桶无关），浓度演化使层0 memo 反复 miss 时此段
  仍命中（数千次子集扫描降为查表）；
- **speciation 嵌套闭包外提**：estimate_state/_buffer_titration 每调用
  重建 4 个函数对象的纯开销消除。

0.3.5 的性能手段（配平数值路径与 sympy bit 级一致，差分验证 33401 次）：
- **手写免分数整数消元零空间**取代 sympy `Matrix.nullspace()`；
- **元素字典只读共享**（构建期 260 万次 `elements_of` 调用）；
- **守恒断言按物种遍历**其元素表。

0.3.4 的性能手段（数值路径与 0.3.3 bit 级一致）：
- **慢标注采样节奏记忆** / **爬行收敛窗口几何外推** / **热路径查表合并**。

## 模块文档

### 设计公理（六条，全部代码都以它们为准）

引擎不是"反应事实表 + 规则库"，而是**从热力学数据推出的平衡求解器**。
六条公理决定了模块划分与每个模块的不变量：

1. **候选 = 平衡**。凡数据表能表达的物种组合，都由 E0 / pKa / Ksp / logβ
   经 Hess 定律派生出 logK，成为候选反应；没有"必须手写"的反应。
   `overrides.json`（27 条）只留原理覆盖不到的硬事实（白磷热碱歧化等）。
2. **连续 pH，单一带符号质子账本**。游离质子不设物种，只记 `H_excess`
   （正 = 残余游离强酸 mol，负 = 残余游离强碱 mol）；方程一律写成
   H⁺ 正则形（不写 OH⁻），OH⁻ 只在人类可读层还原。
3. **驱动力 S = logK − logQ**，平衡求解 = 让 S 归零（`solve_extent` 二分）。
4. **极限记账**：账本只记"还剩下多少"，中间体的净和自然抵消。
5. **守恒契约**：任何执行步都必须保持**非 H/O 元素 + 电荷**守恒
   （H/O 交给溶剂水）；这是 `equations.py` 呈现层单闸门与
   `engine._swap_conserves` 的共同基准。
6. **热力学为体、动力学为例外**：热力学解释不了的硬事实（钝化、光/催化
   依赖、动力学冻结）以 `couples.json` 的 kinetics 字段显式标注，与热力学
   分离；`kinetics: False` 可回到纯热力学基线做对照。

### 依赖链与模块

```
core ← data ← candidates ← normalize ← speciation / acidbase / templates
     ← engine ← system ← interfaces          （joint / converg 为 engine 的协作者）
```

| 模块 | 行数 | 职责 | 关键对象 / 不变量 |
|---|---|---|---|
| `chemkit.interfaces` | 60 | **用户侧唯一门面**，内部实现不从此导出 | `Engine` 与全部公开对象 |
| `chemkit.core` | 336 | 化学式解析、**手写免分数整数消元零空间**配平、温度函数 | `elements_of` / `charge_of` / `balance` / `FormulaError` |
| `chemkit.data` | 414 | 数据表加载、kinetics 字段展开、晶格阴离子电对与 van't Hoff 焓的 Hess 派生、Henry 温度修正 | `load_tables()` / `Tables` |
| `chemkit.candidates` | 889 | 引擎地基：全局常数唯一事实源、`Cand` 与 `logK(T)` 三通道、半反应配平、拆盐兜底、Hess 派生候选 | `Cand` / `logK_T()` / `build_derived()`；派生候选的 logK **必须**等于基候选的线性组合（`tools/hess_audit.py` 常驻核验，H46 病灶即此） |
| `chemkit.normalize` | 272 | 投料 → 引擎账本（电离 / 中和 / 浓酸分子态曲线 / 电离映射） | `normalize()`；`Σz·n + He + c_H = 0` 的建立处 |
| `chemkit.acidbase` | 383 | **精确酸碱求解地基**：质子化族（pKa 连通块）+ 逐级 pKa 严格多级分布 + 电荷平衡闭式解 pH | `build_families()` / `dist_charge()` / `charge_pH()`；13/13 教科书对照（NaAc 8.882 / NH₄Cl 5.125 / NaHCO₃ 8.350 …）。**诊断/校验用，未接入主路径**（替换判据见 architecture §7 F） |
| `chemkit.speciation` | 576 | pH 机器：缓冲滴定记账（强度序堆 + 增量移动）、四分支 pH 估计、多级形态分布、分子态强酸再电离 | `estimate_state()` / `estimate_pH()`；族总量守恒是滴定记账的基本不变量（`tools/titr.py` 逐调用对账） |
| `chemkit.templates` | 961 | 候选宇宙：红氧模板四级缓存 + 三层枚举 memo | `enumerate_candidates()` |
| `chemkit.engine` | 1786 | 判定主循环：S 评估、二分程度求解、驱动力排序 walk（震荡冻结 / 膜封锁 / 气体扫气 / 爬行外推 / 联立跳步） | `judge()`；`_exec` 是 walk 步与联立跳步共用的唯一记账入口 |
| `chemkit.joint` | 385 | 联立平衡求解器：耦合平衡集作为非线性方程组（变量 = 各平衡净程度），阻尼 Newton + 数值 Jacobian + 回溯行搜索 + 非负投影，奇异走 Tikhonov 最小二乘 | `joint_solve()` / `_solve_ph()`；三态返回（内点跳步 / boundary 冻结 / fail 无操作） |
| `chemkit.equations` | 1449 | 方程式装配：步骤聚合、瞬态中间体折叠、H/O 配平、系数有理化、净离子方程式 | `_build_equations()`；化简只是**候选生成器**，唯一裁决是守恒闸门 `_drop_budget_ok` |
| `chemkit.system` | 390 | 高层 API：`Engine` 对象、`Reaction` 包装、`System` 连续投料、共用求解管线 | `Engine` / `Reaction` / `System` / `react()`（三层 API 直落同一管线，无包装链） |
| `chemkit.thermo` | 300 | 独立温度模块：反应热事后分析与绝热耦合（平衡 ⇌ 能量外层温度不动点） | `analyze()` / `coupled()`；与平衡路径完全解耦 |
| `chemkit.converg` | 243 | 收敛质量基准与全量差分：残差画像 dump / diff / top | `resid_live` 口径（排除 frozen / slow / blocked / trace 四类），动 walk 语义前的差分基线制度 |
| `chemkit.testsuit` | 651 | 测试套件（默认不随包加载） | `python -m chemkit.testsuit [用例库] [--out 结果.json]` |

**审计工具（`tools/`，非包内代码）**：`eqcheck.py`（期望方程式精确有理守恒
核验）、`hess_audit.py`（派生候选 Hess 一致性）、`hydrate_audit.py`（水合/
无水相歧义 + 边界敏感性）、`data_audit.py`（数据库机械审计）、`charge_audit.py`
（账本电荷守恒）、`escape_audit.py`（逸出气体质量守恒）、`titr.py`（滴定族
总量对账）、`phlog.py`（pH 估计轨迹）、`osc.py`（走步反复逼近画像）、
`phdiag.py`（pH 机器 vs 电荷平衡对照）、`cliff.py`（pH 连续性探针）、
`topres.py` / `perf.py` / `snap.py` / `xver.py`（残差、性能、快照、跨版本
对照）、`quick.py`（子集测试）、`case.py` / `cands.py` / `extent.py`（单例深探）。

### 数据表（`chemkit/data/`）

| 文件 | 内容 | 规模 |
|---|---|---|
| `couples.json` | 氧化还原电对（ox/red/E0/n/dH/kinetics 闸门；含配合物形态电对） | 163 条 |
| `pka.json` | 酸碱解离常数（强酸 pKa ≤ 0 在规范化即拆解） | 89 条 |
| `ksp.json` | 溶度积（pair 阳阴离子对 + `slight` 微溶标记 + **`phase` 相声明**） | 230 条 |
| `beta.json` | 配合物累积稳定常数 logβ | 158 条 |
| `substance_ex.json` | 物质附加属性（solid/gas/ions 形态、浓酸 conc_forms） | 201 条 |
| `thermo.json` | 物质 ΔHf°（单一事实源：Hess 派生 + 热分析；裸名 = 账本态，`X(g)` = 逸出气相态） | 320 条 |
| `overrides.json` | OVERRIDE 逃生舱（原理覆盖不到的显式注册反应） | 27 条 |
| `tests.json` | 回归测试库 | 1173 例 |

元素覆盖：**78 种元素**有反应数据（H、碱金属/碱土全族、稀土全 17 种、
锕系 9 种、过渡金属全族、铂族全族、后过渡、类金属、卤素/氧族/氮族/碳族
常见非金属）；未覆盖元素（稀有气体、At/Fr/Pa/Po 等）经拆盐兜底仍可正确
电离入账。

**相声明纪律**（`ksp.json` 的 `phase`）：只做**声明**，不改数值——文献
Ksp 的测定相与账本写法可能不一致（无水盐 vs 水合物），换相是化学决策
（要看温度、过饱和度、陈化时间），必须逐条定，不许为了让某条用例变绿而
选值。

### 动力学字段（`couples.json` 唯一书写形式）

| 字段 | 含义 |
|---|---|
| `ox_closed` / `red_closed` | 该电对作氧化剂/还原剂方向封闭 |
| `below_T` / `below_T_only_red` | 低于该温度（K）冻结 / 仅对指定还原剂冻结 |
| `closed_with_red` / `closed_with_ox` / `closed_except_red` | 只对指定搭档封闭 / 除指定还原剂外封闭 |
| `red_pH_min` / `ox_pH_max` | pH 形态闸门（方向感知） |
| `h2o_red_oh_min` | 氧化水的浓碱解锁阈值（如浓碱制锰酸钾） |
| `rev_gate` | 逆向驱动闸门 |

`judge()` 的 conditions 传 `"kinetics": False`（或 react/System 的同名参数）
可回到纯热力学基线（无限时间）做对照。

## 测试

```bash
python -m chemkit.testsuit                  # 全量：1173 例 / 1294 条断言（含酸碱锚点）+ 环闭合
python -m chemkit.testsuit my.json          # 自定义用例库
python -m chemkit.testsuit --out result.json  # 结构化结果（cases/summary/checks）
python tools/quick.py T99 N15 E42           # 子集（按名称前缀/序号），带单例计时
python tools/quick.py --fails               # 复跑上一轮失败项
```

用例断言维度：`changed` / `reacted` / `degree`（int 或 str）/ `has` / `has_not` /
`has_range` / `has_any` / `has_initial` / `has_in_initial_not` / `ph` / `ann` /
`override` / `eq`（净方程）/ `eq_has`（分步方程必须出现者）。

**用例库的写法约定**：断言键"出现即断言、不出现即不管"——想断言"某物种
不存在"要显式写 `null`（如 `"reacted": null`）。`has` 是**下限**断言
（终态 ≥ 给定 mol），`has_range` 是区间，`has_not` 是上限。

**方程式断言是"规范形"比较**（v0.5.0）：`eq`/`eq_has` 先做
**H⁺/OH⁻/H₂O 归一**——`H₂S + OH⁻ → HS⁻ + H₂O` 与 `H₂S → HS⁻ + H⁺` 是同一
化学（只差水的自电离关系），两者都判通过。归一只挪 H/O 的记账位置，
元素/电荷守恒与配平保真：少一个原子、电荷不平、系数不同照样判错；方向
无关（正/逆反应同一）。

**当前基线**：**1290/1294 PASS**（4 例红，均为已知项，见下；1294 = 1202 用例断言 + 92 条酸碱地基断言）。这些数字是
`python -m chemkit.testsuit` 自报的口径，任何改动都必须跑全量复核。

| 长期红例 | 病根（详见 architecture.md §7） |
|---|---|
| `F33 Na[Ga(OH)4]+CO2适量` | 镓酸盐-碳酸盐竞争通道的配比偏差（同铝酸盐通道，铝酸盐 4 例已修） |
| `H89 Ca(HCO3)2+Ca(OH)2` | 分步方程缺一条（通道合并顺序） |
| `P11 PbSO4+Na2CO3->PbCO3` | 分步方程缺一条（通道合并顺序） |
| `V16 SO2+CaCl2溶液 仅溶出不沉淀` | 亚硫酸钙**选相**（无水 vs 二水合物，pKsp 7.2/6.5）——化学决策未定 |

## 性能

**口径**：`python -m chemkit.testsuit` 自报的单遍墙钟（含进程内静态层
构建的一次性摊销之外的热路径），共享沙箱抖动约 ±10%，与机器负载相关；
逐例判读不要用收敛 dump 的 `ms`（含 GC 停顿，单例可被放大数千倍），
要用 `iters` / `steps`。

最近一次全量（1173 例，Windows 沙箱）：

| 分界线 | 用例数 | 占比 |
|---|---|---|
| > 10 ms | 881 | 75.1% |
| > 20 ms | 536 | 45.7% |
| > 50 ms | 186 | 15.9% |
| > 100 ms | 84 | 7.2% |
| > 500 ms | 1（Hg2_1 562 ms） | 0.1% |

- 均值 **34.4 ms**；中位 18.0 ms；P90 80.5 ms；全套 40.4 s。
- 最慢一档：Hg2_1（240 步）、Ni41（159 步）、D32（104 步）、P11（94 步）、
  N06、X10、N30、TS04 —— **慢例全部是"步数多"**，不是单步贵。
- 结构性事实：**22.7% 的用例仍在"反复逼近"**（同一净反应重复执行 ≥3 次，
  见 `tools/osc.py`）。这是逐候选贪心（Gauss-Seidel 坐标松弛）的固有代价，
  也是当前性能与精度的共同天花板——0.5.0 的门槛（算法飞跃）就在这里，
  实验记录见 architecture.md §7 M。
- 每版性能手段的逐年记录见 `changelog.md`（0.3.5 起逐版列出）。
- 静态层：首次判定时构建（电对配对、Hess 派生候选、模板，~3 s），之后
  跨判定驻内存复用；`Engine` 对象长期驻留是推荐用法。

## 已知局限（诚实清单）

1. **6 例断言长期红**（上表）——都不影响主体判定，且各有明确病根。
2. **走步的反复逼近**：22.7% 用例带重复签名；终态对，但收敛靠几何级数，
   精度上限 ~ρⁿ（N15 的 `1.999CO_2` 即此残差的呈现）。
3. **pH 机器是教科书近似**：走步内层（二分探针）用四分支启发式——Ksp /
   气相储库决定的 pH（Cu(OH)₂ 饱和液、浓碱溶铝等）由启发式给出，
   `charge_pH`（精确质子条件）在这些体系上给出的反而是"把账本拉回自洽"
   的假根。**呈现层已按 B3 判据切到精确解**（`exact_proton_pH`：账本
   电中性自洽 ∧ 无固相/水解阳离子储库 ∧ 有可再分配质子化族）——报告的
   `final_pH` 在该档是精确值，B3 档与精确解的偏差例数 **25 → 0**；
   判据、三次接线实测与叙述政策争议见 architecture §7 N/O。
4. **固相选相**：`phase` 只做声明，221 处固相析出量 <20% 投料的体系对
   ±0.5 个 pKsp 敏感（`tools/hydrate_audit.py` 的 ⚠ 清单）。
5. **温度域**：液态水域内（0–100 °C，超临界不建模）；气压非正报错。

## 版本

当前包内版本号 **0.4.4**。版本历史与逐版细节见 **`changelog.md`**；
未完成项、失败实验记录与下一步规格见 **`architecture.md` §7**。

**升级纪律**：0.5.0 的门槛是**求解算法的实质飞跃**（性能、泛用性、鲁棒性
同时提升），不是功能累积；门槛未过就停在 0.4.x——已落地的守恒重构与根因
修复不回退，但不改版本号。0.5.0 阶段已入档的成果见 changelog.md 顶部条目。
