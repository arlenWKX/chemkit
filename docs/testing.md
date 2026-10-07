# 测试与验收

## 1. 怎么跑

```bash
python -m chemkit.testsuit                  # 全量：1176 用例 + 内部自检 + 环闭合
python -m chemkit.testsuit -h               # 用法
python -m chemkit.testsuit my.json          # 自定义用例库
python -m chemkit.testsuit --out result.json  # 结构化结果（cases/summary/checks）
```

**1299 条断言** = 1378 条化学用例 + 123 条内部自检（温度域校验、酸碱地基
断言、高层 API 自检 22 通道、环闭合检查）。

用例库的写法约定：断言键**出现即断言、不出现即不管**；`has` 是下限、
`has_not` 是上限、`has_range` 是区间。想断言"某物种不存在"要显式写 `null`。

## 2. 断言维度

| 键 | 含义 |
|---|---|
| `changed` / `reacted` / `degree` | 判定层（`degree` 支持 `2/1/0` 或 `complete/incomplete/hardly`） |
| `has` / `has_not` / `has_range` / `has_any` | 产物量的下限 / 上限 / 区间 / 任一形态下限 |
| `has_initial` / `has_in_initial_not` | 初态（post-normalize）下限 / 不应含 |
| `ph` | 终态 pH 区间 |
| `ann` / `override` | 标注与逃生舱命中 |
| `eq` | 净方程（`null` = 要求无净方程） |
| `eq_has` | 分步叙述**必须覆盖**的变换（见 §4） |

## 3. 方程式断言是"规范形"比较

`eq` / `eq_has` 比较前先做 **H⁺/OH⁻/H₂O 归一**：`H₂S + OH⁻ → HS⁻ + H₂O` 与
`H₂S → HS⁻ + H⁺` 是同一化学（只差水的自电离关系），两者都判通过；方向无关
（正/逆反应同一）。归一只挪 H/O 的记账位置，**元素/电荷守恒与配平保真**——
少一个原子、电荷不平、系数不同照样判错。

## 4. `eq_has` 是"张成"判据，不是"逐字出现"

把各步与需求都写成净变换向量（精确有理数），只要存在**非负系数**线性组合
把需求**张成**就算覆盖（`testsuit._step_spanned`）。

理由：同一条净变换可以经由不同中间体合法到达（PbSO₄ 既可走 `[Pb(OH)₃]⁻`
也可走 Pb(OH)₂/Pb²⁺；`[Al(OH)₄]⁻` 既可经碳酸盐中间体也可直接碳酸化）。
要求逐字出现等于把**路线**当化学。非负系数才是真约束——负系数意味着要
"倒着走"某一步，那不是合法叙述。

## 5. 断言必须锁化学，不能锁轨迹

判据是 `tools/fragility.py`：把投料做 `(1±1e-11)`／`(1±1e-9)` 相对扰动
（1 mol 只差 1 nmol，物理上无意义）重跑全库断言。**仍翻红的断言就是对数值
噪声敏感、锁了轨迹**，必须改写成区间/下限/守恒式。

同一把尺子也用于放行性能改动：某次二分收敛容差收紧让 42 例走步路径翻转，
但断言 0 例翻红 ⟹ 该改动照常落地。

## 6. 发版前硬门槛

**两套 pKw 约定都必须全绿**——`chemkit.core.pKw_of` 有一个锚定开关
（把 298.15 K 处的 pKw 精确钉在 14.0 的经验式 vs 未锚定式）；两者相差虽小，
却足以让处在"刀口"上的用例翻盆。命令：

```bash
python tools/dev.py anchor off && python tools/dev.py suite
python tools/dev.py anchor on  && python tools/dev.py suite
```

其余固定动作（缺一不可）：

| 动作 | 判据 |
|---|---|
| `python tools/dev.py suite` | 1299/1299（锚定）+ 1299/1299（未锚定） |
| `python tools/dev.py eqcheck` | 全库期望方程式的**精确有理守恒** 0 违规 |
| `python tools/hess_audit.py --T 273.15 --T 363.15 --with-templates` | 派生候选 logK = 基候选线性组合，0 条不自洽 |
| `python tools/data_audit.py` | 数据机械审计硬错误 0 条 |
| `python tools/dev.py perf` | 与 `converg-baseline.json` 逐项对账，位移逐条给出化学理由 |
| `python -m chemkit.testsuit`（环闭合段） | 全部环闭合检查通过 |

## 7. 收敛质量基线

`converg-baseline.json`（`chemkit.converg.dump()` 产出，受版本控制）是当前
口径的**逐例权威快照**：残差画像 + 确定性计数 + `digest_all`。

- 残差口径：`S` 一律在**求解器自己的 pH** 上评（`pH_solver`），呈现 pH 上的
  值单列 `S_pres`——旧口径拿呈现 pH 评残差，等于换一个态去质问求解器。
- 残差口径：**冻结 ≠ 已达平衡**。冻结是防震荡手段，只有 `|S| ≤ 0.1` 的冻结通道
  才豁免；冻结在强驱动上必须计入残差（否则指标全绿而化学错，见 X-38 与
  [`limitations.md`](limitations.md) §6）。
- `digest_all` 只作**辅助信号**，不是门槛：它变了不等于错，不变也不等于对；
  门槛是断言（化学）。
- 逐例判读**不要用 `ms`**（含 GC 停顿，单例可被放大数千倍）——用
  `iters` / `steps` / `sof`。口径详见 [`performance.md`](performance.md)。

## 8. 纪律

改数据、做测量、裁决标准时的工作方式（唯一可信路径、测量点自检、逐例单跑）
见 [`../agents/discipline.md`](../agents/discipline.md)。两条底线：

1. **化学门槛不得为通过而放宽**：标准只锁化学事实（文献 pH、总产量、
   特征步），不锁引擎输出、不锁中间路线；边界用例当哨兵用。
2. **失败实验必须留档**（含实测数字）进 [`../agents/log.md`](../agents/log.md)，
   否则后来者会重复付费。
