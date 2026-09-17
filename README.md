# chemkit

**水溶液反应判定与产物计算引擎**（纯 Python ≥ 3.12，**零第三方依赖**）。

给它一份投料与条件，它回答：会不会反应、反应到什么程度、终态剩下什么、pH 多少、
净离子方程式怎么写、有没有气体逸出、放热多少。

```python
import chemkit

r = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)
r.changed, r.reacted, r.degree        # True True 2（完全反应）
r.net_equation.plain()                # '2H+ + Zn -> H2 + Zn2+'
r.pH                                  # 终态 pH
```

设计哲学：**从热力学数据推出反应行为**——电对 E°、pKa、Ksp、logβ、Henry 常数
经 Hess 定律派生成候选反应，而不是维护一张庞大的"反应事实表"。热力学解释不了
的硬事实（钝化、光/催化依赖、动力学冻结）才以 kinetics 闸门显式标注，与热力学
分离。

| | |
|---|---|
| 化学范围 | 酸碱/缓冲、沉淀溶解、配位、氧化还原、气体逸出、绝热热效应 |
| 数据规模 | 电对 163 · pKa 89 · Ksp 230 · logβ 159 · ΔHf 320 · OVERRIDE 27 |
| 元素覆盖 | 数据表 79 种元素；未收录元素经元素周期表级拆盐兜底仍可电离入账 |
| 回归套件 | 1176 条化学用例 + 内部自检 = **1299 条断言**，两套 pKw 约定均须全绿 |
| 依赖 | 无（含手写精确有理数配平，v0.3.5 起不依赖 sympy） |

当前版本 **0.5.2**；0.5.3 门槛见下方"升级纪律"。

---

## 安装与自检

```bash
git clone https://github.com/arlenWKX/chemkit.git
cd chemkit
python -m chemkit.testsuit        # 全量套件 + 环闭合检查（~40–70 s）
python -m chemkit.testsuit -h     # 用法
```

> **基线如实声明**：全绿只等于"已知问题都被解释过"，不等于没问题。
> 任何改动都必须跑全量复核；"全绿"也不等于"没变"——判定语义的位移要靠
> `chemkit.converg` 的全量差分与 `tools/xver.py` 的逐例对照才看得见。

## 五分钟上手

### 统一入口：`Engine`（持表一次，长期复用）

```python
import chemkit

eng = chemkit.Engine()                              # 加载包内数据表 + 构建缓存
r   = eng.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0) # 一步式
sys = eng.system(V=1.0, T_C=25)                     # 连续投料（共享同一份缓存）
sys.add("NaOH", 0.1)
r2  = sys.add("HCl", 0.15)                          # 按累计投料整体再平衡
sys.feeds, sys.history, sys.result                  # 投料账 / 历次结果 / 最近一次
```

函数式入口 `chemkit.react / chemkit.System` 与 Engine 方法完全同语义
（隐式使用全局默认表）；需要换表或隔离缓存时用 `Engine(tables=...)`。

### 参数

| 参数 | 含义 | 默认 |
|---|---|---|
| `substances` | 投料 `{化学式: mol}` | — |
| `V` | 溶液体积（L） | 1.0 |
| `T` / `T_C` | 温度（K）/（°C），限液态水域 273.15–373.15 K | 298.15 |
| `p` | 外界气压（kPa），影响气体逸出阈值 | 101.3 |
| `isothermal` | 恒温（默认）；`False` = 绝热耦合（温度不动点迭代） | True |
| `kinetics` | 动力学层；`False` = 纯热力学基线（无限时间） | True |
| `gas_escape` | 自产气体逸出；`False` = 闭口体系 | True |
| `tables` | 自定义数据表 | 包内单例 |

### 结果对象：`Reaction`

```python
r.changed / r.reacted / r.degree      # 判定层：有无净变化 / 是否狭义化学反应 / 2·1·0
r.consumption / r.production          # 净消耗 / 净生成 {化学式: mol}
r.initial / r.final                   # 初态（post-normalize）/ 终态组成
r.pH, r.escaped                       # 终态 pH；逸出气相
r.net_equation, r.net_equation_raw    # 总净方程（精编 / 原始），Equation 结构
r.equations, r.steps                  # 分步方程式；走步过程（logK/S/extent/chem）
r.annotations, r.override             # 标注（slow/blocked）；命中的 OVERRIDE 通道
r.heat_kJ, r.dT_K, r.T_final_K        # 热效应（isothermal=False 时）
r.raw, r.consumption_raw, ...         # 引擎 raw dict 与记账层（H_excess 等）
```

方程式是**结构**不是字符串：`str(eq)` / `eq.tex()` 给 TeX（Markdown `$$…$$`
可直接显示），`eq.plain()` 给无标记纯文本。

> **v0.5.2 用户侧变更**：删去了 `judge()`。它与 `react()` 走同一条求解管线、
> 只是少了结果包装，属于重复入口；引擎 raw dict 从 **`r.raw`** 取
> （`chemkit.engine.judge` 仍是内部实现，不承诺跨版本稳定）。

### 化学式写法

类 LaTeX：`H_2SO_4`、`Ca(OH)_2`、`Fe^{3+}`、`SO_4^{2-}`（一价可写 `Cl^-`）、
`[Cu(NH_3)_4]^{2+}`、`K_3[Fe(CN)_6]`、水合 `CuSO_4·5H_2O`。
**离子无需注册即可投料**（如 `Eu^{2+}`），未收录的盐由公式结构 + 电荷平衡 +
常见氧化态表自动电离（`NdCl_3 → Nd³⁺ + 3Cl⁻`）；是否参与反应仍由数据表决定。

---

## 文档

| 文档 | 内容 |
|---|---|
| [`docs/index.md`](docs/index.md) | 文档地图：想了解什么该读哪一篇 |
| [`docs/api.md`](docs/api.md) | 用户侧 API 详解（结果对象逐项、参数、模式开关、错误行为） |
| [`docs/architecture.md`](docs/architecture.md) | 架构与模块职责（设计公理、候选宇宙、走步、pH 机器） |
| [`docs/data.md`](docs/data.md) | 数据表规格、字段与选值纪律 |
| [`docs/testing.md`](docs/testing.md) | 测试与验收：断言维度、发版前固定动作 |
| [`docs/performance.md`](docs/performance.md) | 性能口径与实测指标 |
| [`docs/limitations.md`](docs/limitations.md) | 已知局限（诚实清单） |
| [`docs/changelog.md`](docs/changelog.md) | 版本历史 |
| `agents/` | AI 侧工作记忆（轮次协议、测量纪律、逐轮工程日志）——不是使用文档 |

## 升级纪律

版本号的门槛是**求解算法的实质飞跃**（性能、泛用性、鲁棒性同时提升），
不是功能累积；门槛未过就停在原版本——已落地的根因修复不回退，但不改号。
当前 0.5.2 已入档：停滞计数口径修正、OH⁻ 模板 pKw 双折算根治、自缓冲禁令、
速率/时间尺度维度。0.5.3 的门槛是**羟基络合物族补全 + 水解真式**（详见
[`docs/architecture.md`](docs/architecture.md) §7 X-32）。

## 许可

GPL-3.0（见 `LICENSE`）。
