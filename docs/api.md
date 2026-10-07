# 用户侧 API

`import chemkit` 后可用的一切都来自 `chemkit.interfaces`（v0.3.3 起的唯一门面）；
`chemkit.core / data / candidates / engine / …` 均为实现细节，不承诺跨版本稳定。

```python
chemkit.__all__ == ["Tables", "load_tables", "Engine", "Reaction", "System",
                    "react", "default_tables", "TABLES", "FormulaError",
                    "balance", "__version__"]
```

> **v0.5.2 变更**：删去用户侧 `judge()`。它与 `react()` 走同一条求解管线，
> 只少了结果包装（重复入口）。引擎 raw dict 从 **`Reaction.raw`** 取；
> `chemkit.engine.judge` 仍是内部实现，不承诺稳定。

---

## 1. 三个入口，一条管线

| 入口 | 形态 | 用途 |
|---|---|---|
| `chemkit.Engine()` | 对象 | 进程内建一次、长期复用（持表 + 全部缓存） |
| `chemkit.react(subs, ...)` | 函数 | 一步式判定 → `Reaction` |
| `chemkit.System(subs, ...)` | 类 | 连续投料：每次 `add` 按累计投料整体再平衡 |

```python
eng = chemkit.Engine()                                # 默认包内数据表
r   = eng.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)   # Engine.react
sys = eng.system(V=1.0, T_C=25)                       # Engine.system（共享缓存）
sys.add("NaOH", 0.1); sys.add("HCl", 0.15)
```

`Engine` 的 `react/system` 与模块级 `chemkit.react/System` **完全同语义**，
区别只在数据表来源：前者用 `Engine(tables=...)` 绑定的表，后者用全局默认表。

**`react` / `System` 参数**（关键字专用，`*` 之后）：

| 参数 | 类型 | 默认 | 含义 |
|---|---|---|---|
| `substances` | `dict[str, float]` | —（`System` 可省） | 投料 `{化学式: mol}`，同名累加 |
| `V` | `float` | `1.0` | 溶液体积（L） |
| `T` | `float` | `298.15` | 温度（K） |
| `T_C` | `float \| None` | `None` | 温度（°C）；给了就覆盖 `T` |
| `p` | `float` | `101.3` | 外界气压（kPa），决定泡点 `c_sat = H(T)·p_ext`（气液相界） |
| `isothermal` | `bool` | `True` | `False` = 绝热耦合（温度不动点迭代） |
| `kinetics` | `bool` | `True` | `False` = 纯热力学基线（无限时间） |
| `gas_escape` | `bool` | `True` | `False` = 闭口体系（不建气相库，自产气体全部留在溶液） |
| `tables` | `Tables \| None` | `None` | 临时换表（`System` 亦可） |

**`System` 追加**：

| 成员 | 含义 |
|---|---|
| `sys.add(name, mol)` | 投料并整体再平衡，返回本次 `Reaction` |
| `sys.feeds` | 累计投料 `{化学式: mol}` |
| `sys.history` | 历次 `Reaction` 列表 |
| `sys.result` | 最近一次 `Reaction`（`history[-1]`） |
| `sys.V_L / T_K / p_kpa` | 体系固定条件 |

**数据表入口**：`load_tables(path=None)` 构建/加载表（进程内单例）；
`default_tables()` / `TABLES` 取全局单例；`Tables` 是表对象类型。
`import chemkit` 即预加载包内 `chemkit/data/`。

---

## 2. `Reaction` 结果对象

一次投料平衡后的完整快照。分两层：**人类可读层**（化学习惯，含 H⁺/OH⁻/H₂O）
与**引擎记账层**（`*_raw`：H₂O 为溶剂不入账，H⁺/OH⁻ 合记为带符号质子账本
`H_excess`，正 = 残余游离强酸 mol，负 = 残余游离强碱 mol）。

### 判定层

| 属性 | 类型 | 含义 |
|---|---|---|
| `changed` | `bool` | 体系是否发生显著净变化（宽口径：含纯溶解、电离、水解等形态变化） |
| `reacted` | `bool` | 是否发生**狭义化学反应**：氧化还原 / 中和 / 跨投料的沉淀与配位等；纯溶解、电离、水解为 `False` |
| `degree` | `int` | `2` 完全反应 / `1` 可逆（部分）反应 / `0` 难反应或未反应 |

`bool(r)` 等价于 `r.changed`。

### 组成层

| 属性 | 类型 | 含义 |
|---|---|---|
| `consumption` / `production` | `dict[str, float]` | 净消耗 / 净生成 `{化学式: mol}` |
| `initial` | `dict[str, float]` | 初态组成（post-normalize：强电解质已电离、SO₃ 等已水合、酸碱中和已记账、H₂O 溶剂不入） |
| `final` | `dict[str, float]` | 终态组成（H₂O 不入；`H_excess` 已还原为 H⁺ 或 OH⁻；**总量口径 = 溶解态 + 气相**） |
| `pH` | `float \| None` | 终态 pH（OVERRIDE 路径为 `None`） |
| `gas` | `dict[str, float]` | 气相（第二相）存量 `{化学式: mol}`——恒压气相库。有气相在场时溶解态钉在泡点 `H(T)·p_ext`；进出双向可逆（超泡点鼓泡进入、低于泡点回溶补充），**不销毁物质**。`final`/`production`/`consumption` 均为**总量口径 = 溶解态 + 气相**（v0.6.0 起取代旧的 `escaped`） |

### 方程式层

| 属性 | 类型 | 含义 |
|---|---|---|
| `net_equation` | `Equation \| None` | 总净离子方程式（**精编版**：痕量副过程不叙述）；无显著反应为 `None` |
| `net_equation_raw` | `Equation \| None` | 同一净差的**原始版**；与 `net_equation` 不一致 ⟺ 该体系只发生了痕量副过程 |
| `equations` | `list[Equation]` | 分步离子方程式（按贡献降序） |
| `steps` | `list[dict]` | 走步过程：`kind / equation / logK / S / extent / conversion / chem`。`chem=True` ＝该步是狭义化学反应，`False` ＝单纯形态变化（解离/水解/配位再分布） |

`Equation` 是**结构**不是字符串（`left`/`right` 系数 + `reversible` 布尔 + `kind`），
渲染发生在读取点：

```python
str(eq)          # TeX，等同 eq.tex()：可直接塞进 Markdown 的 $$…$$
eq.plain()       # 无标记纯文本 fallback：'2H+ + Zn -> H2 + Zn2+'
eq.reversible    # True 渲染 <=>（平衡过程），False 渲染 ->
```

### 标注层

| 属性 | 含义 |
|---|---|
| `annotations` | `list[str]`：`"slow"`（存在显著但动力学缓慢的通道）、`"blocked"`（表面膜封锁）等 |
| `override` | 命中的 OVERRIDE 通道 id（如 `"ko2_water"`），未命中为 `None` |

### 热效应层（`isothermal=False` 时）

| 属性 | 类型 | 含义 |
|---|---|---|
| `heat_kJ` | `float \| None` | 放热为正（= −ΔH）；`None` = ΔHf 数据不足（宁缺毋假） |
| `dT_K` | `float \| None` | 溶液温升（水比热容、`V×1000 g` 计；负 = 降温） |
| `T_final_K` | `float` | 终温（绝热耦合收敛值） |
| `thermal` | `dict` | `water_mol / missing / flags / reason`；耦合模式另含 `trace`（逐轮温度/热量记录）与 `converged` |

```python
r = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0, isothermal=False)
r.heat_kJ, r.dT_K, r.T_final_K        # 5.583  1.334  299.48
r.thermal["converged"]                # True（外层温度不动点 2 轮收敛）
```

默认 `isothermal=True`：单遍求解不计热效应，`thermal == {"mode": "isothermal"}`，
`heat_kJ is None`。`isothermal=False` 时温度反馈**真实进入求解**：独立温度模块在
化学平衡与能量平衡间迭代 `T_f = T₀ + q(T_f)/(m·c_p)` 的不动点，化学结果在自洽
终温下给出（K(T) 经 van't Hoff 随温度移动）——是状态函数意义上的绝热终态，
而不是"初温算完再补一个 ΔT 数字"。

### 引擎记账层（raw）

| 属性 | 含义 |
|---|---|
| `raw` | 引擎原始 dict（`steps` / `H_excess` / `cond` / `net_exact` / `override` …） |
| `consumption_raw` / `production_raw` | 记账口径的消耗 / 生成（H₂O 不入账） |
| `initial_raw` / `final_raw` | 记账口径的初态 / 终态（H⁺/OH⁻ 合为 `H_excess`） |
| `H_excess_raw` | 终态带符号质子账本（正 = 残余游离强酸 mol） |

> **v0.5.2 迁移**：原 `judge(subs, cond, T)` 的返回值就是 `r.raw`；
> 原条件键 `c_H / c_OH / pH`（初始强酸/强碱/指定 pH）在用户侧改由**投料**表达
> （要"初始 pH 2"就投对应酸，见 `docs/limitations.md` 的说明）。

---

## 3. 化学式写法

类 LaTeX 语法，解析器在 `chemkit.core`：

- 下标：`H_2SO_4`、`Ca_3(PO_4)_2`、`Cu_2(OH)_2CO_3`
- 电荷：`Fe^{3+}`、`SO_4^{2-}`；一价可简写 `Cl^-` / `Na^+`
- 配合物：`[Fe(CN)_6]^{4-}`、`[Cu(NH_3)_4]^{2+}`、`K_3[Fe(CN)_6]`
- 水合：`CuSO_4·5H_2O`（按元素总量解析）

**离子无需注册即可投料**：任何元素合法的带电物种直接入账（如 `Eu^{2+}`）；
未收录进数据表的盐由公式结构 + 电荷平衡 + 常见氧化态表自动电离
（`NdCl_3 → Nd³⁺ + 3Cl⁻`、`Tb(NO_3)_4 → Tb⁴⁺ + 4NO₃⁻`），保证电荷与质量守恒。
**离子是否参与反应仍由数据表决定**——没有数据的离子就是旁观离子。

---

## 4. 错误行为（宁可报错，不静默）

| 情形 | 行为 |
|---|---|
| 化学式语法错（括号未闭合、`_` 后无数字等） | `chemkit.FormulaError` |
| 温度超出液态水域 | `ValueError`（水的存在形式/活度约定在域外全部失效） |
| 气压非正 | `ValueError` |
| 无法配平 | `chemkit.balance(...)` 返回 `None` |
| 数据不足 | 相应字段给 `None`（如 `heat_kJ`）或如实标注，不猜 |

```python
chemkit.balance(["H_2", "O_2"], ["H_2O"])
# {'reactants': {'H_2': 2, 'O_2': 1}, 'products': {'H_2O': 2}}
```

`balance(reactants, products, free=None)` 用**手写免分数整数消元零空间**
（v0.3.5 起取代 sympy），`free` 可指定允许出现的额外物种。
