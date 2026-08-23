# chemkit

水溶液反应判定引擎（离子基，高中~竞赛水平）。给它一批物质和条件，它从
热力学数据（E0、pKa、Ksp、logβ）推出：反不反应、反应到什么程度、生成
什么、终态 pH、净离子方程式、哪些气体逸出。

```python
import chemkit

r = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)
print(r.reacted)     # True
print(r.equation)    # Zn + 2H^+ -> H_2 + Zn^{2+}
print(r.pH)          # 终态 pH
```

## 安装与运行

纯 Python ≥3.12，唯一第三方依赖是 sympy（配平用）。克隆本仓库后把仓库
根目录加入 `PYTHONPATH` 即可 `import chemkit`。

```bash
git clone https://github.com/arlenWKX/chemkit.git
cd chemkit
python -m chemkit.testsuit     # 运行全部 833 个回归用例 + 数据一致性检查
```

## 快速上手

### 一步式：react()

```python
import chemkit

r = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0)
r.equation      # 'H^+ + OH^- -> H2O' 的净离子式（H+ 正则形显示）
r.reacted       # True
r.degree        # 'complete' / 'incomplete' / 'hardly' / 'none'
```

参数：`substances`（{化学式: mol}）、`V`（体积 L，默认 1.0）、
`T`（K）/ `T_C`（°C，覆盖 T）、`p`（外界气压 kPa，默认 101.3，影响气体
逸出阈值）。

### 连续投料：System

```python
sys = chemkit.System(V=1.0, T_C=25)
sys.add("NaOH", 0.1)          # 纯水 + NaOH
r = sys.add("HCl", 0.15)      # 再投 HCl，体系按累计投料重新平衡
print(r.pH, r.equation)
print(sys.feeds)              # {'NaOH': 0.1, 'HCl': 0.15}
print(sys.history)            # 历次 Reaction
```

### 结果对象：Reaction

人类可读层（化学习惯，含显式 H⁺/OH⁻）：

| 属性 | 含义 |
|---|---|
| `reacted` / `degree` | 是否反应 / 程度（complete · incomplete · hardly · none） |
| `consumed` / `produced` | 净消耗 / 净生成 {化学式: mol} |
| `final` | 终态组成（H₂O 为溶剂不入账） |
| `pH` | 终态 pH |
| `equation` | 净离子方程式（无显著反应为 None） |
| `equations` | 多步离子方程式列表（按贡献降序） |
| `escaped` | 逸出气相 {化学式: mol}（超过亨利溶解上限的自产气体） |
| `annotations` | 标注（如 `"slow"`：该反应动力学缓慢） |

引擎记账层（raw，H₂O 不入账、H⁺/OH⁻ 合记为带符号质子账本）：

`consumed_raw` · `produced_raw` · `final_raw` · `H_excess_raw`（正=残余
游离强酸 mol，负=游离强碱）· `steps`（逐步过程：kind/equation/logK/S/
extent/conversion）· `raw`（judge 原始 dict）。

```python
r = chemkit.react({"CaCO_3": 0.5, "HCl": 1.0}, V=1.0)
r.equations     # 多步离子方程式
r.escaped       # {'CO_2': 0.466...}  逸出的 CO2
r.final         # 溶液里只剩溶解态 CO2
```

## 底层 API

```python
from chemkit import load_tables, judge

T = load_tables()           # 加载包内数据表（可复用，线程内只建一次）
r = judge([{"name": "HCl", "mol": 1.0}, {"name": "Fe", "mol": 0.5}],
          {"V_L": 1.0, "T_K": 298.15, "p_kpa": 101.3}, T)
```

`conditions` 键：`V_L`、`T_K` / `T_C`、`p_kpa`、`c_H` / `c_OH` / `pH`
（指定初始介质）。

辅助：`chemkit.balance(["H_2","O_2"], ["H_2O"])` 配平方程式；
`chemkit.FormulaError` 化学式解析错误。

## 化学式写法

内部表示用 ASCII：`H_2SO_4`、`Ca(OH)_2`、`Fe^{3+}`、`[Fe(CN)_6]^{4-}`、
`NO_3^-`。投料用常见名称/化学式均可（`"NaCl"`、`"HCl"`、`"CO_2"`）。

## 行为边界（读结果前要知道）

- **热力学为体，动力学为例外**：反应性由平衡推出；只有热力学无法解释的
  硬事实（钝化、歧化封锁、温度门槛、单质硫室温惰性等）才以动力学标记
  关闭，命中时结果带 `annotations: ["slow"]` 或直接不反应。
  `chemkit.engine.KINETICS = False` 可回到纯热力学基线。
- **慢的两档语义**：恒慢通道（如 O₂ 氧化酸性 Mn²⁺）不执行、只标注；
  让位档通道（如硫化矿的氧化溶解）慢于离子交换但可发生——有更快
  竞争通道时让位，无竞争时照常执行并标注 slow（ZnS 悬浊液遇 CuSO₄
  发生沉淀转化生成 CuS 而非氧化还原，PbS 则被 H₂O₂ 氧化为 PbSO₄）。
- **温度域**：液态水 273.15–373.15 K；部分通道有温度门（如单质硫的
  氧化/歧化在约 97 °C 以上解锁）。
- **气体**：外界恒压惰性环境；自产气体超过 H(T)·p_ext 溶解上限即逸出
  消失（`escaped`），投料气体视为持续供给不逸出。
- **hardly 也是反应**：degree 按限量试剂转化率分档，trace 量反应
  （如微溶盐溶解再沉淀）会给 `hardly`/`incomplete` 而非 `none`。
- **纯溶解/电离算反应**：NaCl 溶解、弱酸电离目前会给出 reacted=True
  （这符合"体系发生了变化"的事实；如需区分"化学反应"与"形态变化"，见
  architecture.md 的候选类别，proton/dissolve 类步骤即形态变化）。
- ** OVERRIDE**：极少数原理覆盖不到的反应走显式注册（结果带
  `override` 字段），这是逃生舱而非常态。

## 文档

- `architecture.md` — 引擎设计详述（六条公理、候选体系、晶格电对派生、
  动力学三层语义、气体模型、性能结构、数据纪律）
- `data/tests.json` — 818 个带离子方程式断言的回归用例（测试器合计
  833 条通道），是最好的行为示例库
