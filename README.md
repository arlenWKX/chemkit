# chemkit

水溶液反应判定和产物计算引擎。

## 安装

Python>=3.12

```bash
git clone https://github.com/arlenWKX/chemkit.git
cd chemkit
pip install sympy
python -m chemkit.testsuit     # 运行全部 833 个回归用例 + 数据一致性检查
```

## 使用

### 单次反应: react()

```python
import chemkit

r = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0)
r.equation      # 'H^+ + OH^- -> H2O' 的净离子式（H+ 正则形显示）
r.reacted       # True
r.degree        # 'complete' / 'incomplete' / 'hardly' / 'none'
```

参数：
- `substances`（{化学式: mol}）
- `V`（L, 默认1.0L）
- `T`（K, 默认298.15K, 范围273.15–373.15 K）/ `T_C`（°C, 默认25°C, 范围0-100°C）
- `p`（kPa，默认101.3kPa）

### 连续投料: System

```python
sys = chemkit.System(V=1.0, T_C=25)
sys.add("NaOH", 0.1)          # 纯水 + NaOH
r = sys.add("HCl", 0.15)      # 再投 HCl，体系按累计投料重新平衡
print(r.pH, r.equation)
print(sys.feeds)              # {'NaOH': 0.1, 'HCl': 0.15}
print(sys.history)            # 历次 Reaction
```

### 结果对象：Reaction

| 属性 | 含义 |
|---|---|
| `reacted` | 是否反应 |
| `degree` | 反应程度（complete/incomplete/hardly/none） |
| `consumed` / `produced` | 净消耗 / 净生成 {化学式: mol} |
| `final` | 终态组成（H₂O 为溶剂不计入） |
| `pH` | 终态 pH |
| `equation` | 净离子方程式（无显著反应为 None） |
| `equations` | 多步离子方程式列表（按贡献降序） |
| `escaped` | 逸出气相 {化学式: mol}（超过亨利溶解上限的自产气体） |
| `annotations` | 标注（如 `"slow"`：该反应动力学缓慢） |

## 化学式写法

采用类LaTeX语法，上下标严格：`H_2SO_4`、`Ca(OH)_2`、`Fe^{3+}`、`[Fe(CN)_6]^{4-}`、
`NO_3^-`、`"NaCl"`、`"HCl"`、`"CO_2"`

## 注意点

- **热力学为体，动力学为例外**：反应性由平衡推出，热力学无法解释的硬事实以动力学标记。`chemkit.engine.KINETICS = False` 可回到纯热力学基线
- **气体**：外界恒压惰性环境；自产气体超过 H(T)·p_ext 溶解上限即逸出消失（`escaped`），投料气体视为持续供给不逸出。
- **纯溶解/电离算反应**：NaCl 溶解、弱酸电离等目前会给出 reacted=True
- **OVERRIDE**：极少数原理覆盖不到的反应进行显式注册
