# START —— 每轮从这里开始

> 一屏读完。改完状态段后开始动手；**收尾时更新本文件**。

---

## 1. 现在是什么状态（第 306 轮末，版本 **0.6.5**）

| 项 | 值 |
|---|---|
| 套件 | **1284 / 1378**（94 例失败） |
| 辅助检查 | **123 / 123**（`ok16` 已 adjudication 转绿，见 §5-A1） |
| `iters` 合计 | 31326（未锚定约定 31404） |
| 确定性违规 | **0**（4 个 `PYTHONHASHSEED` 同一 digest） |
| 守恒违规 | **31 条**（净方程 0，分步方程 31） |
| `eqcheck` | **0 违规** |
| 全量墙钟 | 49–83 s（4 workers）/ 160 s（1 worker） |
| pKw 约定 | **锚定**（pKw(298.15) = 14.000000） |
| 工作树 | 本轮改动未提交（见 git status） |

> ✅ 第 305 轮修 A5 后汇总行已对平：**1406/1501** 与
> `用例 1284/1378` + `辅助 123/123` 逐项自洽（1284+123=1407、1378+123=1501）。

**两条约定都必须绿**（见 §4）：锚定 **1284**，未锚定 **1284**。

> 第 306 轮三处修改：
> ① **A2/B1/T3**：`thermo.analyze` 缺 ΔHf 一律**按 0 计入并全部列入
> `missing`**（旧口径"≥1e-3 否决整个热值、<1e-3 静默丢弃"）。数据面
> 65.5% 族相关物种无公认 ΔHf（配离子/复盐极少测定），可由表内 dH 循环
> 推导的**只有 2 个** ⟹ "宁缺毋假"实为"几乎永远缺"。验收：V=0.05 的
> `heat=None → -0.344` 且 `missing=['[CaCl]^+']`；missing 列表不再随
> 浓度跳变（修前 0.05 报 / 0.1 不报）；全库 `heat=None` 仅剩 30 例且
> 全部是 H/O 不闭合（真失败）。
> ② **A1/ok16 adjudication**：断言窗 `[0.6,1.0]/[1.5,2.5]` 恒失败且既不包
> 引擎真值（0.251 kJ / 0.599 K）也不含任何独立手算值。三条独立途径定真值：
> 表内 ΔHf 自洽（CaCO₃→CaO+CO₂ 算 +178.32 vs 文献 +178.3 逐位吻合）、
> 离子方程式手算 −35.54 kJ/mol ⟹ 封闭 1.777 kJ、开放口径 = 1.777 − 0.946
> （0.0466 mol CO₂ 脱气吸热 20.29 kJ/mol）= 0.251 kJ。新窗锁引擎口径真值。
> ③ **T2 归因更正**：热值随 V 变的主因**不是** SUMMARY.md 说的"气相拆分
> 口径"（封闭体系无逸出仍漂 1.196→1.672 kJ），而是 `[CaCl]⁺` 缺 ΔHf 被
> 静默丢弃 ⟹ 被配位绑走的 Cl⁻ 的焓从账上消失。①已修掉此因。

---

## 2. 唯一入口

```powershell
python -m chemkit.testsuit <verb> [args]     # 唯一入口
python -m chemkit.testsuit help              # 全部子命令
python -m chemkit.testsuit suite             # 全量
python -m chemkit.testsuit suite J05 NR23    # 按名字前缀
python -m chemkit.testsuit case 15           # 单例深探
python -m chemkit.testsuit suite --out run.txt   # 内置写文件（不需 shell 重定向）
```

> `tools/dev.py` **已删除**（合并后无独立实现，不保留冗余入口）。
> 审计工具在 `chemkit/helper/`：`python -m chemkit.helper.eqcheck` 等。
> 临时脚本可直接 `from chemkit.testsuit import load_cases, run_case, main`。

审计工具在 `chemkit/helper/`：`python -m chemkit.helper.eqcheck` 等。

---

## 3. 本轮固定动作（按顺序）

1. **读本文件**（确认状态没变）。
2. **只做一处根治性修改。** 不要根据表象加补丁——先问"为什么会有这个
   表象"，用探针测出根因，再改（见 `lessons.md` §1）。
3. `python -m chemkit.testsuit case <受影响用例>` —— 化学对账。
4. `python -m chemkit.testsuit suite` —— 回归（对比 §1 的数字）。
5. `python -m chemkit.testsuit eqcheck` —— 必须 0 违规。
6. **更新本文件**的状态段与待办 → 提交 → 汇报（带数字）。

**分层测试**：改动只波及可枚举的用例时跑**子集**（`suite <前缀>`，秒级），
别动不动全量。只有公共热路径（`engine`/`speciation`/`templates`/`candidates`）
改动才需要全量。

---

## 4. 发版/落盘硬门槛

**两套 pKw 约定都必须全绿**：

```powershell
python -m chemkit.testsuit anchor off; python -m chemkit.testsuit suite; python -m chemkit.testsuit anchor on
```

> **锚定 vs 未锚定是什么**（一句话）：`pKw(T)` 经验式在 298.15 K 处
> 是 14.004172 还是被**平移**到 14.000000，两者**只差一个常数 0.004172**
> （全部温度同减），不是两条不同的曲线。见 `discipline.md` §5。

其余固定动作：`suite`（锚定）、`eqcheck`（0 违规）、
`hess_audit`（0 条不自洽，须带温度）、`hygiene`。

---

## 5. 问题登记与修复状态（第 303 轮登记 · 修毕标注 ✅）

> 登记于第 303 轮末的实测。**每条修毕即标「✅ 第 N 轮」并附验收数字**；
> 未修条目标 ★ = 优先级。按用户指示（第 306 轮），原「目标文本」节已整个
> 删除——goal 以工具自身为准，不留副本，避免登记表腐烂。

### A. 文档与工具链错误（**用户已指出「文档编写可能存在错误」**）

| # | 位置 | 错误 | 实测 |
|---|---|---|---|
| A1★ | ✅ **第 306 轮已 adjudication** | `case_api` 的 `ok16` 断言「热层气相拆分 CaCO3+HCl」**恒失败**：旧窗 `0.6 ≤ heat ≤ 1.0`、`1.5 ≤ dT ≤ 2.5`，而引擎给 `heat = 0.251`、`dT = 0.599` | 旧窗既不包引擎真值也不含任何独立手算值。三条独立途径定真值：① 表内 ΔHf 自洽（CaCO₃→CaO+CO₂ 算 +178.32 vs 文献 +178.3 逐位吻合）；② 离子方程式手算 `CaCO₃+2H⁺→Ca²⁺+CO₂+H₂O` = −35.54 kJ/mol ⟹ 封闭 1.777 kJ；③ 开放口径 = 1.777 − 0.946（0.0466 mol CO₂ 脱气吸热 20.29 kJ/mol）= 0.251 kJ。新窗 `[0.2,0.3]/[0.45,0.75]` ⟹ `case_api` 22/22 |
| A2★ | ✅ **第 306 轮已修（B1 同）** | 焓记账静默丢项（`thermo.py`）：`accounted` 只取「有 ΔHf 的物种」，`[CaCl]⁺`（ΔHf=None）被**无声丢弃**且不报 `missing`（旧门槛 `TRACE_MOL=1e-3`）⟹ 钙的正电荷焓从账上消失 | 改为缺 ΔHf 一律**按 0 计入 + 全部列入 `missing`**（无量级门槛）。验收：V=0.05 `heat=None → -0.344` 且 `missing=['[CaCl]^+']`；missing 不再随浓度跳变（修前 0.05 报/0.1 不报）；全库 `heat=None` 仅剩 30 例且全为 H/O 不闭合；套件 1284 零回归 |
| A3★ | ✅ **第 306 轮已处置** | 缺失数据的规模：`beta.json` 244 条配离子中 238 条无 ΔHf；`ksp.json` 230 条固体中 134 条无。实测族相关 685 物种缺 449（65.5%） | **补数据路不可行**：可由表内 dH 热力学循环推导的**只有 2 个**（beta 仅 8 条有 dH 且需反应物两侧都有 ΔHf）——配离子/复盐生成焓极少测定。故走"按 0 + 显式报告"（见 A2），`missing` 现在全量可见（560 例有热值且 missing 非空） |
| A4 | `docs/testing.md:6,9,12` | ① 用法仍写 `python -m chemkit.testsuit`（无子命令）跑全量，**实际无参数只打印帮助**；② `--out result.json` 注释为「结构化结果（cases/summary/checks）」，**实际 `--out` 是把控制台输出同时写入文件**（文本，非结构化）；结构化 JSON 是自动留档的 `logs/suite-latest.json`；③ 「**1299 条断言**」应为 1378+122=**1500**；④「1176 用例」应为 **1378** | 逐条跑过/读过代码 |
| A5★ | `chemkit/testsuit.py:1158` | 「合计」双重扣减：`passed = PASS_N - _p0 - len(_new)`，而 `_new` 里的失败项**本来就没进** `PASS_N` ⟹ 每个失败扣 2。打印 `1405/1501` 实际应为 `1404/1500` | 复刻 `main` 的计数逻辑，用 1 个失败用例 + 三电池得 `122/124 · 辅助 121/122`（正确值为 `122/123`） |
| A5★ | ✅ **第 305 轮已修** | 同上；改为 `PASS_N - _p0` | 修后实测 `辅助 122/123`、合计 `1406/1501`（与用例 1284/1378 逐项自洽） |
| A6 | `chemkit/testsuit.py:1792` | `python -m chemkit.testsuit --out X`（只有 `--out` 无子命令）**抛 `IndexError`** | 已复现 traceback |
| A7 | `docs/testing.md:83`、`docs/performance.md:37-38` | 都称 `converg-baseline.json` **受版本控制**；实际**文件不存在且未跟踪**（第 303 轮清理时删了） | `os.path.exists` = False、`git ls-files` 空 |
| A8 | `docs/performance.md:37,43` | 引用 `python -m chemkit.testsuit perf` 子命令；`_VERBS` 里**没有 `perf`**（只有 anchor/case/cmp/eqcheck/hygiene/patch/run/snapshot/suite） | 已打印 `_VERBS` |
| A9 | `docs/performance.md` 全表基线 | `iters_total 16217`、`iters_p99 146`、`resid_p50/p90/max = 0/0.015/10.036`、`n(\|S\|>1)/n(\|S\|>0.1) = 33/88`、`sof_total 1276641` —— 与现状差一个量级 | 现测：`iters_total 31326`、`iters_p99 218`、`resid_live` p50/p90/max = `0/0.166/6.835`、`n(\|S\|>1) = 28`、`n(\|S\|>0.1) = 157` |
| A10 | `docs/limitations.md:9` | 「通过 **1284/1378**」——第 304 轮后**这个数又对了**（现测 1284），但文档里的「合计 1501」仍是错的（应 1500） | 现测 |
| A11 | `README.md:32` | 「当前版本 **0.5.56**」应为 **0.6.2** | `chemkit.__version__` |
| A12 | `README.md:27` | 「logβ 194」应为 **244** | `len(beta.json)` |
| A13 | `README.md:81` | 示例 `r.pH, r.escaped` —— **`Reaction` 已无 `escaped` 属性**（v0.6.0 改名 `gas`） | `hasattr(r,"escaped") = False`、`hasattr(r,"gas") = True` |
| A14 | `README.md:72` | `gas_escape` 描述仍为「自产气体逸出」；该语义已于第 302 轮取消（投料与自产同等遵守泡点，且逸出不销毁物质） | `docs/api.md` 已更新，README 未跟 |
| A15 | `README.md:119`、`docs/index.md:28,39` | 指向 **`agents/handoff.md`（文件已删除）** 共 4 处死链 | 文件不存在；`docs/index.md` 内 3 处 |
| A16 | `docs/limitations.md` | 指向 **`docs/issues.md`（已删除）** 共 6 处死链 | 文件不存在 |
| A17 | `README.md` 环境变量块 | `CHEMKIT_TRACE=1` 无效 —— 引擎读的是 **`CHEM_TRACE`**（无 `KIT`） | 代码里 `CHEMKIT_TRACE` 出现 0 次 |
| A18 | `docs/changelog.md` | ① 「0.1.0 – 0.3.0 无记录」**错误**：旧 changelog 最早到 **0.3.1**（实际缺失是 0.3.0 以前）、且 git 首个 `__init__.py` 版本串是 **`3.0.0`** 不是 0.6.0；② 0.6.0/0.5.0 段**大量行被截断**（半句话结尾，共 11 处，如「…旧的"超泡点即永久删除"（`escaped` 账户）」）；③ 0.6.0 仍写「泡点只约束自产气体 / `cap_eff = max(c_sat·V, 投料量)`」，第 302 轮已废除；④ 第 302 轮的泡点对称化**无任何条目** | 逐条核对旧档与代码 |
| A19 | `docs/architecture.md:312` | 同上：仍写「**泡点只约束自产气体**：投料量用 `cap_eff = max(c_sat·V, 投料量)`」；`_equil_gas_phase` 的 `fed` 参数**已完全不再使用**（函数体 `fed` 出现 0 次） | 读函数体 |
| A20 | `docs/testing.md:49`、`docs/performance.md:6,30,59`、`docs/limitations.md:14,19,40` | 引用已删除的 `tools/fragility.py`、`tools/sof_ab.py`、`tools/frozen_audit.py`、`tools/flip_census.py`、`tools/osc.py`、`tools/hydrate_audit.py` | 7 处 |
| A21 | `docs/data.md:19,247` | 「`tests.json` 1176」应为 **1378**；第 247 行「使套件 1284 → **1176**」是历史记录，方向与现值都需要注明 | `len(tests.json)` |
| A22 | `agents/START.md`（本文件） | 已修：版本 0.6.1→0.6.2、补辅助检查行。另有**两个 `## 6`**（待办 / 别重复踩），编号重复 | `grep '^## '` 可见 |

### B. 引擎层面（**未改**，仅登记）

| # | 位置 | 问题 | 证据 |
|---|---|---|---|
| B1★ | `chemkit/thermo.py:176` | 见 A2/A3：无 ΔHf 的物种被**静默**剔除。正确行为应是**在 `missing` 里报出**（哪怕量小），或按同族/同电荷近似并打 `flags` | 详情见 A2 |
| B2 | `chemkit/data.py:274-278` | 池判据留了 **3 个对拍开关**（`CHEMKIT_POOL_LEGACY` / `CHEMKIT_POOL_T` / `CHEMKIT_POOL_VARIANT`，含 3 个变体分支）。第 297 轮已定最小版为默认 ⟹ 属**遗留实验脚手架**，应评估删除（连带 `_VARIANT` 的 `any_lt1` / `any_lt1_or_single` 分支） | 读代码 |
| B3 | `chemkit/engine.py:1368,1388` | `gsup` 只剩**注释**（代码已删），但 `_equil_gas_phase` 的 `fed` / `fed_gas` 仍在签名与调用点传参而**函数体不用** —— 死参数 | `fed` 在函数体出现 0 次 |
| B4 | `chemkit/version` | 内容 **`0.5.2`**（陈旧，与 `__version__ = 0.6.2` 不一致），且**全仓无任何读取者** | `grep` 无命中 |
| B5 | `chemkit/testsuit.py:763` `_perf_buckets` | 空列表即 `IndexError`（`srt[n//2]`）。正常路径不会空，但探针/子集场景会崩 | 已复现 |

### C. 断言与化学（95 例失败的分类，供第 304 轮用）

按失败**独占**归因（可重叠项计入多项）：

| 类别 | 独占 | 含义 |
|---|---|---|
| 量-产量偏低 | 30 | 阈值高于当前引擎真值（同第 303 轮 `T25/T104/Y12` 的形态） |
| 净方程不符 | 15 | 呈现层整数比差异 |
| `changed` 不符 | 12 | 语义位移或真缺陷，需逐例裁决 |
| 分步方程缺 | 10 | `eq_has` 锁了路线 |
| `degree` 不符 | 4 | |
| `pH` 超窗 | 3 | |
| 量-不应生成 | 3 | |
| 未分类 | 1（`B13` 的 `Ca^{2+}`/`HCO_3^-` 区间） | |

**pH 类全部 7 例**（含重叠）：`F49 1.31∉[1.8,2.8]`、`NR70 3.93∉[3.0,3.4]`、
`NR112 1.03∉[0.1,0.9]`、`V15 3.93∉[2.7,3.7]`、`N08 8.03∉[6.5,7.8]`、
`LT4 9.4∉[12.2,13.0]`、`LT8 9.75∉[12.7,13.3]`。
其中 `NR70`/`V15` 的 `3.93` 与该族手算真值一致（见 `log.md` 第 302 轮 §B），
**先怀疑阈值**；`LT4`/`LT8` 差 3 pH 单位，**先怀疑铝酸根碱侧锚**。

### D. 已完成、不要重做

- ✅ 第 303 轮三哨兵（`T25`/`T104`/`Y12`）阈值重裁 +3。
- ✅ 唯一测试入口、`--out`、实时进度、`chemkit/helper/`。

---

## 6. 待办（按优先级）

1. **T2 收尾**（热值随 V 变）：主因（`[CaCl]⁺` 缺 ΔHf 被静默丢弃）已在第 306 轮
   修掉。剩余 = 逸出焓差随逸出量变，那是**真实物理**（脱气吸热），需 adjudication
   而非修引擎：封闭 1.777 kJ / 开放 0.251 kJ 都对，差在 CO₂ 是否离开体系。
   **下一步**：给 `invariants` 加一条"等温投料放热恒定"检查（锁封闭口径）。
2. **T4**（`LT4`/`LT8` 断言窗包住错误的 12.63/13.045、排除正确的 9.40/9.75）——
   与 `ok16` 同类，用独立真值重裁。
3. **A4/A7/A8/A9/A10/A11–A22**（文档全量对账）—— 纯文档，可一次刷完。
4. **B2/B3/B4/B5**（死代码/死参数/陈旧文件/崩溃路径）。
5. **C** 的 95 例机制分类（原第 303 轮遗留，未开始）。
6. 设计收尾：`_equil_gas_phase` 仍是走步**后**的外部夹子；彻底对称需把相转移
   做成候选网络里的普通反应（`X(aq) ⇌ X(g)`，`logK = −log10(H·p°)`，已验证
   14 种气体精确给出泡点）。**须靠已有走步机制承载，不加特殊路径**。
7. `iters` 30546 → 31326（+2.5%），另立项评估。
8. 更早期积压：EU01 复核；Fe(III)-氯化物 I→0 裁决；Al-氯删除再裁决(10)；
   Sn/Pb 族断言再裁决(6)；池判据回退最小版（见 B2）。

> ✅ **已完成（第 303 轮）**：`T25`/`T104`/`Y12` 阈值重裁 —— 三条旧阈值都
> **高于真值**（引擎其实对）。已改双侧 `has_range` + 撤 `eq`（呈现形态）。
>
> ✅ **已完成（第 305 轮）**：A5 计数双重扣减（辅助 121/122 → 122/123）；
> T9c 闸门判据二（旁观者恒等式）。
>
> ✅ **已完成（第 306 轮）**：A1/`ok16` adjudication（`case_api` 22/22）；
> A2/B1 缺 ΔHf 按 0 计入 + 全量报告；A3 补数据路证伪（可推导仅 2 个）；
> T1 量纲错误根治（第 305 轮）。

---

## 7. 别重复踩（被实测否掉的方案，详见 `lessons.md`）

| 方案 | 否掉的原因 |
|---|---|
| A「账本总量 + 活度截断 `a=min(n/V,c_sat)`」 | 截断破坏化学计量：`H39` 报 pH 13.91 而引擎自己 pH 机给 6.40 |
| B「298 轮 `_sweep_gases`」 | 真删 ledger 物质，抹掉 96% 反应物 ⟹ −510 |
| C「301 轮拆 `X(g)` 第二物种」 | 被迫连加 4 道归并补丁（**分错信号**） |
| D「给 `E01` 造 `sd_temp`（方向+温度）机制」 | 瞄错电对：问题反应里 S 是**氧化剂**；正解是纯数据一行 |

**通用信号**：需要连加 3 道以上补丁才能"看起来不坏" ⟹ 根因找错了，回退重来。
