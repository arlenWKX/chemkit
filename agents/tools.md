# 开发工具链与"一轮"的固定协议

> 本文件属于 **AI 侧工作记忆**（`agents/`），不是给使用者看的文档。
> 目的：把每轮都要做的事**固化成命令**，让工具调用不再靠现拼脚本。
> 本文件里的教训全部来自实际踩坑（每条都浪费过时间/Token），照做即可避免。
> 使用者的文档从 [`../docs/index.md`](../docs/index.md) 进入。

## 0. 一条命令的入口

```powershell
python tools/dev.py guide             # 打印下面的协议（忘记时先跑这个）
```

所有子命令（`tools/dev.py` 实现，纯标准库）：

| 命令 | 作用 | 关键默认行为 |
|---|---|---|
| `suite [前缀...]` | 全量/子集套件 | 只打摘要 + 失败明细；1200+ 行 `[ok]` 丢弃 |
| `perf [基线.json] [--write]` | 确定性指标 vs `HEAD` 基线 | **不覆盖**受控基线（写 `.tmp_dev_converg.json`） |
| `case 21 T89 [--all]` | 单例深探 | 步表(kind/logK/S/conv/ext) + 账本净差 + 两版净方程 + 套件断言判定 |
| `snapshot OUT [前缀...]` | 全库关键输出快照（含净方程/步表） | 前后对比用，替代 `git stash` |
| `cmp A B` | 只打差异 | 差异为 0 时退出码 0 |
| `anchor on\|off\|status` | pKw 锚定切换 | 幂等；保留文件换行风格 |
| `eqcheck` | 全库 `eq`/`eq_has` 精确有理守恒 | 违规反映到退出码 |
| `hygiene [--fix]` | 行尾噪声、临时文件 | `--fix` 一键 `git checkout` 还原仅换行差异的文件 |
| `patch spec.py [--check]` | 声明式补丁 | 先全量校验（计数断言）再原子落盘；`--check` 只校验 |
| `run script.py [args]` | 任意脚本在 UTF-8 控制台跑 | 消除"探针里一个 `H⁺` 就 UnicodeEncodeError"这一类失误 |

**每个子命令首行都打印 `pKw` 约定**（`锚定 14.0` / `未锚定 14.0042`）——
"这轮跑的到底是哪个约定"是本项目最容易误读的一件事。

## 1. 一轮的固定动作

```
1) python tools/dev.py suite                 # 起点必须绿，记下 N/N
2) 只做一处根治性修改（优先 `edit` 工具；结构性改动写 spec 走 patch）
3) python tools/dev.py case <受影响用例>      # 化学对账：步表/账本/两版净方程/断言
4) python tools/dev.py perf                  # 确定性指标（iters/sof/resid/digest）vs HEAD
5) python tools/dev.py suite                 # 全量回归
6) python tools/dev.py eqcheck && python tools/dev.py hygiene --fix
7) 文档（[`log.md`](log.md) 新增小节 + [`../docs/changelog.md`](../docs/changelog.md)）
   → 提交推送 → 汇报性能表
```

汇报性能表固定四类量：**墙钟**（mean/p50/p90/max，声明为机器噪声）、
**确定性**（`iters_total`/`sof_total`）、**收敛质量**（`resid_p50/p90/max`、
`n(|S|>1)`）、**判定层**（`N/N PASS`，锚定与未锚定两套约定各一行）。

## 2. 禁令（每条都真踩过）

1. **不用 shell heredoc**（PowerShell 不支持 `<<'PY'`，会整段解析失败）；
   不用内联 `python -c "…多行…"`（引号/Unicode/GBK 三重坑）。→ 写文件再跑，
   或用 `dev.py patch`。
2. **不用 shell 的 `grep`/`find`/`type`**：用 `read`/`grep`/`glob` 工具。
3. **控制台只打 ASCII**：方程一律 `.plain()`（`->`/`<=>`）；TeX 只进 Markdown
   `$$…$$` 与 README/architecture 文档。
4. **长输出必须截断**（`Select-Object -Last N`，或 `dev.py` 自带的上限）。
5. **判据看结构，不看渲染串**：v0.5.0 的 `_max_coef(TeX)` 恒 0 事故
   （§7 X-20）就是"回读渲染串"造成的静默失效——凡判据一律吃 dict/list。
6. **断言是化学事实的门槛**：先裁决化学（用 `case` 看步表/账本），再决定
   改引擎还是改标准。**不许**为了让某个约定通过而把标准改松。
7. **前后对比用 `snapshot`/`cmp`**，不要 `git stash`（工作区常有 5+ 个改动，
   stash 失败或 pop 冲突会丢工作）。
8. **`converg.dump()` 会覆盖受控基线**：只跑 `dev.py perf`（默认写临时文件）。
9. 一次只改一个点；改完立刻 `suite` + `perf`，不要攒着一起测。
9b. **大 JSON 数据表（tests.json/beta.json/couples.json）不要用"索引算术改一个
    可变字符串"**：`s[:i] + new + s[j:]` 这类写法在一处失败后位置全错——本轮实测
    把 beta.json 写成非法 JSON（整包 import 失败，只能 `git checkout` 复原）。
    正确做法：先 `grep`/`read` 拿到**原文与行号**，再用 `edit` 工具按唯一片段改；
    或写"按键名定位 → 断言片段唯一 → 替换 → **先 `json.loads` 复验再落盘**"的脚本。
10. 提交前 `hygiene --fix`：Windows 上 `write` 工具会把 CRLF 文件写成 LF，
    `git status` 显示 `M` 但 `git diff` 为空——这类噪声不该进提交。

## 3. 补丁 spec 的写法（结构性改动）

`spec.py`：

```python
# 每条 = (相对路径, 旧串, 新串[, 期望出现次数])
PATCHES = [
    ("chemkit/equations.py", "old exact text", "new exact text"),
    ("chemkit/core.py", "A", "B", 2),        # 期望出现 2 次
]
```

```powershell
python tools/dev.py patch spec.py --check    # 先校验（不写盘）
python tools/dev.py patch spec.py            # 原子落盘（任一条不满足则整体不写）
```

规则：旧串必须是**唯一**片段（或显式给次数）；`patch` 保留原换行风格
（spec 是 LF、Windows 检出的源文件常是 CRLF——匹配/计数/替换都在**归一化
文本**上做，落盘还原原风格，故不会出现"旧串出现 0 次"这类与内容无关的失败）；
校验失败时**一个文件都不写**，因此可以放心一次改多处。

## 4. 已有审计工具的定位（`dev.py` 之外的深挖）

`eqcheck.py`（守恒）、`hess_audit.py`（**Hess 自洽，必须带温度**：派生候选 +
静态模板，`--T 273.15` / `--T 363.15` 是发版前固定动作，各 ~4 s）、
`hydrate_audit.py`
（相歧义/边界敏感性）、`data_audit.py`、`charge_audit.py`、`escape_audit.py`、
`titr.py`、`phlog.py`、`osc.py`、`phdiag.py`、`cliff.py`、`topres.py`、
`perf.py`、`snap.py`、`xver.py`、`fragility.py`（±1e-11/1e-9 投料扰动脆弱性）、
`tension.py`（`--census` 张力普查 / `--class` 分类）、`roots.py`（求根审计）、
`cand_audit.py`（候选不变量）、`eq_semantics.py`（changed/reacted/方程语义矩阵）、
`gas_audit.py`（**气体活度标准态对账**：`--scan` 给全部产气用例的
`a(引擎) vs a(p/p°)` 与 Δlog，是"标准态统一"改动的爆炸半径清单）、
`case.py` / `cands.py` / `extent.py`（单例深探，`dev.py case` 已覆盖常用部分）、
`selfbuf.py`（**自缓冲步普查**：判据"产物配离子的总吸收容量 ≥ 本步释出 H⁺"，
`He_res ≡ 0` 与步长无关的退化步；全库 **133/1173 例（11.3%）** 命中，
主体是四羟基阴离子族 Al/Pb/Cr/Zn/Ga/Sn/Be —— §7 X-33。纯普查，零行为变更）。

**候选 logK 的参考态纪律（§7 X-31 的机械化守卫）**：任何构造 `Cand` 的地方，
`logK` 必须是 **298.15 K 参考值**，OH⁻→H⁺ 之类的折算一律交给 `pkw_coeff`
通道在 `logK_T` 里做**一次**。把运行温度的 `pKw_of(T_K)` 烘进基值同时挂
`pkw_coeff` 会在 T≠298.15 时折算两次，而 298.15 K 恰好恒等（pKw ≡ 14.0）
⟹ 套件全绿、单温 Hess 审计也全绿。`hess_audit.py --T 363.15` 是这条纪律
唯一的机械化守卫（模板族用 `--templates`，含派生用 `--with-templates`）。

## 3. 性能借支审计（X-39）：开关与两列判据

**判据必须两列**：省了多少（`iters`/`sof`）**与**藏了多少（`resid_*`/断言/`digest`）。
反例是"微步排空"——看起来是纯性能手段，实测是**正确性的一部分**（关掉更慢且更错）。

常驻诊断开关（默认关 = 现行为，逐位不变；置 1 关掉该手段）：

```powershell
$env:CHEM_NO_MICRO_FAST="1"; python tools/dev.py perf   # 微步二分早停（清白：+18.5% sof、零真相变化）
$env:CHEM_NO_DRAIN="1";      python tools/dev.py perf   # 微步排空（承重：+48% iters 且残差变差）
```

**爬行窗口几何外推**（原 `CHEM_NO_CRAWL`）已按同一判据**删除**（第 133 轮）：
逐例 A/B 实测全库 1176 例里只影响 2 例（E35 省 17 迭代、Ni41 反而多 5 迭代且
残差更差），合计 +0.1% 却让轨迹依赖窗口缓存 ⟹ 57 行机器（含状态变量与开关）
一次删净，删除处留了记录。

其余诊断开关：`CHEM_TRACE=1`（pick/冻结/微步逐步轨迹）、`CHEM_TRACE_WINDOWS=1`
（逐 32 步窗口的 drift/turnover 标定）、`CHEM_TRACE_JOINT=1`（联立逐迭代残差）。
判据与实测数字见 [`log.md`](log.md) X-39。

## 4. 顺序依赖审计（`tools/order_audit.py`）

**结果是否依赖用例执行顺序**（静态层/缓存/预热路径）——这类缺陷不会在单例复跑里
显形（同进程内逐例独立、同调用重复稳定），只在跨进程换序时冒出来。历史锚点：
J14 `PbCl₂@363K` 曾冷启动 pH 3.0 / Pb²⁺ 0.015631、预热后 pH 6.21 / 0.027474
（后随 X-31 消失；见 X-39 结案）。

```bash
python tools/order_audit.py            # 自跑 normal/reverse/shuffle 三个子进程并比较
python tools/order_audit.py --dump reverse .tmp_r.json
python tools/order_audit.py --cmp A B  # 退出码 0=无差异、1=有差异
```

指纹 = 逐例 `(pH, degree, changed, reacted, 净方程, 步数, H_excess, 终态物种表)`。
**改动静态层、缓存，或任何按 `T_K`/物种缓存的东西之后跑一次**。
