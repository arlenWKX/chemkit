# 工具与协议（命令手册）

> 规则条文见 [`discipline.md`](discipline.md)，经验索引见 [`lessons.md`](lessons.md)。
> 本文件只回答"**敲什么命令、什么时候敲**"。

## 0. 入口

```powershell
python tools/dev.py guide      # 打印协议速查
python tools/readback.py fails # 读上一次套件的结构化留档（不重跑）
```

## 1. 分层测试协议（**默认 T1，禁止无理由全量**）

| 层 | 适用 | 命令 | 成本 |
|---|---|---|---|
| **T0 不跑** | 只改文档/日志/注释；或作用域可证为空 | — | 0 |
| **T1 定向** | 改动只波及**可枚举**的用例（某族化学、某前缀、某数据条目） | `dev.py suite <前缀…>` + `dev.py case <前缀>` | 秒级 |
| **T2 全量（锚定）** | 引擎**公共热路径/公共判据**改动（`speciation`/`engine`/`templates`/`candidates`） | `dev.py suite` | ~45 s |
| **T3 落盘验收** | 提交前 | T2 + 未锚定约定 + `perf` + `eqcheck` + `hygiene --fix` | ~3 min |
| **T3+ 发版** | 版本号变更/发版 | T3 + `hess_audit --with-templates --T 273.15 363.15` + `data_audit` + `db_matrix` | ~5 min |

**跑 T2 之前先在心里写一句"为什么必须全量"**（作用域真的不可枚举？）；
写不出来就退回 T1。改数据/标准时 T1 是主力：先用 `dev.py case` 看化学，
再用 `dev.py suite <前缀>` 看同族是否连带。

## 2. 结果留档与读取（**一次跑，多次读**）

**每条子命令都自动留档**（`_logged` 装饰器，无需人工管道重定向）：
`logs/<cmd>-<时间戳>.log`（完整输出，含全部失败明细）+ `logs/<cmd>-latest.log`。

| 产物 | 内容 | 读法 |
|---|---|---|
| `logs/suite-latest.json` | 套件逐例结果 + `summary`（cases/checks 分开计数）+ `checks.batteries` | `tools/readback.py fails [--grep X] [--full]`、`case <前缀>`、`stats` |
| `logs/perf-latest.log` / `.tmp_dev_converg.json` | `perf` 确定性指标（逐例） | `tools/readback.py converg --top N` |
| `.tmp_dev_before/after.json` | `snapshot` / `cmp` 的前后对比 | `dev.py cmp A B`（差异为 0 时退出码 0） |
| `converg-baseline.json` | 受控基线（只在落盘后 `perf --write` 刷新） | `dev.py perf` |

**纪律**：不要在生成端截断（打印只留摘要/前 12 条失败）；要看更多就**读留档**，
不要重跑。摘要口径形如 `用例 1176/1176 · 辅助检查 123/123 · 合计 1299/1299 PASS`
——用例与辅助检查**分开计数**（历史坑：电池失败进 `FAILS` 却不进 `RESULTS`，
摘要与留档必然打架）。

## 3. 数据读写：`tools/jsondb.py`（**不要再用文本锚串**）

```powershell
python tools/jsondb.py check   chemkit/data/beta.json
python tools/jsondb.py get     chemkit/data/tests.json 12 ph
python tools/jsondb.py set     chemkit/data/tests.json 12 ph "[3.8, 5.5]"
python tools/jsondb.py normalize chemkit/data/*.json     # 统一口径（幂等）
```
批量改动走**一个脚本 + `JsonDoc` API**（一次解析、多处改、一次落盘；
`case_index("用例名")` 定位 tests.json 条目），别连续调 CLI。
数据文件统一口径：`ensure_ascii=False, indent=1` + 末尾换行。


## 3. `dev.py` 子命令

| 命令 | 作用 | 关键默认行为 |
|---|---|---|
| `suite [前缀...]` | 全量/子集套件 | 打摘要 + 失败明细，**每次运行都写 `.tmp_dev_results.json`** |
| `perf [基线] [--write]` | 确定性指标 vs 基线 | **不覆盖**受控基线（写临时文件）；`--write` 才刷新 |
| `case <前缀> [--all]` | 单例深探 | 步表 + 账本净差 + 两版净方程 + 断言判定 |
| `snapshot OUT [前缀...]` | 全库关键输出快照 | 前后对比用，替代 `git stash` |
| `cmp A B` | 只打差异 | 差异为 0 时退出码 0 |
| `anchor on\|off\|status` | pKw 锚定切换 | 幂等；保留换行风格 |
| `eqcheck` | 全库 `eq`/`eq_has` 精确守恒 | 违规反映到退出码 |
| `hygiene [--fix]` | 行尾噪声、临时文件 | `--fix` 一键还原仅换行差异的文件 |
| `patch spec.py [--check]` | 声明式补丁 | **先全量校验（计数断言）再原子落盘**；spec 里只放 `PATCHES` |
| `run script.py [args]` | 任意脚本在 UTF-8 控制台跑 | 消除"一个 `H⁺` 就 UnicodeEncodeError" |

每个子命令首行都打印 **pKw 约定**——"这轮跑的到底是哪个约定"最容易误读。

## 4. 工具创建规范（**先扩后建**）

1. **先查目录**（本文件第 3/5 节 + `ls tools/`）：能用现有工具/子命令解决的，不新建。
2. **能扩不建**：给现有工具加子命令/参数，优于新写脚本；一次性探针**扩展已有
   `_probe_*.py`**（例如 pH 帧扫描、S(x) 剖面、残差画像各留一个通用入口）。
3. **耐久工具放 `tools/`，一次性探针放仓库根 `_probe_*.py`**（gitignore）。
4. **输出结构化留档**（JSON），人读摘要另打；**不在生成端截断信息**。
5. **不得与 stdlib 同名**（`tools/` 会成为 `sys.path[0]`：`inspect.py` 曾让 `argparse` 崩）。
6. **只读优先**：诊断工具不改状态；要改就走 `dev.py patch`。
7. 新增工具后**在 `tools.md` 与 `agents/tools.md`（本文件）登记一行**，否则下次没人找得到。

## 5. 审计工具定位

`readback.py`（**读留档**：失败清单/单例详情/残差榜）、`db_matrix.py`（**跨表完备性**：
pKa 无 thermo ⟹ 无 van't Hoff；Ksp 阳离子无 `nu=1` OH⁻ β ⟹ 仍走 Kh 复合近似；
化合物固相缺 ksp ⟹ 有氧化通道却无沉淀平衡；输出按**用例触达数**排序，是数据库扩充
路线图）、`data_audit.py`（**表内**自洽：守恒/电子数/重复/冲突）、`eqcheck.py`（守恒）、
`hess_audit.py`（**Hess 自洽，必须带温度**，发版前固定动作）、`hydrate_audit.py`、
`charge_audit.py`、`escape_audit.py`、`titr.py`、`phlog.py`、`osc.py`、`phdiag.py`、
`cliff.py`、`topres.py`、`perf.py`、`snap.py`、`xver.py`、`fragility.py`（±1e-11/1e-9
投料扰动脆弱性）、`tension.py`（张力普查/分类）、`roots.py`（求根审计）、
`cand_audit.py`（候选不变量）、`eq_semantics.py`（changed/reacted/方程语义矩阵）、
`gas_audit.py`（气体标准态对账）、`selfbuf.py`（自缓冲步普查）、
`case.py`/`cands.py`/`extent.py`（单例深探）。

**候选 logK 的参考态纪律**：任何构造 `Cand` 处，`logK` 必须是 **298.15 K 参考值**，
OH⁻→H⁺ 折算交给 `pkw_coeff` 通道在 `logK_T` 里做**一次**；`hess_audit --T 363.15`
是这条纪律唯一的机械化守卫。

## 6. 诊断开关（默认关 = 现行为，逐位不变）

```powershell
$env:CHEM_NO_MICRO_FAST="1"; python tools/dev.py perf   # 微步二分早停（清白）
$env:CHEM_NO_DRAIN="1";      python tools/dev.py perf   # 微步排空（承重：+48% iters）
$env:CHEM_TRACE="1"                                     # pick/冻结/微步轨迹
$env:CHEM_TRACE_WINDOWS="1"; $env:CHEM_TRACE_JOINT="1"
```

**性能借支审计（X-39）判据必须两列**：省了多少（iters/sof）**与**藏了多少
（resid/断言/digest）。反例是"微步排空"——看着像纯性能手段，实测是正确性的一部分。

## 7. 顺序依赖审计（`tools/order_audit.py`）

```bash
python tools/order_audit.py            # normal/reverse/shuffle 三个子进程比对
python tools/order_audit.py --cmp A B  # 退出码 0=无差异
```

指纹 = 逐例 `(pH, degree, changed, reacted, 净方程, 步数, H_excess, 终态物种表)`。
**改动静态层、缓存，或任何按 `T_K`/物种缓存的东西之后跑一次**。
