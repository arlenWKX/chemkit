# 工具与命令手册

> 规则条文见 [`discipline.md`](discipline.md)，经验索引见 [`lessons.md`](lessons.md)。
> 本文件只回答"**敲什么命令、什么时候敲**"。

---

## 1. 唯一入口

```powershell
python -m chemkit.testsuit <verb> [args]
python -m chemkit.testsuit help          # 打印全部子命令
python -m chemkit.testsuit <verb> [args]        # 等价旧写法（薄转发，无独立实现）
```

| verb | 作用 |
|---|---|
| `suite [前缀...] [--jobs=N] [--shard=i/N]` | 跑套件（默认全量+并行） |
| `case <前缀> [--all]` | 单例深探：步表/账本净差/两版净方程/断言判定 |
| `snapshot <OUT.json> [前缀...]` | 全库关键输出快照（前后对比用，替代 `git stash`） |
| `cmp <A.json> <B.json>` | 只打差异（无差异则退出码 0） |
| `anchor on\|off\|status` | pKw 锚定切换（幂等） |
| `eqcheck` | 全库 `eq`/`eq_has` 精确守恒（须 0 违规） |
| `patch <spec.py> [--check]` | 声明式补丁（先全量校验 + 计数断言，再原子落盘 + JSON 复验） |
| `hygiene [--fix]` | 行尾噪声 / 临时文件卫生 |
| `run <脚本.py> [参数...]` | 在 UTF-8 控制台跑任意脚本（免受 GBK 编码错误影响） |

**每个子命令首行都打印当前 pKw 约定**——"这轮跑的到底是哪个约定"最容易误读。

---

## 2. 分层测试协议（**默认 T1，禁止无理由全量**）

| 层 | 适用 | 命令 | 成本 |
|---|---|---|---|
| **T0 不跑** | 只改文档/日志/注释；或作用域可证为空 | — | 0 |
| **T1 定向** | 改动只波及**可枚举**的用例（某族化学、某前缀、某数据条目） | `suite <前缀…>` + `case <前缀>` | 秒级 |
| **T2 全量（锚定）** | 引擎**公共热路径/公共判据**改动（`speciation`/`engine`/`templates`/`candidates`） | `suite` | ~70 s |
| **T3 落盘验收** | 提交前 | T2 + 未锚定约定 + `eqcheck` + `hygiene` | ~3 min |
| **T3+ 发版** | 版本号变更 | T3 + `hess_audit`（带温度）+ `data_audit` + `db_matrix` | ~5 min |

**跑 T2 之前先写一句"为什么必须全量"**（作用域真的不可枚举？）；
写不出来就退回 T1。

---

## 3. 结果留档（**一次跑，多次读**）

每条子命令都自动留档：`logs/<cmd>-<时间戳>.log`（完整输出，含全部失败明细）
+ `logs/<cmd>-latest.log`。全量套件的机读产物是 `logs/suite-latest.json`。

```powershell
python -m chemkit.helper.readback fails [--grep 关键词] [--full]
python -m chemkit.helper.readback case <前缀>
python -m chemkit.helper.readback stats
```

**纪律**：不要在生成端截断（打印只留摘要/前 12 条失败）；要看更多就**读留档**，
不要重跑。摘要口径形如
`用例 1280/1378 · 辅助检查 121/122 · 合计 1402/1501 PASS`
——用例与辅助检查**分开计数**。

**长跑时看实时进度**：`logs/suite-progress.json` 持续刷新（原子写），
`pending` 字段直接给出**还在跑的用例**——卡住时用它定位嫌疑名单。

---

## 4. 并行与分片

* 默认 `jobs = min(4, 物理核)`；`--jobs=N` 覆盖（本机实测 4w 55s / 6w 50s /
  8w 51s ⟹ **8 并不快于 6**，超订只增争用）。
* `--jobs=0` 强制串行（供外层并行分片时片内串行）。
* `--shard=i/N`：把用例切成 N 份只跑第 i 份；各片写
  `logs/suite-shard<i>-of<n>.json`，用 `helper/suite_merge.py` 合并。

---

## 5. 数据读写

```powershell
python -m chemkit.helper.jsondb check   chemkit/data/beta.json
python -m chemkit.helper.jsondb get     chemkit/data/tests.json 12 ph
python -m chemkit.helper.jsondb set     chemkit/data/tests.json 12 ph "[3.8, 5.5]"
python -m chemkit.helper.jsondb normalize chemkit/data/*.json   # 统一口径（幂等）
```

批量改动走**一个脚本 + `JsonDoc` API**（一次解析、多处改、一次落盘），
别连续调 CLI。数据文件统一口径：`ensure_ascii=False, indent=1` + 末尾换行。

---

## 6. 审计工具（`chemkit/helper/`）

| 工具 | 用途 |
|---|---|
| `readback` | **读留档**：失败清单 / 单例详情 / 残差榜 |
| `eqcheck` | 全库方程**精确守恒**（发版前固定动作） |
| `hess_audit` | **Hess 自洽**，**必须带温度**（`--T 273.15 --T 363.15 --with-templates`） |
| `data_audit` | 表内自洽：守恒/电子数/重复/冲突 |
| `db_matrix` | 跨表完备性（按用例触达数排序，是数据库扩充路线图） |
| `perf_diff` | 两份留档的**同用例**性能/通过性差分 |
| `order_audit` | 顺序依赖审计（normal/reverse/shuffle 比对） |
| `specdist` | 单例三层视图并列（账本 / 净差 / 呈现） |
| `gasrole` | 气体角色探针（账本量 / `H(T)` / `c_sat`） |
| `suite_parallel` / `suite_merge` | 并行跑套件 / 合并分片 |
| `quick` / `jsondb` | 按序号跑清单 / JSON 读写 |

---

## 7. 诊断开关（默认关 = 现行为）

```powershell
$env:CHEM_NO_MICRO_FAST="1"; python -m chemkit.testsuit suite   # 微步二分早停
$env:CHEM_NO_DRAIN="1";      python -m chemkit.testsuit suite   # 微步排空
$env:CHEM_TRACE="1"                                             # pick/冻结/微步轨迹
$env:CHEMKIT_PROGRESS_S="0"  python -m chemkit.testsuit suite   # 进度每例都写（调试用）
```

**性能借支审计判据必须两列**：省了多少（iters/墙钟）**与**藏了多少
（resid/断言）。反例是"微步排空"——看着像纯性能手段，实测是正确性的一部分。

---

## 8. 工具创建规范（**先扩后建**）

1. **先查本文件第 6 节**：能用现有工具解决的，不新建。
2. **能扩不建**：给现有工具加参数，优于新写脚本。
3. 耐久工具放 `chemkit/helper/`，一次性探针放仓库根 `_probe_*.py`（gitignore）。
4. **输出结构化留档**（JSON），人读摘要另打；**不在生成端截断信息**。
5. **不得与 stdlib 同名**（会被放进 `sys.path[0]`）。
6. **只读优先**；要改数据就走 `patch`。
7. 新增工具后**在 `tools.md` 第 6 节登记一行**，否则下次没人找得到。
