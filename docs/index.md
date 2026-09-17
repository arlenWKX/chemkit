# chemkit 文档地图

仓库根 [`README.md`](../README.md) 只做粗略介绍与五分钟上手。

文档分两侧：

- **`docs/`（给人看）**：使用与架构说明。
- **`agents/`（AI 侧工作记忆）**：轮次协议、测量纪律、逐轮工程日志——
  这些不是使用文档，是"动手前先查、别重复付费"的记录。

## `docs/`：使用与架构

| 我想…… | 读 |
|---|---|
| 知道它能不能解决我的问题 | [`../README.md`](../README.md) |
| 查 API：字段含义、参数怎么传、错误行为 | [`api.md`](api.md) |
| 懂引擎怎么算的（公理 / 模块 / 候选宇宙 / 走步 / pH 机器） | [`architecture.md`](architecture.md) |
| 加或改数据（电对、pKa、Ksp、logβ、ΔHf、OVERRIDE） | [`data.md`](data.md) |
| 跑测试、看懂"全绿"意味着什么、发版前做什么 | [`testing.md`](testing.md) |
| 知道性能口径与实测指标 | [`performance.md`](performance.md) |
| 知道哪里不行（诚实清单） | [`limitations.md`](limitations.md) |
| 知道某版本改了什么 | [`changelog.md`](changelog.md) |

## `agents/`：AI 侧工作记忆

| 我想…… | 读 |
|---|---|
| 改数据或做测量（唯一可信路径、测量点自检、逐例单跑） | [`../agents/discipline.md`](../agents/discipline.md) |
| 用审计工具与"一轮"的固定协议 | [`../agents/tools.md`](../agents/tools.md) |
| 查逐轮实验记录、负结果、根因诊断与交接单 | [`../agents/log.md`](../agents/log.md) |

## 维护约定

1. **README 不做细节**：新增细节写进这里的对应文件，README 最多加一行指针。
2. **一个事实只写一处**：数据规模 → `data.md`；性能口径与指标 → `performance.md`；
   基线与发版门槛 → `testing.md`；版本历史 → `changelog.md`。
3. **数字必须可复现**：文档里的规模/指标都从代码或数据本身现采（有对应的探针），
   不凭记忆写；采集方式写在 `tools.md`。
4. **失败实验必须留档**：证伪过的路线（含实测数字）写进 `architecture.md` §7，
   否则后来者会重复付费。
