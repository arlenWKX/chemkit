# 纪律条文（改数据 / 做测量 / 裁决标准 / 发版）

> **规则式、短**。为什么这样规定见 [`lessons.md`](lessons.md)；命令怎么写见
> [`tools.md`](tools.md)。违反以下条文得到的测量**看起来像"该改动无影响"，
> 实为探针自身失效**。

## 一、数据表增改只能走「patch + 新进程」

`dev.py` 模块级已 `import chemkit` ⟹ **`dev.py run` 中脚本内的数据改写
（自写 JSON 或内存 `T.beta.append`）对加载器不可见**。可信路径只有一条：

1. `python tools/dev.py patch <spec.py>`（唯一合法写盘方式，带 JSON 复验）；
2. 用**新进程**测量：`dev.py case <前缀>` / `suite` / `perf`；
3. 临时红树测完必须复原：`git checkout -- <文件>` 并核对 `git status` 干净。

**测量点必须自检**（一行）：确认数据真的进了表（`assert … in {…} for …`）。

## 二、测量纪律

* 用例清单取自 `run_case` 的**返回值**；**错误文本逐例单跑**核对——
  读 `RESULTS[-1]` 一类共享索引会错位（曾据此把口径问题误判为模型退化）。
* 包装探针（`S_of`/`enumerate_candidates`）要留意**缓存与预热顺序**；
  结论须换路径或重复跑确认。
* **一次运行全量结构化留档**后反复读取；不在生成端截断、不靠重跑换信息。

## 三、断言纪律（化学是门槛）

* 标准只锁**化学事实**（文献 pH、总产量、特征步），不锁引擎输出、不锁中间路线。
* **语义位移量（`changed`/`degree`/呈现形态）不作门槛**。
* **标准答案可能有错**：先判"化学事实 vs 某次引擎输出的快照"。
* **数据库缺口与引擎缺陷分开记账**；不许用改松标准或堆数据互相掩盖。
* 边界用例当哨兵（如 `R05b`：0.1 M FeCl₃ pH ∈ [1.5,1.9]，文献 1.7~1.8），
  防止重裁时把窗改松。

## 四、发版/落盘硬门槛

**两套 pKw 约定都必须全绿**（锚定 14.0 / 未锚定 14.0042）：

```
python tools/dev.py anchor off && python tools/dev.py suite && python tools/dev.py anchor on
```

其余固定动作：`dev.py suite`（锚定）、`dev.py eqcheck`（0 违规）、
`hess_audit.py --T 273.15 --T 363.15 --with-templates`（0 条不自洽）、
`data_audit.py` + `db_matrix.py`（数据侧）、`dev.py perf`（与基线逐位对账）、
`hygiene --fix`（行尾/临时文件）。**只有在 T3 层级才跑这一整套**（见 `tools.md`）。

## 五、命名与编号

* `tools/` 下的工具名**不得与 stdlib 模块同名**（`tools/` 会成为 `sys.path[0]`：
  `inspect.py` 曾直接让 `argparse` 崩溃）。
* 一次性探针命名为 `_probe_*.py` 放仓库根（gitignore），**优先扩展已有探针**。
* 测试用例编号：计划统一为**连续纯数字**（去掉字母+数字前缀）并删除重复用例
  ——未完成，见 `log.md` 第 146 轮待办。
