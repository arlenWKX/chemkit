# START —— 每轮从这里开始

> 一屏读完。改完状态段后开始动手；**收尾时更新本文件**。

---

## 1. 现在是什么状态（第 303 轮末，提交 `cc2063b`）

| 项 | 值 |
|---|---|
| 套件 | **1280 / 1378**（98 例失败） |
| `iters` 合计 | 31326 |
| `eqcheck` | **0 违规** |
| 全量墙钟 | ~70 s（4 workers） |
| pKw 约定 | **锚定**（pKw(298.15) = 14.000000） |
| 工作树 | 干净 |

**两条约定都必须绿**（见 §4）：锚定 1280，未锚定 1280。

---

## 2. 唯一入口

```powershell
python -m chemkit.testsuit <verb> [args]     # 唯一入口
python -m chemkit.testsuit help              # 全部子命令
python -m chemkit.testsuit suite             # 全量
python -m chemkit.testsuit suite J05 NR23    # 按名字前缀
python -m chemkit.testsuit case 15           # 单例深探
python tools/dev.py suite                    # 等价的旧写法（薄转发，无独立实现）
```

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

## 5. 待办（按优先级）

1. **失败用例分类**（98 例）—— 尚未做机制层面的归类。
2. `T25` / `T104` / `Y12` 产量阈值 —— **须先手算真值**再定阈，不可凭空设阈值。
3. 设计收尾：`_equil_gas_phase` 仍是走步**后**的外部夹子；彻底对称需把相转移
   做成候选网络里的普通反应（`X(aq) ⇌ X(g)`，`logK = −log10(H·p°)`，已验证
   14 种气体精确给出泡点）。**须靠已有走步机制承载，不加特殊路径**。
4. iters 30546 → 31326（+2.5%），另立项评估。
5. 更早期积压：EU01 复核；Fe(III)-氯化物 I→0 裁决；Al-氯删除再裁决(10)；
   Sn/Pb 族断言再裁决(6)；池判据回退最小版。

---

## 6. 别重复踩（被实测否掉的方案，详见 `lessons.md`）

| 方案 | 否掉的原因 |
|---|---|
| A「账本总量 + 活度截断 `a=min(n/V,c_sat)`」 | 截断破坏化学计量：`H39` 报 pH 13.91 而引擎自己 pH 机给 6.40 |
| B「298 轮 `_sweep_gases`」 | 真删 ledger 物质，抹掉 96% 反应物 ⟹ −510 |
| C「301 轮拆 `X(g)` 第二物种」 | 被迫连加 4 道归并补丁（**分错信号**） |
| D「给 `E01` 造 `sd_temp`（方向+温度）机制」 | 瞄错电对：问题反应里 S 是**氧化剂**；正解是纯数据一行 |

**通用信号**：需要连加 3 道以上补丁才能"看起来不坏" ⟹ 根因找错了，回退重来。
