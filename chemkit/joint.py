"""chemkit.joint：联立平衡求解器（v0.3.8 深水区主菜）。

walk 是逐平衡 Gauss-Seidel 式坐标松弛：每步 pick 单一候选执行到
（其他物种冻结下的）Q=K。多平衡耦合体系（Ag32 银氨-氧化银 876 次
乒乓、J06 PbCl2-氯配-水解链）以 ρ→1 的微步爬行，数百至近两千步；
且微步禁用（disabled）可把仍有驱动的平衡提前锁死——欠收敛出场
（J06 以 Q=0.26·Ksp 出场，解析不动点 ~0.025）。

本模块把耦合平衡集作为联立非线性方程组一次解到位：
    变量   x_j（每活性平衡的净程度，mol；正向=该平衡当前驱动方向）
    状态   n_s(x) = n_s⁰ + Σ_j ν_{j,s}·x_j —— 线性恒等，元素/电荷
           守恒对任意 x 自动成立（与 walk 步同构，无守恒漂移风险）
    方程   F_j(x) = S_j(state(x)) = logK_j − logQ_j(x) = 0
    解法   阻尼 Newton：数值 Jacobian（前向差分）+ 回溯行搜索 +
           非负投影（在场物种不得穿零；H_excess 有符号自由）；
           Jacobian 奇异/欠定（派生候选是基平衡的 Hess 线性组合、
           净键去重后仍可能线性相关）→ Tikhonov 正规方程最小二乘。

安全边界（不触碰 walk 锁定语义，失败=无操作）：
  · 仅收"两侧在场、未冻结、非 slow/deferred、未被膜封锁、当前
    |S| ≥ 0.02（未平衡）"的平衡——冻结/动力学/膜是判定语义的
    既定部分（Hg2_1 50 冻结、Zn+水膜抑制），联立解无权翻案；
  · 解不出（不收敛/行搜索失败/越界）→ 返回 False，walk 原地继续；
  · 解出 → 走步同款记账（_exec：账本/steps/chem_net/origins/hist）
    逐反应落实，disabled.clear() 后由 walk 全量重评估自校验
    （既有自校正语义：S<0 的反向步自动纠回）。

F 求值与 solve_extent 的口径对齐：redox 用 estimate_state 的
滴定后虚拟账本评 S；含 H+ 的非 redox 用 estimate_pH；其余 pH 无关
（传 7.0 常数）。pH/形态在每次试探状态上只算一遍（纯函数复用）。
"""
from __future__ import annotations

import os as _os

from .speciation import estimate_pH, estimate_state

# 联立求解诊断（默认 None，零行为影响）：CHEM_TRACE_JOINT=1 时逐迭代打印
# (resid, λ, x, pH)——"联立为什么不收敛"只能靠迭代史回答（NR57 锚点实测：
# `[joint-ph-fail] m=2 resid=3.063`，需要知道是停滞、行搜索失败还是病态）。
_JTRACE = bool(_os.environ.get("CHEM_TRACE_JOINT"))

# 联立解参数（改动前先跑 converg dump 基线差分）
#
# **`JOINT_TOL` 必须 ≪ `FREEZE_MAX_S`（engine.py，0.5）**（第 297 轮根因修正）：
# 两者是同一根轴上的两个判据——`|S| < JOINT_TOL` = "**已平衡**"、
# `|S| > FREEZE_MAX_S` = "**仍强驱动**"。原值 0.05 与之只差 **10 倍**，中间
# 形成一条"容差带"：联立收敛到带内（`|S| ≈ 0.03`）就返回 `ok` 交回走步，
# 而走步一动 `|S|` 就涨过 0.5 ⟹ `_try_revive` 判定"仍强驱动"⟹ 解冻 ⟹ 再联立
# ⟹ **冻↔解冻活锁**。实测（`15 CaCO3+HCl` / `44 Al2(SO4)3+Na2S`）：
# `joint_ok ≈ freeze_events ≈ revive ≈ 1000`、iters 3002/3053；
# 收紧到 0.01 后 `44` 3053→154、全库 iters **−14.6%**、通过 1280→1283（零翻红）。
# 取 `FREEZE_MAX_S/50` 留足间隔：**联立判"已平衡"的状态，走步绝不能立刻
# 又判"强驱动"**。
JOINT_TOL = float(_os.environ.get("CHEMKIT_JOINT_TOL", "0.01"))

# **联立解的规模出口**（第 299 轮，正式诊断设施，不再是临时插桩）：
# `joint_solve` 返回的 `x` 有多大量级、覆盖几条腿、`status` 是什么，是判断
# "这次联立到底解决了什么"的唯一依据。此前该信息无处可取，导致第 298 轮
# 反复临时插桩且插桩在回退中丢失。`JOINT_AUDIT is None` 时零开销。
JOINT_AUDIT: list | None = None


def _audit(x: list, status: str, resid: float, m: int = 0,
           touched=None, scale: float = 0.0) -> None:
    """记录一次联立的规模（诊断用，生产路径 `JOINT_AUDIT is None` 零成本）。"""
    if JOINT_AUDIT is None:
        return
    JOINT_AUDIT.append({
        "m": m or len(x), "status": status, "resid": resid,
        "max_x": max((abs(v) for v in x), default=0.0),
        "n_nonzero": sum(1 for v in x if abs(v) > 1e-9),
        "scale": scale,
        "touched": sorted(touched) if touched else [],
    })


# 收敛阈：max|F| < 0.01（log 单位；J06 欠收敛 0.59）
# **第 239 轮：两个参数改为可用环境变量覆盖**，用于回答"B 组是欠迭代
# 还是根本不可解"——`joint_solve` 在 B 组停滞态上返回 `fail, resid≈12`，
# 而迭代上限只有 14 次。先量"给足迭代/放宽阈值能否解开"，再决定
# 是扩参数还是需要重写求解器（初值/缩放/阻尼/固相显式化）。
#
# ⚠️ **第 297 轮实测：`JOINT_MAX_ITER` 在本轮这批 churn 例上不参与**
# （14 → 30 → 60 三个档位的 iters **一个数字都没变**）——因为联立在
# `stagn>=3` 或达标时就退出了，**根本没跑满 14 步**。故本次**不动它**；
# 若将来要调，必须先用 `CHEMKIT_JOINT_MAX_ITER` 做全量扫描拿数据。
JOINT_MAX_ITER = int(_os.environ.get("CHEMKIT_JOINT_MAX_ITER", "14"))
JOINT_MAX_M = 12        # 联立维度上限（超出取 |S| 最大者；数值 Jacobian O(m²) 求值）
JOINT_MIN_M = 2         # 维度下限（单平衡 solve_extent 已一步到位）
JOINT_EPS = 1e-9        # 程度显著阈（|x| 全在此下 = 原地，按无操作处理）
# **步长判据**（第 297 轮）：平衡的数学定义是"F = 0 **且 x 不再变**"，
# 只判残差是半个判据。联立解出 x 后，若**每步增量**也已落到可忽略档
# （`X_MIN = 1e-6` mol，与引擎的"可忽略程度"同一把尺子），则解落在
# 平台上，报 `ok` 与报"原地"等价；但若 `|x|` 不小而**本轮 Δx 仍大**，
# 说明还在移动，不该被当成"到位"。
JOINT_DX_MIN = float(_os.environ.get("CHEMKIT_JOINT_DX_MIN", "1e-6"))


def _state(ledger: dict, H_excess: float, actives: list, x: list,
           H_ION: str, WATER: str) -> tuple[dict, float]:
    """试探状态：账本 + He 沿联合方向量 x 线性外推（守恒恒等）。"""
    led = dict(ledger)
    He = H_excess
    for (c, d), xj in zip(actives, x):
        if xj == 0.0:
            continue
        rr = c.r if d > 0 else c.pr
        pp = c.pr if d > 0 else c.r
        for s, nu in rr.items():
            if s == H_ION:
                He -= nu * xj
            elif s != WATER:
                led[s] = led.get(s, 0.0) - nu * xj
        for s, nu in pp.items():
            if s == H_ION:
                He += nu * xj
            elif s != WATER:
                led[s] = led.get(s, 0.0) + nu * xj
    return led, He


def _solve_linear(J: list, b: list) -> list | None:
    """高斯消元（部分主元）解 J·Δ = b；奇异返回 None。"""
    m = len(b)
    A = [row[:] + [b[i]] for i, row in enumerate(J)]
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(A[r][col]))
        if abs(A[piv][col]) < 1e-12:
            return None
        A[col], A[piv] = A[piv], A[col]
        inv = 1.0 / A[col][col]
        for r in range(m):
            if r == col:
                continue
            f = A[r][col] * inv
            if f == 0.0:
                continue
            for cc in range(col, m + 1):
                A[r][cc] -= f * A[col][cc]
    return [A[i][m] / A[i][i] for i in range(m)]


def _solve_ls(J: list, b: list) -> list:
    """Tikhonov 正规方程最小二乘（最小范数解）：(JᵀJ + μI)Δ = Jᵀb。"""
    m = len(b)
    if not m:
        return []
    JtJ = [[sum(J[k][i] * J[k][j] for k in range(m)) for j in range(m)]
           for i in range(m)]
    Jtb = [sum(J[k][i] * b[k] for k in range(m)) for i in range(m)]
    tr = sum(JtJ[i][i] for i in range(m))
    mu = 1e-10 * (tr / m if tr > 0 else 1.0) + 1e-14
    for i in range(m):
        JtJ[i][i] += mu
    return _solve_linear(JtJ, Jtb) or [0.0] * m


def joint_solve(ledger: dict, H_excess: float, actives: list, V: float,
                T_K: float, T, p_ext_kpa: float, gas_escape: bool,
                S_of, H_ION: str, WATER: str) -> tuple[str, list, float]:
    """联立解耦合平衡集。actives = [(Cand, direction), ...]（已去重过滤）。

    返回 (status, x, resid)：
      "ok"       —— 内点不动点达成（resid < JOINT_TOL），x 为各平衡净
                    程度（mol，沿 direction 定向，可为负=净反向）；调用方
                    走步同款记账落实。
      "boundary" —— Newton 沿联合方向推进中某在场物种被钉在零（正价
                    投影顶死）且残差已显著收敛（< 0.5×初值）但仍过阈：
                    联立不动点在物理域外——典型如数据张力型爬行
                    （Ag32 配位-氧化银：R1 的 Hess 组合与配合/溶度数据
                    冲突，无内点解；走步晚期 freeze 是正确仲裁，此处
                    只是提前检测）。x 为钉零时的最远点，仅供参考。
      "fail"     —— 不收敛/行搜索失败/奇异无进展：严格无操作。
    """
    m = len(actives)
    if m < JOINT_MIN_M:
        _audit([], "fail-m0", 0.0, m)
        return "fail", [], 0.0

    def _F(led: dict, He: float) -> list:
        """全部平衡的 S 向量（pH/形态估计每试探状态只算一次）。"""
        pH_x = estimate_pH(led, He, V, T, T_K)
        led_v = None
        out = []
        for c, d in actives:
            if c.kind == "redox":
                if led_v is None:
                    _, led_v, _ = estimate_state(led, He, V, T, T_K)
                out.append(S_of(c, led_v, V, pH_x, T_K, T,
                                p_ext_kpa, gas_escape))
            elif H_ION in c.r or H_ION in c.pr:
                out.append(S_of(c, led, V, pH_x, T_K, T,
                                p_ext_kpa, gas_escape))
            else:
                out.append(S_of(c, led, V, 7.0, T_K, T,
                                p_ext_kpa, gas_escape))
        return out

    # 变化物种清单（非负投影与穿零保险用）
    touched: set = set()
    for c, d in actives:
        for s in list(c.r) + list(c.pr):
            if s not in (H_ION, WATER):
                touched.add(s)
    # 体系物料规模（mol）：诊断用相对基准（取非水非 H⁺ 物种的最大存量，
    # 代表"这个体系里有多少东西可以动"）。
    _scale = max((v for k, v in ledger.items()
                  if k not in (WATER, H_ION) and isinstance(v, (int, float))),
                 default=0.0)

    x = [0.0] * m
    _x_prev = x[:]          # 上一次接受的解向量（步长判据 `_dx` 用）
    led, He = _state(ledger, H_excess, actives, x, H_ION, WATER)
    F = _F(led, He)
    resid0 = max((abs(f) for f in F), default=0.0)
    if resid0 < JOINT_TOL:
        _audit(x, "fail", resid0, m, touched, _scale)
        return "fail", x, resid0   # 已在阈值内（不该发生——驱动者必有 |S|≥0.02）
    resid_best = resid0
    stagn = 0   # 停滞计数：残差 3 迭代无 30% 改善即中止（失败体系的
    # Newton 常在对流/不连续面附近缓慢爬行——早止损，交回走步）

    for _ in range(JOINT_MAX_ITER):
        # 数值 Jacobian：J[i][j] = ∂F_i/∂x_j（前向差分，逐列扰动）
        J = [[0.0] * m for _ in range(m)]
        for j in range(m):
            h = 1e-6 * max(1.0, abs(x[j])) + 1e-9
            x2 = x[:]
            x2[j] += h
            led2, He2 = _state(ledger, H_excess, actives, x2, H_ION, WATER)
            F2 = _F(led2, He2)
            for i in range(m):
                J[i][j] = (F2[i] - F[i]) / h
        delta = _solve_linear(J, [-f for f in F])
        if delta is None:
            delta = _solve_ls(J, [-f for f in F])
        if delta is None or all(dd == 0.0 for dd in delta):
            _audit(x, "fail", max((abs(f) for f in F), default=0.0), m, touched, _scale)
            return "fail", x, max((abs(f) for f in F), default=0.0)

        # 非负投影：在场物种沿 x+λΔ 不得穿零（状态线性 ⇒ λ 上界可解析）
        n0 = led
        n1, _ = _state(ledger, H_excess, actives,
                       [x[j] + delta[j] for j in range(m)], H_ION, WATER)
        lam_cap = 1.0
        for s in touched:
            v0 = n0.get(s, 0.0)
            dv = n1.get(s, v0) - v0
            if dv < -1e-15 and v0 > 0.0:
                lam_s = 0.95 * v0 / (-dv)
                if lam_s < lam_cap:
                    lam_cap = lam_s
        # 回溯行搜索（充分下降 + 投影钳制 + 穿零浮点保险）
        lam = min(1.0, lam_cap)
        F0n = max((abs(f) for f in F), default=0.0)
        ok_step = False
        for _bt in range(8):
            x_new = [x[j] + lam * delta[j] for j in range(m)]
            led_n, He_n = _state(ledger, H_excess, actives, x_new, H_ION, WATER)
            if any(led_n.get(s, 0.0) < -1e-12 for s in touched):
                lam *= 0.5
                continue
            F_n = _F(led_n, He_n)
            Fn = max((abs(f) for f in F_n), default=0.0)
            if Fn < F0n * (1.0 - 1e-4 * lam) or Fn < JOINT_TOL:
                _x_prev = x[:]          # 位移前快照（双判据用，见下方 `_dx`）
                x, F, led = x_new, F_n, led_n
                ok_step = True
                break
            lam *= 0.5
        if not ok_step:
            # 行搜索失败：判边界（物种被投影钉到耗尽 + 残差已半收）——数据
            # 张力型爬行的联立特征（无内点不动点），交调用方提前 freeze。
            # 钉零判据用相对量（λ 钳制下物种按 0.95 几何趋零，非精确为零）
            led_n, _ = _state(ledger, H_excess, actives, x, H_ION, WATER)
            pinned = any(
                ledger.get(s, 0.0) > 1e-4
                and led_n.get(s, 0.0) <= max(1e-9, 0.02 * ledger[s])
                for s in touched)
            resid = max((abs(f) for f in F), default=0.0)
            if pinned and resid < 0.5 * resid0:
                _audit(x, "boundary", resid, m, touched, _scale)
                return "boundary", x, resid
            _audit(x, "fail", resid, m, touched, _scale)
            return "fail", x, resid
        resid = max((abs(f) for f in F), default=0.0)
        if resid < 0.7 * resid_best:
            resid_best = resid
            stagn = 0
        else:
            stagn += 1
            if stagn >= 3:
                _audit(x, "fail", resid, m, touched, _scale)
                return "fail", x, resid
        if resid < JOINT_TOL:
            # **双判据**（第 297 轮）：平衡 = `F = 0` **且** `x` 不再变。
            # `_dx` 取**本次迭代实际接受的位移** `|x_new − x_prev|`
            # （= `lam·delta`，含行搜索的 λ 缩减），落到 `JOINT_DX_MIN`
            # (1e-6 mol，与引擎"可忽略程度"同一把尺子) 以下 ⟹ 解在平台上，
            # 与"原地无操作"等价：返回 `fail` 让调用方走既有"无操作"分支。
            _dx = max((abs(_xn - _xo) for _xn, _xo in zip(x, _x_prev)),
                      default=0.0)
            if max((abs(xj) for xj in x), default=0.0) < JOINT_EPS:
                _audit(x, "fail-eps", resid, m, touched, _scale)
                return "fail", x, resid   # 解 = 原地（无操作）
            if _dx < JOINT_DX_MIN:
                _audit(x, "fail-plateau", resid, m, touched, _scale)
                return "fail", x, resid   # 残差已达标但本次几乎没动：平台点
            # ⚠️ **已试并否证（第 299 轮）：逐腿自洽复核**。
            # 假说：`resid = max|F|` 在多维里不是充分条件，联立靠其它腿
            # "分摊"残差 ⟹ 单腿仍远离平衡。做法是在此重算落实后的各腿 S，
            # 若 `max|F| ≥ JOINT_TOL` 则报 `fail`。
            # **实测否证**：加入后 `15` 仍是 `996 次 ok`、iters 仍 3002
            # ⟹ 落实后的 12 维联合残差**确实 < 0.01**，没有"分摊"。
            # 该假说不再重试。
            _audit(x, "ok", resid, m, touched, _scale)
            return "ok", x, resid
    _audit(x, "fail", max((abs(f) for f in F), default=0.0), m, touched, _scale)
    return "fail", x, max((abs(f) for f in F), default=0.0)


# ========================================================== pH 一致化联立（v0.4.2 迭代 D）

JOINT_PH_PQ = 0.2      # pH 闭合行软容差（pH 单位；Henderson 近似精度量级）


def _solve_ph(ledger: dict, H_excess: float, actives: list, V: float,
              T_K: float, T, p_ext_kpa: float, gas_escape: bool,
              S_of, H_ION: str, WATER: str) -> tuple[str, list, float]:
    """pH 一致化联立（社区级"一步数学"）。

    与 joint_solve 的差异：pH 提升为第 m+1 个联立变量——
      · 平衡行的 H⁺ 活度通道直接用变量 pH（S_of 对 pH 仿射 ±ν·pH，
        快照 pH 机器的分支跳变不再进入平衡行 Jacobian——He 跨 1e-3
        的悬崖失稳在联立域内被平滑化）；
      · 附加自洽闭合行 F_m+1 = estimate_pH(led(x), He(x)) − pH：
        pH 机器（缓冲滴定/Henderson/弱酸二次式）锚定联立 pH——滴定
        语义原样进入联立，联立解的 pH 与走步收敛后的 pH 机器读数
        一致。
    适用域（调用方族闸）：社区质子交换族数 ≤1（单族缓冲/两性——
    Henderson 语义可信）；≥2 族时 pH 机器闭包只是启发式（E35 型
    假不动点一枪跳 0.589 mol 的翻案教训），由族闸硬失败拦截，
    绝不 fall-through 到 legacy joint_solve（对同样行做同样的坏跳）。

    返回 (status, x, resid)：x 为前 m 位平衡净程度（pH 变量在末位，
    调用方不落实）；status 口径与 joint_solve 一致——行搜索失败时
    的 boundary 判定（钉零 + 残差半收）必须保留（Ag32 教训：张力
    族提前冻结仲裁 ~1786→101 步依赖此判定，丢失会白付 Newton 成本
    还放走正确仲裁）。
    """
    m = len(actives)
    if m < JOINT_MIN_M:
        return "fail", [], 0.0
    n = m + 1   # 末位 = pH 变量

    def _F(led: dict, He: float, ph: float) -> list:
        """平衡行（pH 用联立变量）+ pH 自洽闭合行。"""
        led_v = None
        out = []
        for c, d in actives:
            if c.kind == "redox":
                if led_v is None:
                    _, led_v, _ = estimate_state(led, He, V, T, T_K)
                out.append(S_of(c, led_v, V, ph, T_K, T,
                                p_ext_kpa, gas_escape))
            elif H_ION in c.r or H_ION in c.pr:
                out.append(S_of(c, led, V, ph, T_K, T,
                                p_ext_kpa, gas_escape))
            else:
                out.append(S_of(c, led, V, 7.0, T_K, T,
                                p_ext_kpa, gas_escape))
        # ---- pH 闭合行（**保持启发式**，第 220/237 轮两次验证后决定）----
        # `estimate_pH(led, He) − ph`：把联立 pH 锚在 pH 机器的读数上。
        # **曾试改为守恒方程 `charge_pH(led) − ph`**（第 220 轮单测：安全，
        # 8/1147 例变、零通过性翻转；第 237 轮与"停滞期联立"合用再测）：
        # 结果 `resid_p90` 不变、`n(|S|>1)` 反而 58→59、且 **B 组
        # （`F31`/`H45`/`H43`/`T52`）完全未动** ⟹ **换闭合行无收益**，
        # 故**不落盘**（保留启发式；改动写法见 `handoff` §2.12/§2.16）。
        out.append(estimate_pH(led, He, V, T, T_K) - ph)
        return out

    def _resid(F: list) -> float:
        """联合残差（pH 闭合行给软容差——pH 机器的 Henderson 近似
        本身 ~0.1-0.2 pH 精度，硬 0.05 阈会把可解体系误判失败）。"""
        r = max((abs(f) for f in F[:m]), default=0.0)
        r_ph = abs(F[m])
        return max(r, r_ph - JOINT_PH_PQ)

    touched: set = set()
    for c, d in actives:
        for s in list(c.r) + list(c.pr):
            if s not in (H_ION, WATER):
                touched.add(s)

    x = [0.0] * m
    ph0 = estimate_pH(ledger, H_excess, V, T, T_K)
    led, He = _state(ledger, H_excess, actives, x, H_ION, WATER)
    F = _F(led, He, ph0)
    resid0 = _resid(F)
    if resid0 < JOINT_TOL:
        return "fail", x, resid0
    resid_best = resid0
    stagn = 0
    if _JTRACE:
        print(f"  [joint-ph-start] m={m} resid0={resid0:.4g} ph0={ph0:.4g} "
              f"eqs={[getattr(c, 'kind', '?') for c, _d in actives]}")

    for _ in range(JOINT_MAX_ITER):
        # 数值 Jacobian：(m+1)×(m+1)——前 m 列平衡程度（S 行仿射/闭合
        # 行经 pH 机器）、末列 pH（S 行的 H⁺ 活度通道仿射）
        J = [[0.0] * n for _ in range(n)]
        for j in range(n):
            if j < m:
                h = 1e-6 * max(1.0, abs(x[j])) + 1e-9
                x2, ph2 = x[:], ph0
                x2[j] += h
            else:
                h = 1e-6 * max(1.0, abs(ph0)) + 1e-9
                x2, ph2 = x[:], ph0 + h
            led2, He2 = _state(ledger, H_excess, actives, x2, H_ION, WATER)
            F2 = _F(led2, He2, ph2)
            for i in range(n):
                J[i][j] = (F2[i] - F[i]) / h
        delta = _solve_linear(J, [-f for f in F])
        if delta is None:
            delta = _solve_ls(J, [-f for f in F])
        if delta is None or all(dd == 0.0 for dd in delta):
            return "fail", x, _resid(F)

        # 非负投影（pH 无约束；物种沿 x+λΔ 不得穿零）
        n1, _ = _state(ledger, H_excess, actives,
                       [x[j] + delta[j] for j in range(m)], H_ION, WATER)
        lam_cap = 1.0
        for s in touched:
            v0 = led.get(s, 0.0)
            dv = n1.get(s, v0) - v0
            if dv < -1e-15 and v0 > 0.0:
                lam_s = 0.95 * v0 / (-dv)
                if lam_s < lam_cap:
                    lam_cap = lam_s
        lam = min(1.0, lam_cap)
        F0n = _resid(F)
        ok_step = False
        for _bt in range(8):
            x_new = [x[j] + lam * delta[j] for j in range(m)]
            ph_new = ph0 + lam * delta[m]
            led_n, He_n = _state(ledger, H_excess, actives, x_new,
                                 H_ION, WATER)
            if any(led_n.get(s, 0.0) < -1e-12 for s in touched):
                lam *= 0.5
                continue
            F_n = _F(led_n, He_n, ph_new)
            Fn = _resid(F_n)
            if Fn < F0n * (1.0 - 1e-4 * lam) or Fn < JOINT_TOL:
                x, F, led, ph0 = x_new, F_n, led_n, ph_new
                ok_step = True
                break
            lam *= 0.5
        if not ok_step:
            # boundary 判定与 legacy 同口径（Ag32 教训：pH 模式失败必须
            # 能报 boundary，张力族提前冻结仲裁依赖它）
            led_n, _ = _state(ledger, H_excess, actives, x, H_ION, WATER)
            pinned = any(
                ledger.get(s, 0.0) > 1e-4
                and led_n.get(s, 0.0) <= max(1e-9, 0.02 * ledger[s])
                for s in touched)
            resid = _resid(F)
            if pinned and resid < 0.5 * resid0:
                return "boundary", x, resid
            return "fail", x, resid
        resid = _resid(F)
        if _JTRACE:
            print(f"    [joint-ph-it] resid={resid:.4g} lam={lam:.3g} "
                  f"x=[{' '.join(f'{v:+.4g}' for v in x)}] ph={ph0:.4g} "
                  f"stagn={stagn}")
        if resid < 0.7 * resid_best:
            resid_best = resid
            stagn = 0
        else:
            stagn += 1
            if stagn >= 3:
                return "fail", x, resid
        if resid < JOINT_TOL:
            if max((abs(xj) for xj in x), default=0.0) < JOINT_EPS:
                return "fail", x, resid   # 解 = 原地（无操作）
            return "ok", x, resid
    return "fail", x, _resid(F)
