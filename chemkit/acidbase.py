"""chemkit.acidbase：精确酸碱求解地基（质子化族 + 电荷平衡 pH）。

## 动机

0.4.x 的 pH 由 `speciation.estimate_state` 的**启发式分支机**给出：
强酸直读 / 缓冲滴定 Henderson / 两性中点 / 缓冲对加权 / 弱酸二次式。
三个已知硬伤都与它相关：

  · **J06**（PbCl₂@363K）"酸侧悬崖"：|He| 跨 1e-3 mol/L 时跳分支，
    二分探针内 f 不可微；
  · 多级梯被"最强一级"截断（Z31 硅酸梯只走一级 ⟹ pH 虚高）；
  · 分支 4 取"各来源贡献最大者"，pH 不连续（TS04/H46 的残差行张力）。

本模块给出同一问题的**单一连续表述**，作为 0.5.x 求解器重构的地基：

```
质子化族（pKa 共轭连通块，如 CO₃²⁻/HCO₃⁻/CO₂、NH₃/NH₄⁺）
    · 族总量对任何族内质子转移不变；
    · 任一 pH 下按**逐级 pKa 的严格多级分布**重分配（不截断在最强一级）。
电荷平衡（精确定律，引擎账本恒满足）
    [H⁺] − [OH⁻] = −Σ_{非质子物种} z·c
    ⟹ 方程 F(pH) = Σ_{全部物种} z·c(pH) = 0，在 pH 上严格单调，
       二分必收敛到唯一根。
```

与"质子条件"朴素写法的差别：后者需要 H₂O/H 的原子库存账（水库存
55.6 mol/L 的体相量级会压垮尺度）——0.4.x 的三次试探（Co32/33、
触发点②、pH 钳制）都卡在这里。电荷平衡形式只需要**族总量**（可由
账本一次求和得到）与 pKw，无体相项。

## 用法与边界（诚实记录）

`charge_pH` 是**诊断/估计**工具：它给出的 pH 使**当前账本**（族总量
固定）电荷自洽。若账本本身与某个 pH 不自洽（走步欠收敛留下的混合
状态），`charge_pH` 会给出将其拉回自洽的 pH——这正是可用之处
（0.4.1 的 `_presentation_He` 做的是同类事，但只覆盖强碱分支）。

**它不是完整解**，两类已实测的反例（tools/phdiag.py 全量对照）：
  · **强酸条件 `c_H`**：按"只加质子、不带共轭阴离子"记账，须用 `c_H=`
    参数补回，否则解到碱侧假根；
  · **Ksp/气相储库控 pH**：账本里的溶解量由储库（Cu(OH)₂ 等）决定，
    本身与电荷平衡不自洽（L08 的 Cu²⁺ 0.0005 满足 Ksp 的 pH 5.805，
    却要求 [OH⁻]=0.001）——此时机器的 Ksp 估计比电荷平衡更接近真值。

真正的热力学态需要"pH ⇌ 族分布 ⇌ 走步"三者联立（pH 变化会改变族分布、
族分布改变驱动、驱动改变账本）。0.5.x 的求解器重构以此为目标，本模块
提供其中两块的精确实现。

依赖：core（pKw_of/_vant/charge_of）、candidates（WATER/H_ION）——
与 speciation 同级底层，不反向依赖 engine。
"""
from __future__ import annotations

from math import log10, sqrt

from .core import pKw_of, _vant, charge_of, elements_of
from .candidates import WATER, H_ION

OH_ION = "OH^-"
PH_LO = -1.5                    # pH 搜索域下界（浓酸可到 −1）


def _neg(s: str) -> str:
    """排序取反（用于 max 时选字典序最小者）。"""
    return "".join(chr(0x10FFFF - ord(c)) for c in s)


PH_MARGIN = 2.0                 # 上界 = pKw + 余量


def build_families(T) -> dict:
    """pKa 表 → 质子化族静态表（按 Tables 缓存，`T._ab_families`）。

    返回 {species: (family_root, order, idx, pka, dh)}：
      · `order` 为族内从**最去质子态**到最质子态的排序；
      · `pka[k]` 为 order[k] → order[k+1] 的逐级 pKa；
      · `dh[k]` 为其反应焓（van't Hoff 修正）。

    只用 **n == 1 的逐级条目**建图：n≥2 的合并条目（如 CO₂→CO₃²⁻ 的
    16.7 = 6.4 + 10.3）是同一对端点的冗余边，进图会产生重复边、使
    "最长链"退化（实测曾把 HCO₃⁻→CO₂ 的 pKa 误算成 8.35）。逐级条目
    的 Hess 加和与合并条目一致，故只用逐级条目不丢信息。
    """
    cache = getattr(T, "_ab_families", None)
    if cache is not None:
        return cache
    # 无向邻接（求连通块用）；pKa/焓记在物种对键上（方向无关）
    adj: dict[str, list] = {}
    bond: dict[tuple, tuple] = {}
    # 边的**方向**（acid → base）单独留档：pKa 条目的定义是
    # `A ⇌ B + H⁺`，即 `[B]/[A] = Ka/[H⁺] = 10^(pH − pKa)`。链定序按
    # "最去质子端在前"，多数条目链序与标签序一致（链上 k→k+1 是 base→acid），
    # 但**标记酸恰是少质子一侧**的条目（H₃BO₃/[B(OH)₄]⁻、Tl⁺/TlOH、CO₂/HCO₃⁻…）
    # 会被定序成 acid→base——此时逐级权重比必须取倒数，否则整个族的分布
    # 倒过来（实测硼酸 pH 7 下 [B(OH)₄]⁻/[H₃BO₃] = 174，真值 0.0058；
    # T01/T02 的 TlCl 溶液被算成 pH 1.87）。`bond` 保持原样供 _longest_chain
    # 用（它只需要 pKa 数值）。
    bond_dir: dict[tuple, tuple] = {}
    for e in T.pka:
        if e.get("n", 1) != 1:
            continue
        a, b = e["acid"], e["base"]
        if a == WATER or b == WATER:
            continue
        dh = e["dH"] if "dH" in e else None
        adj.setdefault(a, []).append((b, e["pka"], dh))
        adj.setdefault(b, []).append((a, e["pka"], dh))
        key = (a, b) if a <= b else (b, a)
        bond.setdefault(key, (e["pka"], dh))
        bond_dir.setdefault(key, (a, b))
    seen: set = set()
    out: dict = {}
    for root in sorted(adj):
        if root in seen:
            continue
        block = [root]
        seen.add(root)
        stack = [root]
        while stack:
            cur = stack.pop()
            for nxt, _p, _d in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    block.append(nxt)
                    stack.append(nxt)
        nb: dict[str, list] = {s: [] for s in block}
        for s in block:
            for nxt, pka, dh in adj[s]:
                if nxt in nb:
                    nb[s].append((nxt, pka, dh))
        order = _longest_chain(block, nb, bond)
        pka_from: list[float] = []
        dh_from: list = []
        sgn_from: list[int] = []
        for i in range(len(order) - 1):
            a, b = order[i], order[i + 1]
            key = (a, b) if a <= b else (b, a)
            v = bond.get(key)
            pka_from.append(v[0] if v is not None else 7.0)
            dh_from.append(v[1] if v is not None else None)
            # 逐级权重比的符号：链步 base→acid（与条目标签同向）取 +1
            # （w[acid]/w[base] = 10^(pKa−pH)）；链步 acid→base（标签反向）
            # 取 −1（w[base]/w[acid] = 10^(pH−pKa) = 10^−(pKa−pH)）。
            d = bond_dir.get(key)
            sgn_from.append(1 if (d is not None and d[0] == b) else -1)
        for i, s in enumerate(order):
            out[s] = (min(block), tuple(order), i,
                      tuple(pka_from), tuple(dh_from), tuple(sgn_from))
    T._ab_families = out
    return out


def _longest_chain(block: list, nb: dict, bond: dict | None = None) -> list:
    """族内链定序 → **最去质子端在前**（dist_charge 的基准约定）。

    族由 pKa 共轭对构成，真实数据全为链（H_nA → … → A^{n-}）。定序取
    **最长路**（链上距离最大的点对即两端，与取代基方向无关）：

      · 端点：从任一节点出发两次 BFS 得两个远端；并列时分别取字典序
        最小/最大各探一次，合并候选（多支族取其一，保守）；
      · 方向：**氢原子数少的一端作 order[0]**（最去质子态）；
        氢数相等时（CO₂ vs CO₃²⁻）退回字典序保证跨进程确定。

    教训入档：曾用"沿 pKa 条目的 acid/base 方向链走"定序——`acid` 是
    少一个质子的一侧，但 CO₂ 的 H 数为 0 而 CO₃²⁻ 也为 0，字典序把
    CO₂ 排到 HCO₃⁻ 之后，整条碳酸梯的 pKa 错配一级（HCO₃⁻→CO₂ 被算成
    10.3）。最长路 + 氢计数与取代基方向无关，是链拓扑的稳健判据。
    """
    def _dist(src):
        d = {src: 0}
        st = [src]
        while st:
            cur = st.pop()
            for n, _p, _h in nb[cur]:
                if n not in d:
                    d[n] = d[cur] + 1
                    st.append(n)
        return d

    def _far(src, pick_min: bool):
        d = _dist(src)
        best = max(d.values())
        cands = sorted(s for s in d if d[s] == best)
        return cands[0] if pick_min else cands[-1]

    seed = min(block)
    ends = sorted({_far(seed, True), _far(seed, False),
                   _far(max(block), True), _far(max(block), False)})
    _ = ends                            # 保留供多支族诊断

    def _hcount(s: str) -> int:
        return sum(c for el, c in elements_of(s).items() if el == "H")

    # 链的最去质子端：H 数最少且**邻居 H 数严格更多**的唯一节点。
    # 这一判据对"CO₂ 与 CO₃²⁻ 的 H 数并列最少"的碳酸梯也成立——
    # 只需看邻居：CO₃²⁻ 的邻居 HCO₃⁻ 有 1 个 H，而 CO₂ 的邻居同样有
    # 1 个 H，两者都满足"邻居更多"？不——二元邻居判据需配合"链上
    # 单侧"：用最长路取全链，链上 H 数单调，故唯一极小端是链首。
    # 实现：先取最长路（链拓扑），再按链上 H 数单调性定首端。
    seed = min(block)
    e1 = _far(seed, True)
    e2 = _far(e1, False)
    if e2 == e1:
        return list(sorted(block))
    parent = {e1: None}
    queue = [e1]
    while queue:
        cur = queue.pop(0)
        if cur == e2:
            break
        for n, _p, _h in nb[cur]:
            if n not in parent:
                parent[n] = cur
                queue.append(n)
    path = []
    cur = e2
    while cur is not None:
        path.append(cur)
        cur = parent[cur]
    hs = [_hcount(s) for s in path]
    _ = hs
    # 定向：**逐级 pKa 递减方向即"少质子在前"**。多元酸逐级解离常数
    # 必然递减（K_a1 > K_a2 > …，静电与统计效应），故 pKa 序列单调减的
    # 方向就是最去质子端在前的方向——这是数据自身的化学约束，比氢计数
    # 更可靠（CO₂ 与 CO₃²⁻ 的 H 数并列，氢计数无法定向）。
    rev = list(reversed(path))
    def _pka_seq(chain):
        out = []
        for i in range(len(chain) - 1):
            a, b = chain[i], chain[i + 1]
            k = (a, b) if a <= b else (b, a)
            out.append(bond.get(k, (7.0, None))[0])
        return out

    fwd_seq, rev_seq = _pka_seq(path), _pka_seq(rev)
    fwd_dec = all(fwd_seq[i] >= fwd_seq[i + 1] for i in range(len(fwd_seq) - 1))
    rev_dec = all(rev_seq[i] >= rev_seq[i + 1] for i in range(len(rev_seq) - 1))
    if rev_dec and not fwd_dec:
        path = rev
    elif fwd_dec == rev_dec:
        # 两者皆单调或皆非单调（单级族/异常数据）：退回氢计数 + 名序
        if (tuple(_hcount(s) for s in rev), tuple(rev)) < (
                tuple(_hcount(s) for s in path), tuple(path)):
            path = rev
    pset = set(path)
    # 支链/断链残留：字典序补尾（保证不丢物种；其边 pKa 保守取 7.0）
    return path + [s for s in sorted(block) if s not in pset]


def eff_pka(pka: tuple, dh: tuple, T_K: float) -> tuple:
    """van't Hoff 修正后的逐级有效 pKa（与 speciation._pkapp 同口径）。"""
    if any(d is not None for d in dh):
        return tuple(pka[k] - (_vant(dh[k], T_K) if dh[k] is not None else 0.0)
                     for k in range(len(pka)))
    return pka


def dist_charge(pka: tuple, q: tuple, M: float, V: float, pH: float,
                sgn: tuple | None = None) -> float:
    """族在 pH 下的总电荷（mol，对 V 升）。q[k] = order[k] 的电荷数。

    严格多级分布：以 order[0]（最去质子态）为基准的权重
        w_0 = 1
        w_i = Π_{k<i} 10^(sgn_k·(pKa_k − pH))
    `sgn_k` 由链步与条目标签方向的关系定（见 `build_families`）：同向取 +1
    （w_{k+1}/w_k = [H⁺]/Ka），反向取 −1（w_{k+1}/w_k = Ka/[H⁺]）。
    log 域累加避免上下溢。
    校验：Ac⁻/HAc（order = [Ac⁻, HAc]，pKa 4.76）——pH 7 > pKa 时
    去质子态占优，电荷 → −0.1；NH₄⁺/NH₃（order = [NH₃, NH₄⁺]）——
    pH 7 < pKa 9.25 时质子化态占优，电荷 → +0.1。
    `sgn=None` 等价全 +1（旧口径；仅当调用方确知链序与标签同向时可用）。
    """
    n = len(pka) + 1
    logw = [0.0]
    acc = 0.0
    for k in range(n - 1):
        acc += (pka[k] - pH) * (1 if sgn is None else sgn[k])
        logw.append(acc)
    mx = max(logw)
    sw = 0.0
    sz = 0.0
    for i in range(n):
        wi = 10.0 ** (logw[i] - mx)
        sw += wi
        sz += q[i] * wi
    return M * V * (sz / sw) if V != 1.0 else M * (sz / sw)


def ledger_charge(ledger: dict) -> float:
    """账本净电荷 Σz·n（mol；跳过 `__` 元数据与游离 H⁺/OH⁻ 条目）。

    引擎不变量：`Σz·n + H_excess + c_H = 0`（`c_H` 为 cond 的强酸条件——
    它按"只加质子、不带共轭阴离子"记账，故必须显式补回）。审计发现该
    不变量在 1218 例中大部分成立，偏离即账本自身的守恒缺口。
    """
    return sum(charge_of(s) * m for s, m in ledger.items()
               if not s.startswith("__") and s not in (H_ION, OH_ION))


def charge_pH(ledger: dict, V: float, T, T_K: float,
              lo: float = PH_LO, hi: float | None = None,
              c_H: float = 0.0, tol: float = 0.0, fast: bool = False,
              pinned: tuple = ()):
    """电荷平衡求解 pH；返回 float 或 None（无括号）。

    账本按两轴解读：
      · **质子化族成员**（pKa 连通块）：成员量之和守恒，pH 变化时按逐级
        pKa 的严格分布重分配（多级梯一次到位）；
      · **其他物种**：量固定（强电解质离子、配合物、固相、气相……），
        以电荷数直接进入方程——配合物的酸碱再分布由走步的候选另行处理。
    H⁺/OH⁻ 由 pH 与 pKw 给出，**不来自账本**——账本里若残留 H_ION/OH_ION
    条目它们是被 pH 取代的陈旧副本，计入会与 V·h/V·oh 双重记账。

    `c_H`：cond 的强酸条件量（mol）。`c_H` 按"只加质子、不带共轭阴离子"
    记账 ⟹ 账本净电荷比真实溶液多 `+c_H`，须从固定电荷中扣除，否则
    强酸条件用例（E51/E02/E19/N09/E18/E20/35…）会解出 pH 14+ 的碱侧假根。
    实测：E51（c_H=9）不加此项 → pH 14.96（假）；加 → pH 0.0（真）。

    2026-02 审计作废项：本函数曾有 `H_excess` 形参，把它**加进固定电荷**。
    那是双重记账——`H_excess` 正是 `V·h − V·oh` 应等于的量，加进 fixed
    等于把要求解的未知量减掉。J06 因此报 pH 6.216（真值 3.000），TS04 报
    7.002（真值 1.849）。**账本自身的电荷平衡已足够定 pH，无需 He 输入**；
    需要校验时用 `ledger_charge` 对照。

    边界（诚实记录）：本函数只解**账本内部**的电荷自洽，不重解固相/气相
    储库。若 pH 由 Ksp 储库（Cu(OH)₂、Fe(OH)₃ 等）或气相逸度决定，账本里
    的溶解量本身就与电荷平衡不自洽（L08：Cu²⁺ 0.0005 对 pH 5.805 满足
    Ksp，却要求 [OH⁻]=0.001）——此时本函数给出的是"把账本拉回自洽"的 pH，
    而机器的 Ksp 估计更接近真值。两类反例的清单见 tools/phdiag.py。
    **`pinned` 正是为补这一条而生**（v0.5.0，见下）：把储库自由度写进
    同一个方程，就不必在"账本自洽"与"Ksp 估计"之间二选一了。

    `pinned`：**储库钉住项**。元素 `(z, pKsp, x, y)` 表示"账本里有
    `M_x(OH)_y` 固相在场，自由阳离子 M 的浓度被溶度积钉住"：

        [M] = 10^((y·(pKw − pH) − pKsp)/x)      （Ksp = [M]^x·[OH⁻]^y）

    电荷项 `z·V·[M]` 随 pH **连续**变化（pH 降 → 溶解度升），于是
    "固相在场"从**分支选择**（旧 pH 机器的 Ksp 档：一换档 pH 就跳几个
    单位，§7 X）变成**同一个方程的一项**——固相出现/消失的瞬间，账本里
    的自由离子本就等于饱和浓度，两侧解连续。调用方须先把该阳离子的账本
    条目删掉（钉住值取代它），否则双重记账。

    F(pH) = Σ_{全部物种} z·c(pH) 在 pH 上严格单调减，二分必收敛。

    `tol`：区间宽度收敛阈（0 = 跑到浮点分辨率，90 次上限）。引擎内层
    调用（speciation 的 B3 档）用 `tol=1e-9` 配紧括号——探针是热路径，
    浮点分辨率的最后几十次迭代纯属白跑。

    `fast`：改用**阻尼 Newton**（`_FD` 给出闭式导数：族分布是 pH 的指数族，
    导数 = −ln10·Cov_w(q, c)），典型 4-6 次求值 vs 二分的 ~34 次。用于
    引擎热路径（每判定上万次调用）与全量对照工具。数值上与二分同根
    （同一 F、同一括号），失败时回退二分。
    """
    fams = build_families(T)
    pKw = pKw_of(T_K)
    if hi is None:
        hi = pKw + PH_MARGIN
    fam: dict = {}
    fixed = 0.0
    for s, m in ledger.items():
        if m <= 0.0 or s.startswith("__") or s == WATER:
            continue
        if s == H_ION or s == OH_ION:
            continue          # 游离质子由 pH/pKw 给出，不取账本量
        info = fams.get(s)
        if info is None:
            fixed += charge_of(s) * m
            continue
        fid, order, _idx, pka, dh, sgn = info
        rec = fam.get(fid)
        if rec is None:
            rec = fam[fid] = [eff_pka(pka, dh, T_K),
                              tuple(charge_of(x) for x in order), 0.0, sgn]
        rec[2] += m
    # 强酸条件的无阴离子记账：从固定电荷中扣除（见 docstring）
    if c_H:
        fixed -= c_H
    if not fam:
        if not pinned:
            return _closed_no_family(fixed, V, pKw)

    def _F(pH: float) -> float:
        """F = Σ_{全部物种} z·c（mol）。**随 pH 单调递减**：
        [H⁺] 降、[OH⁻] 升、族平均电荷降（去质子化）、钉住阳离子升
        （溶解随碱升）——故根在 F 由正变负处，二分用 b=mid（F>0 时）收上界。
        """
        h = 10.0 ** (-pH)
        oh = 10.0 ** (pH - pKw)
        tot = V * h - V * oh + fixed
        for ek, q, M, sg in fam.values():
            tot += dist_charge(ek, q, M, V, pH, sg)
        for z, pKsp, xx, y in pinned:
            tot += z * V * 10.0 ** ((y * (pKw - pH) - pKsp) / xx)
        return tot

    def _FD(pH: float) -> tuple[float, float]:
        """(F, dF/dpH)。族分布的 pH 导数有闭式：权重是 pH 的指数族，
        `dlogw_i/dpH = −Σ_{k<i} sgn_k ≡ −c_i`（常数）⟹
        `d⟨q⟩/dpH = −ln10·Cov_w(q, c)`（一次扫描），水项导数 `−ln10·V·(h+oh)`，
        钉住项 `d[M]/dpH = −ln10·(y/x)·[M]`。F 单调递减 ⟹ 导数恒负，
        Newton 步长方向天然正确。
        """
        ln10 = 2.302585092994046
        h = 10.0 ** (-pH)
        oh = 10.0 ** (pH - pKw)
        tot = V * h - V * oh + fixed
        der = -ln10 * V * (h + oh)
        for z, pKsp, xx, y in pinned:
            cm = 10.0 ** ((y * (pKw - pH) - pKsp) / xx)
            tot += z * V * cm
            der -= ln10 * (y / xx) * z * V * cm
        for ek, q, M, sg in fam.values():
            n = len(ek) + 1
            logw = [0.0]
            cs = [0.0]
            acc = 0.0
            cacc = 0.0
            for k in range(n - 1):
                acc += (ek[k] - pH) * sg[k]
                logw.append(acc)
                cacc += sg[k]
                cs.append(cacc)
            mx = max(logw)
            ws = [10.0 ** (v - mx) for v in logw]
            sw = sum(ws)
            mq = sum(q[i] * ws[i] for i in range(n)) / sw
            mc = sum(cs[i] * ws[i] for i in range(n)) / sw
            tot += M * V * mq
            der -= ln10 * M * V * (sum(q[i] * cs[i] * ws[i] for i in range(n))
                                   / sw - mq * mc)
        return tot, der

    if _F(lo) < 0.0:
        return lo                     # 全域为负：下界即解
    if _F(hi) > 0.0:
        return None                   # 全域为正：无根（超强酸/浓碱域外）
    if fast:
        # 阻尼 Newton（导数恒负、F 单调 ⟹ 有根时步长必落在括号内，越界即
        # 折半回退）——典型 4-6 次求值收敛到 1e-12，替代二分的 ~34 次。
        # 引擎若要在热路径用精确质子条件，这一步是前提（性能实测见 §7 O）。
        x = 0.5 * (lo + hi)
        for _ in range(40):
            fx, dx = _FD(x)
            if fx > 0.0:
                lo = x
            else:
                hi = x
            if dx >= 0.0:
                break                 # 导数非负 = 数值异常，交回二分
            nx = x - fx / dx
            if not (lo < nx < hi):
                nx = 0.5 * (x + (hi if fx > 0.0 else lo))
            if abs(nx - x) <= (tol if tol > 0.0 else 1e-12):
                return nx
            x = nx
        # 收敛慢/异常：用已收窄的括号跑二分（此时区间通常已很小）
    a, b = lo, hi
    for _ in range(90):
        if tol > 0.0 and b - a <= tol:
            break
        mid = 0.5 * (a + b)
        if mid <= a or mid >= b:
            break
        if _F(mid) > 0.0:
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


def _closed_no_family(fixed: float, V: float, pKw: float):
    """无质子化族时的闭式解：Σz·c = 0 ⟹ [H⁺] − [OH⁻] = −q。

    q = Σ_{其他物种} z·m/V：强酸给 q < 0（阴离子过剩）⟹ 需游离 H⁺ 补足，
    [H⁺] ≈ −q；强碱给 q > 0 ⟹ [OH⁻] ≈ q。取**量级较小**的根（另一根
    ~|q| 要求阳离子浓度凭空高出投料，非物理）。
    """
    q = fixed / V
    Kw = 10.0 ** (-pKw)
    if -1e-7 < q < 1e-7:
        h = sqrt(Kw)
    elif q < 0.0:
        h = -q + Kw / (-q)
    else:
        h = Kw / (q + Kw / q)
    return -log10(h) if h > 0.0 else None
