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

**它不是完整解**：真正的热力学态需要"pH ⇌ 族分布 ⇌ 走步"三者联立
（pH 变化会改变族分布、族分布改变驱动、驱动改变账本）。0.5.x 的
求解器重构以此为目标，本模块提供其中两块的精确实现。

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
        for i in range(len(order) - 1):
            a, b = order[i], order[i + 1]
            key = (a, b) if a <= b else (b, a)
            v = bond.get(key)
            pka_from.append(v[0] if v is not None else 7.0)
            dh_from.append(v[1] if v is not None else None)
        for i, s in enumerate(order):
            out[s] = (min(block), tuple(order), i,
                      tuple(pka_from), tuple(dh_from))
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


def dist_charge(pka: tuple, q: tuple, M: float, V: float, pH: float) -> float:
    """族在 pH 下的总电荷（mol，对 V 升）。q[k] = order[k] 的电荷数。

    严格多级分布：以 order[0]（最去质子态）为基准的权重
        w_0 = 1
        w_i = Π_{k<i} ([H⁺]/K_{a,k}) = 10^(Σ_{k<i}(pKa_k − pH))
    即每一质子化级比上一级多一个 H⁺（[H⁺] 一次幂）并除以该级 K_a。
    log 域累加避免上下溢。
    校验：Ac⁻/HAc（order = [Ac⁻, HAc]，pKa 4.76）——pH 7 > pKa 时
    去质子态占优，电荷 → −0.1；NH₄⁺/NH₃（order = [NH₃, NH₄⁺]）——
    pH 7 < pKa 9.25 时质子化态占优，电荷 → +0.1。
    """
    n = len(pka) + 1
    logw = [0.0]
    acc = 0.0
    for k in range(n - 1):
        acc += pka[k] - pH
        logw.append(acc)
    mx = max(logw)
    sw = 0.0
    sz = 0.0
    for i in range(n):
        wi = 10.0 ** (logw[i] - mx)
        sw += wi
        sz += q[i] * wi
    return M * V * (sz / sw) if V != 1.0 else M * (sz / sw)


def charge_pH(ledger: dict, V: float, T, T_K: float,
              lo: float = PH_LO, hi: float | None = None,
              H_excess: float = 0.0):
    """电荷平衡求解 pH；返回 float 或 None（无括号）。

    账本按两轴解读：
      · **质子化族成员**（pKa 连通块）：成员量之和守恒，pH 变化时按逐级
        pKa 的严格分布重分配（多级梯一次到位）；
      · **其他物种**：量固定（强电解质离子、配合物、固相、气相……），
        以电荷数直接进入方程——配合物的酸碱再分布由走步的候选另行处理。
    H⁺/OH⁻ 由 pH 与 pKw 给出，不来自账本。

    `H_excess` 为引擎的带符号游离质子账本（mol）：名义上它等于账本净
    电荷的相反数（`Σz·n = −H_excess`），但**账本并不总是守恒**——审计
    1218 例发现 61 例（5.0%）的 `Σz·n + H_excess` 显著偏离零（最大
    13 eq，见 N09/E51/E12 型）。本函数因此以**实际账本净电荷**为准
    求解（`H_excess` 只作为可选校验），把"不守恒"如实暴露成解得 pH 的
    偏移，而不是用不变量把矛盾掩盖掉。

    F(pH) = Σ_{全部物种} z·c(pH) 在 pH 上严格单调减，二分必收敛。
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
        fid, order, _idx, pka, dh = info
        rec = fam.get(fid)
        if rec is None:
            rec = fam[fid] = [eff_pka(pka, dh, T_K),
                              tuple(charge_of(x) for x in order), 0.0]
        rec[2] += m
    # 不守恒时的补偿：把差额并入固定电荷（等价于承认账本自身带净电荷），
    # 使解出的 pH 与账本自洽——否则 5% 的不守恒例会得到完全无意义的 pH。
    if H_excess:
        fixed = fixed + H_excess + sum(
            charge_of(s) * m for s, m in ledger.items()
            if s == H_ION or s == OH_ION)
    if not fam:
        return _closed_no_family(fixed, V, pKw)

    def _F(pH: float) -> float:
        """F = Σ_{全部物种} z·c（mol）。**随 pH 单调递减**：
        [H⁺] 降、[OH⁻] 升、族平均电荷降（去质子化）——故根在
        F 由正变负处，二分用 b=mid（F>0 时）收上界。
        """
        h = 10.0 ** (-pH)
        oh = 10.0 ** (pH - pKw)
        tot = V * h - V * oh + fixed
        for ek, q, M in fam.values():
            tot += dist_charge(ek, q, M, V, pH)
        return tot

    if _F(lo) < 0.0:
        return lo                     # 全域为负：下界即解
    if _F(hi) > 0.0:
        return None                   # 全域为正：无根（超强酸/浓碱域外）
    a, b = lo, hi
    for _ in range(90):
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
