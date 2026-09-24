"""chemkit.speciation：pH 机器（教科书近似的质子条件估计器）。

职能：
    · _buffer_titration：游离强酸/强碱被在账弱碱/弱酸储备（含配离子
      beta_pka 储备）按强度顺序吸收的滴定记账（虚拟账本）
    · estimate_state：四分支 pH 估计（直读/滴定 Henderson/两性/缓冲对
      加权/弱酸二次式/金属水解 Kh）+ 滴定后虚拟账本 + 残余 He
    · _full_speciation：多级形态分布（沿质子化梯走到底）
    · _respeciate_strong_acids：分子态强酸 ⇌ 离子的逐轮再平衡
      （HNO3/H2SO4 的浓溶液分子分数曲线）

被 solve_extent 的二分 f 与 judge 主循环消费；依赖 candidates/
normalize/core/data。数值路径与 v0.3.3 bit 级一致（缓存只做查表合并）。
"""
from __future__ import annotations
import heapq
import os
from math import log10, sqrt
from .core import pKw_of, _vant, charge_of, PKW_298
from .data import Tables
from .candidates import (Cand, logK_T, build_derived, WATER, H_ION, X_MIN,
                         ANN_MIN_EXTENT,
                         _STRONG_ACID, _ksp_xy, STRONG_MOLECULAR_ACIDS)
from .normalize import _mol_fraction
from .acidbase import build_families, charge_pH



# ========================================================== 离子强度修正（第 180 轮）

_SIT_A = 0.51          # Debye-Hückel 常数（25 °C, 水）
# 静态表（滴定储备）用的参考 I：储备强度对 I 的敏感度远低于 S_of 的驱动力，
# 且该表按 T_K 缓存、无法逐态取 I，故取 1 M MCl₃ 体系的量级作近似。
_SIT_REF_I = 3.0
# **离子强度层口径**（第 188 轮，使用者裁定"不为迎合标准而舍弃更本质的做法"）：
# 默认 **1 = 启用**：对**所有**平衡施加 SIT —— Debye-Hückel 主项（Δz² 由反应式
# 自动算）＋有 ε 数据时的特定离子作用项。这是**真实介质**的做法；
# 全部平衡口径一致，不存在"氯络合按 I 修正、酸碱按 I→0"的混合口径。
# 设 CHEMKIT_SIT=0 可回到全 I→0 的旧口径（仅供对拍，不是推荐默认）。
SIT_ALL = os.environ.get("CHEMKIT_SIT", "1") not in ("0", "false", "")

# 第 198 轮 pinned 接管的**对拍开关**（诊断用，待该轮结案后拆除）：
# CHEMKIT_PINNED=0 时分支 4 退化区不做"固相钉住"精确解接管。
PINNED_TAKEOVER = os.environ.get("CHEMKIT_PINNED", "1") not in ("0", "false", "")


_Z_CACHE: dict = {}


def _charge_cached(sp: str) -> float:
    """电荷查表（`charge_of` 要解析化学式，热路径上百万次调用）。"""
    z = _Z_CACHE.get(sp)
    if z is None:
        z = _Z_CACHE[sp] = charge_of(sp)
    return z


def ionic_strength(ledger: dict, V: float) -> float:
    """由账本算离子强度 I = ½·Σ cᵢzᵢ²（只计带电、非固相物种）。

    **数值性能（第 191 轮）**：这是 `S_of` 每次求值的必经点，而 `S_of` 在一例里
    被调用十万次量级 ⟹ 两处优化：① 电荷查表（原 `charge_of` 每次解析化学式）；
    ② 小体系早退（离子项 < 3 个时直接求和，省掉字典遍历的开销）。
    """
    s = 0.0
    for sp, m in ledger.items():
        if m <= 0.0 or sp == WATER or sp.startswith("__"):
            continue
        z = _charge_cached(sp)
        if z:
            s += m * z * z

    # **I 量化**（第 191 轮）：0.05 的桶。I 在二分中连续变化会让
    # `logK_T` 的 (T_K, I) 缓存全部失效；量化后同一状态的多次求值命中同一
    # 键，而 0.05 分辨率对 DH/ε 项的影响 < 0.004 log 单位（远小于数据本身
    # 的不确定度）⟹ 纯数值分辨率选择，不改物理模型。
    return round(0.5 * s / V * 20.0) / 20.0


def ionic_strength_reaction(ledger: dict, V: float, species) -> float:
    """**反应相关离子强度** `I_eff`（第 192 轮）：只由**参与该反应的物种**贡献。

    切断"旁观强电解质 ⟹ 总 I ⟹ 本反应 logK ⟹ 本反应推进"的**正反馈环**
    （第 191 轮实测：总 I 口径下 iters ×2.2、resid_p90 ×15、单例 62 s），
    同时保留"介质越浓、活度修正越强"的真实效应——因为反应物种的浓度本身
    就随介质（同离子/配体过量）而变。
    """
    s = 0.0
    for sp in species:
        m = ledger.get(sp, 0.0)
        if m <= 0.0:
            continue
        z = _charge_cached(sp)
        if z:
            s += m * z * z
    if s <= 0.0:
        return 0.0
    return round(0.5 * s / V * 20.0) / 20.0


def sit_fixpoint(ledger: dict, V: float, species, I0: float,
                 iters: int = 3, relax: float = 0.5) -> float:
    """阻尼不动点：I 与"参与物种浓度"互相依赖时解到自洽（≤ iters 次）。

    **不做缓存**（第 193 轮教训）：曾按 `id(ledger)` + 物种做缓存，但**账本对象
    在走步中被原地修改**（同一 id、内容已变）⟹ 缓存返回陈旧值，iters 30276 →
    69686、单例 117 s、残差全面变差。**按对象标识缓存"可变对象"是错的**——
    要么按内容快照做键（成本高），要么把缓存放到"账本确实冻结"的那一层
    （如单次二分的调用方），不能在通用函数里做。
    """
    I = I0
    for _ in range(iters):
        In = ionic_strength_reaction(ledger, V, species)
        if abs(In - I) < 1e-3:
            return round(In * 10.0) / 10.0
        I = I + relax * (In - I)
    # **量化到 0.1**：让 `logK_T` 的 (T_K, I) 缓存能命中（同一账本状态下
    # 多次求值落在同一桶）。0.1 的 I 分辨率对 DH/ε 项的影响 ≤ 0.01 log 单位，
    # 远小于数据不确定度 ⟹ 纯数值分辨率选择。
    return round(I * 10.0) / 10.0


def sit_logK(logK0: float, dz2: float, eps_rxn, I: float) -> float:
    """SIT 修正：logK(I) = logK° + Δz²·A√I/(1+1.5√I) + ε_rxn·I。

    `eps_rxn is None` 时只施加 Debye-Hückel 主项（缺特定相互作用系数）。
    """
    if I <= 0.0:
        return logK0
    r = I ** 0.5
    v = logK0 + dz2 * _SIT_A * r / (1.0 + 1.5 * r)
    if eps_rxn is not None:
        v += eps_rxn * I
    return v


def _beta_dz2(b, T) -> float:
    """M^{n+} + nu·Cl⁻ ⇌ 络合物 的 Δz²（正比于反应式两侧电荷平方差）。"""
    zc = charge_of(b["complex"])
    zn = charge_of(b["center"])
    nu = b.get("nu", 1)
    return zc * zc - zn * zn - nu * 1.0


# ========================================================== pH 估计器（§3.2，教科书近似）

def _buffer_holds(pH: float, a: float, b: float, pKw: float) -> bool:
    """Henderson 缓冲式在该 pH 上是否成立（X-40）。

    **成立前提**：缓冲剂浓度 ≫ 隐含自由质子（或氢氧根）浓度。否则该酸实际
    已基本完全解离，"比值定 pH"的近似失去意义——极端例子：0.005 M 的
    HSCN/SCN⁻ 1:1（pKa = −1.85）Henderson 给 pH −1.85（[H⁺] = 71 M！），
    而真实解是 pH ≈ 2.3。

    **实测病根**（D32）：走步出口 `He = +9e−09 mol`（纳摩尔级游离酸，物理上
    无意义）触发 plateau，Henderson 用账本里**既有的** 4.5e−5/5.0e−3 比值给出
    pH 0.196 ⟹ 隐含 [H⁺] = 0.64 M，超过该对总量两个数量级、也与账本质子账本
    （He ≈ 0）差 2.2 个 pH 单位：**微量残余把 1 mol 级体系的 pH 钉死**。
    拒绝后落分支 4（在滴定后的完整分布上评估）⟹ D32 得 NH₄⁺ 弱酸的 pH 4.6。

    第 1 版曾用"对总量 > 其他酸碱物种总量"，实测过严（打红 H88 醋酸铵双缓冲、
    T34/T25/H77/E35/N39/F21 的多对体系——它们本就有"对 vs 对"的合法竞争）。
    """
    h = 10.0 ** (-pH)
    oh = 10.0 ** (pH - pKw)
    return max(h, oh) <= max(a, b)


def _inert_acid(T) -> frozenset:
    """**惰性强酸储备**：第一级 pKa ≤ 0 的酸物种（按 T 缓存一次）。

    这些物种被 `_buffer_titration` 的酸堆**按设计排除**（见 `_acids_map`
    构建处的注释"强酸已由 He 直读处理"），也因此**不参与 Henderson 帧**。
    它们出现在账本里时，滴定给出的 pH 就不是该账本的解——见 `_buffer_titration`
    顶部的惰性储备闸（第 150 轮 §7 X-44）。
    """
    s = getattr(T, "_inert_acid", None)
    if s is None:
        s = T._inert_acid = frozenset(
            entries[0]["acid"] for entries in T.pka_acid.values()
            if entries[0]["pka"] <= 0)
    return s


def _frame_ok(ledger: dict, pH: float, V: float, T) -> bool:
    """Henderson 帧（§7 X-44）对该账本是否**自洽**：分子态强酸（pKa ≤ 0，
    被酸堆按设计排除）要能存在，pH 必须低到 h ≳ c；否则它早该电离掉。

    Se33 的 0.5 M H₂SeO₄ 配 pH 4.71（h = 0.0031 ≪ c）⟹ 帧自相矛盾，
    真解是直读 0.301；Z03 的 HClO₃ 0.0116 M 配 pH 0.18（h = 0.657 ≫ c）
    ⟹ 分子态本就合理，不得干预（实测按"存在/成量"判会打乱 Z03 全程轨迹）。
    """
    inert = _inert_acid(T)
    if not inert:
        return True
    m = 0.0
    for sp, mm in ledger.items():
        if sp in inert:
            m += mm
    return m <= 0.0 or 10.0 * m <= 10.0 ** (-pH) * V


def _buffer_titration(ledger: dict, H_excess: float, V: float, T, pKw: float,
                      multilevel: bool = False, T_K: float = 298.15,
                      cache: dict | None = None,
                      touch: frozenset | None = None,
                      no_base: frozenset | None = None,
                      no_acid: frozenset | None = None) -> tuple:
    """He>0：强酸被在账弱碱（Kb 大者先）吸收 B+H+→HB；He<0：强碱被在账弱酸
    （Ka 大者先）吸收 HA+OH-→A-+H2O。全吸收 → Henderson 定 pH（返回）；
    残余超过 1e-3 mol/L → None（交回直读分支）；无储备 → None。

    **守恒纪律（0.5.0 重构，architecture §7 L）**：堆只负责"谁先被滴定"的
    强度定序，**量一律从活账本读**，每一步都是增量移动（−take/+take）。
    旧实现在弹出时做绝对写 `ledger2[base] = m − take`，而 m 只是该条目入堆
    时的局部量；多级模式下产物按部分量重新入堆，再次弹出就把账本里同名
    物种的其余量整块抹掉（FeCl₃+Na₂CO₃：碳 3.000000 → 2.975767，pH 3.09
    虚低到 2.1）。族总量是质子化梯的基本不变量，重构后任何时刻可断言
    （tools/titr.py 逐调用对账）。

    pH 同步改为**从滴定后的账本导出**（该对的 ledger2[碱]/ledger2[酸]），
    不再读循环局部量：分布由守恒决定、pH 由分布决定，两者不再共用一个
    变量——只改记账不改 pH 的两次修复尝试都因此失败（L-2 实测表）。

    返回 (pH|None, 残余He, 虚拟账本)。滴定在虚拟账本上真实记账（base→acid
    或 acid→base 转化），全吸收后分支 4 必须用虚拟账本评估残余酸碱性——
    否则强酸恰好中和全部弱碱时，原账本里的弱碱会虚报碱性（pH 11 假象）。
    pH 非 None：缓冲对 Henderson 定 pH；pH=None 且残余≈0：落分支4（用虚拟账本）。

    **惰性储备闸（第 150 轮，§7 X-44）**：账本里存着**没参与滴定的强酸**
    （pKa ≤ 0，被酸堆按设计排除）时，Henderson 帧不成立——它只描述被吸收的
    那一对，对这块强酸储备一无所知。Se33 实测：账本含 0.5 M H₂SeO₄ 而
    He = −1e−12 时，堆跳过 H₂SeO₄，Henderson 拿 SeO₂/HSeO₃⁻ 的比值给出
    pH 4.710；而同一账本的直读解是 0.301（0.5 M 强酸）。跨阈值 4.4 个 pH
    单位的不连续 ⟹ solve_extent 的 f(x) 在 x ≈ 1e−12 处翻号 ⟹ 二分根落到
    微步线以下 ⟹ 通道被判"零推进"禁用（D2 口径不一致的一个新面貌）。
    处置：放弃 Henderson 帧、返回 None 交回直读/分支 4——那两条路本来就
    正确地处理分子态强酸（分支 4 的强酸哨兵给出 h_c = c）。全库实测
    1176/1176 逐位不变（闸只落在 HSCN 类缓冲体系上，那两侧答案本就相同）。

    **自抵消禁令（§7 X-33）**：`no_base`/`no_acid` 是本步**产物**里"能吞下
    本步全部释出质子"的配离子集，由 `solve_extent` 按容量匹配判定后传入。
    判定式：Σ(产物系数 × νH⁺) ≥ 本步净释出 H⁺。命中时这些产物不得充当
    碱储备——否则 x mol 水解产物恰能吸收 x mol 共生 H⁺，He_res ≡ 0 与 x
    无关 ⟹ pH 与步长解耦（R06 实测：整条二分线上 pH 恒 5.000、S 到 50%
    转化才归零）。化学图像：释出的 H⁺ 留在溶液里（电荷平衡 ⟹ h = 水解量）。"""
    if abs(H_excess) < 1e-12:
        return None, H_excess, ledger
    # **惰性储备闸**见 `_frame_ok`（§7 X-44）：两个 Henderson 出口都过一遍。
    _nb = no_base or frozenset()
    _na = no_acid or frozenset()
    # 静态预计算（每数据表一次）：碱储备列表、酸储备列表、beta_pka 配离子储备。
    # 原实现每次调用全表扫描并做 max/min/列表解析与 startswith 过滤，是最大热点
    if getattr(T, "_titr_static", None) is None:
        # pKa 静态量按"每质子"预存，van't Hoff 修正 dH/n 在调用时按 T_K 施加
        # （缓存与温度无关：dH 本身是常数）
        bases = []
        for base, entries in T.pka_base.items():
            if base in T.solids or base == WATER:
                continue
            e1s = [e for e in entries if e["n"] == 1] or entries
            e_max = max(e1s, key=lambda e: e["pka"] / e["n"])
            pka = e_max["pka"] / e_max["n"]
            dH_pp = e_max["dH"] / e_max["n"] if "dH" in e_max else None
            acid = e1s[0]["acid"] if e1s[0]["acid"] != WATER else None
            if acid is not None:
                bases.append((pka, dH_pp, base, acid))
        bases.sort(key=lambda x: -x[0])
        acids = []
        for acid, entries in T.pka_acid.items():
            if acid in T.solids or acid == WATER:
                continue
            e1 = min(entries, key=lambda e: e["pka"])
            if e1["pka"] <= 0:
                continue   # 强酸已由 He 直读处理
            dH_pp = e1["dH"] / e1["n"] if "dH" in e1 else None
            acids.append((e1["pka"], dH_pp, acid, e1["base"]))
        acids.sort(key=lambda x: x[0])
        beta_pka = []
        # **OH 梯配离子的储备必须"逐级（终止于固相）"，不是"一步脱光"**
        # （第 276 轮）。`beta_pka:` 派生条目一律是
        #   `[M(OH)_k] + k·H⁺ -> M^{z+} + k·H₂O`（k = 梯级、z = 中心电荷），
        # 于是滴定器**无法在中途停在氢氧化物固相上**。实测
        # `H45 Na[Al(OH)4]+HCl 半量`（`tools/ph_path.py`）：传入 He = +0.3747，
        # `_buffer_titration` 把 `[Al(OH)₄]⁻ −0.09368 / Al³⁺ +0.09368`
        # —— 0.3747 / 0.09368 = **4.000** ⟹ 按 4 个 H⁺/铝酸根吸收，`He_res = 0`。
        # 而化学上第一步是 **1 个 H⁺**：`[Al(OH)₄]⁻ + H⁺ -> Al(OH)₃(s) + H₂O`
        # ⟹ 0.3747 mol H⁺ 把 0.3747 mol 铝酸根变固相 ⟹ 铝酸根 0.5 / 固相 0.5
        # （与手算电荷平衡 + `K = β₄·Ksp` **逐位一致**）。
        # 错法还把账本留成"0.78 M 铝酸根 + 0.094 M Al³⁺"这种**不可能共存**的
        # 混合态，随后分支 4/pinned 给出 pH 12.185，让那 0.375 mol 游离强酸
        # 在走步眼里**彻底隐形**（引擎呈现 pH 12.19，而账本自洽解是 **0.426**）。
        #
        # 逐级反应与 logK 全部由库内数据推出（`tools/ladder_reserve_census.py`
        # 普查 **11 个**配离子）：
        #   `[M(OH)_k] + (k−z)·H⁺ -> M(OH)_z(s) + (k−z)·H₂O`
        #   logK = (k−z)·pKw + pKsp − logβ_k
        # 校验 `[Al(OH)_4]^-`：1×14 + 33 − 34.5 = **12.5**
        # ⟹ `Al(OH)₃(s) + OH⁻ -> [Al(OH)₄]⁻` 的 logK = 14 − 12.5 = **+1.5**
        # （与该 beta 条目自己的 `calibrated` 注记逐位一致）。
        _ksp_oh = {}
        for _e in T.ksp:
            _c0, _a0 = _e["pair"]
            if _a0 == "OH^-" and _c0 not in _ksp_oh:
                _ksp_oh[_c0] = _e
        _step = {}
        for _b in T.beta:
            if _b["ligand"] != "OH^-" or _b.get("m", 1) != 1:
                continue
            _k = _b.get("nu", 1)
            _z = charge_of(_b["center"])
            _e = _ksp_oh.get(_b["center"])
            if _k <= _z or _e is None:
                continue
            _n = _k - _z
            _lg = _n * PKW_298 + _e["pKsp"] - _b["logb"]
            # dH：逐级式是本文件上方 `ksp_beta`（solid + ligand -> complex）的
            # **逆向**，故取其反号（`_hess_dH(k·b.dH, m·e.dH)` = b.dH − e.dH）。
            # 水电离那一份由 `pkw_coeff = n` 通道承载，不并入 dH。
            _dh = (_e["dH"] - _b["dH"]
                   if ("dH" in _e and "dH" in _b) else None)
            _step[_b["complex"]] = Cand(
                "derived", {_b["complex"]: 1, H_ION: _n},
                {_e["solid"]: 1, WATER: _n}, _lg, float(_n), dH=_dh,
                meta={"src": f"beta_step:{_e['solid']}/{_b['complex']}"})
        for dc in build_derived(T):
            if not dc.meta.get("src", "").startswith("beta_pka:"):
                continue
            nu_h = dc.r.get(H_ION, 0)
            if nu_h <= 0:
                continue
            comps = [s for s in dc.r if s not in (H_ION, WATER)]
            if len(comps) != 1:
                continue
            _rep = _step.get(comps[0])
            if _rep is not None:
                beta_pka.append((comps[0], _rep, float(_rep.r[H_ION])))
            else:
                beta_pka.append((comps[0], dc, nu_h))
        # 堆构建角色表：物种 -> ('b', 共轭酸) 或 ('c', dc, νH+)；
        # 碱储备优先（与原分支顺序一致：先查 bases 再查 complexes）
        _heap_role = {}
        for b, info in {b: (p, d, a) for p, d, b, a in bases}.items():
            _heap_role[b] = ('b', info[2])
        for c, dc, nh in beta_pka:
            if c not in _heap_role:
                _heap_role[c] = ('c', dc, nh)
        T._titr_static = (bases, acids, beta_pka,
                          {b: (p, d, a) for p, d, b, a in bases},
                          {a: (p, d, b) for p, d, a, b in acids},
                          _heap_role)
    _bases, _acids, _beta_pka, _bases_map, _acids_map, _heap_role = T._titr_static

    # 有效 pKa/配离子储备强度只依赖 T_K——按温度缓存，避免每次调用对
    # 全部储备条目重算 _pka_eff/logK_T（_buffer_titration 是最大热点）。
    # b_entries/a_entries 是堆条目的静态模板 {物种: (堆键, 角色)}：
    # 角色与键都只依赖 T_K，逐次调用重建纯属白开销（重构前每次全账扫描
    # 现算 key/role）。
    eff_cache = getattr(T, "_titr_eff", None)
    if eff_cache is None:
        eff_cache = T._titr_eff = {}
    eff = eff_cache.get(T_K)
    if eff is None:
        _beff = {b: _pka_eff(p, d, T_K) for p, d, b, a in _bases}
        _aeff = {a: _pka_eff(p, d, T_K) for p, d, a, b in _acids}
        # Cl⁻ 系派生条目带 meta['eps'] ⟹ 需 I；本表按 T_K 静态缓存，
        # 无法逐态取 I，故取参考 I（见 _SIT_REF_I）。
        _ceff = {c: logK_T(dc, T_K, _SIT_REF_I
                           if dc.meta.get("eps") is not None else 0.0) / nh
                 for c, dc, nh in _beta_pka}
        _cmap = {c: (dc, nh) for c, dc, nh in _beta_pka}
        # 碱分支：pKa(共轭酸) 越大 Kb 越大 → 堆键取负；配离子按 νH+ 折算容量
        b_entries = {}
        for sp, role in _heap_role.items():
            if role[0] == 'b':
                b_entries[sp] = (-_beff[sp], role[1])
            else:
                b_entries[sp] = (-_ceff[sp], ("__complex__", role[1], role[2]))
        a_entries = {sp: (_aeff[sp], info[2]) for sp, info in _acids_map.items()}
        eff = (_beff, _aeff, _ceff, _cmap, b_entries, a_entries)
        eff_cache[T_K] = eff
    _beff, _aeff, _ceff, _cmap, _b_entries, _a_entries = eff
    # 延迟拷贝：仅当 heap 非空、确实需要修改账本时才 dict(ledger)。
    # heap 为空（强酸/强碱+盐等无弱组分场景）直接返回原账本——judge 中
    # `vled is not ledger` 身份检查据此跳过 _virt_redox_gain（无滴定=无虚拟增益）
    ledger2: dict | None = None
    he = H_excess
    # 缓冲对的下限：两侧都必须在**痕量线以上**才算"定 pH 的对"。产物侧
    # 只查 > 0 会被收尾的碎屑弹出骗到——滴定恰好吸收完时 he 残 1e-17，
    # 最后一个储备条目仍满足 take < avail，若拿它的微量产物算 Henderson，
    # pH 直接崩到 −1（实测 AB05：Na2HPO4 区被算成 pH −1.0）。碎屑对必须
    # 交给分支 4（在滴定后的完整分布上评估）。
    _floor = max(X_MIN, 1e-9 * V)
    # 多级（多元）滴定用堆：转化产物若本身仍可继续质子化/去质子化
    # （HVO3→VO2+、H2PO4-→HPO4^2-…）则按自身 pKa 重新入堆，一次调用沿
    # 质子化梯走到底——快照式单级实现会把中间形态（如 HVO3）滞留为假象，
    # 后续氧化还原竞争因此看不到真实自由形态（VO2+）
    if he > 0:   # 弱碱吸收：pKa(共轭酸) 越大 Kb 越大，先中和
        # 遍历在账物种（通常 ~15 个）而非全储备表（~50 条），热点降载；
        # 配离子作为碱储备（beta_pka 派生：complex + νH+ -> center + ν共轭酸）：
        # 沉淀等步骤释放的 H+ 实际由配离子解离吸收（如 [Cu(NH3)4]2+），
        # 不纳入会使 solve_extent 内部 pH 崩塌、反应假停滞（Cu2+ + 少量氨水）
        heap, cnt = _titration_heap(ledger, _b_entries, cache, "b", touch)
        if _nb:
            heap = [e for e in heap
                    if not (e[2] in _nb and isinstance(e[3], tuple))]
            heapq.heapify(heap)
        if not heap:
            return None, he, ledger   # 无弱碱储备：直接返回原账本（避免无谓拷贝）
        ledger2 = dict(ledger)
        heappush = heapq.heappush
        plateau = None
        # 多级模式下产物（共轭酸）可能仍可继续质子化，按自身强度重新入堆；
        # 同物种只留一个条目（量从活账本读，重复条目只是白弹堆）
        pushed = {e[2] for e in heap} if multilevel else None
        while heap and he > 0.0:
            neg_pka, _cnt, base, acid = heapq.heappop(heap)
            avail = ledger2.get(base, 0.0)
            if avail <= 0.0:
                continue
            if isinstance(acid, tuple):
                # 配离子储备：按派生方程转化，不提供 Henderson 对、不再入堆
                _, dc, nu_h = acid
                take = min(he, avail * nu_h)
                if take <= 0.0:
                    continue
                he -= take
                dx = take / nu_h
                ledger2[base] = avail - dx
                for sp2, nu2 in dc.pr.items():
                    if sp2 != WATER:
                        ledger2[sp2] = ledger2.get(sp2, 0.0) + nu2 * dx
                continue
            take = min(he, avail)
            he -= take
            ledger2[base] = avail - take
            ledger2[acid] = ledger2.get(acid, 0.0) + take
            if take < avail:
                # 未滴完 ⟹ he 已归零，本对即"定 pH 的缓冲对"（部分滴定至多
                # 一次，必为最后一次），留到循环后按账本量算 Henderson
                plateau = (-neg_pka, base, acid)
            elif multilevel and take > 0.0 and acid not in pushed:
                nxt = _bases_map.get(acid)
                if nxt is not None and acid not in _nb:
                    pushed.add(acid)
                    heappush(heap, (-_beff[acid], cnt, acid, nxt[2]))
                    cnt += 1
        if plateau is not None:
            pka, base, acid = plateau
            b_rest = ledger2.get(base, 0.0)
            hb = ledger2.get(acid, 0.0)
            if b_rest > _floor and hb > _floor:
                pH = pka + log10(b_rest / hb)
                # X-40：Henderson 只在"缓冲剂 ≫ 隐含 [H⁺]/[OH⁻]"时成立
                if _buffer_holds(pH, b_rest, hb, pKw) and _frame_ok(
                        ledger, pH, V, T):
                    return min(max(pH, -1.0), pKw + 1.0), he, ledger2
        return None, he, ledger2   # 全吸收（he≈0）→ 分支4；有残余 → 直读
    else:        # 弱酸吸收强碱：pKa 越小 Ka 越大，先中和
        he = -he
        heap, cnt = _titration_heap(ledger, _a_entries, cache, "a", touch,
                                    nominal=pKw + 2)
        if _na:
            heap = [e for e in heap if e[2] not in _na]
            heapq.heapify(heap)
        if not heap:
            return None, -he, ledger
        ledger2 = dict(ledger)
        heappush = heapq.heappush
        plateau = None
        last_full = None      # 最后一个被定量喂满的弱酸对（平衡修正用）
        pushed = {e[2] for e in heap} if multilevel else None
        while heap and he > 0.0:
            pka, _cnt, acid, base = heapq.heappop(heap)
            avail = ledger2.get(acid, 0.0)
            if avail <= 0.0:
                continue
            take = min(he, avail)
            he -= take
            ledger2[acid] = avail - take
            ledger2[base] = ledger2.get(base, 0.0) + take
            if take < avail:
                plateau = (pka, acid, base)
            else:
                if take > 0.0:
                    # **最后一个被定量喂满的对**：碱过量时它是回水解的主体，
                    # 平衡修正（见下）要用它
                    last_full = (pka, acid, base)
                if multilevel and take > 0.0 and base not in pushed:
                    # 产物仍是酸（可再去质子化）→ 重新入堆（仅多级模式）
                    nxt = _acids_map.get(base)
                    if nxt is not None and base not in _na:
                        pushed.add(base)
                        heappush(heap, (_aeff[base], cnt, base, nxt[2]))
                        cnt += 1
        if plateau is not None:
            pka, acid, base = plateau
            a_rest = ledger2.get(acid, 0.0)
            b = ledger2.get(base, 0.0)
            if a_rest > _floor and b > _floor:
                pH = pka + log10(b / a_rest)
                # X-40：同弱碱吸收分支（对称）
                if _buffer_holds(pH, a_rest, b, pKw) and _frame_ok(
                        ledger, pH, V, T):
                    return min(max(pH, -1.0), pKw + 1.0), -he, ledger2
        if he > 0.0 and last_full is not None:
            # **平衡分布修正**（第 143 轮，§7 X-8 第一块）——见文件头。
            # 强碱把最后一个弱酸**定量喂满**后仍有残余游离碱时，定量模型把
            # 该对整块记成共轭碱（HA ≡ 0）并把残余碱记成 B − T；而共轭碱会
            # 回水解（A⁻ + H₂O ⇌ HA + OH⁻，Kb = Kw/Ka），Ka 小者可达百分之
            # 几十。这会把**反应物形态整个抹掉**（D47：H₂O₂ 0.0888 → 0），
            # 按分子态配平的 redox 通道于是在虚拟账本上"反应物不存在"、
            # x* = 0 停摆——而走步 pick 用真实账本判它强驱动（S = +10.24）。
            # 荷衡 [A⁻] + [OH⁻] = B（B = 该对共轭碱 + 残余游离碱）、
            # 物料 [HA] + [A⁻] = T、[HA] = [A⁻]·h/Ka
            # ⟹ f = B − T·Ka·f/(Ka·f + Kw)，右端在 [0, B] 上单调，二分即得。
            _pka_f, _acid_f, _base_f = last_full
            _T = ledger2.get(_acid_f, 0.0) + ledger2.get(_base_f, 0.0)
            if _T > _floor:
                _Ka = 10.0 ** (-_pka_f)
                _Kw = 10.0 ** (-pKw)
                _Tc = _T / V
                _B = (he + _T) / V
                _lo, _hi = 0.0, _B
                for _ in range(60):
                    _mid = 0.5 * (_lo + _hi)
                    _a = _Tc * _Ka * _mid / (_Ka * _mid + _Kw)
                    if _B - _a - _mid > 0.0:
                        _lo = _mid
                    else:
                        _hi = _mid
                _f = 0.5 * (_lo + _hi)
                _a = _Tc * _Ka * _f / (_Ka * _f + _Kw)
                ledger2[_base_f] = _a * V
                ledger2[_acid_f] = (_Tc - _a) * V
                he = _f * V
        return None, -he, ledger2


def _titration_heap(ledger: dict, entries: dict, cache: dict | None, ckey: str,
                    touch: frozenset | None,
                    nominal: float | None = None) -> tuple:
    """滴定强度堆 → `(heap, cnt)`；条目 `(堆键, cnt, 物种, 角色)`。

    **量不入条目**：滴定量从活账本读（守恒纪律，见 `_buffer_titration`）。
    堆只回答"谁先被滴定"——`entries` 是 {物种: (键, 角色)} 的静态模板
    （每 T_K 构建一次），键序即强度序。

    `nominal`：酸分支的"名义酸"阈值——pKa_eff > pKw+2 的物种水溶液中不可能
    给出质子（NH3 pKa≈105），不入堆。

    缓存协议（solve_extent 二分内复用，与逐次重建等价）：
      · 击中 = 账本键序一致 + touch 物种的在场性（> X_MIN）未翻转；
      · 命中按预排序表 heapify 重建——弹堆序由 (键, cnt) 唯一决定
        （cnt 全局唯一 ⟹ 无同键并列），与逐次 push 完全一致；
      · touch 物种跨阈值翻转即作废全量重建：生成型物种从无到有，漏掉它
        等于漏掉一个滴定储备；cnt 恒为首次构建序，两侧等价。"""
    ent = cache.get(ckey) if cache is not None else None
    if ent is not None and ent[0] == tuple(ledger) and all(
            (ledger.get(sp, 0.0) > X_MIN) == ent[1][sp] for sp in ent[2]):
        heap = list(ent[3])
        heapq.heapify(heap)
        return heap, ent[4]
    rows = []
    presence: dict = {}
    cnt = 0
    for sp, m in ledger.items():
        e = entries.get(sp)
        if e is None:
            continue
        presence[sp] = pres = m > X_MIN
        if nominal is not None and e[0] > nominal:
            continue
        if pres:
            rows.append((e[0], cnt, sp, e[1]))
            cnt += 1
    if cache is not None:
        cache[ckey] = (tuple(ledger), presence,
                       tuple(sp for sp in (touch or ()) if sp in presence),
                       tuple(rows), cnt)
    # 两条路径都必须返回**真堆**：弹堆序由 (键, cnt) 唯一决定，与旧实现
    # 逐次 heappush 完全一致；漏了 heapify 就退化成账本序——强度序失效，
    # 滴定会先动弱储备（实测 F35：H2PO4- 抢在 H3PO4 前被滴定，pH 4.65→1.23，
    # 走步据此把 H2PO4- + H+ → H3PO4 跑到 0.5 mol 的幻影碱态）
    heapq.heapify(rows)
    return rows, cnt


def complex_capacity(T) -> dict:
    """beta_pka 配离子储备的**每摩尔吸收质子数** {complex: νH⁺}（每表一次）。

    §7 X-33 用它判定"本步是不是自缓冲步"：候选产物里 β_pka 配离子的总吸收
    容量 ≥ 本步净释出的 H⁺ 时，滴定会把这批质子原地退回给产物自己
    ⟹ He_res ≡ 0 与步长无关、pH 与步长解耦（J14 铅酸根、X-32 铝酸根、
    EU01 锌酸根同机制）。当前仅 `tools/selfbuf.py` 普查规模用（生产路径不调用）。
    """
    cap = getattr(T, "_cx_cap", None)
    if cap is None:
        cap = {}
        for dc in build_derived(T):
            if not dc.meta.get("src", "").startswith("beta_pka:"):
                continue
            nu_h = dc.r.get(H_ION, 0)
            if nu_h <= 0:
                continue
            comps = [s for s in dc.r if s not in (H_ION, WATER)]
            if len(comps) == 1:
                cap[comps[0]] = nu_h
        T._cx_cap = cap
    return cap


def estimate_pH(ledger: dict, H_excess: float, V: float, T, T_K: float,
                cache: dict | None = None,
                touch: frozenset | None = None,
                no_base: frozenset | None = None,
                no_acid: frozenset | None = None,
                pin_mode: bool | None = None) -> float:
    return estimate_state(ledger, H_excess, V, T, T_K, cache, touch,
                          no_base, no_acid, pin_mode)[0]


# 分支标注（**仅审计**：`tools/roots.py` 之外恒为 None，生产路径零成本）。
# §7 X 的 pH 跳变要归因到"换的哪一条分支"，否则只能看到跳变本身。
PH_TAGS: list | None = None

# 分支 4 的 h_c/o_c **来源**（**仅审计**，同 PH_TAGS：生产路径恒 None）。
# 记录 (物种, 贡献类型, 该物种给出的自由度浓度)——pH 冻结在 pKw/2 附近时，
# 必须先看清"赢下 max 的是谁"，否则只能对着 6.2 猜（§7 X-31）。
PH_SRC: list | None = None


def _tag(why: str) -> None:
    if PH_TAGS is not None:
        PH_TAGS.append(why)


def _role_union(T) -> frozenset:
    """分支 4 可能给出酸/碱来源的物种**超集**（按 T 缓存一次）。

    超集是安全方向：命中时只是多跑一次分支 4 扫描；漏判才会让直读早退
    在"账本里还有弱组分"时错误地独占 pH（§7 X-43）。来源与分支 4 的
    四张表同源：pKa 酸侧 / pKa 碱侧 / Ksp-OH 阳离子 / `nu=1` OH⁻ β 的中心。
    """
    u = getattr(T, "_role_union", None)
    if u is None:
        cats = {e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-"}
        first_k = {b["center"] for b in T.beta
                   if b["ligand"] == "OH^-" and b.get("nu") == 1
                   and b.get("m", 1) == 1}
        u = T._role_union = frozenset(
            set(T.pka_acid) | set(T.pka_base) | cats | first_k)
    return u


def _pin_ladders(T) -> dict:
    """各阳离子的**羟合梯**（按 T 缓存一次）：cat -> ((z, logβ, nu, 物种名), …)。

    钉住固相储库时，该金属的全部 OH⁻ 配合物浓度从属于钉住的游离浓度
    （`[M(OH)_nu] = β·[M]·[OH⁻]^nu`），不能冻结在账本量——第 198 轮实测：
    只删游离 Al³⁺、冻结 0.75 mol 铝酸根，电荷平衡靠"凭空造 0.25 M 游离
    Al³⁺"配平出化学上不可能的 pH 3.2（账本固相仅 1e-4 mol 量级）。
    logβ 取 298 K 值（与分支 4 的 `_first_k` 同一口径）。
    """
    cache = getattr(T, "_pin_ladders", None)
    if cache is None:
        cache = {}
        for b in T.beta:
            if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
                continue
            cache.setdefault(b["center"], []).append(
                (charge_of(b["complex"]), b["logb"], b["nu"], b["complex"]))
        T._pin_ladders = cache
    return cache


def _pin_frozen(T) -> dict:
    """各阳离子的**非 OH⁻ 配合物**（按 T 缓存一次）：
    cat -> ((z, logβ, nu, 配体名, 物种名), …)，仅 m==1 单体配合物。

    钉住固相储库时这些配合物同样从属于钉住的游离 [M]
    （`[ML_nu] = β·[M]·cL^nu`，cL 冻结取**配体返还后**的账本值）——只从属
    羟合梯而冻结它们是同一病根的另一半：第 199 轮实测 case 16 账本
    [AlCl]²⁺ 0.22 mol 被冻结成 +0.67 幻影正电荷，pinned 解抬到 pH 12.33；
    N22 的 Zn-Cl 配合物 + 被绑 Cl⁻ 未返还 ⟹ pH 13.26 把锌酸根稳定住，
    沉淀停在 0.18 < 0.7。条目从账本删除时配体按 ν **返还账本**
    （配体守恒：金属由储库吸收，配体回溶液；与 OH⁻ 不同——OH⁻ 来自水储库）。
    """
    cache = getattr(T, "_pin_frozen", None)
    if cache is None:
        cache = {}
        for b in T.beta:
            if b["ligand"] == "OH^-" or b.get("m", 1) != 1:
                continue
            cache.setdefault(b["center"], []).append(
                (charge_of(b["complex"]), b["logb"], b["nu"],
                 b["ligand"], b["complex"]))
        T._pin_frozen = cache
    return cache


def _exact_ok(ledger: dict, He_res: float, T, V: float) -> bool:
    """`charge_pH` 在**本账本**上是否可证可信（第 263 轮抽出，原为 `_deg`
    路径里的内联块——抽出后两处调用点共用同一把尺子，不再各写一遍）。

    三条同时成立才可信：

      ① **无固相在场** —— `charge_pH` 只解账本内部的电荷自洽，不解 Ksp
         储库（见其 docstring 的边界声明）；有固相时账本溶解量本身就可能
         与电荷平衡不自洽（L08 型）。⚠️ 第 263 轮实测 N15：**光靠这条不够**
         —— 走步的探针会在"固相尚未析出"的状态上问 pH，此时闸放行，
         但轨迹已被改变（详见 §5「家族再分配抵消走步」）。
      ② **至少 2 个质子化族成员** —— 族再分配才有对象。注：
         `build_families` 以**物种**为键，故同一族的两个成员（如 H₂S/HS⁻）
         就满足此条；"只有一个成员"（纯两性盐）不满足。
      ③ **账本电荷与残余强酸记账一致**（`|Σz·n + He_res| ≤ 1e-6`）——
         这条是"账本自洽"的直接判据，也是 `charge_pH` 结论可用的前提：
         该函数把 `V·h − V·oh` 与账本电荷配平，若账本自身就不满足
         `Σz·n = −He_res`，它解出的是"把账本拉回自洽"的 pH 而非真值。
    """
    _els = frozenset(T.solids)
    for sp, m in ledger.items():
        if m > X_MIN and sp != WATER and sp in _els:
            return False
    _fam = frozenset(build_families(T))
    if sum(1 for sp, m in ledger.items() if m > 0.0 and sp in _fam) < 2:
        return False
    _net = 0.0
    for sp, m in ledger.items():
        if m > 0.0 and sp != WATER and not sp.startswith("__"):
            _net += charge_of(sp) * m
    return abs(_net + He_res) <= 1e-6


def estimate_state(ledger: dict, H_excess: float, V: float, T, T_K: float,
                   cache: dict | None = None,
                   touch: frozenset | None = None,
                   no_base: frozenset | None = None,
                   no_acid: frozenset | None = None,
                   pin_mode: bool | None = None) -> tuple[float, dict, float]:
    """返回 (pH, 滴定后的虚拟账本, 残余He)。虚拟账本是 pH 一致的自由形态分布；
    残余He是弱酸/弱碱储备吸收后仍未中和的游离强酸/强碱（酸碱平衡后的真实 He）。

    pin_mode：pinned 口径的**整步冻结**开关（D14，第 201 轮）——None（默认）
    = 现状逐点自判定（h_c≈o_c 退化区闸）；False = 本步冻结"弃用"（pinned 块
    整体跳过）；True = 本步冻结"采用"（绕过退化区闸仍尝试 pinned，块内结构性
    前提仍逐点评估）。逐点判定会在二分探针上制造 pH 接缝悬崖被当成根。"""
    pKw = pKw_of(T_K)    # 1) 强酸连续形态分布后，游离 H+ 全部由 He 记账（分子分数不贡献游离 H+，
    #    不再有"分子态浓酸"直读分支——pH 即 -log10(自由 H+) 的自然结果）
    # 2) 缓冲滴定（质子条件近似）：游离强酸/强碱先被在账弱碱/弱酸储备按强度
    #    顺序吸收；被全吸收则由最后缓冲对的 Henderson 式定 pH；
    #    全吸收但无有效缓冲对 → 用残余 He（≈0）落分支 4，而非原 He 直读
    tit, He_res, ledger = _buffer_titration(ledger, H_excess, V, T, pKw,
                                            T_K=T_K, cache=cache, touch=touch,
                                            no_base=no_base, no_acid=no_acid)
    if tit is not None:
        _tag("滴定/Henderson")
        return tit, ledger, He_res
    He = He_res / V
    # **直读必须与分支 4 的来源竞争**（第 146 轮，§7 X-43）：残余强酸/碱是
    # 分支 4 的一侧来源，只有当账本里**没有任何**分支 4 角色物种时，直读
    # 才是精确的（那种账本下分支 4 只会给 10^(−pKw/2)，必输给 |He| ≥ 1e-3）。
    # 否则早退会**独占** pH：M01 实测 He = −1e-3 直读 11.000，而同一账本的
    # Al³⁺ 水解给 1.85e-3（pH 2.73）；x 再走 1e-7 就翻到 2.73 ⟹ 7.5 单位的
    # 悬崖、S 由 +23.5 翻 −1.3 ⟹ solve_extent 无过零点、通道被判"零推进"。
    # **"显著"而不是"在场"**（第 276 轮）：原判据是"账本里**出现**任何角色
    # 物种"，而痕量也会命中。实测 `H45 Na[Al(OH)4]+HCl 半量` 终态——
    #   `Na⁺ 1 / Cl⁻ 0.5 / [Al(OH)₄]⁻ 0.874721 / Al(OH)₃ 0.125279 / He = +0.374721`
    # 电荷平衡 `1 + 0.374721 = 1.374721 = 0.5 + 0.874721` **逐位成立**
    # ⟹ 自洽 pH = −log10(0.374721) = **0.426**。但账本里还有**痕量** `Al³⁺`
    # ⟹ 旧判据 `_role_free = False` ⟹ 直读早退被挡 ⟹ 落到 pinned 精确解得
    # **12.185**（`tools/ph_path.py` 实测 tags = `['电荷平衡精确解(pinned)']`）
    # —— 那 0.375 mol 游离强酸在走步眼里**彻底隐形**，于是 `[Al(OH)₄]⁻ +
    # 4H⁺ -> Al³⁺` 再也走不动、账本停在"0.87 M 铝酸根 + 0.37 M 强酸"这种
    # **化学上不可能共存**的态上。
    # 显著性用引擎**同一把尺子** `ANN_MIN_EXTENT`（mol；`resid_live_ok`/
    # 慢标注同用），不引入新常数。
    _role_free = not any(sp in _role_union(T) and m > ANN_MIN_EXTENT
                         for sp, m in ledger.items())
    if _role_free:
        if He >= 1e-3:
            _tag("H⁺直读")
            return max(-1.0, -log10(He)), ledger, He_res
        if He <= -1e-3:
            _tag("OH⁻直读")
            return min(pKw + 1.0, pKw + log10(-He)), ledger, He_res
    # 4) 缓冲/弱酸弱碱区：取各来源贡献最大者（在滴定后的虚拟账本上评估）
    # **阈值连续化**（第 146 轮，§7 X-43）：残余强酸/碱并入同侧下界。
    # 原先 |He| < 1e-3 时残余被**整块丢弃**，阈值两侧是两套模型：N20 实测
    # He = −1e-3 给 pH 11.000，x 再走 1e-7（He = −0.0009997）立刻掉到分支 4
    # 的 3.501 —— 7.5 个 pH 单位的悬崖，使 pick（在 11.0 上评 S = +22.01）与
    # solve_extent 的 f(x)（x > 0 全在 3.5 上评，S ≈ −0.49）**符号相反** ⟹
    # f 无过零点 ⟹ 通道被判"零推进"禁用、走步带 |S| = 23.5 退出。
    # 残余强酸/碱本就是"同侧最强的一个来源"（与分支 4 的 max 语义一致），
    # 取 max 后两支在阈值极限上给出同一个数；|He| ≥ 1e-3 的直读快路径保留
    # （数值等价），真实强碱体系（He 大）完全不走这一支。
    h_c = 10.0 ** (-pKw / 2)
    o_c = h_c
    # 残余强酸/碱这一侧的**来源标签**（§7 X-45）：只有它是"游离质子池"的
    # 记账，而账本物种给出的 h_c/o_c 是"同一池子"的另一套记账。两者数值
    # 重合时是**同一个量的两种记法**，不能互相竞争（下一段判别）。
    _free_side = 0
    if He > 0.0:
        _free_side = 1
        h_c = max(h_c, He)
    elif He < 0.0:
        _free_side = -1
        o_c = max(o_c, -He)

    # （原 _pka1/_pkapp/_pksp 嵌套定义处——已外提至模块级，T_K 作参数）
    # 分支 4 的静态量（两性资格、各酸第一级 Ka、各碱最强共轭酸 pKa、
    # 金属水解 Kh、两性 pH）只依赖数据表与 T_K——按 T_K 缓存，避免每次
    # 调用全表重算（estimate_state 是 solve_extent 二分的最大热点）。
    # 酸/碱/水解/两性/配离子储备的“有效值”函数全部外提至模块级（T_K
    # 作参数传入），消除每调用 4 个嵌套函数对象重建的白开销。
    # acids_map/bases_map/hyd_map 用 dict 而非 list，热路径改为遍历在账
    # 物种（~15 项）而非全表（~75 项），减少无物种命中的空转 get/除法。
    est_cache = getattr(T, "_est_static", None)
    if est_cache is None:
        est_cache = T._est_static = {}
    sc = est_cache.get(T_K)
    if sc is None:
        amph_eligible = set()
        amph_pH = {}
        for sp in set(T.pka_acid) & set(T.pka_base):
            e_a = min(T.pka_acid[sp], key=lambda e: e["pka"])
            e_bs = [e for e in T.pka_base[sp] if e["n"] == 1] or T.pka_base[sp]
            pka_b = max(_pkapp(e, T_K) for e in e_bs)
            if _pka1(e_a, T_K) <= pKw + 2 and pka_b <= pKw + 2:
                amph_eligible.add(sp)   # 如 HCO3-；NH3(pKa105)/HS-(pKa19) 不算
                amph_pH[sp] = (_pka1(e_a, T_K) + pka_b) / 2
        # **显式一级水解**：beta 表里该阳离子存在 nu = 1 的 OH⁻ 配离子
        # （M^{n+} + H₂O ⇌ M(OH)^{(n−1)+} + H⁺）⟹ Ka = β₁·Kw 是实测值。
        # 这些阳离子在分支 4 里按**一元弱酸**处理（见 acids_map 追加与
        # hyd_map 的跳过），不再叠加 Ksp 派生的"水解到底"复合式：后者把
        # 一级水解产物当固相，与显式数据是两套模型（§7 X-32 capstone）。
        # **只收单核一级**（第 170 轮）：`logb` 对多核（m>1）或高配位（nu>1）
        # 是**高级累积常数**，代入会把"第 n 级"当成"一级" ⟹ pH 严重失真。
        # 当前库内该集合的取值不受影响（190 条全 m=1、一级条目全 nu=1）。
        _first_k = {b["center"]: 10.0 ** (b["logb"] - pKw) for b in T.beta
                    if b["ligand"] == "OH^-" and b.get("nu") == 1
                    and b.get("m", 1) == 1}
        acids_map: dict = {}   # acid -> Ka 或 _STRONG_ACID 哨兵
        for acid, entries in T.pka_acid.items():
            if acid in amph_eligible:
                continue
            e1 = min(entries, key=lambda e: e["pka"])   # 第一级
            acids_map[acid] = _STRONG_ACID if e1["pka"] <= 0 else 10.0 ** (-_pka1(e1, T_K))
        for _cat, _ka in _first_k.items():
            if _cat not in amph_eligible:
                acids_map.setdefault(_cat, _ka)
        bases_map: dict = {}   # base -> Kb
        for base, entries in T.pka_base.items():
            if base in T.solids or base == WATER or base in amph_eligible:
                continue
            e1s = [e for e in entries if e["n"] == 1] or entries
            pka = max(_pkapp(e, T_K) for e in e1s)           # 最强一级共轭酸
            bases_map[base] = 10.0 ** (pka - pKw)
        # **两性金属的含氧酸根 → 碱侧角色**（第 276 轮，**在第 275 轮的否证 A
        # 之上重做**）。Kb 完全由库内数据推出，无新数据、无阈值：
        #   `[M(OH)_k] ⇌ [M(OH)_{k-1}] + OH⁻` ⟹ Kb = β_{k-1}/β_k（β₀ ≡ 1）
        # 中间级缺数据（`Al³⁺` 恰缺 ν=3）时**终止于固相**：
        #   `[M(OH)_k] ⇌ M(OH)_z(s) + (k−z)·OH⁻` ⟹ Kb = 10^(pKsp − logβ_k)/(k−z)…
        # 取值见 `tools/ladder_base_census.py`（11 个物种，`[Al(OH)_4]^-`
        # 得 `Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻` 的 logK = **+1.5**，与该 beta 条目
        # 自己的 `calibrated` 注记逐位一致）。
        #
        # **为什么第 275 轮的否证 A 现在可以重做**：那次失败的直接原因是
        # `_buffer_titration` 用"一步脱光"（4 H⁺/铝酸根）把账本留成
        # **0.78 M 铝酸根 + 0.094 M Al³⁺** 这种不可能共存的混合态，再给它
        # 一个碱角色，`o_c` 就被那个虚高的浓度撑爆（教科书例
        # `16 AlCl3+3NaOH` 从 Al(OH)₃ 变成纯铝酸根）。本轮先把储备改成
        # **逐级（终止于固相）**（见 `_buffer_titration` 的 `_step`），
        # 滴定后的铝酸根回到它应有的量，角色才有意义。
        _ksp_oh = {}
        for _e in T.ksp:
            _c0, _a0 = _e["pair"]
            if _a0 == "OH^-" and _c0 not in _ksp_oh:
                _ksp_oh[_c0] = _e
        _lad: dict = {}
        for b in T.beta:
            if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
                continue
            cx = b.get("complex")
            if not cx or cx in bases_map or cx in T.solids:
                continue
            if charge_of(cx) >= 0:
                continue                      # 只收含氧酸根（净负电）
            _lad[(b["center"], b.get("nu", 1))] = (cx, float(b["logb"]))
        for (_ctr, _k), (cx, _logbk) in _lad.items():
            _e = _ksp_oh.get(_ctr)
            _z = charge_of(_ctr)
            if _e is None:
                continue                      # 无 Ksp ⟹ 这条梯没有固相终点
            # **一律以固相为终点**（与滴定储备同一口径）：`k − z` 个 OH⁻
            _n = _k - _z
            if _n <= 0:
                continue
            kb = 10.0 ** ((_pksp(_e, T_K) - _logbk) / _n)
            if kb > 0.0:
                bases_map[cx] = kb

        # 实测（`tools/branch4_src.py H43`）：`H43 AlCl3+NaOH 1:3.5` 的**呈现
        # pH = 3.26**，而同一账本的电荷平衡精确解（引擎自己的 `pH_solver`）
        # 是 **12.037**，与手算 12.03 逐位吻合。成因：分支 4 的 `o_c` 扫的是
        # `rolemap`，而账本里 **0.845 M 的 `[Al(OH)_4]^-`** 在其中**没有角色**
        # ⟹ `o_c` 看不见它 ⟹ 只剩痕量 `Al^{3+}`（0.031 M）撑起 `h_c = 3.2e-3`
        # ⟹ 启发式以为这是酸液，选酸侧给 3.26。
        #
        # **第一版修法（把含氧酸根补进 `bases_map`）被全量实测否掉**：
        # Kb 本身可由库内数据精确推出（`[M(OH)_k] ⇌ [M(OH)_{k-1}] + OH⁻`
        # ⟹ Kb = β_{k-1}/β_k；中间级缺数据时退回固相 Kb = 1/(β_k·Ksp)，
        # 校验 `[Al(OH)_4]^-` 得 Kb = 0.0316 ⟹ `Al(OH)₃ + OH⁻ ->
        # [Al(OH)₄]⁻` 的 logK = **+1.5**，与该 beta 条目自己的 `calibrated`
        # 注记逐位一致；普查见 `tools/ladder_base_census.py`，11 个物种）。
        # 但 `o_c` 是"游离 [OH⁻]"的**代理量**，弱碱式给它 0.148（0.845 M ×
        # Kb=0.0316），而该态的电荷平衡真值只要 0.011 ⟹ **高估 13 倍**：
        #   通过 1208 → **1179（翻红 29 / 翻绿 0）**、`n(|S|>1)` 40 → 54；
        #   受害例形态高度一致 —— `16 AlCl3+3NaOH`（**恰好 3 当量**）给出
        #   pH 13.08 + 净方程 `4OH⁻ + Al³⁺ -> [Al(OH)₄]⁻`，而教科书产物是
        #   **Al(OH)₃**（要第 4 个 OH⁻ 才溶成铝酸根）；`19 ZnCl2+2NaOH`、
        #   `E55 明矾+适量NaOH`、`E41/E43/N18/N19/N20` 同形。
        # ⟹ **"角色表缺谁"与"该用多大浓度"是两件事**：补角色会把启发式的
        # 代理量算错。正确的方向见下面 `_blind` 的**前提判据**。
        # ⚠️ 第 277 轮删掉了此处**重复的一份** `bases_map` 定义与一份过期的
        # `_ladder_cx`：第 275/276 轮两次插入把"pKa 碱储备 + 含氧酸根碱侧角色"
        # 写成了**两遍**，后一遍 `bases_map: dict = {}` 把前一遍**整块清掉**
        # ⟹ 含氧酸根的角色一条也没生效（`tools/branch4_src.py` 实测
        # `[Al(OH)_4]^-` 不在 182 项的 `rolemap` 里 ⟹ 分支 4 给 H45 pH 7.0）。
        # 教训：**同一名字的 `x = {}` 出现两次就是静默清空**，见 lessons。
        hyd_map: dict = {}    # cat -> (Kh 每阳离子, qc)
        for e in T.ksp:
            cat, an = e["pair"]
            if an != "OH^-" or cat in _first_k:
                # 一级水解已有实测 β ⟹ 酸侧的一元弱酸路线接管（上面）；
                # 这里若再放 Ksp 复合式就是同一物种两套模型并存。
                continue
            _x, _y = _ksp_xy(e)
            # Kh = c/Kw 通道常数；第二元 = **每阳离子释放的质子数** n = y/x
            # （与 charge_of(cat) 在全表上逐条相等，此处直接取 y/x 以免
            #  将来出现 M₂O 型条目时指数悄悄错位）
            hyd_map[cat] = (10.0 ** ((_pksp(e, T_K) - _y * pKw) / _x),
                            _y / _x)
        # **两性物种 → (共轭酸, 共轭碱)**（第 263 轮）：判定账本里的两性物种
        # 是"纯两性盐"还是"**缓冲对的一员**"。中点式 `(pKa1+pKa2)/2` 的推导
        # 前提是 `[H₂A] = [A²⁻]` **由歧化自身产生**（即该物种是唯一质子条件
        # 物种）；共轭伙伴同时在账时，`[H₂A]`（或 `[A²⁻]`）由投料定，
        # 前提整个不成立 ⟹ 那是缓冲体系，不是纯两性盐。
        # 实测（`tools/amph_exact.py`，独立闭式解为判据）：
        #   H₂S 半中和 中点 10.500 / 精确 7.0000；CO₂+NaOH 1:1 中点 8.350 / 6.3999。
        # 而**纯**两性盐（NaHCO₃ 8.350/8.3456、纯 NaHS 10.500/9.4972）里，
        # 前者的残差来自"丢了 h/oh"的近似误差（≤1.0），不是模型错——
        # 第 263 轮先只在**模型错**的那一侧接管（见 §5 的记账）。
        amph_partner: dict = {}
        for sp in amph_eligible:
            _eca = None
            for _e in T.pka_base.get(sp, ()):
                if _e.get("n", 1) == 1:
                    _eca = _e.get("acid")
                    break
            _ecb = None
            for _e in T.pka_acid.get(sp, ()):
                if _e.get("n", 1) == 1:
                    _ecb = _e.get("base")
                    break
            amph_partner[sp] = (_eca, _ecb)
        # 共轭酸碱对映射（缓冲对识别）：base -> (acid, pKa)，取最强一级
        conj: dict[str, tuple[str, float]] = {}
        for e in T.pka:
            if e.get("n", 1) != 1:
                continue
            b, a = e["base"], e["acid"]
            if b not in conj or _pkapp(e, T_K) > conj[b][1]:
                conj[b] = (a, _pkapp(e, T_K))
        # 物种角色表（酸/碱/水解/两性一次解析）：热路径每物种 1 次
        # dict.get 取代 5-6 次独立查表（原 54k 次调用/慢例的最大空转）。
        # 角色互斥性与原分支顺序一致：amph 从 acids/bases 排除；Kb
        # 分支 continue 跳过水解/两性；Ka 分支不 continue（可叠加水解）
        rolemap: dict = {}
        for sp in set(acids_map) | set(bases_map) | set(hyd_map) | set(amph_pH):
            rolemap[sp] = (acids_map.get(sp), bases_map.get(sp),
                           conj.get(sp) if sp in bases_map else None,
                           hyd_map.get(sp),
                           amph_pH.get(sp) if sp not in T.solids else None)
        # **只收净负电的含氧酸根**（铝酸根/锌酸根/锡酸根…）：它们是"完全水解
        # 的产物"，在账本里可以成为**主要物种**却对角色表完全隐形。阳离子型
        # 部分水解产物（`[Ca(OH)]⁺`、`[Al(OH)]^{2+}` 等）**不收**——`hyd_map`/
        # `_first_k` 已经覆盖它们，且实测把 `[Ca(OH)]⁺` 算进来会打坏
        # `Y18 CaO+水 放熟石灰`（`Ca(OH)₂` 0.9999 → 0.0202、残差 0 → 4.556）。
        _ladder_cx = {b["complex"] for b in T.beta
                      if b["ligand"] == "OH^-" and b.get("m", 1) == 1
                      and b.get("complex") and charge_of(b["complex"]) < 0}
        # 配离子 -> 其中心的氢氧化物固相：**第 280 轮已删除**（`_pin_pair`
        # 判据连同它的 `_cx_solid` 表一起回退，见下方 `_blind_dom` 处的否证记录）
        sc = (amph_eligible, amph_pH, acids_map, bases_map, hyd_map, conj,
              rolemap, amph_partner, _ladder_cx)
        est_cache[T_K] = sc
    (amph_eligible, amph_pH, acids_map, bases_map, hyd_map, conj, rolemap,
     amph_partner, _ladder_cx) = sc
    # 单遍扫描在账物种：合并原 acids_l/bases_l/hyd_l/amph/buf 五个独立循环。
    # max/sum 可换序，h_c/o_c/amph/buf 的最终值与原实现等价。
    amph: list = []
    buf: list = []
    _floor_V = 1e-12 * V
    _role_get = rolemap.get
    _src = PH_SRC      # 诊断开关：一次全局查找，循环内只判 None
    # **被角色表忽略的羟合梯配合物**（第 275 轮）：账本里存在、但 `rolemap`
    # 没有它的酸碱角色 ⟹ 分支 4 的 `h_c`/`o_c` 对它**完全盲**。
    _blind = 0.0
    for sp, m in ledger.items():
        if m <= _floor_V:
            continue
        c = m / V
        role = _role_get(sp)
        if role is None:
            if sp in _ladder_cx and c > _blind:
                _blind = c
            continue
        Ka, Kb, conj_pair, Kh_qc, amph_v = role
        if Ka is not None:
            if Ka is _STRONG_ACID:
                h_c = max(h_c, c)
                if _src is not None:
                    _src.append((sp, "强酸", c))
            else:
                _v = (-Ka + sqrt(Ka * Ka + 4 * Ka * c)) / 2
                h_c = max(h_c, _v)
                if _src is not None:
                    _src.append((sp, "Ka", _v))
        if Kb is not None:
            if Kb >= 1.0:                                # 水解近完全（S2-、C2^2- 等）
                o_c = max(o_c, c)
                if _src is not None:
                    _src.append((sp, "Kb≥1", c))
            else:
                _v = (-Kb + sqrt(Kb * Kb + 4 * Kb * c)) / 2
                o_c = max(o_c, _v)
                if _src is not None:
                    _src.append((sp, "Kb", _v))
            # 共轭缓冲对（仅碱在账时检查其共轭酸是否也在账）
            if conj_pair is not None:
                acid_conj, pka_c = conj_pair
                ca = ledger.get(acid_conj, 0.0) / V
                if ca > 1e-12:
                    buf.append((pka_c + log10(c / ca), min(c, ca)))
            continue
        if Kh_qc is not None:
            Kh, npp = Kh_qc
            if npp <= 1.0:
                # 一价金属水解到底即纯固相（活度 1），无累积共轭碱：
                # cat + H2O → ½M2O + H+ 给出 h = Kh·c；套弱酸二次式
                # 会把 Ag+ 类高估 ~1/sqrt(Kh·c) 倍（D26：Ag+ 被估成 pH 3
                # 的酸，驱动铬酸根幻影质子化死循环）。
                _v = Kh * c
                h_c = max(h_c, _v)
                if _src is not None:
                    _src.append((sp, "Kh(1)", _v))
            else:
                # **n 质子真式**（§7 X-32 capstone，第 140 轮落地）：
                # `h^n = Kh·(c − h/n)`。旧实现是"一元弱酸"形态的二次式
                # `h² = Kh·c`——它在 n≥3 时把 pH 抬得过高（Al³⁺ 1 M 给
                # 4.50、真值 3.00；Fe³⁺ 0.1 M 给 2.51、真值 ~1.7）。
                # 之所以能换：Fe/Al 两轴的显式羟基络合物已入库（第 137/139
                # 轮），旧式原本**补偿**缺失中间体的那部分职责已由数据承担。
                _v = _nth_root_h(Kh, c, npp)
                h_c = max(h_c, _v)
                if _src is not None:
                    _src.append((sp, f"Kh({npp:g})", _v))
            continue
        if amph_v is not None and c > 1e-6:
            # 共轭伙伴是否**同时在账**（第 263 轮）：伙伴在场 ⟹ 缓冲对，
            # 非纯两性盐 ⟹ 中点式的推导前提不成立。
            _pa, _pb = amph_partner.get(sp, (None, None))
            _has_pair = bool(
                (_pa is not None and ledger.get(_pa, 0.0) > _floor_V)
                or (_pb is not None and ledger.get(_pb, 0.0) > _floor_V))
            amph.append((c, amph_v, _has_pair))
    # **X-45 的两种判据均被实测否掉**（第 151/152 轮），当前**不修**，只留标签
    # `_free_side`（标记 h_c/o_c 的哪一侧来自游离质子池 `He_res/V`）：
    #   · 判据 A（按数值相当去重，1e-2 相对容差）：D47 悬崖消除、iters −0.33%，
    #     但打红**未锚定约定**的 H88（醋酸铵：o_c 来自 NH₃ 的 Kb 式，账本真来源）；
    #   · 判据 B（按来源标签去重：游离质子池侧与账本对侧相差 10 倍以内即视为
    #     同一池子、弃用游离侧）：D47 悬崖**完全消除**（pH 全程 1.2318），
    #     但**锚定约定 12 例翻红**（16/E55/N20/N34/T51/T52/H42/M01/B26/X04/
    #     U01/Amp14，全是 Al³⁺/Fe³⁺ 碱量不足体系）——那些体系的游离质子池与
    #     账本侧数值同样接近，**但两者是独立来源**（外加 OH⁻ vs 水解产物），
    #     弃用游离侧 ⟹ Al(OH)₃ 产量 0.545 < 0.9。
    # ⟹ "数值接近"与"来源标签"都不足以判别"同一量的两种记法"。真判据必须是
    # **代数同一性**（把账本代回电荷平衡，看该项是否恒等于 He_res/V），
    # 而不是任何形式的数值比较。全部否证数字见 log §7 X-45。
    # 两性物种（HCO3-、HS-、H2PO4- 等）：pH ≈ (pKa_酸 + pKa_共轭酸)/2，
    # 其浓度远大于其他酸碱贡献时以两性平衡为准（NaHCO3 溶液 pH≈8.3）
    if amph:
        c_a, pH_a, _pair = max(amph, key=lambda t: t[0])
        if c_a > 100.0 * max(h_c, o_c):
            # **精确质子条件优先——仅当模型结构性错误时**（第 263 轮根治）：
            # 中点式 `(pKa1+pKa2)/2` 不是独立的化学模型，它是**质子条件的
            # 近似**——由 `[H₂A] = [A²⁻]`（忽略 h 与 oh）推出，并默认
            # "该两性物种是**唯一**的质子条件物种"。**两个前提都会破**，
            # 但性质不同：
            #   · **模型错**（`_pair`）：其共轭伙伴同时在账 ⟹ `[H₂A]` 由投料
            #     定、不由歧化定 ⟹ 它是**缓冲对的一员**，中点式的前提整个
            #     不成立。H₂S 半中和 中点 10.500 / 精确 7.0000 (= pKa1)；
            #     CO₂+NaOH 1:1 中点 8.350 / 精确 6.3999。**没有任何近似能修**，
            #     只能换成方程 ⟹ 本支路接管。
            #   · **近似误差**（纯两性盐、无伙伴）：中点式丢了 h/oh，
            #     纯 NaHS 0.01 M 中点 10.500 / 精确 9.4972（+1.0）。
            #     第 263 轮**不接管**这一侧——不是不能（`charge_pH` 在
            #     NaHS 上与独立闭式解差 0.0000），而是实测接管会把 N15
            #     （硅酸体系）推进走步层极限环（见 §5「家族再分配抵消走步」），
            #     收益/风险不成比例 ⟹ 记账待后续，不在本轮赌。
            # 触发面实测：中点支路仅占全部 `estimate_pH` 调用的 **0.71%**
            # （1,314,020 次中 9,362 次，`tools/amph_survey.py`）。
            # 判据与正/负对照见 `tools/amph_exact.py`（引擎无关的闭式解）。
            if _pair and _exact_ok(ledger, He_res, T, V):
                _ex = charge_pH(ledger, V, T, T_K, fast=True)
                if _ex is not None:
                    _tag("电荷平衡精确解(两性)")
                    return min(max(_ex, -1.0), pKw + 1.0), ledger, He_res
            _tag("两性中点")
            return min(max(pH_a, -1.0), pKw + 1.0), ledger, He_res
    # 共轭缓冲对：弱碱与其共轭酸（或反之）同时在账时，pH 由
    # Henderson-Hasselbalch 决定（pKa + log(c_b/c_a)），sqrt(Kb·c) 的
    # 单物种估计会把缓冲体系误判为强碱性（NH4Ac 曾被估到 pH 10.4——
    # NH3 浓度稍高即触发；真实是 NH4+/NH3 与 Ac-/HAc 双缓冲 ≈7）。
    # 多对并存按弱组分浓度加权平均（双水解盐 → (pKa+pKa')/2 语义）
    if buf:
        w_sum = sum(w for _, w in buf)
        if w_sum > max(h_c, o_c):
            pH_buf = sum(p * w for p, w in buf) / w_sum
            _tag("缓冲对")
            return min(max(pH_buf, -1.0), pKw + 1.0), ledger, He_res
    # **X-45 已诊断、未落地**（第 151 轮）：h_c 与 o_c 可能是**同一个质子失衡的
    # 两种记法**——D47 实测酸侧 `Fe³⁺ 弱酸式` 与碱侧 `−He_res/V` 在
    # He ≈ −0.0586459 处给出**同一个 0.0586**，两者相对差随 He 平滑穿过 0，
    # 谁大谁小由浮点末位决定 ⟹ pH 在 1.2318 / 12.7682 之间跳 11.5 个单位，
    # 走步恰好骑在刃上：该沉淀的不沉淀、该溶解的不溶解（退出残差 S = +34.086）。
    # 想过的"两侧相当就去掉碱侧那一项"**实测不可行**：1e-2 相对容差会打红
    # 未锚定约定的 H88（醋酸铵真缓冲：h_c = 1.8e-5 与 o_c = 1.3e-5 差 32%，
    # 但账本里没有"游离强碱"这一项，o_c 来自 NH₃ 的 Kb 式，是真来源不能删）。
    # ⟹ 判据必须能区分"同一量的两种记法"与"两个独立来源的巧合接近"，
    # 下一轮从**来源标签**（哪一项是 He_res/V、哪一项是账本物种）入手。
    # **精确质子条件接管**（第 197 轮，X-45 的根本修法）：
    # 分支 4 的 h_c/o_c 都是**启发式估计**；当两侧**量级相当**（比值 ≤ 10）时，
    # 谁大谁小由浮点末位/分支条件决定 ⟹ 会**稳定选错侧**（第 196 轮 E55 实测：
    # 引擎停在碱侧 pH 10.396，而同账本的电荷平衡自洽解是 3.604，差 6.8 单位，
    # 沉淀通道因此永远零推进）。
    # 此时改用 `charge_pH`（账本电荷平衡的**精确**解，已验证 12/12 教科书锚点
    # + 80/80 条目自洽）：它把"选哪一侧"变成"解一个方程"，即第 161 轮判定的
    # **代数同一性**。只在"量级相当"这个**退化区**触发 ⟹ 热路径成本可控
    # （`fast=True` 为阻尼 Newton，4–6 次求值；非退化区仍走原启发式）。
    # **角色盲区判据**（第 275 轮·§1.12）：上面那句"启发式只在两侧量级相当
    # （退化区）时才让位给精确解"隐含一个**前提**——`h_c`/`o_c` 确实看到了
    # 账本里的酸碱内容。前提不成立时这个"量级相当"判断本身没有意义。
    #
    # 实测（`tools/branch4_src.py H43`）：`H43 AlCl3+NaOH 1:3.5` 的账本里有
    # **0.845 M `[Al(OH)_4]^-`**，而角色表里没有它 ⟹ `o_c` 只看到纯水的
    # 1e-7、`h_c` 只看到痕量 `Al³⁺` 给的 3.2e-3 ⟹ 比值 3 万，"不退化"
    # ⟹ 启发式选酸侧给 **pH 3.26**；而同账本的电荷平衡精确解（引擎自己的
    # `pH_solver`）是 **12.037**，与手算 12.03 逐位吻合 —— **差 8.8 个
    # pH 单位**，`resid_max` 被算成 24.511。
    #
    # 判据只用"读账本"、不含新阈值：**被忽略的羟合梯配合物比启发式实际
    # 用到的任何量都大** ⟹ 启发式无权选支，改用精确解。尺度直接取
    # `max(h_c, o_c)`（启发式自己的量级），不引入新常数。
    #
    # ⚠️ **不要改成"给含氧酸根补 Kb 角色"**（第 275 轮已否，数字见
    # `bases_map` 上方注释）：Kb 可精确推出，但 `o_c` 是"游离 [OH⁻]"的
    # 代理量，补角色会把它算成 0.148 而真值 0.011 ⟹ 通过 1208→1179、
    # 教科书例 `16 AlCl3+3NaOH`（恰好 3 当量）从 Al(OH)₃ 变成纯铝酸根。
    #
    # ⚠️ **只在启发式用的是"软估计"时才推翻它**（第 276 轮补的第二个前提）：
    # `h_c`/`o_c` 有两个来源 —— **残余游离强酸/碱池**（`He_res/V`，
    # **精确记账**，不是估计）与**角色表派生的弱酸/弱碱估计**。前者是硬量，
    # 由它选出的支**可信**；只有后者才可能因为"看不见某个物种"而选错。
    # 实测 `H45 Na[Al(OH)4]+HCl 半量`：终态账本 `Na⁺ 1 / Cl⁻ 0.5 /
    # [Al(OH)₄]⁻ 0.874721 / Al(OH)₃ 0.125279 / He = 0.374721`
    # —— 电荷平衡 `1 + 0.374721 = 1.374721 = 0.5 + 0.874721` **逐位成立**
    # ⟹ 自洽 pH = −log10(0.374721) = **0.426**；而 `_blind_dom` 不设此闸时
    # 会推翻它、改走精确解给 **12.19**（差 11.8 个单位，且让那 0.375 mol
    # 游离强酸在走步眼里彻底隐形）。
    _blind_dom = _blind > max(h_c, o_c) and _free_side == 0
    # ⛔ **第 278–280 轮否证并已回退：`_pin_pair`（固相 + 其含氧酸根同时在账
    # ⟹ 交 pinned 精确解）**。理由本来很干净 —— 分支 4 的弱碱式
    # `[OH⁻] = (−Kb+√(Kb²+4Kb·c))/2` 假设"该碱是溶液里唯一的碱"，
    # 而固相在场时自由配离子由 `[M(OH)_k] = β_k·Ksp·[OH⁻]^{k−z}` 钉死
    # （`H45` 的滴定后态：铝酸根 0.5 + 固相 0.5，启发式 13.045、
    # 精确解 **12.186**）。但**实测三段全负**：
    #   · 第 278 轮：只把 `_pin_pair` 加进**外层**闸 ⟹ 全量**逐位中性**
    #     （因为内层两条路都还没含它）；
    #   · 第 279 轮：给内层"游离阳离子 > 0"闸加豁免 ⟹ `H45` **没修好**
    #     （仍 13.045 / 0 步），`H43` 12.185 → 13.045、残差 0.147 → **0.984**；
    #   · 第 280 轮：查清真正的卡点 —— 内层**两条路**（pinned 与
    #     非 pinned 精确解）的条件里都没有 `_pin_pair`；给它们都补上后
    #     `_pin_pair` **终于真的生效**，结果：`16 AlCl3+3NaOH`
    #     pH 6.26 → **7.0**（教科书产物 Al(OH)₃ 被毁）、`H43` 同样
    #     **13.045 / 残差 0.984**，而 `H45` **依旧 13.045 / 0 步**。
    # ⟹ "有固相就把该态交给 pinned 精确解"是**错的**：pinned 口径会
    # 把正确的沉淀体系（`16`）算坏，也救不了 `H45`。
    # **下一轮的正确方向**（本轮已算清，见 handoff §1.19）：
    # 不再借道 pinned 机器，而是**单独**解"阴离子被固相钉住"的那一个未知量
    # `[A] = β_k·Ksp·[OH⁻]^{k−z}` 与账本电荷平衡的联立（`H45` 一行可解：
    # `a + oh = 0.5` 且 `a = 31.6·oh` ⟹ `oh = 0.01534` ⟹ pH **12.186**）。
    # 需先查 `_pin_ladders(T)` 对 `Al³⁺` 的梯（ν=1,2,**4**，ν=3 是已知数据
    # 缺口）是否按"相邻 ν"处理 —— 这很可能是 pinned 口径在此失准的原因。
    _deg = (h_c > 0.0 and o_c > 0.0 and 0.1 * o_c <= h_c <= 10.0 * o_c)
    if _deg or pin_mode or _blind_dom:
        # **固相在场时用 `pinned` 把储库自由度写进同一个方程**（第 198/199 轮）：
        # 第 197 轮查明，E55 这类病灶的 `h_c ≈ o_c`（触发条件本就满足），
        # 被挡是因为"排除固相储库"那条闸——而 `charge_pH` 的 `pinned` 参数
        # 正是为储库而生：`[M] = 10^((y(pKw−pH) − pKsp)/x)` 作为电荷项的
        # **一项**随 pH 连续变化，于是"固相是否在场"不再是分支选择。
        # 钉住后须删掉该阳离子**及其全部配合物**的账本条目（电荷由
        # 钉住项+羟合梯+非 OH⁻ 从属项表达），否则双重记账/配平幻觉（第 198 轮
        # 教训：只删游离离子 ⟹ 0.75 mol 铝酸根被冻结，方程凭空造 0.25 M 游离
        # Al³⁺ 配平，pH 冻结在 3.2 且走步失去梯度；第 199 轮残余：非 OH⁻
        # 配合物（[AlCl]²⁺、Zn-Cl）仍冻结 ⟹ 幻影正电荷把 pinned 解抬到
        # 12.3/13.26——第 200 轮落地：一并删除 + 配体按 ν 返还账本）。
        _pin = []
        _inv = []
        _led2 = None
        # pin_mode（第 201 轮，D14）：True=整步冻结的"采用"口径——绕过退化区闸
        # （闸随 f(x) 的 x 穿界 = 接缝悬崖，B26 实测 pH 13.73→1.52 跳变被二分
        # 当根）；False=冻结"弃用"；None=逐点自判定。结构性前提仍逐点评估。
        if (PINNED_TAKEOVER and pin_mode is not False
                and (_deg or pin_mode or _blind_dom)):
            _lad_all = _pin_ladders(T)
            _frz_all = _pin_frozen(T)
            _seen = set()
            for _e in T.ksp:
                _cat, _an = _e["pair"]
                if _an != "OH^-" or _cat in _seen:
                    continue
                # **第 278 轮：不再要求该阳离子在账本里"是键"**。
                # 第 267 轮已把"游离阳离子量 > 0"这一条放宽（`pin_mode is True`
                # 时不受限），但保留了一道更外层的闸 —— `_cat not in ledger`
                # （键都不在就直接跳过）。实测 `H45 Na[Al(OH)4]+HCl 半量`：
                # 滴定把 0.5 mol H⁺ 吸收成 `[Al(OH)₄]⁻ 0.5 / Al(OH)₃(s) 0.5`，
                # 虚拟账本里**没有 `Al^{3+}` 这个键**（它的量由 Ksp 给出）
                # ⟹ 整个 pinned 块从不进入 ⟹ 退回分支 4 的弱碱启发式
                # **13.045**，而该态电荷平衡 + `a/[OH⁻] = 10^1.5` 的精确解是
                # **12.19**。**固相在场才是钉住的前提**（下一道闸已把关），
                # 键在不在不是前提。
                # 固相必须真的在场（账本里有该固相且量显著）
                _solid = _e["solid"]
                _ns = ledger.get(_solid, 0.0)
                if _ns <= X_MIN:
                    continue
                # **第 267 轮：冻结采用 pinned 时，账本里的游离阳离子量不作数。**
                # 钉住口径下 `[M] = 10^((y(pKw−pH) − pKsp)/x)` 由 **Ksp 给出**，
                # 账本里那点游离量是走步的记账、不是约束。而某一步可能把该阳离子
                # **整步消耗到 0**——`H45 Na[Al(OH)4]+HCl 半量` 实测：残差通道
                # `Al³⁺ -> [Al(OH)₄]⁻ + 4H⁺` 的 `x_max` **恰等于**账本 `Al³⁺`
                # (0.016999)，于是 `x_max` 端 `Al³⁺ = 0` ⟹ 原闸跳过 ⟹ `_pin` 空
                # ⟹ 即使 `pin_mode=True` 也会在端点退回逐点自判定 ⟹ 退回
                # `酸侧max`（pH 12.185 → 3.70，跳 −8.48）⟹ 二分落在跳上
                # ⟹ `x*=0`、残差 24.991。
                # 故：**整步冻结为 True 时不受此闸限制**（前提仍由固相在场把关）；
                # 逐点自判定（`pin_mode is None`）与显式弃用（`False`）时
                # 行为**逐位不变**。
                # ⛔ **第 279 轮否证**：再加一个豁免口 `not _pin_pair`
                # （理由是"固相 + 其含氧酸根同时在账 ⟹ 固相就是储库"）——
                # 实测 `H45` **仍是 13.045 / 0 步**（没修好），而 `H43`
                # 12.185 → **13.045**、残差 0.147 → **0.984**（打坏了）
                # ⟹ 已回退。`_pin_pair` 进不去 pinned 块**不是**卡在这一行。
                if ledger.get(_cat, 0.0) <= 0.0 and pin_mode is not True:
                    continue
                _seen.add(_cat)
                _x, _y = _ksp_xy(_e)
                _lad = _lad_all.get(_cat, ())
                if _led2 is None:
                    _led2 = dict(ledger)
                _led2.pop(_cat, None)
                _n_inv = ledger[_cat] + _x * _ns     # 金属库存（回检用）
                for _zc, _lb, _nu, _cx in _lad:
                    if _led2.pop(_cx, 0.0) > 0.0:
                        _n_inv += ledger[_cx]
                # 非 OH⁻ 配合物：先删除并**返还配体**（ν×m 回账本），
                # 再按返还后的账本值冻结 cL，coef = β·cL^ν 折进钉住项。
                _frz_c = []
                for _zc, _lb, _nu, _lig, _cx in _frz_all.get(_cat, ()):
                    _mc = _led2.pop(_cx, 0.0)
                    if _mc > 0.0:
                        _n_inv += _mc
                        _led2[_lig] = _led2.get(_lig, 0.0) + _nu * _mc
                        _frz_c.append((_zc, _lb, _nu, _lig))
                _frz = tuple((_zc, (10.0 ** _lb)
                              * (_led2.get(_lig, 0.0) / V) ** _nu)
                             for _zc, _lb, _nu, _lig in _frz_c)
                _pin.append((charge_of(_cat), _pksp(_e, T_K), _x, _y,
                             _lad, _frz))
                _inv.append(_n_inv)
        if _pin:
            _ex = charge_pH(_led2, V, T, T_K, fast=True, pinned=tuple(_pin))
            if _ex is not None:
                # **储库耗尽回检**：钉住隐含的溶解总量超过库存 ⟹ "固相在场"
                # 的前提不成立（固相本该溶完），退回原闸——自洽性判据。
                _oh = 10.0 ** (_ex - pKw)
                _ok = True
                for _p, _n_inv in zip(_pin, _inv):
                    _z, _pk, _x2, _y2, _lad = _p[:5]
                    _frz = _p[5] if len(_p) > 5 else ()
                    _cm = 10.0 ** ((_y2 * (pKw - _ex) - _pk) / _x2)
                    _d = V * _cm * (1.0 + sum(10.0 ** _lb * _oh ** _nu
                                              for _zc, _lb, _nu, _cx in _lad)
                                    + sum(_cf for _zc, _cf in _frz))
                    if _d > _n_inv:
                        _ok = False
                        break
                if _ok:
                    _tag("电荷平衡精确解(pinned)")
                    return min(max(_ex, -1.0), pKw + 1.0), ledger, He_res
                _tag("pinned耗尽回绝")
        # 非 pinned 精确解路径：退化区（_deg）**或角色盲区（_blind_dom）**走。
        # pin_mode=True 绕过退化区闸只为 pinned 接管，不扩大此路径的触发面。
        # `_blind_dom`（第 275 轮）是**前提判据**：角色表看不见账本里的羟合梯
        # 配合物时，`h_c/o_c` 的"量级相当"比较本身无意义，交回精确解。
        if _deg or _blind_dom:
            # 判据已抽为 `_exact_ok`（第 263 轮）：与本函数上方"两性支路"
            # 共用同一把尺子。此处仍取 `min_fams=2`（行为逐位不变）。
            if _exact_ok(ledger, He_res, T, V):
                _ex = charge_pH(ledger, V, T, T_K, fast=True)
                if _ex is not None:
                    _tag("电荷平衡精确解")
                    return min(max(_ex, -1.0), pKw + 1.0), ledger, He_res
    pH = -log10(h_c) if h_c >= o_c else pKw + log10(o_c)
    _tag("酸侧max" if h_c >= o_c else "碱侧max")
    if _src is not None:
        # 审计（第 152 轮）：标出**哪一侧来自游离质子池**——"同一量的两种记法"
        # 的判别要靠它，光看数值分不出（§7 X-45）。
        _src.append(("__branch4__", "h_c" if h_c >= o_c else "o_c",
                     h_c if h_c >= o_c else o_c))
        _src.append((f"__free_side={_free_side}__", "He_res/V",
                     abs(He_res) / V))
    return min(max(pH, -1.0), pKw + 1.0), ledger, He_res


# 分子态强酸的再平衡参数：酸 -> (共轭阴离子, 每分子释出质子数)
_RESPECIATE_ACIDS = {"HNO_3": ("NO_3^-", 1), "H_2SO_4": ("SO_4^{2-}", 2)}


def _pka1(e: dict, T_K: float) -> float:
    """该 pKa 条目在 T_K 的有效值（pKa=−logKa，van't Hoff 变号）。
    模块级函数（原为 estimate_state 内嵌套定义——每次调用都重建函数
    对象，estimate_state 是二分探针的最大热点，38k 次/慢例的白开销）。"""
    return e["pka"] - (_vant(e["dH"], T_K) if "dH" in e else 0.0)


def _pkapp(e: dict, T_K: float) -> float:
    """每质子有效 pKa（多元酸/n>1 条目按 1/n 折算；对 pKa 变号）。"""
    return e["pka"] / e["n"] - (_vant(e["dH"] / e["n"], T_K) if "dH" in e else 0.0)


def _pksp(e: dict, T_K: float) -> float:
    """pKsp(T)：pKsp = −logKsp，van't Hoff 变号。"""
    return e["pKsp"] - (_vant(e["dH"], T_K) if "dH" in e else 0.0)


def _pka_eff(pka: float, dH_pp, T_K: float) -> float:
    """有效 pKa（_buffer_titration 用；原内嵌套定义，同上外提）。"""
    # pKa = −logKa：van't Hoff 修正对 pKa 变号（吸热电离 T 升 pKa 降）
    return pka - (_vant(dH_pp, T_K) if dH_pp is not None else 0.0)


def _nth_root_h(Kh: float, c: float, n: float) -> float:
    """解 **n 质子水解真式** `h^n = Kh·(c − h/n)`（h ≥ 0）。

    `f(h) = h^n − Kh·(c − h/n)` 在 [0, n·c] 上单调递增（h^n 与 +Kh·h/n 同向），
    二分 60 次即达浮点精度。n = 1 时退化为 `h = Kh·c/(1+Kh)`（与旧的一元式
    `h = Kh·c` 在 Kh≪1 时同值，故 n≤1 仍走原路，行为不变）。
    """
    if c <= 0.0 or Kh <= 0.0:
        return 0.0
    lo, hi = 0.0, n * c
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mid ** n - Kh * (c - mid / n) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _respeciate_strong_acids(ledger: dict, H_excess: float, V: float, T) -> float:
    """分子态强酸 ⇌ 离子的逐轮再平衡：目标分子池 = f(c_tot)×库存，
    c_tot 由影子库存（"__tot_" 键，仅记本酸，盐类同离子不污染）给出。
    反应消耗游离 H+ 后分子池必须再电离，否则质子被锁死（Cu+4M HNO3
    停滞于 88%、Ba(OH)2+H2SO4 等量中和后残碱）；氧化剂通道消耗分子池/
    阴离子池后库存自然缩水，不会被无限回补。H2SO4 按 1:2 释出质子
    （全电离约定记账的表观近似；HSO4- 中间态依裁决不注册）。
    反应生成的酸（如 NO2 溶于水产生的 HNO3）无影子库存，不做分子化——
    稀区分子分数本来就 <1%，与全电离约定一致。返回修正后的 H_excess。"""
    He = H_excess
    for acid, (anion, n_H) in _RESPECIATE_ACIDS.items():
        shadow = ledger.get(f"__tot_{acid}", 0.0)
        if shadow <= 0.0:
            continue
        m_mol = ledger.get(acid, 0.0)
        m_an = ledger.get(anion, 0.0)
        # 游离酸库存：阴离子被金属离子聘为反离子（Cu(NO3)2）后不算"游离酸"——
        # 分子化曲线只对游离酸浓度有意义（否则 T78 停点时曲线反而把残余
        # 质子抽回分子池）。阴离子池被消耗时库存同步缩水。
        free_an = min(m_an, max(He, 0.0))
        avail = m_mol + min(free_an, max(shadow - m_mol, 0.0))
        tot = min(shadow, avail)
        ledger[f"__tot_{acid}"] = tot
        if tot <= 0.0:
            continue
        target = _mol_fraction(acid, tot / V, T) * tot
        d = target - m_mol          # >0 缔合（抽阴离子+质子）；<0 再电离
        if d > 0.0:
            d = min(d, m_an, He / n_H if He > 0.0 else 0.0)
        else:
            d = -min(-d, m_mol)
        if abs(d) <= X_MIN:
            continue
        ledger[acid] = m_mol + d
        ledger[anion] = m_an - d
        He -= n_H * d                 # 缔合(d>0)吸走质子，再电离(d<0)释放
    return He


# ========================================================== 闭式 pH 快路径（v0.4.2 迭代 A）

def weak_species_set(T) -> frozenset:
    """所有可能影响 pH 的物种超集（T 级静态，一次构建）：
    pKa 酸/碱两侧 + OH^- 型 Ksp 阳离子（金属水解）+ beta_pka 派生配
    离子的反应物侧（滴定储备）。超集方向保守——闭式判据"在账物种与
    本集合交为空"时，完整路径的两侧滴定堆必空、分支 4 角色扫描必无
    命中（rolemap 键集 ⊆ 本超集），pH 退化为纯闭式解。"""
    s = getattr(T, "_weak_set", None)
    if s is None:
        names = set(T.pka_acid) | set(T.pka_base)
        for e in T.ksp:
            cat, an = e["pair"]
            if an == "OH^-":
                names.add(cat)
        for dc in build_derived(T):
            if not dc.meta.get("src", "").startswith("beta_pka:"):
                continue
            for sp in dc.r:
                if sp not in (H_ION, WATER):
                    names.add(sp)
        s = T._weak_set = frozenset(names)
    return s


def closed_pH(ledger: dict, H_excess: float, V: float, T, T_K: float):
    """无弱组分体系的闭式 pH：返回 (pH, vled, He_res) 或 None（含弱
    组分——走 estimate_state 完整路径）。bit 级等价依据：
      · 两侧滴定堆均空（在账弱组分物种集为空 ⊇ heap 判据）→
        _buffer_titration 返回 (None, H_excess, ledger 原身份)；
      · 分支 4 扫描无角色命中 → h_c=o_c=10^(-pKw/2) 不被推进，
        amph/buf 空，收尾退化为纯水；
      · 收尾/直读公式逐字符复刻（含 -log10(10.0**(-pKw/2)) 的浮点
        路径——与 pKw/2 直写可能有 1 ulp 差，不可化简）。
    量口径取分支 4 的 floor（1e-12·V，两侧堆的 X_MIN 更大）：漏判
    方向安全（该物种在完整路径同样被跳过）。"""
    weak = weak_species_set(T)
    floor_V = 1e-12 * V
    for sp, m in ledger.items():
        if m > floor_V and sp in weak:
            return None
    pKw = pKw_of(T_K)
    He = H_excess / V
    if He >= 1e-3:
        return max(-1.0, -log10(He)), ledger, H_excess
    if He <= -1e-3:
        return min(pKw + 1.0, pKw + log10(-He)), ledger, H_excess
    pH = -log10(10.0 ** (-pKw / 2))
    return pH, ledger, H_excess


def exact_proton_pH(ledger: dict, H_excess: float, V: float, T,
                    T_K: float) -> float | None:
    """B3 档的**精确质子条件 pH**；非 B3 档返回 None（呈现层冷路径专用）。

    判据（architecture §7 N-6 三条前提，全库分类验证过）：
      ① 账本**电中性自洽**：`|Σz·n + He| ≤ 1e-6`（投料本身不电中性的
         裸离子体系、以及带 c_H 的强酸条件体系都落在这一条之外）；
      ② **无储库物种**在场：无固相、无 Ksp-OH 水解阳离子——这些体系的
         pH 由 Ksp/逸度（账本之外的自由度）决定，账本内部的质子条件
         看不到它们（L08 型：机器 5.80 由 Cu(OH)₂ 的 Ksp 定，质子条件
         只会把残余 +0.001 电荷拉平到 pH 11，错的）；
      ③ 至少有一个**可再分配**的质子化族成员。
    满足三条时 `charge_pH` 是权威解（族逐级严格分布 + 水自电离 +
    固定离子电荷平衡，12/12 教科书锚点 + 80/80 条目自洽）。

    **只许在冷路径调用**（`final_pH` / 收敛探针 / 质量口径各一次）：
    分支 4 的两条启发式在 B3 档有实测误差（两性中点式 AB03 10.50 vs 9.84、
    Q05 1:1 缓冲对 10.50 vs pKa₁ 7.00、缓冲对加权 H88 7.246 vs 7.005），
    但精确解进二分探针要 +16% 墙钟（§7 O-2/O-3），故热路径维持启发式。
    """
    b3 = getattr(T, "_b3_static", None)
    if b3 is None:
        b3 = T._b3_static = (
            frozenset(T.solids) | frozenset(e["pair"][0] for e in T.ksp
                                            if e["pair"][1] == "OH^-"),
            frozenset(build_families(T)))
    res_set, fam_set = b3
    net = 0.0
    nfam = 0
    for sp, m in ledger.items():
        if m <= 0.0 or sp == WATER or sp.startswith("__"):
            continue
        if m > X_MIN and sp in res_set:
            # 痕量（≤ X_MIN）固相不构成储库：引擎全局把 X_MIN 以下的量视为
            # 噪声（D33 的 Al(OH)₃ 1.3e-7 曾因此挡住精确解，而 1e-7 的固相
            # 定不了 pH）。储库判据必须与引擎自己的痕量线一致。
            return None
        net += charge_of(sp) * m
        if sp in fam_set:
            nfam += 1
    if not nfam or abs(net + H_excess) > 1e-6:
        return None
    return charge_pH(ledger, V, T, T_K, fast=True)


def _full_speciation(ledger: dict, H_excess: float, V: float, T, T_K: float) -> tuple[dict, float]:
    """多级酸碱全形态分布 + 残余 He（仅供"惰性实现"判定）。与 estimate_state
    相同的分支结构，但滴定沿质子化梯走到底（VO3-→HVO3→VO2+ 一次调用完成），
    用于发现账本上不存在、却在当前酸度下真实存在的氧化还原物种。
    pH 估计不走此路——快照语义是 466 用例验证过的行为，两处各司其职。"""
    pKw = pKw_of(T_K)
    for acid in STRONG_MOLECULAR_ACIDS:
        # 分子态酸 ≥1M 才视为浓酸区制、跳过多级形态重排（连续形态分布下
        # 稀酸也有小量分子形态，不能一见分子就跳过——那是旧二元世界的语义）
        if ledger.get(acid, 0.0) / V >= 1.0:
            return ledger, H_excess
    _, He_res, vled = _buffer_titration(ledger, H_excess, V, T, pKw, multilevel=True,
                                        T_K=T_K)
    return vled, He_res

