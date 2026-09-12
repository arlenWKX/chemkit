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
from math import log10, sqrt
from .core import pKw_of, _vant, charge_of
from .data import Tables
from .candidates import (Cand, logK_T, build_derived, WATER, H_ION, X_MIN,
                         _STRONG_ACID, _ksp_xy, STRONG_MOLECULAR_ACIDS)
from .normalize import _mol_fraction

# ========================================================== pH 估计器（§3.2，教科书近似）

def _buffer_titration(ledger: dict, H_excess: float, V: float, T, pKw: float,
                      multilevel: bool = False, T_K: float = 298.15,
                      cache: dict | None = None,
                      touch: frozenset | None = None) -> tuple:
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
    pH 非 None：缓冲对 Henderson 定 pH；pH=None 且残余≈0：落分支4（用虚拟账本）。"""
    if abs(H_excess) < 1e-12:
        return None, H_excess, ledger
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
        for dc in build_derived(T):
            if not dc.meta.get("src", "").startswith("beta_pka:"):
                continue
            nu_h = dc.r.get(H_ION, 0)
            if nu_h <= 0:
                continue
            comps = [s for s in dc.r if s not in (H_ION, WATER)]
            if len(comps) != 1:
                continue
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
        _ceff = {c: logK_T(dc, T_K) / nh for c, dc, nh in _beta_pka}
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
                if nxt is not None:
                    pushed.add(acid)
                    heappush(heap, (-_beff[acid], cnt, acid, nxt[2]))
                    cnt += 1
        if plateau is not None:
            pka, base, acid = plateau
            b_rest = ledger2.get(base, 0.0)
            hb = ledger2.get(acid, 0.0)
            if b_rest > _floor and hb > _floor:
                pH = pka + log10(b_rest / hb)
                return min(max(pH, -1.0), pKw + 1.0), he, ledger2
        return None, he, ledger2   # 全吸收（he≈0）→ 分支4；有残余 → 直读
    else:        # 弱酸吸收强碱：pKa 越小 Ka 越大，先中和
        he = -he
        heap, cnt = _titration_heap(ledger, _a_entries, cache, "a", touch,
                                    nominal=pKw + 2)
        if not heap:
            return None, -he, ledger
        ledger2 = dict(ledger)
        heappush = heapq.heappush
        plateau = None
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
            elif multilevel and take > 0.0 and base not in pushed:
                # 产物仍是酸（可再去质子化）→ 重新入堆（仅多级模式）
                nxt = _acids_map.get(base)
                if nxt is not None:
                    pushed.add(base)
                    heappush(heap, (_aeff[base], cnt, base, nxt[2]))
                    cnt += 1
        if plateau is not None:
            pka, acid, base = plateau
            a_rest = ledger2.get(acid, 0.0)
            b = ledger2.get(base, 0.0)
            if a_rest > _floor and b > _floor:
                pH = pka + log10(b / a_rest)
                return min(max(pH, -1.0), pKw + 1.0), -he, ledger2
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


def estimate_pH(ledger: dict, H_excess: float, V: float, T, T_K: float,
                cache: dict | None = None,
                touch: frozenset | None = None) -> float:
    return estimate_state(ledger, H_excess, V, T, T_K, cache, touch)[0]


def estimate_state(ledger: dict, H_excess: float, V: float, T, T_K: float,
                   cache: dict | None = None,
                   touch: frozenset | None = None) -> tuple[float, dict, float]:
    """返回 (pH, 滴定后的虚拟账本, 残余He)。虚拟账本是 pH 一致的自由形态分布；
    残余He是弱酸/弱碱储备吸收后仍未中和的游离强酸/强碱（酸碱平衡后的真实 He）。"""
    pKw = pKw_of(T_K)
    # 1) 强酸连续形态分布后，游离 H+ 全部由 He 记账（分子分数不贡献游离 H+，
    #    不再有"分子态浓酸"直读分支——pH 即 -log10(自由 H+) 的自然结果）
    # 2) 缓冲滴定（质子条件近似）：游离强酸/强碱先被在账弱碱/弱酸储备按强度
    #    顺序吸收；被全吸收则由最后缓冲对的 Henderson 式定 pH；
    #    全吸收但无有效缓冲对 → 用残余 He（≈0）落分支 4，而非原 He 直读
    tit, He_res, ledger = _buffer_titration(ledger, H_excess, V, T, pKw,
                                            T_K=T_K, cache=cache, touch=touch)
    if tit is not None:
        return tit, ledger, He_res
    He = He_res / V
    if He >= 1e-3:
        return max(-1.0, -log10(He)), ledger, He_res
    if He <= -1e-3:
        return min(pKw + 1.0, pKw + log10(-He)), ledger, He_res
    # 4) 缓冲/弱酸弱碱区：取各来源贡献最大者（在滴定后的虚拟账本上评估）
    h_c = 10.0 ** (-pKw / 2)
    o_c = h_c

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
        acids_map: dict = {}   # acid -> Ka 或 _STRONG_ACID 哨兵
        for acid, entries in T.pka_acid.items():
            if acid in amph_eligible:
                continue
            e1 = min(entries, key=lambda e: e["pka"])   # 第一级
            acids_map[acid] = _STRONG_ACID if e1["pka"] <= 0 else 10.0 ** (-_pka1(e1, T_K))
        bases_map: dict = {}   # base -> Kb
        for base, entries in T.pka_base.items():
            if base in T.solids or base == WATER or base in amph_eligible:
                continue
            e1s = [e for e in entries if e["n"] == 1] or entries
            pka = max(_pkapp(e, T_K) for e in e1s)           # 最强一级共轭酸
            bases_map[base] = 10.0 ** (pka - pKw)
        hyd_map: dict = {}    # cat -> (Kh 每阳离子, qc)
        for e in T.ksp:
            cat, an = e["pair"]
            if an != "OH^-":
                continue
            _x, _y = _ksp_xy(e)
            hyd_map[cat] = (10.0 ** ((_pksp(e, T_K) - _y * pKw) / _x),
                            charge_of(cat))
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
        sc = (amph_eligible, amph_pH, acids_map, bases_map, hyd_map, conj,
              rolemap)
        est_cache[T_K] = sc
    amph_eligible, amph_pH, acids_map, bases_map, hyd_map, conj, rolemap = sc
    # 单遍扫描在账物种：合并原 acids_l/bases_l/hyd_l/amph/buf 五个独立循环。
    # max/sum 可换序，h_c/o_c/amph/buf 的最终值与原实现等价。
    amph: list = []
    buf: list = []
    _floor_V = 1e-12 * V
    _role_get = rolemap.get
    for sp, m in ledger.items():
        if m <= _floor_V:
            continue
        c = m / V
        role = _role_get(sp)
        if role is None:
            continue
        Ka, Kb, conj_pair, Kh_qc, amph_v = role
        if Ka is not None:
            if Ka is _STRONG_ACID:
                h_c = max(h_c, c)
            else:
                h_c = max(h_c, (-Ka + sqrt(Ka * Ka + 4 * Ka * c)) / 2)
        if Kb is not None:
            if Kb >= 1.0:                                # 水解近完全（S2-、C2^2- 等）
                o_c = max(o_c, c)
            else:
                o_c = max(o_c, (-Kb + sqrt(Kb * Kb + 4 * Kb * c)) / 2)
            # 共轭缓冲对（仅碱在账时检查其共轭酸是否也在账）
            if conj_pair is not None:
                acid_conj, pka_c = conj_pair
                ca = ledger.get(acid_conj, 0.0) / V
                if ca > 1e-12:
                    buf.append((pka_c + log10(c / ca), min(c, ca)))
            continue
        if Kh_qc is not None:
            Kh, qc = Kh_qc
            if qc == 1:
                # 一价金属水解到底即纯固相（活度 1），无累积共轭碱：
                # cat + H2O → ½M2O + H+ 给出 h = Kh·c；套弱酸二次式
                # 会把 Ag+ 类高估 ~1/sqrt(Kh·c) 倍（D26：Ag+ 被估成 pH 3
                # 的酸，驱动铬酸根幻影质子化死循环）。多价金属分步水解
                # 经可溶羟基中间体，实测行为近弱酸二次式，保持不变。
                h_c = max(h_c, Kh * c)
            else:
                h_c = max(h_c, (-Kh + sqrt(Kh * Kh + 4 * Kh * c)) / 2)
            continue
        if amph_v is not None and c > 1e-6:
            amph.append((c, amph_v))
    # 两性物种（HCO3-、HS-、H2PO4- 等）：pH ≈ (pKa_酸 + pKa_共轭酸)/2，
    # 其浓度远大于其他酸碱贡献时以两性平衡为准（NaHCO3 溶液 pH≈8.3）
    if amph:
        c_a, pH_a = max(amph, key=lambda t: t[0])
        if c_a > 100.0 * max(h_c, o_c):
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
            return min(max(pH_buf, -1.0), pKw + 1.0), ledger, He_res
    pH = -log10(h_c) if h_c >= o_c else pKw + log10(o_c)
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

