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
            if b_rest > max(X_MIN, 1e-9 * V) and hb > 0.0:
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
            if a_rest > max(X_MIN, 1e-9 * V) and b > 0.0:
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
    return list(rows), cnt
