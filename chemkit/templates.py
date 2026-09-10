"""chemkit.templates：候选宇宙（模板构建与枚举）。

职能（四级缓存体系，architecture.md §性能）：
    · _redox_pair_static：温度无关静态层（配平/电子数/dH/vs 闸门/
      恒慢标记预解析）——每数据表一次
    · _redox_templates：温度过滤层（T_min/慢速温度规则/halate/pKw
      动态阈值）——按 (T_K, kinetics) 缓存；93.5% 无温度规则电对
      预组装双 kinetics 成品
    · _build_static_cands：静态候选（配位/解离/沉淀/溶解），按
      (T, T_K) 缓存，预构建 Cand 挂模板
    · enumerate_candidates：存在性/pH/浓度闸门过滤 + 三层枚举 memo
      （完整结果 / idxs 收集 / gate·形态·慢方向）

依赖 candidates/core/data。
"""
from __future__ import annotations
from functools import reduce
from math import gcd, log10
from .core import (elements_of, charge_of, K_NERNST_298, PKW_298,
                   k_nernst, _vant, pKw_of)
from .data import Tables, henry_of, _half_balance
from .candidates import (Cand, logK_T, build_derived, _bal, _bal_fast,
                         _half_pairs, _redox_mix_ok, _half_scale, _ksp_xy,
                         WATER, H_ION, X_MIN, ACT_FLOOR, SAT_SKIP,
                         _NONMETAL_ELEMS, WATER_FIRST_METALS,
                         NONMETAL_SOLIDS, HALATE_DISP_T, HALATE_BASE_PH,
                         P_EXT_KPA)

# ========================================================== 候选枚举

def _redox_pair_static(T):
    """电对配对 (a,b) 的**温度无关**静态层，缓存于 T._redox_pair_static。

    配平（_bal_fast/_bal/配体补充）、电子数折算、dH 组装、same_elem、
    solid_acid、E0 选择性闸门（vs_red_*）、closed_with_red 掩蔽与全部
    "恒慢"动力学标记只依赖数据表——多温度场景（跨温度用例、非恒温 ΔT
    重判、温度扫描）只在冷启动计算一次；T_min/slow_below/halate 歧化
    温度/pKw 动态 pH 阈值等温度相关规则以轻量数据记录，由
    _redox_templates 按 (T_K, kinetics) 过滤产出模板（毫秒级）。

    项：(a, b, variants, dE, same_elem, solid_acid, slow_static,
        deferred_raw, sd_ph_hi, closed_red, gates_only, slow_T,
        halate_rule, h2o_oh_min)
      variants        [(a_ox, r, pr, n, dH_full)]（NO3- 双变体含 HNO3）
      slow_static     恒慢标记（温度无关）
      deferred_raw    b 电对 derived_from（晶格氧化让位档，kinetics 时生效）
      sd_ph_hi        [(方向, pH 上限)] ox_pH_max（kinetics 时生效）
      closed_red      closed_with_red 掩蔽命中（kinetics 时生效）
      gates_only      [(T_min, partner_red, reds)] only_vs_red 闸门
                      （T_K < T_min 且对方 red ∈ reds 时跳过）
      slow_T          [(T_thr, species, except_reds)] 温度阈值慢（T_K < T_thr
                      且对方 red ∉ except_reds 时 slow；species=None 恒命中）
      halate_rule     [(方向, hx)] 卤酸歧化慢方向（T_K < HALATE_DISP_T[hx]
                      时按 HALATE_BASE_PH 落 sd_static/sd_ph）
      h2o_oh_min      a 电对 h2o_red_oh_min（b.ox 为 O_2 时 sd_ph 附加
                      pH 阈值 pKw(T_K)+log10(oh_min)，随温度变）
    """
    hit = getattr(T, "_redox_pair_static", None)
    if hit is not None:
        return hit
    # 固相电对的还原侧 **阳离子/中性** 物种受酸限量约束——但仅限金属自
    # 腐蚀族（电对共享元素为金属：ox 是该金属的氧化膜/单质，酸限量建模
    # 金属-水界面腐蚀受膜保护）。晶格阴离子电对（S/PbS：共享元素为 S）
    # 不是金属腐蚀，其氧化耗 H+ 在中性悬浮液中照常进行（PbS + H2O2 ->
    # PbSO4 是真实反应），不得纳入酸限量。阴离子还原剂（I-/S2- 等）同前。
    solid_metals = {c["red"] for c in T.couples
                    if c["ox"] in T.solids and charge_of(c["red"]) >= 0
                    and ((set(elements_of(c["ox"])) & set(elements_of(c["red"])))
                         - _NONMETAL_ELEMS)}
    # 固相 -> 阳离子映射（水优先金属规则的固相形态判据）
    ksp_cat = {e["solid"]: e["pair"][0] for e in T.ksp}
    # 预计算：每电对的元素集（歧化/归中判定、固相配体提取共用；
    # 免去 26k 对 × 每对重建两个 set 的重复开销）
    _cs = T.couples
    _eas = [set(elements_of(c["ox"])) | set(elements_of(c["red"]))
            for c in _cs]
    # 预计算：固相 -> 晶格 Ksp 单元（[(cat, an), ...]，通常 1-3 项；
    # 配体提取时免去对 ksp_by_pair 226 项全扫描）
    _solid_cells: dict[str, list] = {}
    for (cat, an), cell in T.ksp_by_pair.items():
        _solid_cells.setdefault(cell["solid"], []).append((cat, an))
    pairs: list = []
    for i, a in enumerate(_cs):
        _ea = _eas[i]
        for j, b in enumerate(_cs):
            if j == i:
                continue
            # ---- 静态选择性闸门（vs_red_only/block/E_max 与温度无关，
            # 此处直接过滤；closed_with_red 掩蔽仅记录，kinetics 层判定）
            closed_red = False
            skip = False
            for c, partner in ((a, b), (b, a)):
                g = c.get("gate")
                if g:
                    if "vs_red_only" in g and partner["red"] not in g["vs_red_only"]:
                        skip = True
                        break
                    if "vs_red_block" in g and partner["red"] in g["vs_red_block"]:
                        skip = True
                        break
                    if "vs_red_E_max" in g:
                        # 配对级判定用 partner 电对自身的 E0（物种级
                        # redox_red_E 取最低值会把配合物形态的电对串扰到
                        # 裸金属电对——[Cu(CN)2]^-/Cu -0.90 曾把 red="Cu"
                        # 的最低电位拉破 vs_red_E_max -0.6 闸门，使 N2O
                        # 路径对裸 Cu 误开放、NO 路径被 vs_red_block_E_max
                        # 误拦截）；本对的氧化/还原能力由本对 E0 唯一决定
                        e_red = partner["E0"]
                        if e_red > g["vs_red_E_max"]:
                            skip = True
                            break
                    if "vs_red_block_E_max" in g:
                        e_red = partner["E0"]
                        if e_red <= g["vs_red_block_E_max"]:
                            skip = True
                            break
                kin = c.get("kinetics")
                if kin and partner["red"] in kin.get("closed_with_red", ()):
                    closed_red = True
            if skip:
                continue
            # only_vs_red 温度闸门预解析（T_min 对指定还原剂生效，
            # 如 MnO2 对 Cl- 需加热）
            gates_only = []
            for c, partner in ((a, b), (b, a)):
                g = c.get("gate") or {}
                if "T_min" in g and "only_vs_red" in g:
                    gates_only.append((g["T_min"], partner["red"],
                                       frozenset(g["only_vs_red"])))
            # 共享氧化形或还原形的电对配对：电子 bookkeeping 净零，净反应实为
            # 酸碱/形态转换（如 HClO/Cl- ⊗ ClO-/Cl- 净得 HClO->ClO-+H+），
            # 属 pKa 模块管辖，作为 redox 枚举会产生幻影自发通道
            if a["red"] == b["red"] or a["ox"] == b["ox"]:
                continue
            # 同元素电对族（歧化/归中，如 Cu2+/Cu+ ⊗ Cu+/Cu、
            # MnO4-/MnO4^2- ⊗ MnO4^2-/MnO2）：溶液中快速氧化还原平衡，
            # 允许逆向存在性（a.red + b.ox 在账即可参与）——中间态
            # 的排空通道不因正向反应物耗尽而消失（公理1：候选是平衡）。
            # 跨元素对维持正向存在性过滤，避免逆向闸门语义爆炸。
            _eb = _eas[j]
            # 仅限共享金属元素的电对族（Cu/Hg/Mn/Fe/V/Cr…）：非金属
            # 多价链（N/S/卤素）的逆向通道受浓度/气体逸出/温度强约束，
            # 逆向存在性放开会冲垮正向闸门语义（经验：238 对中 N/S/卤素
            # 占 166 对，全部爆炸来自它们）
            same_elem = bool((_ea & _eb) - _NONMETAL_ELEMS)
            # ---- 恒慢标记（温度无关部分；T 阈值规则拆到 slow_T）
            slow_static = False
            slow_T: list = []
            if b["red"] in a.get("slow_with_red", ()):
                slow_static = True
            if a["red"] in b.get("slow_with_red", ()):
                slow_static = True
            # 该电对作还原剂侧（其 red 被氧化）恒慢：如水被阳极氧化为 H2O2
            if b.get("slow_as_reductant"):
                slow_static = True
            # 指定氧化剂组合恒慢：如 SO3^2- 还原 H+/水析 H2（亚硫酸盐溶液动力学稳定）
            if b["ox"] in a.get("slow_with_ox", ()):
                slow_static = True
            if a["ox"] in b.get("slow_with_ox", ()):
                slow_static = True
            # S(IV)（SO2/亚硫酸）作还原剂把金属离子还原为金属单质：水溶液中动力学
            # 封闭——SO2 还原只到中间价态（如 Cu2+->Cu+），从不在水溶液析出金属
            if (b["red"] == "SO_2" and a["red"] in T.solids
                    and a["red"] not in NONMETAL_SOLIDS):
                slow_static = True
            # H2 作还原剂在水溶液中常温恒慢（需催化剂或高温加热才表现还原性，
            # 如 H2 还原 CuO 需加热；水溶液中 H2 不还原 Cu2+/Fe3+ 等）
            if b["red"] == "H_2":
                slow_static = True
            # 单质硫作还原剂在水溶液中室温恒慢（硬事实：升华硫在水中可无限期
            # 稳定存在；S+热浓碱歧化、S+沸腾浓硝酸氧化为硫酸、S+Na2SO3 煮沸制
            # 硫代硫酸钠、橡胶硫化 ~140°C 均需加热）。逆向（S 作氧化剂被还原，
            # 如 Hg+S 室温研磨汞珠回收）不受影响。例外：H2O2 室温即可氧化硫
            # （浓 H2O2 与硫粉反应剧烈放热，分析化学消化法即用此）。解锁温度
            # 取 ~97°C（沸腾稀溶液仍可稳定 S——稀硝酸溶 CuS 止步于单质硫即此事实）。
            if b["red"] == "S" and a["ox"] != "H_2O_2":
                slow_T.append((370.0, None, ()))
            # 活泼金属（Li/Na/K/Rb/Cs/Ca/Sr/Ba）在水溶液中优先还原水而非金属物种
            # （教材规则：钠投入盐溶液只与水反应，再生成氢氧化物沉淀；Be/Mg 可直接置换）。
            # 覆盖两类氧化剂形态：金属阳离子、含金属阳离子的固相氢氧化物/氧化物
            if (b["red"] in WATER_FIRST_METALS and a["ox"] != H_ION
                    and (charge_of(a["ox"]) > 0
                         or (a["ox"] in ksp_cat
                             and charge_of(ksp_cat[a["ox"]]) > 0))):
                slow_static = True
            # ---- 温度阈值慢（slow_below 及 O2 阳极特例）
            for c, other, is_a in ((a, b, True), (b, a, False)):
                sb = c.get("slow_below")
                if not sb:
                    continue
                if c["ox"] == "O_2" and not is_a:
                    # O2 析出（该电对作 b 方、O2 为产物 = 水被氧化）：4e- 阳极
                    # 过程，温度不解锁（与 O2 作氧化剂的 slow_below 区分开）。
                    # 例外：对方电对声明 h2o_red_oh_min（如浓碱中 MnO4-/MnO4^2-
                    # 氧化水制锰酸根——实验硬事实）。此时跳过本电对全部
                    # slow_below 规则（该规则针对 O2 作氧化剂方向，作还原剂
                    # 侧不适用），改由模板构建期的 pH 动态闸门控制（稀溶液仍慢，
                    # 见 h2o_oh_min 的 sd_ph）。
                    if not other.get("h2o_red_oh_min"):
                        slow_static = True
                    continue
                slow_T.append((sb, other["red"],
                               tuple(c.get("slow_except_red", ()))))
            # ---- 卤素歧化到卤酸根方向感知的慢标记（教材规则：冷稀碱 -> ClO-，
            # 热碱 -> ClO3-；歧化经 XO- 中间体，其进一步歧化室温慢——漂白液
            # 室温稳定即此动力学事实；Br/I 为碱催化解锁，见 HALATE_DISP_T/
            # HALATE_BASE_PH）。同一平衡有两种枚举取向：歧化取向正向慢、归中
            # 取向逆向慢；归中方向（IO3-+I-）与异种氧化（Cl2+I2 -> IO3-）始终快。
            halate_rule: list = []
            for _d, _hx in ((1, b["ox"]), (-1, a["ox"])):
                _src_ok = (a["ox"] == b["red"]) if _d == 1 else (b["ox"] == a["red"])
                if _src_ok and HALATE_DISP_T.get(_hx) is not None:
                    halate_rule.append((_d, _hx))
            # 水作还原剂析 O2 的浓碱解锁闸门（稀溶液动力学封闭——pH <
            # pKw+log10(oh_min) 时正向慢。实验事实：KMnO4 只在浓碱中可观察
            # 地氧化水（锰酸钾制备），稀溶液中分解慢到可忽略）
            h2o_oh_min = (a.get("h2o_red_oh_min")
                          if b["ox"] == "O_2" else None)
            # 电对声明 ox_pH_max：其 ox 被还原的方向（a 侧正向 d=1 / b 侧逆向
            # d=-1）在 pH 高于阈值时动力学封闭。硬事实例：绿色锰酸根在浓碱中
            # 亚稳（热力学预测部分歧化，实验上冻结）——该冻结必须方向感知：
            # 歧化的逆向枚举（归中取向）不经过 a 侧 gate，正向闸门拦不住。
            sd_ph_hi: list = []
            for _c, _d in ((a, 1), (b, -1)):
                _pm = _c.get("ox_pH_max")
                if _pm is not None:
                    sd_ph_hi.append((_d, _pm))
            # 派生晶格电对（derived_from，如 S/CuS）作还原剂 = 硫化物晶格被
            # 氧化：破晶格电子转移，动力学上慢于同固体的离子交换（复分解），
            # 但并非不发生。实验事实：浮选活化 ZnS+Cu2+ -> CuS+Zn2+ 是快速
            # 交换而非氧化（交换优先）；无快通道时氧化照常进行（PbS 被
            # H2O2 氧化、CuS 溶于热硝酸、硫化矿生物浸出）。标 deferred：
            # 有快候选竞争时让位，无快候选时执行——区别于 slow（恒慢，
            # 只标注不执行）。
            deferred_raw = bool(b.get("derived_from"))
            # NO3- 电对双变体：离子态 / 浓硝酸分子态（运行时按存在性选择）
            variants = []
            for a_ox in ([a["ox"], "HNO_3"] if a["ox"] == "NO_3^-" else [a["ox"]]):
                pool = [WATER, H_ION]
                # 快路径：半反应 lcm 组合（same_elem 族零空间欠定，仍走 sympy）
                bal = (None if same_elem else
                       _bal_fast(a_ox, a["red"], b["ox"], b["red"], T))
                if bal is None:
                    bal = _bal([a_ox, b["red"]], [a["red"], b["ox"]], pool)
                if bal is None:
                    # 配合物电对的配体参与配平（[AuCl4]-/Au 经王水溶解：
                    # Au + 4Cl- -> [AuCl4]- + 3e，Cl- 必须可入账）；先按常规
                    # 池配平，失败才补配体——不影响既有模板
                    lig = []
                    for sp in (a["ox"], a["red"], b["ox"], b["red"]):
                        be = T.beta_by_complex.get(sp)
                        if be and be["ligand"] not in lig:
                            lig.append(be["ligand"])
                        # 固相电对的晶格旁观离子参与配平（AgI/Ag：Ag + I- ->
                        # AgI + e，I- 必须可入账；与配合物配体同型，均为电对
                        # 半反应的固有参与物种）。取与电对另一方不共享元素的
                        # 晶格离子：AgI/Ag 共享 Ag → 取 I-（金属氧化还原，
                        # 阴离子旁观）；S/PbS 共享 S → 取 Pb2+（晶格阴离子
                        # 被氧化，阳离子释放）。双向皆共享（Fe3O4 类）时不取。
                        if sp in _solid_cells:
                            _coup = a if sp in (a["ox"], a["red"]) else b
                            _oth = (_coup["red"] if sp == _coup["ox"]
                                    else _coup["ox"])
                            _shared = (set(elements_of(sp))
                                       & set(elements_of(_oth)))
                            for cat, an in _solid_cells[sp]:
                                for _ion in (cat, an):
                                    if any(_e2 not in _shared for _e2
                                           in elements_of(_ion)) \
                                            and _ion not in lig:
                                        lig.append(_ion)
                    if lig:
                        pool = [WATER, H_ION] + lig
                        bal = _bal([a_ox, b["red"]], [a["red"], b["ox"]], pool)
                if bal is None:
                    continue
                # 同元素电对族：欠定零空间中的配平必须能分解为两个半反应
                # 的干净组合（残余仅限非氧化还原旁观），否则为混入第三电对
                # 的幻影（见 _redox_mix_ok）
                if same_elem and not _redox_mix_ok(
                        a, b, a_ox, dict(bal["reactants"]),
                        dict(bal["products"]), pool,
                        (_ea & _eb) - _NONMETAL_ELEMS):
                    continue
                r, pr = dict(bal["reactants"]), dict(bal["products"])
                # 实际转移电子数（配平已约简，lcm 会高估，如 I2 歧化 n_a=2/n_b=5
                # 的 3I2 方程实际转 5e 而非 lcm=10）：按半反应系数折算；
                # 归中体系产物侧物种共享，须从反应物侧折算。
                # 氧化还原+沉淀/酸碱合并的混合候选中系数含旁观计量（如
                # 2Fe(OH)3+3Fe2+ -> 3Fe(OH)2+2Fe3+ 实际只转 2e，1 个 Fe2+ 仅沉淀），
                # 任一单物种折算都是电子数的上界 → 取四向折算的最小正值（最紧上界）
                _opts = []
                if a["red"] != b["ox"]:
                    _opts.append((_half_scale(a), pr.get(a["red"], 0)))
                _opts.append((_half_scale(b), r.get(b["red"], 0)))
                _opts.append((a["n"], r.get(a_ox, 0)))
                _opts.append((b["n"], pr.get(b["ox"], 0)))
                _vals = [sc * nu for sc, nu in _opts if nu]
                n = min(_vals) if _vals else 0
                if n <= 0:
                    continue
                if H_ION in r and H_ION in pr:
                    continue
                # van't Hoff：全反应 ΔH = ν_a·dH_a − ν_b·dH_b（ν=n/dH_n 为半反应
                # 折算系数；H+ 的 ΔHf 约定为 0，H+ 项自动含于半反应 dH）。
                # 任一电对缺 dH → None，回退 Nernst k(T) 缩放（ΔS≈0 近似）。
                dH_full = None
                if "dH" in a and "dH" in b:
                    dH_full = n / a["dH_n"] * a["dH"] - n / b["dH_n"] * b["dH"]
                variants.append((a_ox, r, pr, n, dH_full))
            if variants:
                # 预组装：无温度规则（gates_only/slow_T/halate/h2o_oh/
                # 静态 T_min）的对，任意温度下模板不变——直接缓存两种
                # kinetics 形态的成品模板，温度层 O(1) 引用（哨兵：
                # None=该模式跳过，False=动态处理，tuple=成品）。
                # 进一步：sd_ph_hi 为空时（约 99% 模板）slow_dirs 恒空，
                # meta 全常量——预构建 Cand，枚举时直接 out.append。
                _g = a.get("gate") or {}
                is_dyn = bool(gates_only or slow_T or halate_rule
                              or h2o_oh_min is not None
                              or ("T_min" in _g and "only_vs_red" not in _g))
                _dE = a["E0"] - b["E0"]
                _acid = b["red"] in solid_metals
                if is_dyn:
                    pre_kin = pre_nokin = False
                else:
                    _mk = lambda slow, defer, sd_hi: (
                        None if sd_hi else
                        tuple(Cand("redox", v[1], v[2], 0.0, 0.0,
                                   redox=(v[3], _dE), dH=v[4],
                                   meta={"slow": slow, "slow_dirs": frozenset(),
                                         "ox_couple": a["ox"],
                                         "acid_limited": _acid,
                                         "deferred": defer})
                              for v in variants))
                    pre_kin = (None if (a.get("ox_inert") or closed_red)
                               else (a, b, variants, _dE,
                                     slow_static, frozenset(), [],
                                     _acid, same_elem,
                                     tuple(sd_ph_hi), deferred_raw,
                                     _mk(slow_static, deferred_raw, sd_ph_hi)))
                    pre_nokin = (a, b, variants, _dE,
                                 False, frozenset(), [],
                                 _acid, same_elem,
                                 [], False, _mk(False, False, None))
                pairs.append((a, b, variants, _dE, same_elem,
                              _acid, slow_static,
                              deferred_raw, tuple(sd_ph_hi), closed_red,
                              tuple(gates_only), tuple(slow_T),
                              tuple(halate_rule), h2o_oh_min,
                              pre_kin, pre_nokin))
    T._redox_pair_static = pairs
    return pairs


def _redox_templates(T_K: float, T, kinetics: bool = True) -> list:
    """(T, T_K, kinetics) 下 redox 候选的静态模板，按 (T_K, kinetics)
    字典缓存于 T._redox_tmpl。

    重活（配平/电子数/静态闸门）在 _redox_pair_static 与温度无关地缓存
    ——本函数仅做温度过滤（T_min/only_vs_red/slow_below/halate/pKw 动态
    阈值），单个新温度点构建为毫秒级。

    主循环是不动点迭代：电对配对、配平、电子数折算、慢标记、温度闸门在
    同一 (T, T_K, kinetics) 下是"定值"，只算一次；存在性/pH 闸门/产物
    形态随迭代动态过滤（见 enumerate_candidates §4.1）。模板顺序保持原
    (a, b) 嵌套顺序以维持候选序列（tie-break 依赖枚举顺序）。

    模板项：(a, b, variants, dE, slow, sd_static, sd_ph, acid_limited,
            same_elem, sd_ph_hi, deferred)
    variants：[(a_ox, r, pr, n)]——常规电对单变体；NO3- 电对附加 HNO3
    分子态变体（浓硝酸），运行时按存在性二选一。
    """
    cache = getattr(T, "_redox_tmpl", None)
    if cache is None:
        cache = T._redox_tmpl = {}
    hit = cache.get((T_K, kinetics))
    if hit is not None:
        return hit
    pairs = _redox_pair_static(T)
    tmpls = []
    append = tmpls.append
    for item in pairs:
        # 预组装快路径：无温度规则的对直接引用成品模板
        # （哨兵：None=该模式跳过，False=动态处理）
        pre = item[14] if kinetics else item[15]
        if pre is not False:
            if pre is not None:
                append(pre)
            continue
        (a, b, variants, dE, same_elem, solid_acid, slow_static,
         deferred_raw, sd_ph_hi, closed_red, gates_only, slow_T,
         halate_rule, h2o_oh_min) = item[:14]
        if kinetics:
            # ox_inert：氧化形永不参与电子转移的惰性物种（动力学宣告）
            if a.get("ox_inert"):
                continue
            # T_min（无 only_vs_red）为静态温度闸门
            g = a.get("gate") or {}
            if "T_min" in g and "only_vs_red" not in g and T_K < g["T_min"]:
                continue
            # only_vs_red：该电对的 T_min 闸门只对指定还原剂生效
            # （MnO2 对 Cl- 需加热）
            if any(T_K < tmin and sp in reds
                   for tmin, sp, reds in gates_only):
                continue
            # closed_with_red：配位掩蔽闸门（化学事实：强配位剂/沉淀剂
            # 优先与氧化剂离子形成配合物/沉淀，redox 通道不开放）。典型：
            # Hg2+ + 4I- → [HgI4]2- (logβ4=30) 优先于 Hg2+ 还原 I-；
            # Ag+ + 2CN- → [Ag(CN)2]- (logβ2=21) 优先于 Ag+ 还原 CN-。
            # kinetics=False（纯热力学基线）时让位：redox 与配位按平衡竞争。
            if closed_red:
                continue
            slow = slow_static or any(
                T_K < thr and (sp is None or sp not in exc)
                for thr, sp, exc in slow_T)
            # 派生晶格电对的让位档
            deferred = deferred_raw
            sd_static = frozenset()
            sd_ph: list = []
            # 卤酸歧化慢方向（温度阈值 + 碱催化解锁阈值）
            for _d, _hx in halate_rule:
                if T_K < HALATE_DISP_T[_hx]:
                    _ph = HALATE_BASE_PH.get(_hx)
                    if _ph is None:
                        sd_static |= {_d}
                    else:
                        sd_ph.append((_d, _ph))
            # 浓碱解锁阈值随 pKw(T) 漂移
            if h2o_oh_min is not None:
                sd_ph.append((1, pKw_of(T_K) + log10(h2o_oh_min)))
        else:
            slow = False
            deferred = False
            sd_static = frozenset()
            sd_ph = []
            sd_ph_hi = []
        no_sd = not (sd_static or sd_ph or sd_ph_hi)
        pre_cands = (None if not no_sd else
                     tuple(Cand("redox", v[1], v[2], 0.0, 0.0,
                                redox=(v[3], dE), dH=v[4],
                                meta={"slow": slow, "slow_dirs": frozenset(),
                                      "ox_couple": a["ox"],
                                      "acid_limited": solid_acid,
                                      "deferred": deferred})
                           for v in variants))
        append((a, b, variants, dE, slow, sd_static, sd_ph,
                solid_acid, same_elem, sd_ph_hi, deferred, pre_cands))
    # a_ox -> 模板下标：每迭代只检查氧化形在账的模板（原实现每迭代全表扫
    # 数千模板，反而慢于旧的双重循环）；下标排序后遍历保持原 (a,b) 嵌套序。
    # by_ared：同元素族的 a.red -> 下标，供逆向存在性预取（歧化排空通道）。
    by_aox: dict = {}
    by_ared: dict = {}
    for i, tpl in enumerate(tmpls):
        by_aox.setdefault(tpl[0]["ox"], []).append(i)
        if tpl[8]:
            by_ared.setdefault(tpl[0]["red"], []).append(i)
    T._redox_tmpl[(T_K, kinetics)] = (tmpls, by_aox, by_ared)
    return tmpls, by_aox, by_ared


def _oxide_dissolve_info(sp: str, T):
    """识别 form:solid 的二元金属氧化物 M_xO_y，返回 (阳离子, x, n, y, pKsp)；
    否则 None。要求 x·n = 2y（单一价态），且对应氢氧化物在 Ksp 表中。"""
    ex = T.ex.get(sp)
    if ex is None or ex.get("form") != "solid":
        return None
    elems = elements_of(sp)
    if len(elems) != 2 or "O" not in elems:
        return None
    (metal, x), (_, y) = sorted(elems.items(), key=lambda kv: kv[0] == "O")
    if 2 * y % x:
        return None
    n = 2 * y // x
    cat = f"{metal}^+" if n == 1 else f"{metal}^{{{n}+}}"
    pKsp = getattr(T, "_oh_ksp", {}).get(cat)
    if pKsp is None:
        return None
    return cat, x, n, y, pKsp


def _build_static_cands(T_K: float, T, pKw: float) -> list:
    """构建 §4.2/4.3/4.3b/4.4 静态候选表（质子转移、沉淀/溶解、氧化物酸溶、
    配位/解离），返回 [(Cand, 需求集 frozenset)]，按 (T, T_K) 缓存于
    T._cand_static。需求集排除 WATER/H_ION（二者恒在账），每迭代只做
    req ⊆ present 子集过滤。顺序与原内联实现一致：质子 → 沉淀 → 氧化物 → 配位。
    """
    out: list = []

    # ---- 4.2 质子转移（pKa ≤ 0 强酸条目不参与，规范化已拆解）
    # 配平与 H+ 挂侧检查静态（只依赖数据表），缓存于 T；存在性每迭代过滤
    if getattr(T, "_proton_static", None) is None:
        ps = []
        for e in T.pka:
            if e["pka"] <= 0:
                continue
            acid, base = e["acid"], e["base"]
            pka_e = e["pka"]
            if acid in T.solids:
                # 固相酸的脱质子 = 溶解 + 电离两步之和，表观 pKa 必须含溶解度：
                # pKa_eff = pKsp(H+·酸根) − 下游脱质子链 pKa 之和。
                # 直接用溶液 pKa 会让固体酸越过溶解度无限"再溶解"
                # （H2SiO3 凝胶在中间 pH 被 pKa 9.8 的幻影解离吃光；
                # 正确值 pKa_eff = 24.5 − 12.0 = 12.5，即 10^-12.5 的
                # K = h·[HSiO3-] 恰好给出真实溶度 ~0.03 mol/L @pH 11）
                adj = _solid_acid_pka(base, T)
                if adj is not None:
                    pka_e = adj
            e = dict(e)
            e["pka"] = pka_e
            bd = _bal([acid], [base, H_ION], [WATER])
            diss = (dict(bd["reactants"]), dict(bd["products"])) \
                if bd is not None and bd["products"].get(H_ION, 0) > 0 else None
            bp = _bal([base, H_ION], [acid], [WATER])
            prot = (dict(bp["reactants"]), dict(bp["products"])) \
                if bp is not None and bp["reactants"].get(H_ION, 0) > 0 else None
            ps.append((acid, base, e["pka"], diss, prot, e.get("dH")))
        T._proton_static = ps
    for acid, base, pka, diss, prot, dH_ion in T._proton_static:
        if diss is not None:
            out.append((Cand("proton", diss[0], diss[1], -pka, dH=dH_ion,
                             meta={"dir": "diss"}),
                        frozenset((acid,))))
        if prot is not None:
            out.append((Cand("proton", prot[0], prot[1], pka,
                             dH=-dH_ion if dH_ion is not None else None,
                             meta={"dir": "prot"}),
                        frozenset((base,))))

    # ---- 4.3 沉淀 / 溶解（氢氧化物用 H+ 正则形，与阳离子水解同一平衡）
    for e in T.ksp:
        cat, an = e["pair"]
        solid, pK = e["solid"], e["pKsp"]
        n_cat, n_an = _ksp_xy(e)
        # dH 不含水电离部分（OH- 型固体的水电离项由 pkw_coeff 通道承载），
        # 沉淀/溶解均直接取 ∓溶解焓
        dH_sol = e.get("dH")
        dH_p = -dH_sol if dH_sol is not None else None
        if an == "OH^-":
            out.append((Cand("precip", {cat: n_cat, WATER: n_an},
                             {solid: 1, H_ION: n_an},
                             pK - n_an * pKw, -n_an, dH=dH_p,
                             meta={"solid": solid}),
                        frozenset((cat,))))
            out.append((Cand("dissolve", {solid: 1, H_ION: n_an},
                             {cat: n_cat, WATER: n_an},
                             n_an * pKw - pK, n_an, dH=dH_sol,
                             meta={"solid": solid}),
                        frozenset((solid,))))
        else:
            r: dict = {cat: n_cat}
            if n_an:
                r[an] = n_an
            out.append((Cand("precip", r, {solid: 1}, pK, dH=dH_p,
                             meta={"solid": solid}),
                        frozenset(r)))
            pr: dict = {cat: n_cat}
            if n_an:
                pr[an] = n_an
            out.append((Cand("dissolve", {solid: 1}, pr, -pK, dH=dH_sol,
                             meta={"solid": solid}),
                        frozenset((solid,))))

    # ---- 4.3b 难溶金属氧化物酸溶（form:solid 二元氧化物 M_xO_y）
    # M_xO_y + 2yH+ -> xM^n+ + yH2O；logK 由对应氢氧化物 Ksp 估算。
    # _oh_ksp：阳离子 -> 氢氧化物 pKsp（供 _oxide_dissolve_info 查表）。
    oh = {}
    for e in T.ksp:
        if e["pair"][1] == "OH^-":
            oh[e["pair"][0]] = e["pKsp"]
    T._oh_ksp = oh
    for sp in T.ex:
        info = _oxide_dissolve_info(sp, T)
        if info is None:
            continue
        cat, x, n, y, pKsp = info
        out.append((Cand("dissolve", {sp: 1, H_ION: 2 * y}, {cat: x, WATER: y},
                         x * (n * 14 - pKsp - 1.5), x * n,
                         meta={"solid": sp}),
                    frozenset((sp,))))

    # ---- 4.4 配位 / 解离（OH- 配体正则化）
    for b in T.beta:
        center, lig, comp = b["center"], b["ligand"], b["complex"]
        nu = b["nu"]
        dH_b = b.get("dH")
        dH_nb = -dH_b if dH_b is not None else None
        if lig == "OH^-":
            out.append((Cand("complex", {center: 1, WATER: nu},
                             {comp: 1, H_ION: nu},
                             b["logb"] - nu * pKw, -nu, dH=dH_b),
                        frozenset((center,))))
            out.append((Cand("decomplex", {comp: 1, H_ION: nu},
                             {center: 1, WATER: nu},
                             nu * pKw - b["logb"], nu, dH=dH_nb),
                        frozenset((comp,))))
        else:
            out.append((Cand("complex", {center: 1, lig: nu}, {comp: 1},
                             b["logb"], dH=dH_b),
                        frozenset((center, lig))))
            out.append((Cand("decomplex", {comp: 1}, {center: 1, lig: nu},
                             -b["logb"], dH=dH_nb),
                        frozenset((comp,))))
    return out


def _gate_species(T) -> tuple:
    """浓度闸门（c_min/c_max/ligand_max/浓酸分子态）涉及的账本物种：
    枚举 memo 的浓度桶只对这些物种敏感（其余浓度不影响候选集）。"""
    gs = {"HNO_3"}
    for c in T.couples:
        g = c.get("gate")
        if not g:
            continue
        if "c_min" in g or "c_max" in g:
            acid = g.get("acid", c["ox"])
            gs.add(f"__tot_{acid}" if g.get("c_basis") == "shadow" else acid)
        for sp in g.get("ligand_max", ()):
            gs.add(sp)
    return tuple(sorted(gs))


def enumerate_candidates(ledger: dict, H_excess: float, pH: float, V: float, T_K: float, T,
                         kinetics: bool = True, memo: dict | None = None) -> list[Cand]:
    """枚举当前账本状态下全部可用候选。

    memo：主循环传入的枚举结果缓存（key = (present 集, pH 桶 0.1,
    闸门物种浓度 log 桶 0.25, V, T_K, kinetics)）。不动点迭代问
    present 集与 pH 常常数百步不变（中间体饥饿慢收敛：乒乓/催化循环），
    每次重扫 30k 模板纯属浪费——命中时直接返回候选列表（Cand 只读共享）。
    桶化阈值处（pH 闸门边界 ±0.05、浓度闸门边界 ×10^±0.125）行为可能与
    未缓存时有个别差异，由测试套件全量回归把关。"""
    pKw = pKw_of(T_K)
    aOH = 10.0 ** (pH - pKw)
    # present 含 H_ION（始终在账，由 H_excess 记账）与 WATER（溶剂），
    # 使原 _species_present(sp, present) 简化为 sp in present（4.2M 次调用消除）。
    # subset 判定 req <= present 不受影响：req 已排除 (WATER, H_ION)。
    present = {s for s, m in ledger.items() if m > X_MIN}
    present.add(H_ION)
    out: list[Cand] = []

    mkey = None
    if memo is not None:
        gs = getattr(T, "_gate_species", None)
        if gs is None:
            gs = T._gate_species = _gate_species(T)
        conc = tuple(sorted(
            (sp, round(log10(ledger.get(sp, 0.0) + 1e-30) * 4)) for sp in gs))
        mkey = (frozenset(present), round(pH, 1), conc, V, T_K, kinetics)
        hit = memo.get(mkey)
        if hit is not None:
            return hit

    # ---- 4.1 电子转移
    # 静态部分（电对配对/配平/电子数/慢标记/温度闸门）由 _redox_templates 按
    # (T, T_K) 缓存——主循环是不动点迭代，每迭代重建候选时这些"定值"不再重算；
    # 存在性、pH/浓度闸门、产物形态规则随状态变化，每迭代动态过滤。
    # 按 a_ox 索引预取在账氧化形的模板，下标排序后遍历保持原 (a,b) 嵌套序
    tmpls, by_aox, by_ared = _redox_templates(T_K, T, kinetics)
    # 层1 缓存：idxs 收集 + 存在性过滤只依赖 (present 集, 浓硝酸标志,
    # T_K, kinetics)——与 pH/浓度桶解耦。pH 桶漂移（缓冲体系平衡移动）导致
    # 层0 memo 反复 miss 时，层1 仍命中（沉淀转化类用例的主降载点）。
    hno3_hot = ledger.get("HNO_3", 0.0) / V >= 1.0
    pkey = (frozenset(present), hno3_hot, T_K, kinetics)
    pres_cache = getattr(T, "_enum_pres", None)
    if pres_cache is None:
        pres_cache = T._enum_pres = {}
    rows = pres_cache.get(pkey)
    if rows is None:
        idxs = []
        for k, ii in by_aox.items():
            if k in present:
                idxs.extend(ii)
            # 浓硝酸分子态：NO3- 缺席时以 HNO3 分子充当氧化剂
            elif k == "NO_3^-" and hno3_hot:
                idxs.extend(ii)
        fwd_set = set(idxs)
        # 同元素电对族的逆向存在性：a.red 在账即预取（歧化/归中排空通道，
        # 见 _redox_templates 注释）；跨元素对不做逆向预取
        for k, ii in by_ared.items():
            if k not in present:
                continue
            for i in ii:
                if i not in fwd_set:
                    idxs.append(i)
        idxs.sort()
        rows = []
        for i in idxs:
            (a, b, variants, _dE, _slow, _sds, _sdp, _al,
             same_elem, *_r) = tmpls[i]
            a_ox = a["ox"]
            fwd_ox = True
            if a_ox not in present:
                if a_ox == "NO_3^-" and hno3_hot:
                    a_ox = "HNO_3"
                else:
                    fwd_ox = False
            fwd = fwd_ox and b["red"] in present
            # 逆向可用仅限同元素族：a.red 与 b.ox 在账（如 Cu+ 在账、Cu2+
            # 在账时 2Cu+ ⇌ Cu+Cu2+ 参与竞争，方向由 S 决定）
            rev = same_elem and a["red"] in present and b["ox"] in present
            if fwd or rev:
                rows.append((i, fwd, rev, a_ox))
        if len(pres_cache) < 4096:
            pres_cache[pkey] = rows
    for i, fwd, rev, a_ox in rows:
        (a, b, variants, dE, slow, sd_static, sd_ph, acid_limited,
         same_elem, sd_ph_hi, deferred, pre_cands) = tmpls[i]
        # gate（pH_max/c_min/c_max）描述"该电对作为氧化剂"的可及性，
        # 只对 a（氧化剂）侧生效；b 为还原剂侧，其 ox 是产物，不应被闸门拦截
        # （例：NO2 溶于水被氧化为 NO3- 不要求 pH≤1.5）
        # 内联 _couple_gate_dyn -> _gate_check(a.get("gate"), ...) 消除 1.9M 次调用
        if fwd and not _gate_check(a.get("gate"), a, ledger, pH, V, T):
            fwd = False
        # 固相形态规则仅适用于溶剂/酸背景（a 为 H+）；强氧化剂在场时
        # 金属被氧化为游离离子是正常的（例：Cu+AgNO3 -> Cu2+ + Ag）
        if fwd and a_ox == H_ION and not _product_form_ok(b, aOH, T, T_K):
            fwd = False
        # b 侧还原通道 pH 下限：如 Mn2+ 氧化为 MnO2 实际经 Mn(OH)2，
        # 酸性中该方向动力学/形态上不可达
        if fwd and kinetics and pH < b.get("red_pH_min", -1e9):
            fwd = False
        # b 侧电对被逆向驱动（red 消耗、ox 生成）时的可选反向闸门 rev_gate：
        # 如 NO2 作还原剂被氧化回 NO3-，仅当残余游离酸稀薄——浓介质中生成
        # 的 NO2 以气体逸出，不被 Hg2+ 等边际氧化剂（ΔE≈0.05V）回氧化，
        # 与歧化支路同一"逸出气体"粗粒化（c_basis shadow，1M 阈值）
        rg = b.get("rev_gate") if kinetics else None
        if fwd and rg is not None and not _gate_check(rg, b, ledger, pH, V, T):
            fwd = False
        if not (fwd or rev):
            continue
        # 预构建 Cand 快路径：无 pH 依赖慢方向（约 99% 模板）时 meta 全
        # 常量——直接引用模板构建期的 Cand，省去每迭代的构造/字典开销。
        # 单变体模板（绝大多数）跳过变体搜索。
        n_var = len(variants)
        if n_var == 1:
            v = variants[0]
            if v[0] != a_ox:
                v = None
        else:
            v = next((x for x in variants if x[0] == a_ox), None)
        if v is None:
            continue
        if pre_cands is not None:
            out.append(pre_cands[0] if n_var == 1 else
                       pre_cands[variants.index(v)])
            continue
        _, r, pr, n, dH_full = v
        # 碱催化解锁的卤素歧化慢方向：pH 低于阈值才标记慢（碘量法窗口 pH<8、
        # 强碱中 I2/Br2 歧化至卤酸根为快反应）
        slow_dirs = sd_static
        for d, ph in sd_ph:
            if pH < ph:
                slow_dirs |= {d}
        for d, ph in sd_ph_hi:
            if pH > ph:
                slow_dirs |= {d}
        # 有固相电对的金属受酸限量约束（耗尽后走固相/膜路径）；
        # 无固相电对的金属（Na/Ca/K）水直接氧化，不受限
        out.append(Cand("redox", r, pr, 0.0, 0.0,
                        redox=(n, dE), dH=dH_full,
                        meta={"slow": slow, "slow_dirs": slow_dirs,
                              "ox_couple": a["ox"],
                              "acid_limited": acid_limited,
                              "deferred": deferred}))

    # ---- 4.2/4.3/4.3b/4.4 静态候选（质子/沉淀溶解/氧化物酸溶/配位）：
    # 候选本体（配平、logK、pkw_coeff、dH）只依赖数据表与 T_K，按
    # (T, T_K) 缓存一次构建；每迭代只按存在性（req ⊆ present）过滤，
    # 不再重建 Cand/字典（原实现每迭代全量重建，是 enumerate 的主开销）。
    # Cand 的 r/pr 创建后不可变（key 缓存依赖此），跨迭代复用安全。
    sc = getattr(T, "_cand_static", None)
    if sc is None or sc[0] != T_K:
        sc = (T_K, _build_static_cands(T_K, T, pKw))
        T._cand_static = sc
    # 静态段+派生段只依赖 present 集（与 pH/浓度闸门桶无关）——按
    # present 二级缓存：浓度演化使层0 memo（含 conc 桶）反复 miss 时，
    # 此段仍命中（2k 次子集扫描降为一次查表；Cu+AgNO3 类逐步氧化体系
    # 每 3~4 步就有浓度桶漂移，present 集数十步不变）。列表按原顺序
    # 构建（静态序 + 派生序），缓存命中与现场过滤 bit 级一致。
    stc = getattr(T, "_enum_static2", None)
    if stc is None:
        stc = T._enum_static2 = {}
    pset = frozenset(present)
    st_out = stc.get(pset)
    if st_out is None:
        st_out = [cand for cand, req in sc[1] if req <= pset]
        # 派生段快速过滤（元素位预筛 + 元组结构）：pset miss 时原实现
        # 对 ~13k 派生候选逐个调 _derived_present（函数调用 + 2 次
        # frozenset 子集判定 ≈ 0.3μs/个，Hg2_1 类 47 次 miss 共 ~0.2s）。
        # 预筛原理（按侧独立）：f ⊆ pset ⟹ elements(f) ⊆ elements(pset)
        # （逆否保证逐侧预筛永不误杀）；典型体系 pset 仅含 ≤6 种元素，
        # 13k 候选中逐侧元素相容者 ~1%（元素掩码 ≤72 位 int，单次 AND），
        # frozenset 子集判定只对幸存侧执行。结果序与逐个判定 bit 级一致
        # （2000 组随机 pset 等价性验证通过）。
        dfast = getattr(T, "_derived_fast", None)
        if dfast is None:
            dfast = _build_derived_fast(T)
        spe = T._sp_emask
        pm = 0
        for s in pset:
            m = spe.get(s)
            if m:
                pm |= m
        st_out.extend(cand for ef, er, f, r, cand in dfast
                      if (ef & pm == ef and f <= pset)
                      or (er & pm == er and r <= pset))
        if len(stc) < 512:
            stc[pset] = st_out
    out.extend(st_out)
    if mkey is not None and len(memo) < 4096:
        memo[mkey] = out
    return out


def _build_derived_fast(T):
    """派生候选过滤快速结构（T 级一次构建，与 build_derived 的 T._derived
    同生命周期）：[(fwd 元素掩码, rev 元素掩码, fwd_req, rev_req, cand)]
    + 物种→元素掩码表 T._sp_emask。惰性 req 补算在此一次性完成（与
    _derived_present 的幂等写一致）。元素掩码按侧独立（"任一侧 req ⊆
    pset" 的预筛必须逐侧：f ⊆ pset ⟹ elements(f) ⊆ elements(pset)，
    逆否保证逐侧预筛永不误杀）。"""
    d = build_derived(T)
    for c in d:
        if c.meta.get("fwd_req") is None:
            c.meta["fwd_req"] = frozenset(
                s for s in c.r if s not in (WATER, H_ION))
            c.meta["rev_req"] = frozenset(
                s for s in c.pr if s not in (WATER, H_ION))
    emap: dict = {}          # 元素 -> 位
    spe: dict = {}           # 物种 -> 元素掩码

    def _mask(species: frozenset) -> int:
        m = 0
        for s in species:
            v = spe.get(s)
            if v is None:
                v = 0
                for e in elements_of(s):
                    b = emap.get(e)
                    if b is None:
                        b = emap[e] = 1 << len(emap)
                    v |= b
                spe[s] = v
            m |= v
        return m

    rows = []
    for c in d:
        f, r = c.meta["fwd_req"], c.meta["rev_req"]
        rows.append((_mask(f), _mask(r), f, r, c))
    T._derived_fast = rows
    T._sp_emask = spe
    return rows


def _derived_present(c: Cand, pset: frozenset) -> bool:
    """派生候选在账判定（fwd/rev 任一方向 req ⊆ present）；
    req 的惰性补算写一次（Cand 跨迭代共享，meta 幂等）。"""
    req = c.meta.get("fwd_req")
    if req is None:   # 直接构造（非 add()）的派生候选：惰性补算，写一次
        req = c.meta["fwd_req"] = frozenset(
            s for s in c.r if s not in (WATER, H_ION))
        c.meta["rev_req"] = frozenset(
            s for s in c.pr if s not in (WATER, H_ION))
    return req <= pset or c.meta["rev_req"] <= pset


def _gate_check(g: dict | None, c: dict, ledger: dict, pH: float, V: float, T) -> bool:
    """单个闸门字典的动态检查（pH/浓度随迭代变化；T_min 静态部分在模板
    构建时过滤）。浓度闸门逐轮动态（真实化学：反应消耗酸使介质由浓转稀，
    氧化剂通道随之中途切换——T77/T78 同为 8M HNO3，仅凭 3Cu 耗尽酸库转稀
    才给出 NO）。"c_basis": "shadow" 时改读游离酸影子库存（仅本酸、不含
    盐类同离子）：用于 NO2 的两条回氧支路（歧化、金属离子回氧化）——
    3NO2+H2O->2HNO3+NO 是"NO2 通入水"的吸收平衡，被产物游离硝酸
    Le Chatelier 抑制；残余游离酸达摩尔量级时生成的 NO2 以气体逸出不被
    回吸/回氧（引擎不建气相的粗粒化，实际裕度 4M vs 0M，阈值宽区间不
    敏感）。反应生成的酸无影子库存（读作 0）→ NO2 通入纯水通道始终开。
    g 为 None 或空字典视为无闸门，避免上层 1.9M 次 `c.get('gate') or {}`
    的临时字典分配。"""
    if not g:
        return True
    if "pH_max" in g and pH > g["pH_max"]:
        return False
    if "pH_min" in g and pH < g["pH_min"]:
        return False
    if "c_min" in g or "c_max" in g:
        acid = g.get("acid", c["ox"])
        if g.get("c_basis") == "shadow":
            conc = ledger.get(f"__tot_{acid}", 0.0) / V
        else:
            conc = _acid_conc(acid, ledger, V, T)
        if "c_min" in g and conc < g["c_min"]:
            return False
        if "c_max" in g and conc >= g["c_max"]:
            return False
    # ligand_max: 配位掩蔽闸门——当指定配体浓度超过阈值时，该电对不再作为
    # 氧化剂（化学事实：强配位剂如 S2O3^2-/CN- 把金属阳离子掩蔽为配阴离子，
    # 有效 E0 大幅下降，原本的 redox 不再发生）。例如 Ag+/Ag (E0=0.80) 在
    # S2O3^2- 存在时被掩蔽为 [Ag(S2O3)2]3-/Ag (E0=0.002)，无法被 S2O3 还原。
    ligand_max = g.get("ligand_max")
    if ligand_max:
        for sp, cmax in ligand_max.items():
            if ledger.get(sp, 0.0) / V >= cmax:
                return False
    return True


def _acid_conc(acid: str, ledger: dict, V: float, T) -> float:
    """分子态浓酸的浓度（浓 H2SO4 / 浓 HNO3）；非分子形态返回 0。"""
    ex = T.ex.get(acid)
    if not ex:
        return 0.0
    cf = ex.get("conc_forms")
    if cf and cf.get("concentrated") == "molecule":
        return ledger.get(acid, 0.0) / V
    return 0.0


def _solid_acid_pka(base: str, T) -> float | None:
    """固相酸表观 pKa：沿 base 的逐级脱质子链走到 ksp 记录的 (H+, 酸根)
    端点，返回 pKsp − 链上 pKa 之和；走不到端点返回 None（保持原 pKa）。"""
    for (cat, an), cell in T.ksp_by_pair.items():
        if cat != H_ION:
            continue
        sp, s, seen = base, 0.0, set()
        while sp != an and sp not in seen:
            seen.add(sp)
            nxt = T.pka_acid.get(sp)
            if not nxt:
                break
            e1 = min(nxt, key=lambda e: e["pka"])
            s += e1["pka"]
            sp = e1["base"]
        if sp == an:
            return cell["pKsp"] - s
    return None




def _ksp_asat(cell: dict, aOH: float, T_K: float) -> float:
    """氢氧化物固相在 aOH 下其游离阳离子的饱和活度（Ksp = a_cat^x·aOH^y）。"""
    x, y = _ksp_xy(cell)
    pksp = cell["pKsp"] - (_vant(cell["dH"], T_K) if "dH" in cell else 0.0)
    return (10.0 ** (-pksp) / aOH ** y) ** (1.0 / x)


def _product_form_ok(b: dict, aOH: float, T, T_K: float = 298.15) -> bool:
    """固相形态规则（§4.6）：b.ox 为游离金属离子、存在同金属固相电对、
    且当前 pH 下其氢氧化物饱和活度低于 SAT_SKIP 时，跳过裸离子路径。"""
    ox = b["ox"]
    cell = T.ksp_by_pair.get((ox, "OH^-"))
    if cell is None:
        return True
    solid = cell["solid"]
    if not any(c["ox"] == solid and c["red"] == b["red"] for c in T.couples):
        return True
    a_sat = _ksp_asat(cell, aOH, T_K)
    return a_sat > SAT_SKIP

