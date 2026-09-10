"""chemkit.candidates：候选反应对象与配平层（引擎地基模块）。

职能（自底向上，被 normalize/speciation/templates/engine 依赖）：
    · 全局常数唯一事实源（X_MIN/ACT_FLOOR/动力学阈值/气体活度约定……）
    · Cand 候选反应（kind/r/pr/logK/dH/meta + 热路径缓存字段）
    · logK_T：logK(T) 三通道（redox Nernst 缩放 / van't Hoff / pKw 耦合）
    · 半反应配平（_bal / _bal_fast / _half_pairs / _redox_mix_ok）
    · 兜底电离（元素周期表级拆盐 _split_salt_general：公式结构+
      电荷平衡+氧化态交叉验证，未注册盐电离为旁观离子）
    · 派生候选（build_derived：pKa/pKsp/logβ 的 Hess 精确加和）

自身只依赖 core/data——模块依赖链的最底层。
"""
from __future__ import annotations
import heapq
from dataclasses import dataclass, field
from functools import reduce
from math import gcd, log10, sqrt
from .core import (elements_of, charge_of, balance,
                   K_NERNST_298, PKW_298, k_nernst, _vant, pKw_of)
from .data import Tables, _half_balance

# ---------------------------------------------------------- 全局常数（§7 纪律：唯一、文档化）
# ---- 外界气相模型：总压 101.3 kPa 惰性环境（外界气相不含产物气体）----
# 热力学气体活度 a = p/p°（p°=101.325 kPa 标准压力）。
# 自产气体：惰性环境下外界 p_i≈0，气体不断逸出直至残余分压与扫气速率达到
# 稳态；按残余分压约定 P_RES=1 kPa（开放烧杯自然对流量级）→ a≈1e-2。
# 外加供给气体（投料中的气体试剂）：可溶气体试剂（氨水、氯水等）真实形态
# 是溶质，按溶质活度 a=c/c°（"持续供给"语义由账本浓度承载）——不用
# p/p°≈1，否则 NH3 络合溶解（AgCl+浓氨水）等全部被低估。
# 自产气体按泡点判据（见 S_of 气体分支）。
P_STD_KPA = 101.325             # 热力学标准压力
P_EXT_KPA = 101.3               # 外界气相总压（惰性环境）
P_RES_KPA = 1.0                 # 逸出气体残余分压约定（惰性扫气稳态）
A_GAS = P_RES_KPA / P_STD_KPA   # 逸出气体活度（产物侧）≈ 1e-2，全局唯一

# 动力学层总开关——现为 judge() 的 kinetics 参数（默认 True，经
# conditions["kinetics"]
# 传入，贯穿 enumerate_candidates/_vs_gate_ok/_redox_templates/主循环）。
# False = 纯热力学基线（无限时间）：数据中的 slow/gate/ox_inert/T_min/
# rev_gate/red_pH_min 等动力学标记、卤素歧化慢方向规则、膜封锁（钝化膜/
# 覆盖层）机制一律不生效。用于暴露热力学预测与现实的偏差清单。
# 非金属元素集：同元素电对族逆向存在性（歧化/归中排空）只对共享金属
# 元素的配对开放；非金属多价链的逆向通道由既有正向闸门+动力学标记约束
_NONMETAL_ELEMS = frozenset(
    {"H", "O", "N", "S", "F", "Cl", "Br", "I", "C", "P", "Se", "Si",
     "B", "As", "Sb", "Te"})
ACT_FLOOR = 1e-12               # 溶质活度数值下限（只是种子，程度由零点决定）
X_MIN = 1e-6                    # 可忽略程度（mol）
SOLVENT_FIRST = 10.0            # 溶剂优先阈值（界面动力学规则，全引擎仅此一条）
BLOCKED_EXTENT = 0.02           # 膜封锁时的痕量程度
SAT_SKIP = 0.05                 # 固相形态规则：饱和活度低于此值则跳过裸离子路径
ANN_MIN_EXTENT = 1e-3           # 慢反应标注的显著程度阈值（mol，低于此量级不标注）
DEGREE_COMPLETE = 0.99
DEGREE_PARTIAL = 0.05
MAX_ITER = 3000
# 水溶液中优先与水反应的活泼金属（教材规则）
WATER_FIRST_METALS = {"Li", "Na", "K", "Rb", "Cs", "Ca", "Sr", "Ba"}
# 非金属单质固相：S(IV)-金属规则的例外（SO2 还原含氧酸制非金属单质是工业反应，
# 如 H2SeO3+2SO2+H2O->Se+2H2SO4 用于硒提纯）
NONMETAL_SOLIDS = {"C", "S", "Se", "Si", "P_4", "I_2", "B"}
# 卤酸根集合与歧化解锁温度（教材规则：卤素+冷稀碱 -> XO-，+热碱 -> XO3-；
# 歧化经 XO- 中间体，其进一步歧化室温慢——漂白液室温稳定即此动力学事实）
# 卤素歧化到卤酸根的慢化温度（分卤素，文献动力学）：ClO- -> ClO3- 室温慢、
# 热碱快（教材"冷稀碱 -> ClO-，热碱 -> ClO3-"，漂白液室温稳定）→ T_lim=340。
# BrO-/IO- 歧化快，但直接歧化需经 XO- 中间体、由 OH- 催化：中性水中 X2 歧化
# 动力学封闭（氯水/溴水/碘水皆稳定），强碱中解锁——T_lim 取极大值表示"无温度
# 解锁"，由 HALATE_BASE_PH 给出碱解锁阈值：Br2+2OH- -> BrO- 平衡常数 ~1e11
# （Kh(Br2)=5.8e-9, Ka(HBrO)=2e-9），pH≥9 即显著；I2 体系 K≈4.6e4
# （Kh(I2)=2e-13, Ka(HIO)=2.3e-11），pH≥12 才显著（碘量法 pH 5-8 操作窗口、
# pH 11 碘仍稳定，即 D44 NaClO+KI 终态）。
HALATE_DISP_T = {"ClO_3^-": 340.0, "BrO_3^-": 1e9, "IO_3^-": 1e9}
HALATE_BASE_PH = {"BrO_3^-": 9.0, "IO_3^-": 13.0}
# 固相 -> 阳离子映射（ksp 表），供活泼金属规则识别固相氢氧化物氧化剂
_HALF_CACHE: dict = {}

def _half_scale(c: dict) -> float:
    # 半反应中每个 red 粒子对应的电子数 = n / (red 系数)。
    # red 系数由活性元素（非 O/H）原子守恒推出：IO3-/I2 -> 0.5，I2/I- -> 2
    key = (c["ox"], c["red"])
    if key in _HALF_CACHE:
        return _HALF_CACHE[key]
    eo = elements_of(c["ox"])
    er = elements_of(c["red"])
    nu = 1.0
    for el, cnt in eo.items():
        if el in ("O", "H"):
            continue
        if er.get(el):
            nu = cnt / er[el]
            break
    else:
        # 两侧均只含 O/H（如 H2O2/H2O、O2/OH-）：按 O 原子数折算
        if eo.get("O") and er.get("O"):
            nu = eo["O"] / er["O"]
    scale = c["n"] / nu if nu else float(c["n"])
    _HALF_CACHE[key] = scale
    return scale
WATER = "H_2O"
H_ION = "H^+"
WATER_MOL_PER_L = 55.6
# 浓溶液中以分子态记账的强酸（HCl 任何浓度全电离，不在此列——原条目永不命中，已移除）
STRONG_MOLECULAR_ACIDS = {"H_2SO_4": 2, "HNO_3": 1}
# estimate_state 静态缓存中的强酸标记（pKa<=0，已在 normalize 拆解，此处仅作
# h_c 直读上限贡献）；用对象哨兵避免与 None（不在表中）歧义
_STRONG_ACID = object()

# ========================================================== 兜底电离（元素周期表级）
# 离子注册（T.cations/T.anions）是反应数据的入口（Ksp/pKa/E0/β），不是
# 离子"存在"的许可——拆盐兜底让任何化学上成立的盐都能电离成旁观离子。

# 元素符号全集（118 号）：投料化学式的合法性校验
_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
    "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())

# 单原子阴离子典型价（绝对值）：卤素/氢 -1，氧族 -2，氮族 -3，碳族/硅 -4、硼 -3
_MONO_ANION_Q = {"F": 1, "Cl": 1, "Br": 1, "I": 1, "At": 1, "H": 1,
                 "O": 2, "S": 2, "Se": 2, "Te": 2, "Po": 2,
                 "N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3,
                 "C": 4, "Si": 4, "Ge": 4, "B": 3}

# 金属常见水溶液氧化态（成盐阳离子电荷）。兜底拆盐的电荷交叉验证：
# 阳离子电荷由化学计量 + 阴离子电荷推出（x·qc = y·qa），必须落在该金属
# 的常见价态集内（AlCl4 这类不成立组合被拒绝；Na2O2/KO2 由电荷平衡自然
# 区分为过氧化物/超氧化物）
_METAL_Q = {
    # 碱金属 / 碱土金属
    "Li": (1,), "Na": (1,), "K": (1,), "Rb": (1,), "Cs": (1,), "Fr": (1,),
    "Be": (2,), "Mg": (2,), "Ca": (2,), "Sr": (2,), "Ba": (2,), "Ra": (2,),
    # 第 3 族 + 镧系
    "Sc": (3,), "Y": (3,), "La": (3,), "Ce": (3, 4), "Pr": (3, 4),
    "Nd": (3,), "Pm": (3,), "Sm": (2, 3), "Eu": (2, 3), "Gd": (3,),
    "Tb": (3, 4), "Dy": (3,), "Ho": (3,), "Er": (3,), "Tm": (2, 3),
    "Yb": (2, 3), "Lu": (3,),
    # 锕系（水溶液化学常见价）
    "Th": (4,), "Pa": (4, 5), "U": (3, 4, 6), "Np": (4, 5), "Pu": (3, 4, 6),
    "Am": (3,), "Cm": (3,), "Bk": (3, 4), "Cf": (3,), "Es": (3,), "Fm": (3,),
    "Md": (2, 3), "No": (2, 3), "Lr": (3,),
    # 过渡金属
    "Ti": (2, 3, 4), "Zr": (4,), "Hf": (4,),
    "V": (2, 3, 4, 5), "Nb": (4, 5), "Ta": (4, 5),
    "Cr": (2, 3, 6), "Mo": (3, 4, 5, 6), "W": (4, 5, 6),
    "Mn": (2, 3, 4, 6, 7), "Tc": (4, 7), "Re": (3, 4, 6, 7),
    "Fe": (2, 3), "Ru": (2, 3, 4, 6, 8), "Os": (3, 4, 6, 8),
    "Co": (2, 3), "Rh": (1, 2, 3, 4), "Ir": (2, 3, 4, 6),
    "Ni": (2, 3, 4), "Pd": (2, 4), "Pt": (2, 4, 6),
    "Cu": (1, 2), "Ag": (1, 2, 3), "Au": (1, 3),
    "Zn": (2,), "Cd": (2,), "Hg": (1, 2),
    # 后过渡 / 类金属（成盐价态）
    "Al": (3,), "Ga": (3,), "In": (1, 3), "Sn": (2, 4), "Pb": (2, 4),
    "Tl": (1, 3), "Bi": (3, 5), "Ge": (2, 4), "Sb": (3, 5), "As": (3, 5),
    "Si": (4,), "Po": (2, 4),
}


def _ion_name(el: str, q: int, positive: bool) -> str:
    """单原子离子名：Na^+ / Fe^{3+} / Cl^- / S^{2-}（引擎记法）。"""
    sign = "+" if positive else "-"
    return f"{el}^{sign}" if q == 1 else f"{el}^{{{q}{sign}}}"


_SALT_GENERAL_CACHE: dict[str, dict | None] = {}


def _split_salt_general(name: str, T):
    """兜底拆盐（注册表未覆盖时）：元素周期表级的结构性电离。

    化学事实依据：
      - 阳离子取公式中的金属元素 M（常见水溶液氧化态集见 _METAL_Q）；
      - 阴离子按其余元素组成匹配注册阴离子（组成整倍数），或单原子典型价
        （卤/氧/氮/碳族）；
      - 电荷平衡 x·qc = y·qa 推出阳离子电荷，交叉验证落在 M 的常见价态内。
    无反应数据（Ksp/pKa/电对/β）的离子成为旁观离子——正确电离入账但
    不参与反应，与"未知整盐静默丢弃"相比保留了电荷与质量守恒。
    返回 {离子: 计量} 或 None。"""
    if name in _SALT_GENERAL_CACHE:
        return _SALT_GENERAL_CACHE[name]
    result = None
    try:
        elems = elements_of(name)
        if len(elems) >= 2 and all(e in _ELEMENTS for e in elems):
            # 注册阴离子组成缓存（每数据表一次；排除误入 anions 的中性碱；
            # 排序保证跨进程确定性——hash 随机化不影响候选序）
            comps = getattr(T, "_anion_comps", None)
            if comps is None:
                comps = sorted(
                    ((an, -charge_of(an), elements_of(an))
                     for an in T.anions if charge_of(an) < 0))
                T._anion_comps = [c for c in comps]
            cands = []
            for M in (e for e in elems if e in _METAL_Q):
                x = elems[M]
                rest = {k: v for k, v in elems.items() if k != M}
                # (a) 注册阴离子：rest 必须是阴离子组成的整倍数——
                # 整除性对全部元素强制（仅筛可整除键会把 NO2- 错配到
                # 硝酸盐：O 9%2≠0 被静默跳过，y 误判）
                for an, qa, ce in comps:
                    if set(ce) != set(rest) or not ce:
                        continue
                    if any(rest.get(k, 0) % v for k, v in ce.items()):
                        continue
                    ys = {rest[k] // ce[k] for k in ce}
                    if len(ys) != 1:
                        continue
                    y = ys.pop()
                    if y < 1:
                        continue
                    qc, rem = divmod(y * qa, x)
                    if rem == 0 and qc in _METAL_Q[M]:
                        cands.append((_ion_name(M, qc, True), x, an, y, qc))
                # (b) 单原子阴离子典型价：rest 为单一元素
                if len(rest) == 1:
                    (E, y), = rest.items()
                    qa = _MONO_ANION_Q.get(E, 0)
                    if qa:
                        qc, rem = divmod(y * qa, x)
                        if rem == 0 and qc in _METAL_Q[M]:
                            cands.append((_ion_name(M, qc, True), x,
                                          _ion_name(E, qa, False), y, qc))
            if cands:
                # 多解择优：阳离子电荷最小者（同组成罕见歧义，如混价前体）
                cat, x, an, y, _qc = min(cands, key=lambda c: (c[4], -c[3]))
                result = {cat: x, an: y}
    except Exception:
        result = None
    _SALT_GENERAL_CACHE[name] = result
    return result









@dataclass(slots=True)
class Cand:
    kind: str                 # redox / proton / precip / dissolve / complex / decomplex / derived
    r: dict                   # 反应物（H+ 正则形，可含 "H^+" 与 "H_2O"）
    pr: dict                  # 产物
    logK: float               # 298 K 值；redox 为 n*dE0/0.05916
    pkw_coeff: float = 0.0    # logK(T) = logK + pkw_coeff*(pKw(T)-14)
    redox: tuple | None = None  # (n_e, dE0) —— redox 的 logK 随 T 变
    dH: float | None = None   # 反应焓 kJ/mol（正向书写式，**不含水电离部分**——
                              # 该部分恒由 pkw_coeff 通道承载；None=无数据回退）
    meta: dict = field(default_factory=dict)
    _key: tuple | None = field(default=None, repr=False, compare=False)
    # ---- 热路径缓存（Cand 生命周期内不变量，跨 judge 调用共享） ----
    # logK(T_K) 纯函数值缓存（同一 T_K 直接命中，避免 22 万次/慢例的
    # Nernst/van't Hoff 重算）；执行计划/净键/在场物种元组同理
    _logK_c: dict = field(default_factory=dict, repr=False, compare=False)
    _plan: tuple | None = field(default=None, repr=False, compare=False)
    _plan_tid: int = field(default=-1, repr=False, compare=False)
    _pres: tuple | None = field(default=None, repr=False, compare=False)
    _nk_f: tuple | None = field(default=None, repr=False, compare=False)
    _nk_r: tuple | None = field(default=None, repr=False, compare=False)

    @property
    def key(self) -> tuple:
        """候选唯一键（(sorted_r, sorted_pr)）。r/pr 创建后不可变，缓存避免
        重复 sorted+tuple（原 206k 次调用每代重排，是 judge 主循环的隐性热点）。"""
        k = self._key
        if k is None:
            k = (tuple(sorted(self.r.items())), tuple(sorted(self.pr.items())))
            self._key = k
        return k

    @property
    def pres_specs(self) -> tuple:
        """(反应物物种元组, 产物物种元组)——H2O/H+ 之外的在场检查对象。"""
        p = self._pres
        if p is None:
            p = (tuple(s for s in self.r if s not in (WATER, H_ION)),
                 tuple(s for s in self.pr if s not in (WATER, H_ION)))
            self._pres = p
        return p

    @property
    def netkey_fwd(self) -> tuple:
        """净反应键（忽略 H2O/H+ 的 sorted 物种对），正向。"""
        k = self._nk_f
        if k is None:
            k = (tuple(sorted((s, n) for s, n in self.r.items()
                              if s not in (WATER, H_ION))),
                 tuple(sorted((s, n) for s, n in self.pr.items()
                              if s not in (WATER, H_ION))))
            self._nk_f = k
        return k

    @property
    def netkey_rev(self) -> tuple:
        """净反应键，反向（netkey_fwd 的对偶）。"""
        k = self._nk_r
        if k is None:
            f = self.netkey_fwd
            k = (f[1], f[0])
            self._nk_r = k
        return k


def logK_T(c: Cand, T_K: float) -> float:
    """logK(T) 三通道：redox Nernst 缩放（无 dH 时的回退）/ van't Hoff（有 dH）
    + pkw_coeff 水电离部分。van't Hoff 在 298.15 K 恒等，493 套件零扰动。
    值只依赖 (Cand, T_K)——按 T_K 缓存在 Cand 上（S_of 每次求值的必经点，
    缓存后与逐次重算 bit 级一致）。"""
    k = c._logK_c.get(T_K)
    if k is None:
        if c.redox is not None:
            n, dE0 = c.redox
            if c.dH is not None:
                k = n * dE0 / K_NERNST_298 + _vant(c.dH, T_K)
            else:
                k = n * dE0 / k_nernst(T_K)
        else:
            vant = _vant(c.dH, T_K) if c.dH is not None else 0.0
            k = c.logK + c.pkw_coeff * (pKw_of(T_K) - PKW_298) + vant
        c._logK_c[T_K] = k
    return k


# ========================================================== 配平缓存

_BAL_CACHE: dict = {}


def _bal(reactants, products, free):
    key = (tuple(reactants), tuple(products), tuple(free))
    if key not in _BAL_CACHE:
        _BAL_CACHE[key] = balance(list(reactants), list(products), free=list(free))
    return _BAL_CACHE[key]


# ---------------------------------------------------------- 电对配对配平快路径
# _redox_templates 对 couples² 全对枚举配平（~5 万次/温度点），sympy nullspace
# 每次约 1.2ms 是最大热点。常规电对的配平结果是两个半反应的电子 lcm 组合
# （教科书标准方程），可用整数运算 O(1) 构造：s_a·half_a(还原) +
# s_b·half_b(氧化逆转)，e- 净零，H+/H2O 自然对消，再归一化到最小整数比。
# 仅同元素电对族（same_elem：歧化/归中，零空间欠定需幻影校验）与半反应
# 构造失败的电对（配合物/固相晶格需配体入账）仍走 sympy。

def _half_pairs(T) -> dict:
    """(ox, red) → 半反应整数系数 (k, h, m, w, ne)（缓存于 T._half_pairs）。

    k·ox + h·H+ + ne·e- -> m·red + w·H2O（w 可为负：H2O 在反应物侧）。
    HNO_3 分子态变体（NO3- 电对用）：h 减 1。"""
    hm = getattr(T, "_half_pairs", None)
    if hm is None:
        hm = {}
        for c in T.couples:
            try:
                hb = _half_balance(c["ox"], c["red"])
            except Exception:
                hb = None
            if hb is None:
                continue
            k, h, m, w, ne = hb
            hm.setdefault((c["ox"], c["red"]), (k, h, m, w, ne))
            if c["ox"] == "NO_3^-" and h >= 1:
                hm.setdefault(("HNO_3", c["red"]), (k, h - 1, m, w, ne))
        T._half_pairs = hm
    return hm


def _bal_fast(a_ox: str, a_red: str, b_ox: str, b_red: str, T):
    """电对配对 (a_ox + b_red -> a_red + b_ox) 的半反应组合配平快路径。

    返回与 balance() 同构的 {'reactants': {...}, 'products': {...}}；
    构造不出半反应 / 方向校验 / 元素-电荷守恒失败返回 None（上层回退
    sympy _bal——保证与慢路径结果完全一致的兜底）。

    严格等价约束：balance() 把同时出现在 free 池与指定侧的物种锁定在
    指定侧（H+/H2O 电对作 a 或 b 侧时 H+ 的侧向语义），半反应组合可能
    把它消到另一侧——此时必须回退，否则会引入慢路径不存在的模板。
    故四个端点物种任一属于 {H2O, H+} 时直接走 sympy。"""
    if any(sp in (WATER, H_ION) for sp in (a_ox, a_red, b_ox, b_red)):
        return None
    # 零空间维数 = 1 保证唯一性（6 物种 − (元素行 + 电荷行)）：四个端点
    # 物种的元素与 {H, O} 并集 ≥ 4（跨元素对）时 rank = 5、dim = 1，半反应
    # 组合即唯一解；单元素族（卤素/氮/硫歧化网络，E=3 → dim ≥ 2）回退
    # sympy——欠定零空间中 sympy 挑选的组合与半反应 lcm 可能不同点。
    _els = (set(elements_of(a_ox)) | set(elements_of(a_red))
            | set(elements_of(b_ox)) | set(elements_of(b_red)) | {"H", "O"})
    if len(_els) < 4:
        return None
    hm = _half_pairs(T)
    ha, hb = hm.get((a_ox, a_red)), hm.get((b_ox, b_red))
    if ha is None or hb is None:
        return None
    k_a, h_a, m_a, w_a, ne_a = ha
    k_b, h_b, m_b, w_b, ne_b = hb
    if ne_a <= 0 or ne_b <= 0:
        return None
    g = gcd(ne_a, ne_b)
    s_a, s_b = ne_b // g, ne_a // g   # 电子数 lcm 缩放
    v: dict = {}

    def _acc(sp, nu):
        v[sp] = v.get(sp, 0) + nu

    # a 半反应（还原方向）× s_a
    _acc(a_ox, -k_a * s_a)
    _acc(a_red, m_a * s_a)
    _acc(H_ION, -h_a * s_a)
    _acc(WATER, w_a * s_a)
    # b 半反应（氧化方向 = 还原式取负）× s_b
    _acc(b_red, -m_b * s_b)
    _acc(b_ox, k_b * s_b)
    _acc(H_ION, h_b * s_b)
    _acc(WATER, -w_b * s_b)
    v = {sp: nu for sp, nu in v.items() if nu}
    # 固定物种方向校验（与 balance._try_vec 的侧一致性同口径）；
    # H+/H2O 属自由池物种（如 a 侧电对 ox 即为 H+ 时，H+ 系数可与
    # b 侧氧化释出的 H+ 对消归零，方向不受限定）
    _free = (WATER, H_ION)
    if any(v.get(sp, 0) > 0 for sp in (a_ox, b_red) if sp not in _free):
        return None
    if any(v.get(sp, 0) < 0 for sp in (a_red, b_ox) if sp not in _free):
        return None
    # 归一化到最小整数比
    g2 = reduce(gcd, (abs(nu) for nu in v.values()))
    if g2 > 1:
        v = {sp: nu // g2 for sp, nu in v.items()}
    # 守恒断言（元素 + 电荷；整数精确）：任何失败回退 sympy。
    # 按物种遍历其元素表累积差值（整数和与次序无关，共享只读 dict）。
    diff: dict[str, int] = {}
    for sp, nu in v.items():
        for el, cnt in elements_of(sp).items():
            diff[el] = diff.get(el, 0) + cnt * nu
    if any(diff.values()):
        return None
    if sum(charge_of(sp) * nu for sp, nu in v.items()):
        return None
    return {"reactants": {sp: -nu for sp, nu in v.items() if nu < 0},
            "products": {sp: nu for sp, nu in v.items() if nu > 0}}


def _redox_mix_ok(a: dict, b: dict, a_ox: str, r: dict, pr: dict,
                  pool: list, metals: set) -> bool:
    """同元素电对族配平结果的电子守衡校验。

    同元素多电对族的配平问题欠定（零空间含第三个电对的半反应方向），
    _bal 任取的一点可能实为三个半反应的混合（幻影：5MnO2 + MnO4^2- +
    2H2O -> 4MnOOH + 2MnO4- 混入了 MnO4-/MnO2 电对），其 n 折算与
    logK = n·dE/k 毫无意义。合法候选必须是两个半反应的干净组合
    （R = x·h_a + y·h_b）；残余只允许非氧化还原旁观（沉淀/酸碱——
    共享金属的氧化数不变，如 2Fe(OH)3+3Fe2+ -> 3Fe(OH)2+2Fe3+ 中
    1 个 Fe2+ 仅沉淀）。仅在模板构建期调用（按 (T, T_K) 缓存）。
    """
    def _half_vec(c, ox_name):
        """电对半反应（ox -> red，还原方向）的物种系数向量（产物正、
        反应物负；电子不入向量——_bal 无法配平半反应的电荷）。金属
        （非 O/H）原子守恒定 ox/red 计量比，O 差补 H2O、H 差补 H+。
        含配体/晶格离子等无法按此规则配平的电对返回 None（上层放行）。"""
        eo, er = elements_of(ox_name), elements_of(c["red"])
        nu_r = None
        for el, cnt in eo.items():
            if el in ("O", "H"):
                continue
            if not er.get(el):
                return None
            k = cnt / er[el]
            if nu_r is None:
                nu_r = k
            elif abs(nu_r - k) > 1e-9:
                return None
        if nu_r is None:
            return None
        # 还原形含氧化形没有的金属（晶格电对 S/HgS：半反应 S + Hg2+ + 2e
        # -> HgS 需游离离子入账，此处无法构造）→ None（上层保守接受）
        for el in er:
            if el not in ("O", "H") and el not in eo:
                return None
        v = {ox_name: -1.0, c["red"]: nu_r}
        dO = eo.get("O", 0) - nu_r * er.get("O", 0)
        v[WATER] = float(dO)
        v[H_ION] = -(eo.get("H", 0) - nu_r * er.get("H", 0) + 2.0 * dO)
        return v

    va = _half_vec(a, a_ox)      # a 还原方向
    vb = _half_vec(b, b["ox"])   # b 还原方向
    if va is None or vb is None:
        return True   # 半反应配不出时不比现状更糟，保持原行为

    R: dict = {}
    for sp, nu in r.items():
        R[sp] = R.get(sp, 0.0) - nu
    for sp, nu in pr.items():
        R[sp] = R.get(sp, 0.0) + nu

    def _scale(v_self, v_other, anchors):
        for sp in anchors:
            cs, co = v_self.get(sp, 0.0), v_other.get(sp, 0.0)
            if cs and not co and sp in R:
                return R[sp] / cs
        return None

    # 合法候选须为干净组合 R = x·va − y·vb（a 还原、b 氧化）
    x = _scale(va, vb, (a_ox, a["red"]))
    y = _scale(vb, va, (b["ox"], b["red"]))
    if x is None or y is None:
        return True   # 退化情形无法判定时保持原行为
    y = -y   # R 中 vb 分量以氧化方向（−vb）出现
    resid = {sp: nu - x * va.get(sp, 0.0) + y * vb.get(sp, 0.0)
             for sp, nu in R.items()}
    resid = {sp: nu for sp, nu in resid.items() if abs(nu) > 1e-9}
    # 半反应独有物种在 R 中缺席时 R 缺失对应分量，按 0 处理即可。
    # sorted：set 迭代序随 PYTHONHASHSEED 变化，会改变 resid 的插入序、
    # 进而改变下方金属守衡浮点累加顺序——N08/H25/UO02 三例曾因此跨进程
    # 掷骰子（digest 双态翻转，v0.3.8 差分基线跨进程不可复现的根因）。
    # 排序后插入序确定，累加序固定为字典序（两种历史形态各半出现的
    # 边界例从此固定取其一——1247 在两种形态下均 PASS，语义无损）
    for sp in sorted((set(va) | set(vb)) - set(R)):
        nu = -x * va.get(sp, 0.0) + y * vb.get(sp, 0.0)
        if abs(nu) > 1e-9:
            resid[sp] = resid.get(sp, 0.0) + nu
    if not resid:
        return True   # 干净的两电对组合
    # 残余必须为氧化数守衡的非氧化还原旁观（沉淀/酸碱）；只对
    # 组成元素为 金属+O/H 的物种可可靠判定（O=-2/H=+1 规则），
    # 含配体等复杂物种无法判定时保持原行为
    for m in metals:
        s = 0.0
        for sp, nu in resid.items():
            el = elements_of(sp)
            if m not in el:
                continue
            if any(e2 not in ("O", "H") for e2 in el if e2 != m):
                return True
            on = (charge_of(sp) + 2.0 * el.get("O", 0)
                  - el.get("H", 0)) / el[m]
            s += nu * on * el[m]
        if abs(s) > 1e-6:
            return False
    return True



def _ksp_xy(e: dict) -> tuple[int, int]:
    """Ksp 溶解反应的化学计量 (x 阳离子, y 阴离子)。

    以固相分子式中阳离子元素计数为准（Ag2O → x=2），阴离子数由电中性
    推出（y = x·q_cat/q_an）；纯电荷 gcd 只对 M(OH)n 型 1:n 盐正确，
    对 Ag2O 这类 2:1 固体会配出 Ag2O + H+ -> Ag+ 这种不守恒方程。
    分子式缺失/不整除时退回电荷 gcd。
    """
    cat, an = e["pair"]
    qc, qa = charge_of(cat), -charge_of(an)
    g = gcd(qc, qa)
    el = next(iter(elements_of(cat)))
    x = elements_of(e["solid"]).get(el, 0)
    if x > 0 and (x * qc) % qa == 0:
        return x, x * qc // qa
    return qa // g, qc // g


# ========================================================== 派生候选（§4.5：Hess 精确加和）

def build_derived(T) -> list[Cand]:
    if hasattr(T, "_derived"):
        return T._derived
    out: list[Cand] = []
    def add(kind, r, pr, logK, pkw_coeff=0.0, src="", dH=None):
        r = {k: v for k, v in r.items() if v}
        pr = {k: v for k, v in pr.items() if v}
        if not r or not pr:
            return
        # 静态需求集（除 H2O/H+）：§4.5 存在性过滤每迭代做子集判定，
        # 避免反复 genexpr/all 扫描（派生候选近千条，曾是最大热点之一）
        out.append(Cand(kind, r, pr, logK, pkw_coeff, dH=dH,
                        meta={"src": src,
                              "fwd_req": frozenset(s for s in r if s not in (WATER, H_ION)),
                              "rev_req": frozenset(s for s in pr if s not in (WATER, H_ION))}))

    def _hess_dH(*terms):
        """派生候选的 dH = 组元 dH 的同号 Hess 加和（组元缺数据→None 整体回退；
        水电离部分不进 dH，由 pkw_coeff 通道承载）。"""
        if any(t is None for t in terms):
            return None
        return sum(terms)

    for e in T.ksp:
        cat, an = e["pair"]
        solid, pK = e["solid"], e["pKsp"]
        # Ksp ⊗ β：solid + ν·ligand ⇌ complex + n_an·anion（logK = logβ − pKsp）
        # 计量系数由 _bal 现场配平给出，不预计算 n_cat/n_an
        for b in T.beta:
            if b["center"] != cat:
                continue
            bal = _bal([solid, b["ligand"]], [b["complex"], an], [WATER, H_ION])
            if bal is None:
                continue
            k_complex = bal["products"].get(b["complex"], 0)
            if k_complex <= 0:
                continue
            lg = b["logb"] * k_complex - pK * bal["reactants"][solid]
            r, pr = dict(bal["reactants"]), dict(bal["products"])
            # 总是正则化 OH-：配体为 OH- 或固体为氢氧化物时方程会产生 OH-
            coeff = _canonicalize_oh(r, pr)
            dh = _hess_dH(k_complex * b["dH"] if "dH" in b else None,
                          bal["reactants"][solid] * e["dH"] if "dH" in e else None)
            add("derived", r, pr, lg + coeff * PKW_298, coeff,
                f"ksp_beta:{solid}/{b['complex']}", dH=dh)
        # Ksp ⊗ pKa：solid + n·H+ ⇌ cation + 共轭酸（logK = pKa(和) − pKsp）
        for pe in T.pka_base.get(an, []):
            bal = _bal([solid, H_ION], [cat, pe["acid"]], [WATER])
            if bal is None:
                continue
            k_solid = bal["reactants"].get(solid, 0)
            if k_solid <= 0 or bal["reactants"].get(H_ION, 0) <= 0:
                continue
            n_acid = bal["products"].get(pe["acid"], 0)
            lg = pe["pka"] * n_acid - pK * k_solid
            r2, pr2 = dict(bal["reactants"]), dict(bal["products"])
            coeff2 = _canonicalize_oh(r2, pr2)
            dh = _hess_dH(-n_acid * pe["dH"] if "dH" in pe else None,
                          k_solid * e["dH"] if "dH" in e else None)
            add("derived", r2, pr2, lg + coeff2 * PKW_298, coeff2,
                f"ksp_pka:{solid}/{pe['acid']}", dH=dh)
        # Ksp ⊗ Ksp：solid1 + anion2 ⇌ solid2 + anion1（同阳离子，logK = pKsp1 − pKsp2）
        for e2 in T.ksp:
            if e2 is e or e2["pair"][0] != cat:
                continue
            an2 = e2["pair"][1]
            bal = _bal([solid, an2], [e2["solid"], an], [WATER])
            if bal is None:
                continue
            k1 = bal["reactants"].get(solid, 0)
            k2 = bal["products"].get(e2["solid"], 0)
            if k1 <= 0 or k2 <= 0:
                continue
            r3, pr3 = dict(bal["reactants"]), dict(bal["products"])
            coeff3 = _canonicalize_oh(r3, pr3)
            dh = _hess_dH(k1 * e["dH"] if "dH" in e else None,
                          -k2 * e2["dH"] if "dH" in e2 else None)
            add("derived", r3, pr3,
                -pK * k1 + e2["pKsp"] * k2 + coeff3 * PKW_298, coeff3,
                f"ksp_ksp:{solid}->{e2['solid']}", dH=dh)
        # Ksp ⊗ Ksp 跨阳离子（同阴离子）：solid1 + cat2 ⇌ solid2 + cat1
        # （logK = pKsp2·k2 − pKsp1·k1）。沉淀转化的另一半：闪锌矿遇 Cu2+
        # 变铜蓝（ZnS+Cu2+→CuS+Zn2+，铜蓝次生富集成因）即此通道；
        # 同阳离子版本只覆盖换阴离子（AgCl→AgBr），跨阳离子换阳离子同样
        # 是沉淀转化的教科书主体
        for e2 in T.ksp:
            if e2 is e or e2["pair"][1] != an:
                continue
            cat2 = e2["pair"][0]
            if cat2 == cat:
                continue
            bal = _bal([solid, cat2], [e2["solid"], cat], [WATER, H_ION])
            if bal is None:
                continue
            k1 = bal["reactants"].get(solid, 0)
            k2 = bal["products"].get(e2["solid"], 0)
            if k1 <= 0 or k2 <= 0:
                continue
            r5, pr5 = dict(bal["reactants"]), dict(bal["products"])
            coeff5 = _canonicalize_oh(r5, pr5)
            dh = _hess_dH(k1 * e["dH"] if "dH" in e else None,
                          -k2 * e2["dH"] if "dH" in e2 else None)
            add("derived", r5, pr5,
                -pK * k1 + e2["pKsp"] * k2 + coeff5 * PKW_298, coeff5,
                f"ksp_ksp_xcat:{solid}->{e2['solid']}", dH=dh)

    # Ksp ⊗ β（归中沉淀）：complex + center ⇌ ν·solid（配体即沉淀阴离子时，
    # logK = ν·pKsp − logβ）。游离中心离子与配离子共存必超饱和
    # （[Ag(CN)2]- + Ag+ → 2AgCN↓，Q/Ksp 可达 1e5），独立分步求解需经
    # 痕量自由配体微步乒乓（数量级瓶颈），Hess 合并候选一步直达终态
    for b in T.beta:
        if b["ligand"] == "OH^-":
            # 羟基络合-氢氧化物沉淀族已有 OH- 通道（ksp_beta 正向 + 酸碱
            # 再平衡），叠加归中沉淀候选会改写既有方程路径而无新增化学
            continue
        for e in T.ksp:
            if e["pair"][0] != b["center"] or e["pair"][1] != b["ligand"]:
                continue
            bal = _bal([b["complex"], b["center"]], [e["solid"]], [WATER, H_ION])
            if bal is None:
                continue
            k_complex = bal["reactants"].get(b["complex"], 0)
            k_solid = bal["products"].get(e["solid"], 0)
            if k_complex <= 0 or k_solid <= 0:
                continue
            lg = k_solid * e["pKsp"] - k_complex * b["logb"]
            r4, pr4 = dict(bal["reactants"]), dict(bal["products"])
            coeff4 = _canonicalize_oh(r4, pr4)
            dh = _hess_dH(-k_complex * b["dH"] if "dH" in b else None,
                          k_solid * e["dH"] if "dH" in e else None)
            add("derived", r4, pr4, lg + coeff4 * PKW_298, coeff4,
                f"ksp_beta_disprop:{e['solid']}/{b['complex']}", dH=dh)

    # β ⊗ pKa：complex + ν·H+ ⇌ center + ν·共轭酸（logK = ν·pKa − logβ）
    for b in T.beta:
        lig = b["ligand"]
        if lig == "OH^-":
            lg = b["nu"] * PKW_298 - b["logb"]
            bal = _bal([b["complex"], H_ION], [b["center"], WATER], [WATER])
            if bal is None or bal["reactants"].get(H_ION, 0) <= 0:
                continue
            out.append(Cand("derived", dict(bal["reactants"]), dict(bal["products"]),
                            lg, b["nu"], dH=-b["dH"] if "dH" in b else None,
                            meta={"src": f"beta_pka:{b['complex']}"}))
        else:
            for pe in T.pka_base.get(lig, []):
                if pe["n"] != 1:
                    continue
                lg = b["nu"] * pe["pka"] - b["logb"]
                bal = _bal([b["complex"], H_ION], [b["center"], pe["acid"]], [WATER])
                if bal is None or bal["reactants"].get(H_ION, 0) <= 0:
                    continue
                k_acid = bal["products"].get(pe["acid"], 0)
                kc = bal["reactants"][b["complex"]]
                dh = _hess_dH(-k_acid * pe["dH"] if "dH" in pe else None,
                              -kc * b["dH"] if "dH" in b else None)
                out.append(Cand("derived", dict(bal["reactants"]), dict(bal["products"]),
                                pe["pka"] * k_acid - b["logb"] * kc,
                                0.0, dH=dh, meta={"src": f"beta_pka:{b['complex']}"}))

    # Ksp⊗弱酸（酸为反应物）：solid + HA -> cat + 共轭碱（弱酸溶蚀沉淀，
    # 如 CaCO3 + CO2 + H2O -> Ca2+ + 2HCO3-，钟乳石/暂时硬水；强酸溶解由
    # ksp_pka 处理，此处只收比阴离子共轭酸更弱的 HA）
    for e in T.ksp:
        solid, cat, an, pK = e["solid"], e["pair"][0], e["pair"][1], e["pKsp"]
        if an == "OH^-":
            continue
        for pe_an in T.pka_base.get(an, []):
            if pe_an["n"] != 1:
                continue
            h_an, pka_hi = pe_an["acid"], pe_an["pka"]
            if h_an == WATER:
                continue
            for ha, entries in T.pka_acid.items():
                if ha in T.solids or ha == WATER:
                    continue
                e1 = min(entries, key=lambda x: x["pka"])
                if e1["pka"] <= 0 or e1["pka"] >= pka_hi:
                    continue
                bal = _bal([solid, ha], [cat, h_an, e1["base"]], [WATER, H_ION])
                if bal is None or H_ION in bal["reactants"] or H_ION in bal["products"]:
                    continue
                k_s = bal["reactants"].get(solid, 0)
                k_ha = bal["reactants"].get(ha, 0)
                if k_s <= 0 or k_ha <= 0:
                    continue
                lg = -pK * k_s + (pka_hi - e1["pka"]) * k_ha
                dh = _hess_dH(k_s * e["dH"] if "dH" in e else None,
                              -k_ha * pe_an["dH"] if "dH" in pe_an else None,
                              k_ha * e1["dH"] if "dH" in e1 else None)
                add("derived", dict(bal["reactants"]), dict(bal["products"]), lg, 0.0,
                    f"ksp_acid:{solid}/{ha}", dH=dh)

    # β_pka ⊗ Ksp(氢氧化物)：complex + (ν/n−1)·center + ν·H2O -> (ν/n)·solid + ν·共轭酸
    # （H+ 抵消形派生：沉淀-解配耦合通道，如 [Cu(NH3)4]2+ + Cu2+ + 4H2O ->
    #   2Cu(OH)2 + 4NH4+。两单独候选各自被 pH 反馈锁死，组合通道直接可达平衡）
    for b in T.beta:
        if b["ligand"] == "OH^-":
            continue
        cell = T.ksp_by_pair.get((b["center"], "OH^-"))
        if cell is None:
            continue
        n_OH = charge_of(b["center"])
        for pe in T.pka_base.get(b["ligand"], []):
            if pe["n"] != 1:
                continue
            bal = _bal([b["complex"], b["center"]], [cell["solid"], pe["acid"]],
                       [WATER, H_ION])
            if bal is None:
                continue
            if H_ION in bal["reactants"] or H_ION in bal["products"]:
                continue   # 只保留 H+ 抵消形
            kc = bal["reactants"].get(b["complex"], 0)
            ks = bal["products"].get(cell["solid"], 0)
            ka = bal["products"].get(pe["acid"], 0)
            if kc <= 0 or ks <= 0 or ka <= 0:
                continue
            lg = (kc * (b["nu"] * pe["pka"] - b["logb"])
                  + ks * (cell["pKsp"] - n_OH * PKW_298))
            dh = _hess_dH(-kc * b["nu"] * pe["dH"] if "dH" in pe else None,
                          -kc * b["dH"] if "dH" in b else None,
                          -ks * cell["dH"] if "dH" in cell else None)
            add("derived", dict(bal["reactants"]), dict(bal["products"]),
                lg, -ks * n_OH, f"beta_ksp:{b['complex']}/{cell['solid']}", dH=dh)

    # 去重
    seen: dict = {}
    for c in out:
        seen.setdefault(c.key, c)
    T._derived = list(seen.values())
    return T._derived


def _canonicalize_oh(r: dict, pr: dict) -> float:
    """把方程中的 OH- 改写为 H2O/H+ 正则形；返回 pKw 系数（logK 修正）。"""
    coeff = 0.0
    k = r.pop("OH^-", 0)
    if k:
        r[WATER] = r.get(WATER, 0) + k
        pr[H_ION] = pr.get(H_ION, 0) + k
        coeff -= k          # logK -= k·pKw
    k = pr.pop("OH^-", 0)
    if k:
        pr[WATER] = pr.get(WATER, 0) + k
        r[H_ION] = r.get(H_ION, 0) + k
        coeff += k          # logK += k·pKw
    for d in (r, pr):
        if d.get(WATER, 0) == 0:
            d.pop(WATER, None)
    return coeff

