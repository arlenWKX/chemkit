"""chemkit.system：高层 API（Engine/react/System 与 Reaction 包装）。

v0.3.6 起用户侧统一入口是 **Engine 对象**：持表一次、长期复用，
judge / react / system 三个层级是它的**平级方法**——没有 react()→
System()→judge() 的层层包装（一步式不再为了调 judge 而建一个带
history 的 System；System 也不再经由 react 之类的上层函数中转）。
依赖链：engine（judge）/equations（方程式）/thermo（绝热耦合）。
"""
from __future__ import annotations

from .data import Tables
from .engine import judge
from . import thermo
from .equations import (_parse_equation, _balance_h_o, _restore_oh,
                        _fmt_term, _rationalize, _try_integer_snap,
                        make_equation, Equation,
                        _step_key, _collect_steps,
                        _collapse_transient_intermediates, _step_to_ionic,
                        _build_equations, _balanced_quick, _filter_trace,
                        _SPECIES_NORMALIZE, _TERM_RE, _SIDE_SEP,
                        WATER, H_ION, OH_ION, _TRACE, STEP_MIN,
                        TABLES, default_tables)


class Reaction:
    """一次反应的结果（一次投料平衡后的完整快照）。

    人类可读层（化学习惯，含 H+/OH-/H2O）：
        changed        体系是否发生显著净变化（bool，净账本判定；含纯溶解/
                       弱酸弱碱电离/水解等单投料形态变化——较宽口径）
        reacted        是否发生狭义化学反应（bool：氧化还原/中和/跨投料的
                       沉淀配位等；纯溶解、电离、水解等单投料形态变化为 False）
        degree         程度整数：2=完全反应 / 1=可逆（部分）反应 /
                       0=难反应或未反应
        consumption    净消耗 {化学式: mol}
        production     净生成 {化学式: mol}
        initial        初态组成 {化学式: mol}（post-normalize：强电解质已电离、
                       气体如 SO3 已与水反应为酸、酸碱中和已记账；H2O 为溶剂不入）
        final          终态组成 {化学式: mol}（H2O 为溶剂不入）
        pH             终态 pH（None 表示不适用，如 OVERRIDE 路径）
        net_equation   总净离子反应方程式（**Equation 结构** | None）；无显著
                       反应时为 None。**读取时**才渲染：
                       `str(eq)`/`eq.tex()` = TeX（Markdown `$$…$$` 可显示），
                       `eq.plain()` = 不带上下标标注的纯字符串 fallback
        equations      多步离子方程式列表（Equation 结构，按贡献降序）
        annotations    标注列表（slow / blocked 等）
        override       命中的 OVERRIDE id（或 None）
        escaped        逸出气相 {化学式: mol}（泡点扫气：超过 H(T)·p_ext
                       溶解上限的自产气体，逸出即消失，不再参与反应；
                       final 中同种气体只剩溶解态，production 含逸出部分）

    raw 属性（引擎原始记账，H2O 不入账、H+/OH- 合记为 H_excess）：
        consumption_raw  消耗 {化学式: mol}
        production_raw   生成 {化学式: mol}
        initial_raw      初态（post-normalize 引擎记账，H+/OH- 合为 H_excess_initial）
        final_raw        终态 {化学式: mol}
        H_excess_raw     终态带符号质子账本（正 = 残余游离强酸 mol，
                          负 = 残余游离强碱 mol）
        steps            逐步反应过程（kind/equation/logK/S/extent/conversion）
        raw              judge() 原始 dict（备用，不鼓励直接读取）

    热效应层（独立温度模块 thermo.py，与平衡求解完全解耦）：
        heat_kJ      放热为正（= −ΔH；None = 数据不足或恒温模式）
        dT_K         溶液温升（以水的比热容、V×1000 g 计；负 = 降温；
                     绝热耦合口径 = T_final − T₀）
        T_final_K    终温（绝热耦合：外层温度不动点的收敛值）
        thermal      完整温度分析 dict（heat/dT/T_final/water_mol/mass_g/
                     missing/flags/reason；绝热耦合模式另含 trace 逐轮
                     记录与 converged 收敛标志，详见 thermo.py 模块文档）
        isothermal=True（默认）单遍求解不计热效应，
        thermal={"mode": "isothermal"}，heat/dT 为 None；
        isothermal=False 绝热耦合：化学平衡在自洽终温下求解。
    """
    __slots__ = (
        # 人类可读层
        "changed", "reacted", "degree",
        "consumption", "production", "initial", "final",
        "pH", "net_equation", "net_equation_raw", "equations",
        "annotations", "override", "escaped",
        # 热效应层（独立温度模块）
        "heat_kJ", "dT_K", "T_final_K", "thermal",
        # 引擎记账层（raw）
        "consumption_raw", "production_raw", "initial_raw", "final_raw",
        "H_excess_raw", "raw",
    )

    def __init__(self, r: dict, tables: "Tables | None" = None):
        """包装 judge()/thermo.coupled() 的结果 dict。

        tables：热效应回退分析用的数据表（直接包装 judge 结果且
        isothermal=False 时，Reaction 需要数据表查 ΔHf；未提供则用
        模块级 TABLES 单例）。经 System/react 创建时自动传入。"""
        self.raw = r
        self.changed: bool = r["changed"]
        self.reacted: bool = r["reacted"]
        self.degree: int = r["degree"]
        self.pH: float | None = r["final_pH"]
        self.annotations: list[str] = list(r["annotations"])
        self.override: str | None = r.get("override")
        self.escaped: dict[str, float] = {e["name"]: e["mol"]
                                          for e in r.get("escaped", [])}

        # ---- raw（引擎记账）----
        self.consumption_raw: dict[str, float] = {e["name"]: e["mol"]
                                                 for e in r["consumption"]}
        self.production_raw: dict[str, float] = {e["name"]: e["mol"]
                                                for e in r["production"]}
        self.initial_raw: dict[str, float] = {e["name"]: e["mol"]
                                              for e in r.get("initial", [])}
        self.final_raw: dict[str, float] = {e["name"]: e["mol"] for e in r["final"]}
        self.H_excess_raw: float = r.get("H_excess", 0.0)

        # ---- 人类可读层（化学习惯）----
        # net_equation 是**精编版**（痕量副过程不叙述）；net_equation_raw 是
        # **原始版**（账本净差的格式化，永远给出）。两者不一致 ⟹ 该体系只发生
        # 了痕量副过程（如 1 M CuSO₄ 在自身 Ksp 线上析出 5e-4 mol Cu(OH)₂）。
        # 结构化对象（Equation：left/right/arrow/kind），字符串只是渲染视图
        (self.equations, self.consumption, self.production,
         self.net_equation, self.net_equation_raw) = \
            _build_equations(r["steps"], r)
        # **不做预先转换**（用户口径）：这里只存 Equation 结构；需要文本时
        # 由读取方在读取点转换——`str(eq)`/`eq.tex()` = TeX（Markdown $$ 块可
        # 直接显示），`eq.plain()` = 无标记纯字符串 fallback。

        # OVERRIDE 路径无 steps，退化为 raw 层面的化学式方程式
        if (not self.consumption and not self.production and self.consumption_raw
                and self.override):
            self.consumption = dict(self.consumption_raw)
            self.production = dict(self.production_raw)
            self.net_equation = make_equation(self.consumption, self.production)
            self.net_equation_raw = self.net_equation
            if self.net_equation is not None and not self.equations:
                self.equations = [self.net_equation]

        # 终态组成：raw final + H+/OH-（H2O 为溶剂不入）
        self.final: dict[str, float] = dict(self.final_raw)
        if self.H_excess_raw > 1e-6:
            self.final[H_ION] = self.final.get(H_ION, 0.0) + self.H_excess_raw
        elif self.H_excess_raw < -1e-6:
            self.final[OH_ION] = self.final.get(OH_ION, 0.0) + (-self.H_excess_raw)

        # 初态组成：raw initial + H+/OH-（H2O 为溶剂不入）。raw initial 已含
        # post-normalize 的所有物种（强电解质电离、SO3+H2O→HSO4-+H+ 等），
        # 仅需补出 H+/OH-（来自 H_excess_initial）
        self.initial: dict[str, float] = dict(self.initial_raw)
        He0 = r.get("H_excess_initial", 0.0)
        if He0 > 1e-6:
            self.initial[H_ION] = self.initial.get(H_ION, 0.0) + He0
        elif He0 < -1e-6:
            self.initial[OH_ION] = self.initial.get(OH_ION, 0.0) + (-He0)

        # ---- 热效应层（独立温度模块；不参与上述任何平衡/方程逻辑）----
        # isothermal=True（默认）恒温：不计热效应，单遍求解（判定引擎的
        # 主用途）；False 绝热耦合：热层在外层迭代温度不动点，化学平衡在
        # 自洽终温下求解（见 thermo.coupled）。调用隔离在 thermo 单点，
        # 异常不拖垮 Reaction 其余属性。
        cd = r.get("cond") or {}
        if cd.get("isothermal"):
            self.thermal = {"mode": "isothermal"}
            self.heat_kJ = None
            self.dT_K = None
            self.T_final_K = float(cd.get("T_K", 298.15))
        else:
            th = r.get("thermal")
            if th is None:
                # 直接包装 judge() 结果（未经 thermo.coupled）时回退单遍分析；
                # ΔHf 查表需 Tables（v0.3.5 起热层与数据表单一事实源）
                th = thermo._safe_analyze(
                    r, tables if tables is not None else TABLES)
                th.setdefault("converged", None)
            self.thermal = th
            self.heat_kJ = th.get("heat_kJ")
            self.dT_K = th.get("dT_K")
            self.T_final_K = th.get("T_final_K",
                                    float(cd.get("T_K", 298.15)))

    # ---- 便捷转发 ----
    @property
    def steps(self) -> list[dict]:
        """逐步反应过程（引擎记账）。每项含 kind/equation/logK/S/extent/
        conversion 等字段。"""
        return self.raw["steps"]

    # ---- 魔术方法 ----
    def __bool__(self) -> bool:
        return self.changed

    def __repr__(self) -> str:
        pro = ", ".join(f"{k}×{v:.3g}" for k, v in self.production.items())
        tag = "反应" if self.reacted else ("变化" if self.changed else "无变化")
        return f"<Reaction {tag} degree={self.degree} [{pro}] pH={self.pH}>"


# ============================================================ 共用求解管线

def _solve(subs: list[dict], cond: dict, tables: Tables,
           isothermal: bool) -> dict:
    """唯一求解入口（v0.3.6 起三条 API 层级共用）：
    isothermal → engine.judge 单遍；False → thermo.coupled 绝热耦合
    （外层温度不动点，温度反馈真实进入求解）。
    react / System.add / Engine 三者都直落此处，无中间包装层。"""
    if isothermal:
        return judge(subs, cond, tables)
    return thermo.coupled(subs, cond, tables)


def _cond(V: float, T: float, T_C: float | None, p: float,
          isothermal: bool, kinetics: bool, gas_escape: bool) -> dict:
    """关键字参数风格 → judge conditions（三个模式开关一并写入）。"""
    return {"V_L": float(V),
            "T_K": float(T) if T_C is None else float(T_C) + 273.15,
            "p_kpa": float(p),
            "isothermal": bool(isothermal),
            "kinetics": bool(kinetics),
            "gas_escape": bool(gas_escape)}


def _subs(substances: dict[str, float]) -> list[dict]:
    """投料 dict → judge 的 list[{name, mol}]（顺序稳定，直接展开）。"""
    return [{"name": n, "mol": float(m)} for n, m in substances.items()]


# ============================================================ Engine

class Engine:
    """判定引擎对象：用户侧统一入口（v0.3.6 起）。

    进程内建一次、长期复用：持有数据表与全部构建缓存（电对模板、
    派生候选、枚举缓存均随 Tables 存活，同温度判定命中热路径），
    三个 API 层级是它的平级方法，无层层包装：

        eng = chemkit.Engine()                        # 默认包内数据表
        r = eng.react({"NaOH": 0.1, "HCl": 0.1})      # 一步式
        sys = eng.system(V=1.0, T_C=25)               # 连续投料（共享缓存）
        raw = eng.judge(subs, {"V_L": 1.0})           # 引擎 raw dict 直通

    参数：
        tables  数据表（默认包内 TABLES 单例；自定义表经 load_tables
                构建，同一 Engine 换表需新建）

    模块级 react/System/judge 函数与此对象方法完全同语义——只是
    隐式使用全局默认表；需要换表/隔离缓存时用 Engine 对象。
    """

    def __init__(self, tables: Tables | None = None):
        self.tables: Tables = tables if tables is not None else TABLES

    # ---- 底层直通：judge（纯函数，raw dict）----
    def judge(self, substances: list[dict],
              conditions: dict | None = None) -> dict:
        """投料 + conditions → 引擎 raw dict（无包装；conditions 键
        V_L/T_K/T_C/c_H/c_OH/pH/p_kpa/isothermal/kinetics/gas_escape）。"""
        return judge(substances, conditions, self.tables)

    # ---- 一步式：react ----
    def react(self, substances: dict[str, float], *,
              V: float = 1.0, T: float = 298.15, T_C: float | None = None,
              p: float = 101.3, tables: Tables | None = None,
              isothermal: bool = True, kinetics: bool = True,
              gas_escape: bool = True) -> "Reaction":
        """一步式反应：直接调求解管线并包装 Reaction（不建 System）。
        参数与模块级 react() 一致（tables 可临时换表，默认本表）。"""
        tb = tables if tables is not None else self.tables
        r = _solve(_subs(substances),
                   _cond(V, T, T_C, p, isothermal, kinetics, gas_escape),
                   tb, isothermal)
        return Reaction(r, tb)

    # ---- 连续投料体系：system ----
    def system(self, substances: dict[str, float] | None = None, *,
               V: float = 1.0, T: float = 298.15, T_C: float | None = None,
               p: float = 101.3, tables: Tables | None = None,
               isothermal: bool = True, kinetics: bool = True,
               gas_escape: bool = True) -> "System":
        """建立绑在本引擎上的 System（共享数据表与缓存）。"""
        return System(substances, V=V, T=T, T_C=T_C, p=p,
                     tables=tables if tables is not None else self.tables,
                     isothermal=isothermal, kinetics=kinetics,
                     gas_escape=gas_escape)

    def __repr__(self) -> str:
        return f"<Engine tables={self.tables!r}>"


# ============================================================ System

class System:
    """反应体系：固定体积/温度/气压，支持连续投料。

    参数：
        substances  初始投料 {化学式: mol}，默认空（纯水）
        V           溶液体积（L），默认 1.0
        T           温度（K），默认 298.15（25 °C）
        T_C         温度（°C），若给定则覆盖 T
        p           外界气压（kPa），默认 101.3；影响气体逸出阈值
        tables      自定义数据表（默认用模块级 TABLES）
        isothermal  恒温假设（bool，默认 True）。True 时单遍求解不计热
                    效应（判定引擎主用途）；False 时绝热耦合求解：独立
                    温度模块在平衡与能量平衡间迭代温度不动点，化学结果
                    在自洽终温下给出（反应热按水比热容计温升，见
                    Reaction.heat_kJ/dT_K/T_final_K/thermal.trace）
        kinetics    动力学层开关（bool，默认 True）。False = 纯热力学基线
                    （无限时间）：slow/gate/膜封锁等动力学标记一律不生效
        gas_escape  自产气体逸出开关（bool，默认 True）。False = 闭口体系：
                    反应产生的气体不逸出（保留在溶液账本参与平衡，
                    相当于密闭容器；逸出账户 escaped 为空）

    建立时（若给了 substances）与每次 add() 自动触发反应——按累计投入量
    重新平衡（化学上等价于连续投料的再平衡），返回本次 Reaction。
    全部历史结果保存在 history。

    示例：
        sys = chemkit.System(V=1.0)
        sys.add("NaOH", 0.1)              # 纯水 + 0.1 mol NaOH
        r = sys.add("HCl", 0.15)          # 再投入 HCl，体系重新平衡
        r.net_equation, r.pH
    """

    def __init__(self,
                 substances: dict[str, float] | None = None,
                 *,
                 V: float = 1.0,
                 T: float = 298.15,
                 T_C: float | None = None,
                 p: float = 101.3,
                 tables: Tables | None = None,
                 isothermal: bool = True,
                 kinetics: bool = True,
                 gas_escape: bool = True):
        self.V_L: float = float(V)
        self.T_K: float = float(T) if T_C is None else float(T_C) + 273.15
        self.p_kpa: float = float(p)
        self.isothermal: bool = bool(isothermal)
        self.kinetics: bool = bool(kinetics)
        self.gas_escape: bool = bool(gas_escape)
        self._tables: Tables = tables if tables is not None else TABLES
        self._feeds: dict[str, float] = {}
        self.history: list[Reaction] = []
        self.result: Reaction | None = None
        if substances:
            for name, mol in substances.items():
                self._feeds[name] = self._feeds.get(name, 0.0) + float(mol)
            self._react()

    def add(self, name: str, mol: float) -> Reaction:
        """加入物质（mol），自动触发反应，返回本次 Reaction。"""
        self._feeds[name] = self._feeds.get(name, 0.0) + float(mol)
        return self._react()

    def _react(self) -> Reaction:
        # v0.3.6：System 直落共用求解管线（不经 react/System 包装链）
        cond = {"V_L": self.V_L, "T_K": self.T_K, "p_kpa": self.p_kpa,
                "isothermal": self.isothermal,
                "kinetics": self.kinetics,
                "gas_escape": self.gas_escape}
        subs = [{"name": n, "mol": m} for n, m in self._feeds.items()]
        r = _solve(subs, cond, self._tables, self.isothermal)
        self.result = Reaction(r, self._tables)
        self.history.append(self.result)
        return self.result

    @property
    def feeds(self) -> dict[str, float]:
        """累计投料 {化学式: mol}（副本，外部修改不影响内部状态）。"""
        return dict(self._feeds)

    def __repr__(self) -> str:
        fd = ", ".join(f"{k}×{v:.3g}" for k, v in self._feeds.items())
        return f"<System [{fd}] V={self.V_L}L T={self.T_K}K p={self.p_kpa}kPa>"


# ============================================================ 一步式 API

def react(substances: dict[str, float],
          *,
          V: float = 1.0,
          T: float = 298.15,
          T_C: float | None = None,
          p: float = 101.3,
          tables: Tables | None = None,
          isothermal: bool = True,
          kinetics: bool = True,
          gas_escape: bool = True) -> Reaction:
    """一步式反应：直接调求解管线，返回 Reaction（v0.3.6 起不再
    经由 System 包装链——一次判定零中间对象）。

    参数与 System 一致（substances 必填）；isothermal/kinetics/gas_escape
    三个模式开关的含义见 System 文档。需要换表/长期持有缓存时用
    chemkit.Engine 对象。

    示例：
        r = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)
        print(r.net_equation.plain())   # '2H+ + Zn -> H2 + Zn2+'（纯字符串）
        print(r.net_equation.tex())     # Markdown $$ 块可直接显示
        rt = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0,
                           isothermal=False)
        print(rt.heat_kJ, rt.dT_K)  # 绝热耦合：放热与温升
    """
    tb = tables if tables is not None else TABLES
    r = _solve(_subs(substances),
               _cond(V, T, T_C, p, isothermal, kinetics, gas_escape),
               tb, isothermal)
    return Reaction(r, tb)

