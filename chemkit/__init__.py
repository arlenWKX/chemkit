"""chemkit：水溶液反应判定引擎（离子基、高中~竞赛水平）。

零第三方依赖（v0.3.5 起，含精确配平）；import 即自动加载数据表
（interfaces.TABLES 单例）。

用户侧 API 的唯一门面是 chemkit.interfaces（v0.3.3 起）；v0.3.7 起
统一入口为 **Engine 对象**——react/system 两个层级是它的平级方法，
无层层包装：

    import chemkit
    eng = chemkit.Engine()                              # 持表一次、长期复用
    r = eng.react({"Zn": 1.0, "H_2SO_4": 1.0})          # 一步式
    sys = eng.system(V=1.0)                             # 连续加料（共享缓存）
    r = sys.add("NaOH", 0.1)

函数式便捷入口同语义（隐式使用全局默认表）：

    r = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0)
    sys = chemkit.System(V=1.0, T=298.15, p=101.3)
    r = sys.add("HCl", 0.15)

    r.changed, r.reacted, r.degree, r.pH                 # 判定层
    r.consumption, r.production, r.initial, r.final      # 人类可读层
    r.net_equation, r.equations    # Equation 结构，读取时 .tex()/.plain()
    r.consumption_raw, r.H_excess_raw                    # 引擎记账层

引擎 raw dict（steps / H_excess / cond / override 等）经 **`r.raw`** 取；
用户侧不再单列 `judge()`（v0.5.2 删去与 `react` 重叠的重复入口：
两者走同一条求解管线，raw 层从 Reaction 上取即可）。

辅助 API（化学式与方程式）：
    from chemkit import FormulaError, balance
    balance(["H_2", "O_2"], ["H_2O"])  # -> {"reactants": {...}, "products": {...}}

测试套件为独立模块 chemkit.testsuit（默认不加载）：
    python -m chemkit.testsuit        # 运行 data/tests.json 全部用例 + 环闭合检查
"""
from .interfaces import (
    Tables, load_tables, Engine,
    Reaction, System, react, default_tables, TABLES,
    FormulaError, balance,
)

__version__ = "0.5.22"
__all__ = ["Tables", "load_tables", "Engine",
           "Reaction", "System", "react", "default_tables",
           "TABLES", "FormulaError", "balance", "__version__"]
