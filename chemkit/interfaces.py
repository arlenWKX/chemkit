"""chemkit.interfaces：用户侧 API 的唯一门面（v0.3.3 起）。

本模块回答一个问题：**"作为 chemkit 的使用者，我能调什么？"**——
全部公开对象从这里导入（或直接 `import chemkit` 后使用属性）。
内部模块（system/engine/data/core/thermo）的其余符号均为实现细节，
不承诺跨版本稳定。

    from chemkit.interfaces import Engine, react, System, Reaction

统一入口（v0.3.6 起）：**Engine 对象**——持表一次、长期复用，
三个层级是它的平级方法（无 react()→System()→judge() 层层包装）：

  Engine()                                 进程内建一次、长期运行
      eng.react(substances, ...) -> Reaction     一步式判定
      eng.system(substances, ...) -> System      连续加料体系
      eng.judge(subs, cond) -> dict              引擎 raw dict 直通

函数式便捷入口（与 Engine 方法同语义，隐式使用全局默认表）：

  react(substances, V, T_C, p, ...) -> Reaction
      一步式：直接落求解管线，零中间对象
  System(substances, V, T_C, p, ...)
      反应体系：连续加料、累计再平衡、history 过程保留
      sys.add("NaOH", 0.1) -> Reaction
  Reaction(judge 结果) -> 人类可读对象
      changed / reacted / degree / pH
      consumption / production / initial / final
      net_equation / equations（Equation 结构，读取时渲染） / annotations
      raw（引擎记账层：steps / H_excess / consumption_raw 等）

底层（判定引擎直通）：

  judge(substances, conditions, T) -> dict
      纯函数式：投料+条件 → 引擎 raw dict（无包装）
  load_tables() / Tables / default_tables() / TABLES
      数据表加载（数据与引擎分离；import chemkit 即预加载）

辅助（化学式与方程式）：

  balance(reactants, products, free) -> 配平结果 | None
  FormulaError

模式开关（react/System 与 judge conditions 通用）：
    kinetics    动力学层（slow/gate/T_min/钝化膜等；默认 True）；
                False = 纯热力学基线
    gas_escape  自产气体逸出（开放体系；默认 True）；False = 闭口
    isothermal  恒温假设（**默认 True**）：单遍求解不计热效应，判定
                引擎主用途；False = 绝热耦合：独立温度模块在化学平衡
                与能量平衡间迭代温度不动点，化学结果在自洽终温下给出
                （Reaction.heat_kJ / dT_K / T_final_K / thermal.trace
                含逐轮记录与 converged 标志）
"""
from .data import Tables, load_tables
from .engine import judge
from .system import (Engine, Reaction, System, react, default_tables, TABLES)
from .core import FormulaError, balance

__all__ = ["Tables", "load_tables", "judge", "Engine",
           "Reaction", "System", "react", "default_tables",
           "TABLES", "FormulaError", "balance"]
