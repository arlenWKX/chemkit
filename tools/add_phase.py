"""一次性：给 ksp.json 中水合形态歧义的条目补 `phase` 字段。

背景（V16 的教训）：引擎账本只有**无水式**物种 + 水（溶剂活度 1），同一
固相在不同晶型/水合态下的 Ksp 是不同的实测量。数据里若只写无水式而不声明
用的是**哪一个相**的值，读数据的人无法判断，边界用例（析出量小到 ±0.5 个
pKsp 就翻转）也无法复核。

本脚本只**补充声明**，不改任何 pKsp 数值——换值是化学决策，须逐条定。
用法：python tools/add_phase.py [--write]
"""
from __future__ import annotations

import json
import sys

PATH = "chemkit/data/ksp.json"

PHASE = {
    "CaSO_3": "无水 CaSO3（pKsp 7.2）。**二水 CaSO3·2H2O 为 6.5**——V16 "
              "SO2+CaCl2 正落在此二值之间（实测 pQ=7.207），故该用例的"
              "「是否析出」由选相决定，不选相则无意义",
    "CaSO_4": "二水石膏 CaSO4·2H2O（4.6）；硬石膏 CaSO4 为 4.4–5.0。"
              "账本用无水式书写，Ksp 按**二水相**取",
    "CaC_2O_4": "一水 CaC2O4·H2O（8.6）；无水为 8.0。账本用无水式，"
                "Ksp 按**一水相**取",
    "MgCO_3": "三水 MgCO3·3H2O（≈5.0）；菱镁矿 MgCO3 为 5.2–8。"
              "本表取水合相（溶液中析出的是水合碳酸镁）",
    "CoCO_3": "无水 CoCO3（12.84，CRC）；水合形态约 11。本表取无水值",
    "NiCO_3": "无水 NiCO3（8.0，CRC）；水合形态约 6.9。本表取无水值",
    "MnCO_3": "无水菱锰矿 MnCO3（10.7）",
    "FeCO_3": "无水菱铁矿 FeCO3（10.5）",
    "Fe(OH)_3": "无定形水合氧化铁 Fe(OH)3·xH2O（38.0）——新鲜沉淀的相；"
                "老化成针铁矿/赤铁矿后溶解度低几个数量级，本表取无定形"
                "（与实验「Fe3+ 加碱得红褐胶状沉淀」一致）",
    "Al(OH)_3": "无定形 Al(OH)3·xH2O（33.0）——新鲜沉淀的相；三水铝石"
                "（gibbsite）33.5、勃姆石（boehmite）更高温稳定。"
                "本表取无定形",
    "Zn(OH)_2": "无定形 Zn(OH)2（16.9）；ε-Zn(OH)2 晶态约 15.5。本表取无定形",
    "Mg(OH)_2": "六方 brucite（无水，11.2）——氢氧化镁不形成稳定水合物",
    "BaSO_4": "重晶石（无水，10.0）——硫酸钡无稳定水合物",
    "Li_2CO_3": "无水（2.77）——碳酸锂室温即无水相，590 ℃ 以上才熔化",
    "NaHCO_3": "无水（-0.08，表「极易溶，仅同离子效应下析出」）——"
               "碳酸氢钠无稳定水合物",
}


def _detect_indent(text: str) -> int:
    """探测原文缩进宽度——各数据文件约定不同（ksp.json 是 2，beta/tests 是
    1），写回时若不匹配会把整文件重排版，淹没真实 diff。"""
    for line in text.splitlines():
        if line.strip() == "{":
            return len(line) - len(line.lstrip(" "))
    return 1


def main(write: bool) -> None:
    raw = open(PATH, encoding="utf-8", newline="").read()
    ind = _detect_indent(raw)
    rows = json.loads(raw)
    n = 0
    for e in rows:
        ph = PHASE.get(e["solid"])
        if ph and e.get("phase") != ph:
            e["phase"] = ph
            n += 1
    print(f"已为 {n} 条补 phase 字段（共 {len(rows)} 条；探测缩进 indent={ind}）")
    if not write:
        print("（--write 写回）")
        return
    # newline="" 保留 LF：Python 文本模式默认把 \n 翻成 CRLF，会让整文件
    # 行尾全变（ksp.json 首版即此，4057 行假 diff）
    with open(PATH, "w", encoding="utf-8", newline="") as f:
        f.write(json.dumps(rows, ensure_ascii=False, indent=ind) + "\n")
    print(f"已写回 {PATH}")


if __name__ == "__main__":
    main("--write" in sys.argv)
