# solenoid_separator

基于超导螺线管的同位素分离器（Rb 同位素）项目的计算代码、结果与报告页面。
内容同步自项目共享文件夹 `solenoid-separator/`（快照日期 2026-10-26）。

## 目录

| 目录 | 内容 |
|---|---|
| `sim/` | 早期光学模拟（基线、环形源、空间电荷） |
| `source/` | 环形平行引出离子源设计（v2/v3 几何、引出路径） |
| `beamline/` | 束线计算：空间电荷、误差/公差、真空、散射、收集器 |
| `highcurrent/` | 大流强空间电荷与电荷交换分析 |
| `neutralization/` | 空间电荷中和研究（电子温度 Te 模型） |
| `expert_review/` | 专家评审复核：工作点扫描、`newwp/`（机器 B：2 次交叉、R0=104 mm、30 keV、6 mA）、公差 |
| `integrated/` | 统一的端到端光学代码（`sep.py`、`source.py`、`collector.py`、`mc.py` 等）及 `results/` |
| `review/` | 评审意见与修正记录 |
| `report/` | 报告构建脚本与 HTML 中间/最终产物 |
| `*.html` | 早期独立报告页面 |

结果主要为 JSON / Markdown，运行日志为 `.txt` / `.log` / `.out`。

## 未纳入仓库的大文件

以下文件体积过大（或为可重新生成的缓存），已通过 `.gitignore` 排除，只保存在项目共享文件夹中：

| 文件 | 大小 | 说明 |
|---|---|---|
| `beamline/base2785.npz` | 330 MB | 基线 2.785 T 粒子追踪输出 |
| `beamline/base_out.npz` | 330 MB | 基线粒子追踪输出 |
| `expert_review/newwp/recs_h6_1.912_5.npy` | 20 MB | 机器 B 1.912 T、5 mA 轨迹记录 |
| `expert_review/newwp/recs_h6_1.912_6.npy` | 20 MB | 机器 B 1.912 T、6 mA 轨迹记录 |
| `integrated/cache/sol_L1.73205_..._zmin-1.8.pkl` | 13 MB | 螺线管场图缓存（可由 `sep.py` 重新生成） |
| `beamline/sol_cache_5740354900026072187.pkl` | 6.3 MB | 螺线管场图缓存（可重新生成） |
| `beamline/src/src_0p1mA.pkl` | 6.4 MB | 0.1 mA 源粒子分布 |
| `beamline/src/src_0p3mA.pkl` | 6.4 MB | 0.3 mA 源粒子分布 |
| `source/src_v3.pkl` | 6.4 MB | v3 源粒子分布 |

另外排除了 `__pycache__/` 和 `*.pyc`。较小的 `.pkl`（约 0.6–0.7 MB 的源分布等）已纳入仓库。
