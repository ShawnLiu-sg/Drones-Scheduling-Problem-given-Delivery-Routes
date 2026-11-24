# Drones Scheduling Problem given Delivery Routes

本库提供固定卡车访问顺序与随车配送无人机架次时的共享补给调度算法和实验算例。晚释放包裹须由仓库的共享补给无人机在沿途站点交接给所属车辆。模型决定补给分批、交接站、补给机起飞顺序，以及卡车与配送机的同步时刻。

## 实验与入口

| 实验 | 算例 | 使用方法 |
| --- | --- | --- |
| 1　总返库延迟 | 216 个固定车机协同场景 | 在同一输入上分别运行 `full` 和 `liu` |
| 2　算法组件消融 | 主实验中共享补给机数少于配送车数的 120 个场景 | 在同一子集运行 `base`、`no_capacity`、`no_interval`、`full` |
| 3　纯卡车路线 | 24 个 CVRPLIB 场景 | 在相同的纯卡车输入上分别运行 `full` 和 `liu` |
| 4　全程完工时间 | 216 对场景，432 个配置 | `makespan` 同时考虑最后一辆车与最后一架在用补给机的返库时刻 |
| 5　配送机续航余量 | 216 个场景，各取 0%、5%、15%、30%，共 864 个配置 | `endurance` 可先回放本地 15% 参考排程，再在新余量下调度 |
| 6　释放扰动 | 主实验输入和固定排程，四种分布、三档强度 | `release` 对本地已核验排程计算模拟概率和解析概率 |
| 7　释放与补给运力 | 复用主实验场景的单因素配对与机队三元组 | 对本地计算的上下界进行配对 |

另外收录早期方案 A、方案 D 的网络与定价实现、阈值变体、冲突搜索和小规模枚举程序，以及开发阶段的历史算例。算法入口与依赖见 [ALGORITHMS.md](ALGORITHMS.md)。正式实例字段见 [数据说明](datasets/README.md)，历史算例及提取方法见 [历史数据说明](datasets/legacy/README.md)。

## 本地运行

在仓库根目录执行下列命令。建议使用 Python 3.11。优化实验使用 Gurobi 及其有效许可，部分方法还会调用 IBM ILOG CP Optimizer。将 `cpoptimizer` 加入 PATH 或用环境变量 `CP_OPTIMIZER` 指向其可执行文件。

```sh
python -m pip install -r requirements.txt
python datasets/validate_inputs.py
python run_suite.py main full --tier small --first 1
python run_suite.py main liu --tier small --first 1
python run_suite.py main release --tier small --first 1 --plans local_runs/main/full
```

`run_suite.py` 按清单串行运行，默认将本次计算写入 `local_runs/<数据集>/<方法>/`。去掉 `--first` 可运行所选规模的全部输入；去掉 `--tier` 可运行全部规模。小、中规模每例默认 120 秒，大规模 300 秒。用 `--ablation` 选择实验 2 的固定子集。单例也可用 `python run_solver.py --help` 查看各方法与参数。

```sh
python run_suite.py main full
python run_suite.py cvrp liu
python run_suite.py makespan makespan
python run_suite.py endurance endurance --plans local_runs/main/full
python run_suite.py main release --plans local_runs/main/full
python compare_management.py --runs local_runs/main/full --out local_runs/management
```

实验 6 的概率以所用固定排程为条件。本地重新求解可能得到另一份同目标值排程，概率会随补给起飞时刻和分批而变化。批量入口会列出缺少本地排程的场景。正式四类输入清单共含 1,536 个配置，跨实验共享的主场景仍分别保留在相应清单中。

## 数据来源

主实验固定车机路线来自 Ren、Froger、Jabali 和 Liang 的公开归档。纯卡车对照使用六组 CVRPLIB 路线。几何、站序、配送架次及来源标识见 [路线输入](datasets/source_routes.json)。[Ren 等（2024）](https://doi.org/10.1016/j.ejor.2024.05.031) · [CVRPLIB](https://galgos.inf.puc-rio.br/cvrplib/en/instances)
