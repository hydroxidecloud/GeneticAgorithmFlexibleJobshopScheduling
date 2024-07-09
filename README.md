# 作业车间调度问题优化项目

这是一个基于遗传算法的作业车间调度问题（Job Shop Scheduling Problem, Jm||Cmax）的开源项目。本项目旨在优化调度算法，以最小化最大完工时间（Makespan），适用于标准测试算例和实际案例。

## 项目简介

作业车间调度问题（Job Shop Scheduling Problem, Jm||Cmax）是调度领域的一个重要研究课题。本项目设计并实现了一种基于遗传算法的调度算法，旨在优化最大完工时间（Makespan）。项目包含标准测试算例（如ta01, ta40, ta60）和实际案例的调度优化。

### 系统描述

1. 系统包含𝑚台机器，需要完成𝑛个工作。
2. 每个作业需依次完成𝑚道工序，每道工序需在一台机器上加工，加工时间为常数。
3. 每台机器最多同时加工一个作业。
4. 假设机器的准备时间（如装夹和拆卸等）已包括在加工时间中，无需另外考虑。
5. 假设机器间的运输时间可以忽略不计。

## 项目结构

- `GA_jobshop_makespan.py`: 用于计算标准测试算例（如ta01, ta40, ta60）的最大完工时间。
- `GA_jobshop_realcase.py`: 用于计算实际案例的最大完工时间。
- `sequence_makespan_calculator.py`: 用于计算某一种标准测试算例染色体序列的总工时。
- `sequence_realcase_calculator.py`: 用于计算某一种实际案例染色体序列的总工时。

## 使用说明

### 计算标准测试算例的最大完工时间

存储ta01, ta40, ta60于~/instance/中，运行 `GA_jobshop_makespan.py` 脚本以优化标准测试算例的最大完工时间。
```bash
python GA_jobshop_makespan.py
```

### 计算测试算例的最大完工时间

存储parts.csv于~/instance/中，运行 `GA_jobshop_realcease.py` 脚本以优化标准测试算例的最大完工时间。
```bash
python GA_jobshop_realcase.py
```
