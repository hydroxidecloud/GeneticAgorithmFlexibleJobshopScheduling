# Job Shop Scheduling Problem Optimization Project

This repository contains an open-source project focused on optimizing the Job Shop Scheduling Problem (Jm||Cmax) using Genetic Algorithm. The primary objective is to minimize the Makespan, applicable to both standard test cases and real-world scenarios.

## Project Overview

The Job Shop Scheduling Problem (Jm||Cmax) is a significant research area in scheduling. This project implements a Genetic Algorithm-based scheduling algorithm to optimize the Makespan. It includes optimizations for standard test cases (e.g., ta01, ta40, ta60) and real-world scenarios.

### System Description

1. The system involves 𝑚 machines and 𝑛 jobs to be completed.
2. Each job consists of 𝑚 operations processed sequentially on different machines, each with a fixed processing time.
3. Each machine can handle only one job at a time.
4. Setup times (e.g., mounting and dismounting) are included in the processing times.
5. Transport times between machines are negligible.

## Project Structure

- `GA_jobshop_makespan.py`, `GA_jobshop_makespan.rb`: Calculates the Makespan for standard test cases (e.g., ta01, ta40, ta60).
- `GA_jobshop_realcase.py`, `GA_jobshop_realcase_refactored.py`: Calculates the Makespan for real-world scenarios.
- `sequence_makespan_calculator.py`: Computes the total processing time for a chromosome sequence of a standard test case.
- `sequence_realcase_calculator.py`: Computes the total processing time for a chromosome sequence of a real-world scenario.

## Usage Instructions

### Generate Gantt and Calculate Makespan for Standard Test Cases

Place ta01, ta40, ta60 in the ./instance/ directory. Run the following command to optimize the Makespan for standard test cases:

```bash
python GA_jobshop_makespan.py
```
or

```bash
ruby --jit GA_jobshop_makespan.rb
```

### Generate Gantt and Calculate Makespan for Real-world Scenarios

Place parts.csv in the ./instance/ directory. Run the following command to optimize the Makespan for real-world scenarios:

```bash
python GA_jobshop_realcase.py
```
or

```bash
python GA_jobshop_realcase_refactored.py
```

### Calculate Makespan for Existed Standard Test Case Genetic Sequence

```bash
python sequence_makespan_calculator.py
```

### Calculate Makespan for Existed Real-world Scenarios Genetic Sequence

```bash
python sequence_realcase_calculator.py
```

