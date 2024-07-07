import datetime
import random
# import time
# import copy

import numpy as np
import plotly.figure_factory as ff
# import matplotlib.pyplot as plt
# from tqdm import tqdm
# import pandas as pd

instance_name = str(input('Please input the instance name (default value ta01): ') or 'ta01')

with open("instance/" + instance_name, "r") as file:
    lines = file.readlines()

num_pt = int(lines[0].split()[3])
num_mc = int(lines[0].split()[7])
num_gene = num_mc * num_pt  # number of genes in a chromosome

pt_start = lines.index('Times\n') + 1
pt = []
for i in range(num_pt):
    pt.append(list(map(int, lines[pt_start + i].split())))

ms_start = lines.index('Machines\n') + 1
ms = []
for i in range(num_pt):
    ms.append(list(map(int, lines[ms_start + i].split())))

# 获取命令行输入的字符串
input_string = input("Please input sequence_best (eg.[11, 2, 14,...]): ")
# 去掉字符串两端的方括号并分割成列表
string_list = input_string.strip('[]').split(',')
# 去掉每个元素的多余空格并转换为整数
sequence_best = list(map(lambda x: int(x.strip()), string_list))
# sequence_best = list(input('Please input sequence_best (eg.[11, 2, 14,...]): '))

'''--------plot gantt chart-------'''

m_keys = [j + 1 for j in range(num_mc)]
j_keys = [j for j in range(num_pt)]
key_count = {key: 0 for key in j_keys}
j_count = {key: 0 for key in j_keys}
m_count = {key: 0 for key in m_keys}
j_record = {}
for i in sequence_best:
    gen_t = int(pt[i][key_count[i]])
    gen_m = int(ms[i][key_count[i]])
    j_count[i] = j_count[i] + gen_t
    m_count[gen_m] = m_count[gen_m] + gen_t

    if m_count[gen_m] < j_count[i]:
        m_count[gen_m] = j_count[i]
    elif m_count[gen_m] > j_count[i]:
        j_count[i] = m_count[gen_m]

    start_time = j_count[i] - pt[i][key_count[i]]  # convert seconds to hours, minutes and seconds
    end_time = j_count[i]

    j_record[(i, gen_m, key_count[i] + 1)] = [start_time, end_time]

    key_count[i] = key_count[i] + 1


def generate_random_color():
    color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
    return color


colors = [generate_random_color() for _ in range(num_pt)]

start_date = datetime.datetime(2024, 6, 28)
df = []
for (job, machine, procedure), (start, end) in j_record.items():
    start_time = start_date + datetime.timedelta(minutes=start)
    end_time = start_date + datetime.timedelta(minutes=end)
    df.append(dict(Task='Machine %s' % machine,
                   Start=start_time.strftime('%Y-%m-%d %H:%M:%S'),
                   Finish=end_time.strftime('%Y-%m-%d %H:%M:%S'),
                   Resource='Job %s' % (job + 1)))
df = sorted(df, key=lambda x: int(x['Task'].replace('Machine ', '')))
fig = ff.create_gantt(df,
                      index_col='Resource',
                      colors=colors,
                      show_colorbar=True,
                      group_tasks=True,
                      showgrid_x=True,
                      title='Flexible Job shop Schedule')
fig.show()
# pio.plot(fig, filename='GA_jobshop_realcase.html', auto_open=False)


Tbest = max(j_count.values())
print("Optimal Value: %f" % Tbest)


def write_solution_to_file(solution, filename):
    # solution_t = solution.T
    with open(filename, 'w', encoding='utf-8') as file:
        file.write('Nb of jobs ' + str(num_pt) + ' Nb of Machines ' + str(num_mc) + '\n')
        file.write(str(Tbest) + '\n')
        file.write('Solutions' + '\n')
        for row in solution:
            file.write('\t\t'.join(map(str, row)) + '\n')


def create_solution_matrix(sequence_best, num_pt, num_mc, ms):
    # 初始化Solution矩阵，大小为num_mc * num_pt，初始值为-1表示未分配
    solution = np.full((num_mc, num_pt), -1)

    # 用于记录每个工件的工序次数
    job_operation_count = {i: 0 for i in range(num_pt)}

    for operation in sequence_best:
        job = operation  # 工件编号
        op_num = job_operation_count[job]  # 当前工件的工序编号
        machine = ms[job][op_num] - 1  # 获取工件的工序对应的机器编号并转换为索引
        job_operation_count[job] += 1  # 增加工件的工序次数

        # 将该工件的当前工序放入Solution矩阵对应机器的下一空闲位置
        for idx in range(num_pt):
            if solution[machine][idx] == -1:
                solution[machine][idx] = job
                break
    solution = solution + 1
    return solution


solution_best = create_solution_matrix(sequence_best, num_pt, num_mc, ms)
write_solution_to_file(solution_best, instance_name + '_sol.txt')
