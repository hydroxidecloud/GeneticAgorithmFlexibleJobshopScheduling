import numpy as np
import time
import copy
import random
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import datetime
from tqdm import tqdm
import pandas as pd

# from collections import Counter
# import plotly.offline as pio

''' ================= initialization setting ======================'''

instance_name = 'realcase'

df = pd.read_csv('instance/parts.csv', header=1)

columns = df.columns.tolist()
columns[0] = 'Part'
columns[1] = 'Quantity'
columns[2] = 'Feature'
columns[3] = 'Process'
df.columns = columns
# df.columns = ['Part', 'Quantity', 'Feature', 'Process', 'C1', 'C2', 'C3', 'C4', 'L1', 'PM1', 'NYM1', 'WYM1', 'RCL1',
#               'RCL2', 'WX1', 'LX1', 'LX2', 'LX3', 'Z1', 'Z2', 'DHH1', 'XQG1', 'JC1', 'QG1']

machine_mapping = {col: idx + 1 for idx, col in enumerate(df.columns[4:])}
# machine_mapping = {
#     'C1': 1, 'C2': 2, 'C3': 3, 'C4': 4, 'L1': 5, 'PM1': 6, 'NYM1': 7, 'WYM1': 8, 'RCL1': 9, 'RCL2': 10,
#     'WX1': 11, 'LX1': 12, 'LX2': 13, 'LX3': 14, 'Z1': 15, 'Z2': 16, 'DHH1': 17, 'XQG1': 18, 'JC1': 19, 'QG1': 20
# }
id_to_machine = {v: k for k, v in machine_mapping.items()}

df['Part'] = df['Part'].ffill()

grouped = df.groupby('Part')

# Split each set of data into separate DataFrames and store in a dict
parts_dict = {part: data.reset_index(drop=True) for part, data in grouped}

parts = []

for part, data in grouped:
    part_info = {
        "part_type": part,
        "quantity": int(data['Quantity'].iloc[0]),
        "operations": []
    }

    for _, row in data.iterrows():
        operation_type = int(row['Process'])
        machines = {machine_mapping[machine]: row[machine] for machine in machine_mapping if not pd.isna(row[machine])}
        part_info["operations"].append((operation_type, machines))

    parts.append(part_info)

num_mc = len(machine_mapping)

# num_pt = len(parts)
# num_gene = num_mc * num_pt  # number of genes in a chromosome
num_pt = 0
num_gene = 0
parts_single = []
id_index = 0
for part in parts:
    num_pt = num_pt + part["quantity"]
    num_gene = num_gene + part["quantity"] * len(part["operations"])
    # print(len(part["operations"]))
    # print(num_pt, num_gene)
    # print(part)

    for i in range(part["quantity"]):
        part_single_info = {
            "part_id": id_index,
            "operations": part["operations"],
            "num_mc_in_need": len(part["operations"])
        }
        parts_single.append(part_single_info)
        id_index = id_index + 1

# 获取命令行输入的字符串
input_string = input("Please input sequence_best (eg.[11, 2, 14,...]): ")
# 去掉字符串两端的方括号并分割成列表
string_list = input_string.strip('[]').split(',')
# 去掉每个元素的多余空格并转换为整数
sequence_best = list(map(lambda x: int(x.strip()), string_list))
# sequence_best = list(input('Please input sequence_best (eg.[11, 2, 14,...]): '))

'''--------plot gantt chart-------'''

j_keys = [j for j in range(num_pt)]  # j_keys 是列表，包含所有零件的索引
key_count = {key: 0 for key in j_keys}  # key_count 是字典，用于记录每个零件已经执行了多少操作。初始值为0
j_count = {key: 0 for key in j_keys}  # j_count 是字典，用于记录每个零件的总加工时间
m_keys = [j + 1 for j in range(num_mc)]  # m_keys 是列表，包含所有机器的索引。注意这里的键从1开始
m_count = {key: 0 for key in m_keys}  # m_count 是字典，用于记录每个机器的总加工时间
j_record = {}

for i in sequence_best:  # i是零件的索引号
    # 正在被分配的零件 temp_dict sample {1: 6.0, 2: 6.0, 3: 6.0, 4: 6.0}
    temp_dict = parts_single[i]['operations'][key_count[i]][1]
    # 可分配机器列表 sample [1, 2, 3, 4]
    gen_m_list = [key for key, value in temp_dict.items()]
    gen_m_list = sorted(gen_m_list, key=lambda x: (m_count[x] + temp_dict[x]))

    # gen_m 分配机器
    gen_m = gen_m_list[0]
    # gen_t 工序时间
    gen_t = int(temp_dict[gen_m])

    # j_count[i]用于记录零件i的总加工时间
    j_count[i] = j_count[i] + gen_t
    # m_count[gen_m]用于记录机器gen_m(被分配到的机器)的总加工时间
    m_count[gen_m] = m_count[gen_m] + gen_t

    # m_count[gen_m]和j_count[i]取其中较大的值，赋值
    if m_count[gen_m] < j_count[i]:
        m_count[gen_m] = j_count[i]
    elif m_count[gen_m] > j_count[i]:
        j_count[i] = m_count[gen_m]

    start_time = j_count[i] - gen_t
    end_time = j_count[i]

    j_record[(i, gen_m, key_count[i] + 1)] = [start_time, end_time]

    key_count[i] = key_count[i] + 1


Tbest = max(j_count.values())
print("Optimal Value: %f" % Tbest)

# 生成随机颜色
def generate_random_color():
    color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
    return color


colors = [generate_random_color() for _ in range(len(parts_single))]

start_date = datetime.datetime(2024, 6, 28)
df = []
for (job, machine, procedure), (start, end) in j_record.items():
    start_time = start_date + datetime.timedelta(minutes=start)
    end_time = start_date + datetime.timedelta(minutes=end)
    machine_id = id_to_machine.get(machine, '')
    df.append(dict(Task='Machine %s (%s)' % (machine, machine_id),
                   Start=start_time.strftime('%Y-%m-%d %H:%M:%S'),
                   Finish=end_time.strftime('%Y-%m-%d %H:%M:%S'),
                   Resource='Job %s' % (job + 1)))
df = sorted(df, key=lambda x: int(x['Task'].split()[1].replace('Machine', '')))
fig = ff.create_gantt(df,
                      index_col='Resource',
                      colors=colors,
                      show_colorbar=True,
                      group_tasks=True,
                      showgrid_x=True,
                      title='Flexible Job shop Schedule')
fig.show()


# pio.plot(fig, filename='GA_jobshop_realcase.html', auto_open=False)

# def write_solution_to_file(solution, filename):
#     # solution_t = solution.T
#     with open(filename, 'w') as file:
#         file.write('Nb of jobs ' + str(num_pt) + ' Nb of Machines ' + str(num_mc) + '\n')
#         file.write('\n')
#         file.write('Solutions' + '\n')
#         for row in solution:
#             file.write('\t\t'.join(map(str, row)) + '\n')
#
#
# def create_solution_matrix(sequence_best, num_pt, num_mc, ms):
#     # 初始化Solution矩阵，大小为num_mc * num_pt，初始值为-1表示未分配
#     solution = np.full((num_mc, num_pt), -1)
#
#     # 用于记录每个工件的工序次数
#     job_operation_count = {i: 0 for i in range(num_pt)}
#
#     for operation in sequence_best:
#         job = operation  # 工件编号
#         op_num = job_operation_count[job]  # 当前工件的工序编号
#         machine = ms[job][op_num] - 1  # 获取工件的工序对应的机器编号并转换为索引
#         job_operation_count[job] += 1  # 增加工件的工序次数
#
#         # 将该工件的当前工序放入Solution矩阵对应机器的下一空闲位置
#         for idx in range(num_pt):
#             if solution[machine][idx] == -1:
#                 solution[machine][idx] = job
#                 break
#     solution = solution + 1
#     return solution


def align_columns(data, separator='\t'):
    # Split the data into rows
    rows = [row.split(separator) for row in data.strip().split('\n')]
    # Find the maximum width of each column
    col_widths = [max(len(item) for item in col) for col in zip(*rows)]
    # Create a format string with appropriate width for each column
    format_string = separator.join([f'{{:{width}}}' for width in col_widths])
    # Format each row using the format string
    aligned_rows = [format_string.format(*row) for row in rows]
    # Join the formatted rows into a single string
    return '\n'.join(aligned_rows)


filename = instance_name + '_sol.txt'

with open(filename, 'w', encoding='utf-8') as file:
    file.write('Nb of jobs ' + str(num_pt) + ' Nb of Machines ' + str(num_mc) + '\n')
    file.write(str(Tbest) + '\n')
    file.write('Solutions' + '\n')
    # file.write('设备编号\t零件编号\t工序编号\t开始时间\t结束时间\n')

solution_column_names = ['设备编号', '零件编号', '工序编号', '开始时间', '结束时间']

solution_data = []
for (job, machine, procedure), (start, end) in j_record.items():
    # start_minutes =
    # end_minutes =
    solution_data.append([machine, job + 1, procedure, start, end])

solution_data = sorted(solution_data, key=lambda x: (x[0], x[2]))
# Sample
# solution_data = [
#     [1, 1, 1, '0', '10'],...
# ]

solution_pd = pd.DataFrame(solution_data, columns=solution_column_names)
# 将DataFrame转换为适当格式的字符串
solution_str = solution_pd.to_csv(sep='\t', index=False, header=True)
solution_str = solution_str.replace('\r', '')
solution_str = align_columns(solution_str)
with open(filename, 'a', encoding='utf-8') as file:
    file.write(solution_str)

# solution_best = create_solution_matrix(sequence_best, num_pt, num_mc, ms)
# write_solution_to_file(solution_best, instance_name + '_sol.txt')