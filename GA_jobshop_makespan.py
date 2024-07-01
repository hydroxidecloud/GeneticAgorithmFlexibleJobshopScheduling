import copy
import datetime
import random
import time

import matplotlib.pyplot as plt
import numpy as np
import plotly.figure_factory as ff
from tqdm import tqdm

''' ================= initialization setting ======================'''

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

population_size = int(input('Please input the Size of Population (default value 100): ') or 100)
crossover_rate = float(input('Please input the Size of Crossover Rate (default value 0.8): ') or 0.8)
mutation_rate = float(input('Please input the Size of Mutation Rate (default value 0.1): ') or 0.1)
mutation_selection_rate = float(input('Please input the Mutation Selection Rate (default value 0.05): ') or 0.05)
# roulette_offset = int(input('Please input the Roulette Offset (default value 1500): ') or 1500)
# roulette_expo_base = float(input('Please input the Roulette Exponential f. Base (default value 2): ') or 2)
tournament_size = int(input('Please input the Tournament Size (default value 2): ') or 2)
num_iteration = int(input('Please input the Number of Iteration (default value 3000): ') or 3000)
num_mutation_jobs = round(num_gene * mutation_selection_rate)

if (population_size <= 0 or
        num_iteration <= 0 or
        not (0 <= crossover_rate <= 1 and
             0 <= mutation_rate <= 1 and
             0 <= mutation_selection_rate <= 1)):
    print('超参数不符合要求，程序终止')
    exit()

start_time = time.time()

'''==================== main code ==============================='''
'''----- generate initial population -----'''
Tbest = int(1e10)
best_list, best_obj = [], []
population_list = []
makespan_record = []
iteration_makespan_mean_record = []
iteration_makespan_min_record = []
sequence_now = []
sequence_best = []
for i in range(population_size):
    nxm_random_num = list(np.random.permutation(num_gene))  # generate a random permutation of 0 to num_pt*num_mc-1
    population_list.append(nxm_random_num)  # add to the population_list
    for j in range(num_gene):
        population_list[i][j] = population_list[i][
                                    j] % num_pt  # convert to job number format, every job appears m times

# Add tqdm to the iteration loop
for n in tqdm(range(num_iteration), desc="Iterations"):
    Tbest_now = int(1e10)
    iteration_makespan_list = []

    '''-------- two point crossover --------'''
    parent_list = copy.deepcopy(population_list)
    offspring_list = copy.deepcopy(population_list)
    S = list(np.random.permutation(
        population_size))  # generate a random sequence to select the parent chromosome to crossover

    for m in range(int(population_size / 2)):
        crossover_prob = np.random.rand()
        if crossover_rate >= crossover_prob:
            parent_1 = population_list[S[2 * m]][:]
            parent_2 = population_list[S[2 * m + 1]][:]
            child_1 = parent_1[:]
            child_2 = parent_2[:]
            cutpoint = list(np.random.choice(num_gene, 2, replace=False))
            cutpoint.sort()

            child_1[cutpoint[0]:cutpoint[1]] = parent_2[cutpoint[0]:cutpoint[1]]
            child_2[cutpoint[0]:cutpoint[1]] = parent_1[cutpoint[0]:cutpoint[1]]
            offspring_list[S[2 * m]] = child_1[:]
            offspring_list[S[2 * m + 1]] = child_2[:]

    '''----------repairment-------------'''
    for m in range(population_size):
        job_count = {}
        larger, less = [], []
        for i in range(num_pt):
            if i in offspring_list[m]:
                count = offspring_list[m].count(i)
                pos = offspring_list[m].index(i)
                job_count[i] = [count, pos]  # store the above two values to the job_count dictionary
            else:
                count = 0
                job_count[i] = [count, 0]
            if count > num_mc:
                larger.append(i)
            elif count < num_mc:
                less.append(i)

        for k in range(len(larger)):
            chg_job = larger[k]
            while job_count[chg_job][0] > num_mc:
                for d in range(len(less)):
                    if job_count[less[d]][0] < num_mc:
                        offspring_list[m][job_count[chg_job][1]] = less[d]
                        job_count[chg_job][1] = offspring_list[m].index(chg_job)
                        job_count[chg_job][0] = job_count[chg_job][0] - 1
                        job_count[less[d]][0] = job_count[less[d]][0] + 1
                    if job_count[chg_job][0] == num_mc:
                        break

    '''--------mutation--------'''
    for m in range(len(offspring_list)):
        mutation_prob = np.random.rand()
        if mutation_rate >= mutation_prob:
            m_chg = list(
                np.random.choice(num_gene, num_mutation_jobs, replace=False))  # chooses the position to mutation
            t_value_last = offspring_list[m][m_chg[0]]  # save the value which is on the first mutation position
            for i in range(num_mutation_jobs - 1):
                offspring_list[m][m_chg[i]] = offspring_list[m][m_chg[i + 1]]  # displacement

            offspring_list[m][m_chg[num_mutation_jobs - 1]] = t_value_last

    '''--------fitness value(calculate makespan)-------------'''
    total_chromosome = copy.deepcopy(parent_list) + copy.deepcopy(
        offspring_list)  # parent and offspring chromosomes combination
    chrom_fitness, chrom_fit = [], []
    total_fitness = 0
    for m in range(population_size * 2):
        j_keys = [j for j in range(num_pt)]
        key_count = {key: 0 for key in j_keys}
        j_count = {key: 0 for key in j_keys}
        m_keys = [j + 1 for j in range(num_mc)]
        m_count = {key: 0 for key in m_keys}

        for i in total_chromosome[m]:
            gen_t = int(pt[i][key_count[i]])
            gen_m = int(ms[i][key_count[i]])
            j_count[i] = j_count[i] + gen_t
            m_count[gen_m] = m_count[gen_m] + gen_t

            if m_count[gen_m] < j_count[i]:
                m_count[gen_m] = j_count[i]
            elif m_count[gen_m] > j_count[i]:
                j_count[i] = m_count[gen_m]

            key_count[i] = key_count[i] + 1

        makespan = max(j_count.values())
        chrom_fit.append(makespan)
        # chrom_fitness.append(1 / 2**(makespan-roulette_offset))
        chrom_fitness.append(1 / makespan)
        total_fitness = total_fitness + chrom_fitness[m]
        iteration_makespan_list.append(makespan)

    pass

    '''----------selection(roulette wheel approach)----------'''
    # pk, qk = [], []
    #
    # for i in range(population_size * 2):
    #     pk.append(chrom_fitness[i] / total_fitness)
    # for i in range(population_size * 2):
    #     cumulative = 0
    #     for j in range(0, i + 1):
    #         cumulative = cumulative + pk[j]
    #     qk.append(cumulative)
    #
    # selection_rand = [np.random.rand() for i in range(population_size)]
    #
    # for i in range(population_size):
    #     if selection_rand[i] <= qk[0]:
    #         population_list[i] = copy.deepcopy(total_chromosome[0])
    #     else:
    #         for j in range(0, population_size * 2 - 1):
    #             if qk[j] < selection_rand[i] <= qk[j + 1]:
    #                 population_list[i] = copy.deepcopy(total_chromosome[j + 1])
    #                 break
    '''----------selection(tournament approach)----------'''
    selected_population = []
    selected_crom_fit = []
    for _ in range(population_size):
        # 随机选择 tournament_size 个体
        tournament_indices = np.random.choice(len(total_chromosome), tournament_size, replace=False)  # 一个个体可以被多次选择
        # 从中选择适应度最高的个体
        best_index = tournament_indices[0]
        for index in tournament_indices[1:]:
            if chrom_fit[index] < chrom_fit[best_index]:
                best_index = index
        # 将选择的个体添加到新的种群中并更新chrom_fit列表
        selected_population.append(copy.deepcopy(total_chromosome[best_index]))
        selected_crom_fit.append(chrom_fit[best_index])
        population_list = selected_population

    '''----------comparison----------'''
    for i in range(population_size * 2):
        if chrom_fit[i] < Tbest_now:
            Tbest_now = chrom_fit[i]
            sequence_now = copy.deepcopy(total_chromosome[i])
    if Tbest_now <= Tbest:
        Tbest = Tbest_now
        sequence_best = copy.deepcopy(sequence_now)

    makespan_record.append(Tbest)
    iteration_makespan_mean_record.append(np.mean(iteration_makespan_list))
    iteration_makespan_min_record.append(np.min(iteration_makespan_list))

'''----------result----------'''
print("Optimal Sequence", sequence_best)
print("Optimal Value: %f" % Tbest)
print('Program Processing Time: %s' % (time.time() - start_time))

# %matplotlib inline
plt.plot([i for i in range(len(makespan_record))], makespan_record, 'b')
plt.plot([i for i in range(len(iteration_makespan_mean_record))], iteration_makespan_mean_record, 'r')
plt.plot([i for i in range(len(iteration_makespan_min_record))], iteration_makespan_min_record, 'g')
plt.ylabel('makespan', fontsize=15)
plt.xlabel('generation', fontsize=15)
plt.show()

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


# pio.write_html(fig, file='GA_job_shop_scheduling.html', auto_open=True)

def write_solution_to_file(solution, filename):
    # solution_t = solution.T
    with open(filename, 'w', encoding='utf-8') as file:
        file.write('Nb of jobs ' + str(num_pt) + ' Nb of Machines ' + str(num_mc) + '\n')
        file.write(str(Tbest)+'\n')
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
