import numpy as np
import time
import copy
import random
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import datetime
from tqdm import tqdm
import pandas as pd

'''==================== initial function ==============================='''
def load_and_process_data(file_path):
    global parts, machine_mapping, id_to_machine,num_mc
    df = pd.read_csv(file_path, header=1)
    columns = df.columns.tolist()
    columns[0] = 'Part'
    columns[1] = 'Quantity'
    columns[2] = 'Feature'
    columns[3] = 'Process'
    df.columns = columns

    machine_mapping = {col: idx + 1 for idx, col in enumerate(df.columns[4:])}
    num_mc=len(machine_mapping)
    id_to_machine = {v: k for k, v in machine_mapping.items()}
    df['Part'] = df['Part'].ffill()
    grouped = df.groupby('Part')

    parts = []
    for part, data in grouped:
        part_info = {
            "part_type": part,
            "quantity": int(data['Quantity'].iloc[0]),
            "operations": []
        }
        for _, row in data.iterrows():
            operation_type = int(row['Process'])
            machines = {machine_mapping[machine]: row[machine] for machine in machine_mapping if
                        not pd.isna(row[machine])}
            part_info["operations"].append((operation_type, machines))
        parts.append(part_info)


def initialize_parts():
    global parts_single, num_pt, num_gene
    parts_single = []
    num_pt = sum(part["quantity"] for part in parts)
    num_gene = sum(part["quantity"] * len(part["operations"]) for part in parts)

    id_index = 0
    for part in parts:
        for i in range(part["quantity"]):
            parts_single.append({
                "part_id": id_index,
                "operations": part["operations"]
            })
            id_index += 1


def repairment(population_list):
    target_counts = {part['part_id']: len(part['operations']) for part in parts_single}

    for m in range(len(population_list)):
        job_count = {}
        larger, less = [], []
        for i in range(num_pt):
            count = population_list[m].count(i)
            job_count[i] = [count, population_list[m].index(i) if count > 0 else 0]
            target_count = target_counts.get(i, 0)

            if count > target_count:
                larger.append(i)
            elif count < target_count:
                less.append(i)

        for chg_job in larger:
            while job_count[chg_job][0] > target_counts.get(chg_job, 0):
                for d in less:
                    if job_count[d][0] < target_counts.get(d, 0):
                        population_list[m][job_count[chg_job][1]] = d
                        job_count[chg_job][1] = population_list[m].index(chg_job)
                        job_count[chg_job][0] -= 1
                        job_count[d][0] += 1
                    if job_count[chg_job][0] == target_counts.get(chg_job, 0):
                        break
    return population_list


def initialize_population(population_size):
    population_list = [list(np.random.permutation(num_gene)) for _ in range(population_size)]
    for chrom in population_list:
        for j in range(num_gene):
            chrom[j] = chrom[j] % num_pt
    return population_list


'''==================== genetic process ==============================='''
def two_point_crossover(parent_list, crossover_rate):
    offspring_list = copy.deepcopy(parent_list)
    S = list(np.random.permutation(len(parent_list)))

    for m in range(len(parent_list) // 2):
        if np.random.rand() < crossover_rate:
            parent_1 = parent_list[S[2 * m]]
            parent_2 = parent_list[S[2 * m + 1]]
            child_1, child_2 = parent_1[:], parent_2[:]
            cutpoints = sorted(np.random.choice(num_gene, 2, replace=False))

            child_1[cutpoints[0]:cutpoints[1]] = parent_2[cutpoints[0]:cutpoints[1]]
            child_2[cutpoints[0]:cutpoints[1]] = parent_1[cutpoints[0]:cutpoints[1]]

            offspring_list[S[2 * m]] = child_1
            offspring_list[S[2 * m + 1]] = child_2

    return offspring_list


def mutate(offspring_list, mutation_rate, num_mutation_jobs):
    for chrom in offspring_list:
        if np.random.rand() < mutation_rate:
            m_chg = np.random.choice(num_gene, num_mutation_jobs, replace=False)
            t_value_last = chrom[m_chg[0]]
            for i in range(num_mutation_jobs - 1):
                chrom[m_chg[i]] = chrom[m_chg[i + 1]]
            chrom[m_chg[num_mutation_jobs - 1]] = t_value_last
    return offspring_list


def decode_chromosome(chromosome,return_value):
    j_keys = list(range(num_pt))
    key_count = {key: 0 for key in j_keys}
    j_count = {key: 0 for key in j_keys}
    m_keys = list(range(1, num_mc + 1))
    m_count = {key: 0 for key in m_keys}
    j_record={}

    for i in chromosome:
        temp_dict = parts_single[i]['operations'][key_count[i]][1]
        gen_m = min(temp_dict.keys(), key=lambda x: ( temp_dict[x],m_count[x]))
        gen_t = int(temp_dict[gen_m])

        j_count[i] += gen_t
        m_count[gen_m] += gen_t
        if m_count[gen_m] < j_count[i]:
            m_count[gen_m] = j_count[i]
        elif m_count[gen_m] > j_count[i]:
            j_count[i] = m_count[gen_m]

        key_count[i] += 1
        if return_value==0:
            start_time = j_count[i] - gen_t  # convert seconds to hours, minutes and seconds
            end_time = j_count[i]

            j_record[(i, gen_m,key_count[i])] = [start_time, end_time]
    if return_value:
        return max(j_count.values())
    else:
        return j_record
def calculate_fitness(total_chromosome):


    chrom_fitness, chrom_fit = [], []
    for m in range(len(total_chromosome)):
        makespan = decode_chromosome(total_chromosome[m],1)
        chrom_fit.append(makespan)    # chrom_fit=population makespanlist
        # chrom_fitness.append(1 / 2 ** (makespan - roulette_offset))
        chrom_fitness.append(1 / makespan)     # chrom_fitness=population ms function list

    return chrom_fitness, chrom_fit


def roulette_wheel_selection(chrom_fitness, total_chromosome, population_size):
    pk = [fitness / sum(chrom_fitness) for fitness in chrom_fitness]
    qk = np.cumsum(pk).tolist()
    selection_rand = [np.random.rand() for _ in range(population_size)]
    population_list = []

    for rand in selection_rand:
        for i, q in enumerate(qk):
            if rand <= q:
                population_list.append(copy.deepcopy(total_chromosome[i]))
                break

    return population_list

def tournament_selection(total_chromosome, chrom_fit, population_size, tournament_size=2):
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

    return selected_population,selected_crom_fit

def tabu_search(current_solution, best_fitness,num_iterations, tabu_list_length=5):
    current_solution = current_solution[:]
    best_solution = current_solution[:]
    tabu_list = []
    for _ in range(num_iterations):
        neighborhood = local_search(current_solution)
        best_move = None
        best_move_fitness = float('inf')

        for move in neighborhood:
            move_fitness = decode_chromosome(move,1)  # 可引入藐视准则
            if (move not in tabu_list) or (move_fitness<best_fitness):  # 符合藐视准则的和不在禁忌列表中的
                if move_fitness < best_move_fitness:    # 找到邻域中最优解
                    best_move = move
                    best_move_fitness = move_fitness

        if best_move is not None:
            current_solution = best_move
            tabu_list.append(best_move)
            if len(tabu_list) > tabu_list_length:
                tabu_list.pop(0)

        if decode_chromosome(current_solution,1) < best_fitness:
            best_solution = current_solution[:]
            best_fitness = decode_chromosome(best_solution,1)
    return best_solution,best_fitness

# 示例局部搜索函数：随机交换两个基因位置  （邻域结构应该考虑关键路径
def local_search(current_solution):
    neighborhood = []
    for _ in range(50):
        neighbor = current_solution[:]
        idx1, idx2 = random.sample(range(len(neighbor)), 2)
        neighbor[idx1], neighbor[idx2] = neighbor[idx2], neighbor[idx1]
        neighborhood.append(neighbor)
    return neighborhood


def compare_and_update_best(chrom_fit, population_list):
    Tave = sum(chrom_fit) / len(chrom_fit)
    Tbest_index = np.argmin(chrom_fit)
    sequence_now=[]
    Tbest_now=int(1e10)
    if chrom_fit[Tbest_index] < Tbest_now:
        Tbest_now = chrom_fit[Tbest_index]
        sequence_now = copy.deepcopy(population_list[Tbest_index])

    return Tbest_now, sequence_now, Tave


'''==================== result ==============================='''

def print_results(sequence_best, Tbest, start_time):
    print("Optimal Sequence", sequence_best)
    print("Optimal Value: %f" % Tbest)
    print('Program Processing Time: %s' % (time.time() - start_time))

def plot_results(makespan_record, iteration_makespan_mean_record, iteration_makespan_min_record):
    plt.plot([i for i in range(len(makespan_record))], makespan_record, 'b', label='Best Makespan')
    plt.plot([i for i in range(len(iteration_makespan_mean_record))], iteration_makespan_mean_record, 'r', label='Average Makespan')
    plt.plot([i for i in range(len(iteration_makespan_min_record))], iteration_makespan_min_record, 'g', label='Minimum Makespan')
    plt.ylabel('Makespan', fontsize=15)
    plt.xlabel('Generation', fontsize=15)
    plt.legend()
    plt.show()

def generate_random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def plot_gantt_chart(j_record):

    colors = [generate_random_color() for _ in range(len(parts_single))]
    start_date = datetime.datetime(2024, 6, 28)
    df = []
    for (job, machine,procedure), (start, end) in j_record.items():
        start_time = start_date + datetime.timedelta(minutes=start)
        end_time = start_date + datetime.timedelta(minutes=end)
        machine_id = id_to_machine.get(machine, '')
        df.append(dict(Task='Machine %s (%s)' % (machine, machine_id),
                       Start=start_time.strftime('%Y-%m-%d %H:%M:%S'),
                       Finish=end_time.strftime('%Y-%m-%d %H:%M:%S'),
                       Resource='Job %s' % (job + 1)))
    df = sorted(df, key=lambda x: int(x['Task'].split()[1].replace('Machine', '')))
    fig = ff.create_gantt(df, index_col='Resource', colors=colors, show_colorbar=True, group_tasks=True, showgrid_x=True, title='Flexible Job shop Schedule')
    fig.show()

def align_columns(data, separator='\t'):
    rows = [row.split(separator) for row in data.strip().split('\n')]
    col_widths = [max(len(item) for item in col) for col in zip(*rows)]
    format_string = separator.join([f'{{:{width}}}' for width in col_widths])
    aligned_rows = [format_string.format(*row) for row in rows]
    return '\n'.join(aligned_rows)

def save_results(filename, Tbest, j_record):
    solution_column_names = ['设备编号', '零件编号', '工序编号', '开始时间', '结束时间']
    solution_data = []
    for (job, machine, procedure), (start, end) in j_record.items():
        solution_data.append([machine, job + 1, procedure, start, end])
    solution_data = sorted(solution_data, key=lambda x: (x[0], x[2]))
    solution_pd = pd.DataFrame(solution_data, columns=solution_column_names)
    solution_str = solution_pd.to_csv(sep='\t', index=False, header=True)
    solution_str = solution_str.replace('\r', '')
    solution_str = align_columns(solution_str)

    with open(filename, 'w', encoding='utf-8') as file:
        file.write(f'Nb of jobs {num_pt} Nb of Machines {num_mc}\n')
        file.write(f'{Tbest}\n')
        file.write('Solutions\n')
        file.write(solution_str)

'''==================== run function ==============================='''
def genetic_algorithm(file_path, population_size, crossover_rate, mutation_rate, mutation_selection_rate, num_iteration):

    load_and_process_data(file_path)

    initialize_parts()
    population_list = initialize_population(population_size)
    population_list = repairment(population_list)

    num_mutation_jobs = round(num_gene * mutation_selection_rate)
    Tbest = float('inf')
    sequence_best = []
    makespan_record = []
    iteration_makespan_mean_record = []
    iteration_makespan_min_record = []

    for _ in tqdm(range(num_iteration), desc="Iterations"):
        parent_list = copy.deepcopy(population_list)
        offspring_list = two_point_crossover(parent_list, crossover_rate)
        offspring_list = repairment(offspring_list)
        offspring_list = mutate(offspring_list, mutation_rate, num_mutation_jobs)

        total_chromosome = parent_list + offspring_list

        chrom_fitness, chrom_fit = calculate_fitness(total_chromosome)

        # population_list = roulette_wheel_selection(chrom_fitness, total_chromosome, population_size)
        population_list,chrom_fit = tournament_selection(total_chromosome, chrom_fit, population_size)

        # for i in range(population_size):
        #     population_list[i],chrom_fit[i] = tabu_search(population_list[i],chrom_fit[i], 50)  # 10为禁忌搜索的迭代

        Tbest_now, sequence_now, Tave = compare_and_update_best(chrom_fit, population_list)

        makespan_record.append(Tbest)
        iteration_makespan_mean_record.append(Tave)
        iteration_makespan_min_record.append(Tbest_now)

        if Tbest_now <= Tbest:
            Tbest = Tbest_now
            sequence_best = copy.deepcopy(sequence_now)

    instance_name = 'realcase'
    filename = instance_name + '_sol.txt'
    print_results(sequence_best, Tbest, start_time)
    plot_results(makespan_record, iteration_makespan_mean_record, iteration_makespan_min_record)
    j_record = decode_chromosome(sequence_best, 0)
    plot_gantt_chart(j_record)
    save_results(filename, Tbest, j_record)




if __name__ == "__main__":
    file_path = 'instance/parts.csv'
    population_size = int(input('Please input the Size of Population (default value 100): ') or 100)
    crossover_rate = float(input('Please input the Crossover Rate (default value 0.8): ') or 0.8)
    mutation_rate = float(input('Please input the Mutation Rate (default value 0.1): ') or 0.1)
    mutation_selection_rate = float(input('Please input the Mutation Selection Rate (default value 0.05): ') or 0.05)
    num_iteration = int(input('Please input the Number of Iteration (default value 3000): ') or 3000)

    start_time=time.time()
    if (population_size <= 0 or num_iteration <= 0 or
            not (0 <= crossover_rate <= 1 and 0 <= mutation_rate <= 1 and 0 <= mutation_selection_rate <= 1)):
        print('Hyperparameters not valid, program terminated')
    else:
        genetic_algorithm(file_path, population_size, crossover_rate, mutation_rate, mutation_selection_rate,
                          num_iteration)


#  1500轮 3574  882.52t    换位3377  time:1370   3425.000000   tour:3549