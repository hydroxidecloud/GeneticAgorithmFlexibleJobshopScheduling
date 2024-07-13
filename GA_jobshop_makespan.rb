require 'date'
require 'time'
require 'securerandom'
require 'ruby-progressbar'

# ================ initialization setting ======================

print "Please input the instance name (default value ta01): "
instance_name = gets.chomp
instance_name = 'ta01' if instance_name.empty?

lines = File.readlines("instance/#{instance_name}")

num_pt = lines[0].split[3].to_i
num_mc = lines[0].split[7].to_i
num_gene = num_mc * num_pt  # number of genes in a chromosome

pt_start = lines.index("Times\n") + 1
pt = lines[pt_start, num_pt].map { |line| line.split.map(&:to_i) }

ms_start = lines.index("Machines\n") + 1
ms = lines[ms_start, num_pt].map { |line| line.split.map(&:to_i) }

print "Please input the Size of Population (default value 100): "
population_size = gets.chomp.to_i
population_size = 100 if population_size.zero?

print "Please input the Size of Crossover Rate (default value 0.8): "
crossover_rate = gets.chomp.to_f
crossover_rate = 0.8 if crossover_rate.zero?

print "Please input the Size of Mutation Rate (default value 0.1): "
mutation_rate = gets.chomp.to_f
mutation_rate = 0.1 if mutation_rate.zero?

print "Please input the Mutation Selection Rate (default value 0.05): "
mutation_selection_rate = gets.chomp.to_f
mutation_selection_rate = 0.05 if mutation_selection_rate.zero?

print "Please input the Number of Iteration (default value 3000): "
num_iteration = gets.chomp.to_i
num_iteration = 3000 if num_iteration.zero?

num_mutation_jobs = (num_gene * mutation_selection_rate).round

if population_size <= 0 || num_iteration <= 0 || !(0..1).cover?(crossover_rate) || !(0..1).cover?(mutation_rate) || !(0..1).cover?(mutation_selection_rate)
  puts '超参数不符合要求，程序终止'
  exit
end

start_time = Time.now

# ==================== main code ===============================

# ----- generate initial population -----
Tbest = Float::INFINITY
sequence_best = nil
best_list, best_obj = [], []
population_list = []
makespan_record = []

population_size.times do
  nxm_random_num = (0...num_gene).to_a.shuffle
  population_list << nxm_random_num.map { |num| num % num_pt }
end

progressbar = ProgressBar.create(title: "Iterations", total: num_iteration, format: "%t: |%B| %p%% %e")

num_iteration.times do |n|
  Tbest_now = Float::INFINITY
  sequence_now = nil

  # -------- two point crossover --------
  parent_list = Marshal.load(Marshal.dump(population_list))
  offspring_list = Marshal.load(Marshal.dump(population_list))
  S = (0...population_size).to_a.shuffle

  (population_size / 2).times do |m|
    if crossover_rate >= rand
      parent_1 = parent_list[S[2 * m]]
      parent_2 = parent_list[S[2 * m + 1]]
      child_1 = parent_1.dup
      child_2 = parent_2.dup
      cutpoint = (0...num_gene).to_a.sample(2).sort

      child_1[cutpoint[0]...cutpoint[1]] = parent_2[cutpoint[0]...cutpoint[1]]
      child_2[cutpoint[0]...cutpoint[1]] = parent_1[cutpoint[0]...cutpoint[1]]
      offspring_list[S[2 * m]] = child_1
      offspring_list[S[2 * m + 1]] = child_2
    end
  end

  # ---------- repairment -------------
  population_size.times do |m|
    job_count = Hash.new { |hash, key| hash[key] = [0, 0] }
    larger, less = [], []

    offspring_list[m].each_with_index do |job, pos|
      job_count[job][0] += 1
      job_count[job][1] = pos if job_count[job][1] == 0

      if job_count[job][0] > num_mc
        larger << job
      elsif job_count[job][0] < num_mc
        less << job
      end
    end

    larger.each do |chg_job|
      while job_count[chg_job][0] > num_mc
        less.each do |l_job|
          if job_count[l_job][0] < num_mc
            offspring_list[m][job_count[chg_job][1]] = l_job
            job_count[chg_job][1] = offspring_list[m].index(chg_job)
            job_count[chg_job][0] -= 1
            job_count[l_job][0] += 1
          end
          break if job_count[chg_job][0] == num_mc
        end
      end
    end
  end

  # -------- mutation --------
  offspring_list.each_with_index do |offspring, m|
    if mutation_rate >= rand
      m_chg = (0...num_gene).to_a.sample(num_mutation_jobs)
      t_value_last = offspring[m_chg[0]]

      (num_mutation_jobs - 1).times do |i|
        offspring[m_chg[i]] = offspring[m_chg[i + 1]]
      end

      offspring[m_chg[num_mutation_jobs - 1]] = t_value_last
    end
  end

  # -------- fitness value (calculate makespan) -------------
  total_chromosome = parent_list + offspring_list
  chrom_fitness, chrom_fit = [], []
  total_fitness = 0

  total_chromosome.each do |chrom|
    j_count = Hash.new(0)
    m_count = Hash.new(0)
    key_count = Hash.new(0)

    chrom.each do |job|
      gen_t = pt[job][key_count[job]]
      gen_m = ms[job][key_count[job]]
      j_count[job] += gen_t
      m_count[gen_m] += gen_t

      m_count[gen_m] = [m_count[gen_m], j_count[job]].max
      j_count[job] = m_count[gen_m]

      key_count[job] += 1
    end

    makespan = j_count.values.max
    chrom_fitness << 1.0 / (2**(makespan-Tbest))
    chrom_fit << makespan
    total_fitness += chrom_fitness.last
  end

  # ---------- selection (roulette wheel approach) ----------
  pk = chrom_fitness.map { |fitness| fitness / total_fitness }
  qk = pk.each_with_object([]) { |p, arr| arr << (arr.last || 0) + p }

  selection_rand = Array.new(population_size) { rand }

  population_size.times do |i|
    selected = qk.index { |q| selection_rand[i] <= q } || (population_size * 2 - 1)
    population_list[i] = total_chromosome[selected]
  end

  # ---------- comparison ----------
  population_size.times do |i|
    if chrom_fit[i] < Tbest_now
      Tbest_now = chrom_fit[i]
      sequence_now = total_chromosome[i]
    end
  end

  if Tbest_now <= Tbest
    Tbest = Tbest_now
    sequence_best = sequence_now.dup
  end

  makespan_record << Tbest

  # 更新进度条
  progressbar.increment
end

# ---------- result ----------
puts "Optimal Sequence #{sequence_best}"
puts "Optimal Value: #{Tbest}"
puts "Program Processing Time: #{Time.now - start_time}"


# # Plotting (using a simpler plotting approach in Ruby)
# require 'gnuplot'

# Gnuplot.open do |gp|
#   Gnuplot::Plot.new(gp) do |plot|
#     plot.title 'Genetic Algorithm'
#     plot.xlabel 'Generation'
#     plot.ylabel 'Makespan'

#     x = (0...makespan_record.size).to_a
#     y1 = makespan_record
#     # y2 = iteration_makespan_mean_record
#     # y3 = iteration_makespan_min_record

#     plot.data << Gnuplot::DataSet.new([x, y1]) do |ds|
#       ds.with = 'lines'
#       ds.title = 'Best Makespan'
#     end

#     # plot.data << Gnuplot::DataSet.new([x, y2]) do |ds|
#     #   ds.with = 'lines'
#     #   ds.title = 'Mean Makespan'
#     # end

#     # plot.data << Gnuplot::DataSet.new([x, y3]) do |ds|
#     #   ds.with = 'lines'
#     #   ds.title = 'Min Makespan'
#     # end
#   end
# end
