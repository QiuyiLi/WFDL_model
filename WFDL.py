import numpy as np

### arguments ###
population_size = 1000
total_generation = 10000
# tree_height_in_coalescent_unit = total_generation/population_size

species_tree=[(2000,[(5000,[(8000,[])])])] # caterpillar
# species_tree=[(2000,[(5000,[])]),(8000,[])] # balanced
# species_tree=[(0,[])] # naive species tree of only one leaf

dup_rate_base = 10**(-6)
dup_rate_scale = 1
dup_rate = dup_rate_scale * dup_rate_base
# dup_rate_per_coalescent_unit = dup_rate*population_size

n_rep = 1000000
initial_allele_count = 1
# only one mutant right after the duplication
first_order_dup_count = 0
higher_order_dup_count = 0

### functions ###
# joint probability that no gene is sampled from descendants
def joint_prob(descendant_list):
  prod=1
  for i in descendant_list:
    prod*=(1-i)
  return 1-prod

# when a duplicaion occurs at generation, return the child unilocus tree
def get_dup_subtree(tree, generation):
  dup_subtree=[]
  for item in tree:
    if item[0]>=generation:
      dup_subtree.append(item)
  return dup_subtree

# return 1. allele freq for each descendant species at currect locus 
#        2. new duplications occured at currect locus
def dup_sampling_locus(allele_count,start_generation,total_generation,tree):
  dup_list=[]
  res = []
  for subtree in tree:
    if start_generation == subtree[0]:
      # print('new branch at ',start_generation)
      new_res, new_dup = dup_sampling_locus(allele_count,start_generation,total_generation,subtree[1])
      if(new_res[-1]=='end'):
        # print('direct end')
        return ['end'],[]
      res+=new_res
      dup_list+=new_dup
  for generation in range(start_generation+1, total_generation):
      allele_count = np.random.binomial(population_size, allele_count/population_size)
      # print('allele_count at',generation,'is',allele_count)
      if allele_count == population_size:
        res.append('end')
        return res, dup_list
      elif allele_count > 0:
        if generation != total_generation-1:
          new_dup_number = np.random.binomial(allele_count,dup_rate)
          if new_dup_number > 0:
            for j in range(new_dup_number):
              dup_list.append((generation,get_dup_subtree(tree,generation)))
            # print('new dup',dup_list)
          for subtree in tree:
            if generation == subtree[0]:
              # print('new branch at ',generation)
              new_res, new_dup = dup_sampling_locus(allele_count,generation,total_generation,subtree[1])
              if(new_res[-1]=='end'):
                # print('direct end')
                return ['end'],[]
              res+=new_res
              dup_list+=new_dup
      else:
        # print('end with 0')
        break
        
  res.append(allele_count)
  # print('return res ',res)
  # print('return dup ',dup_list)
  return res, dup_list

# the main funtion, sample all duplications
def dup_sampling(total_generation, initial_allele_count, species_tree):
  isFirst = True # first order duplication
  dup_list = [(0,species_tree)]
  while dup_list:
    allele_count = initial_allele_count
    # print('start from generation',dup_list[0])
    dup_generation, dup_tree = dup_list.pop(0)
    descendant_allele_freq, new_dup = dup_sampling_locus(allele_count,dup_generation,total_generation,dup_tree)
    dup_list += new_dup
    # print('descendant_allele_freq',descendant_allele_freq)
    # print('new_dup',new_dup)
    if (descendant_allele_freq[-1]=='end'):
      # print('direct end')
      if isFirst:
        # print('first_order')
        return 0
      else:
        # print('higher_order')
        return 1
    descendant_allele_frac=[x/population_size for x in descendant_allele_freq]
    if np.random.binomial(1, joint_prob(descendant_allele_frac)) == 1:
      if isFirst:
        # print('joint first_order')
        return 0
      else:
        # print('joint higher_order')
        return 1
    else:
      # print('total end with 0')
      isFirst = False

### main ###
for k in range(n_rep):
  # if k%10000 ==0:
  #   print(k/n_rep)
  result = dup_sampling(total_generation,initial_allele_count,species_tree)
  if result == 0:
    # print('first1',first_order_dup_count)
    first_order_dup_count += 1
  elif result == 1:
    # print('later1',higher_order_dup_count)
    higher_order_dup_count += 1

print('first order dups: ', first_order_dup_count)
print('higher order dups: ', higher_order_dup_count)
print('total dup number: ', first_order_dup_count+higher_order_dup_count)
