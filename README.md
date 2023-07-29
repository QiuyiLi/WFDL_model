# Wright-Firsher with Duplication and Loss

This python3 code simulates the allele frequency according to the WFDL (Wright-Firsher with Duplication and Loss) model. Note that this code does not produce any gene trees, which requires too much memory/runing time to track the genealogy of every individual in the whole population. It only simulates the allele frequency and thereafter the fraction of higher-order duplications, which is ignored in the [MLMSC-II model](https://github.com/QiuyiLi/MLMSC-II). For more details, please refer to [The Effect of Copy Number Hemiplasy on Gene].

##  Usage
The users will need to open WFDL.py and manually change the arguments listed on top of the code:

**population_size** stands for the effective population size in the Wright-Firsher model, e.g.,
```
population_size = 1000
```

**total_generation** stands for the species tree height in number of generations, e.g.,
```
total_generation = 10000
```

**species_tree** new format, detail explanation, e.g.,
```
species_tree=[(2000,[(5000,[(8000,[])])])] # caterpillar
species_tree=[(2000,[(5000,[])]),(8000,[])] # balanced
species_tree=[(0,[])] # naive species tree of only one leaf
```

**dup_rate_base** stands for the unit duplication rate, and **dup_rate_scale** is a scaling constant, e.g.,
```
dup_rate_base = 10**(-6)
dup_rate_scale = 1
dup_rate = dup_rate_scale * dup_rate_base
```

**n_rep** sets the number of repetitions, e.g.,
```
n_rep = 1000000
```
By WF model, the survival probability of a duplication is O(1/N), therefore the number of sampled duplications is supposed to be around n_rep/population_size.

To execute the code:
```
python3 WFDL.py
```

Output:
```
first order dups:  1113
higher order dups:  1
total dup number:  1114
```
