import Bio.Phylo
import scipy.stats
import math
import random

def _iQ(s,birth_rate,death_rate,n_taxa):
   if birth_rate == death_rate:
      return 1.0 / (birth_rate * (s**(-1.0/n_taxa) - 1.0))
   else:
      t1 = 1.0 - death_rate / birth_rate * s**(1.0/n_taxa)
      t2 = 1.0 - s**(1.0/n_taxa)
      return 1.0 / (birth_rate - death_rate) * math.log(t1/t2)


def _iF(s,birth_rate,death_rate,n_taxa,max_height):
   if birth_rate == death_rate:
      return s * max_height / (1.0 + birth_rate * max_height * (1.0 - s))
   else:
      e = math.exp((death_rate - birth_rate) * max_height)
      t1 = birth_rate - death_rate * e
      t2 = (1.0 - e) * s
      t3 = t1 - death_rate * t2
      t4 = t1 - birth_rate * t2
      return 1.0 / (birth_rate - death_rate) * math.log(t3/t4)


def simulate(birth_rate=1.0, death_rate=0.0, n_taxa=100, max_height=None):
   if birth_rate < death_rate:
      raise ValueError("death_rate has to be lower than birth_rate")

   r = scipy.stats.uniform.rvs(size=n_taxa-1)
   r.sort()

   if max_height is None:
      max_height = _iQ(scipy.stats.uniform.rvs(),birth_rate,death_rate,n_taxa)

   tree = Bio.Phylo.BaseTree.Tree()
   tree.clade.branch_length = max_height
   population = [tree.clade]

   for i in xrange(n_taxa-1):
      time = _iF(r[n_taxa-2-i],birth_rate,death_rate,n_taxa,max_height)
      clade = random.choice(population)
      clade.branch_length -= time
      clade.name = ""
      clade.split(n=2,branch_length=time)

      population.remove(clade)
      population += [clade[0],clade[1]]

   for i in xrange(len(population)):
      population[i].name = "S%.3d"%(i+1)
   
   return tree
