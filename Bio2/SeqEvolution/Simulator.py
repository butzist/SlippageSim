import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import Bio.Align
import Bio.Alphabet.IUPAC
import Bio2.Phylo.BDTree
import Bio2.SeqEvolution
import Bio2.SeqEvolution.AncestralMSA
import scipy.stats
import copy
import numpy

class SimulationResult:
   def __init__(self,rootseq=None,tree=None,events=None):
      if rootseq == None:
         length = scipy.stats.geom.rvs(1.0/300)
         rootseq = Bio2.SeqEvolution.random_sequence(length,pi=Bio2.SeqEvolution.WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein)

      if tree == None:
         tree = Bio2.Phylo.BDTree.simulate(3,1,10)

      if events == None:
         events = [Bio2.SeqEvolution.SubstitutionFactory(rate=1.0,Q=Bio2.SeqEvolution.WAG.Q),
               Bio2.SeqEvolution.GeometricInsertionFactory(rate=0.005,alen=3.5),
               Bio2.SeqEvolution.GeometricDeletionFactory(rate=0.005,alen=3.5)]

      _label_tree(tree.clade)
      info = _evolve(rootseq,tree.clade,events)
      info[tree.clade.name] = Bio2.SeqEvolution.mutate_sequence(rootseq,0,events,return_all=True)
      msa = Bio2.SeqEvolution.AncestralMSA.ancestral_msa(tree,info)

      self.rootseq = rootseq
      self.tree = tree
      self.leaves = map(lambda x:x.name,tree.get_terminals())

      self.ancestral_msa = msa
      self.seqs = [Bio.SeqRecord.SeqRecord(i.seq,id=n,name=n,description="simulated sequence") for n,i in info.iteritems()]
      self.leaf_seqs = filter(lambda x:x.id in self.leaves,self.seqs)
      self.info = info

      seqs = filter(lambda x:x.id in self.leaves,self.ancestral_msa)
      self.msa = _msa_remove_gap_cols(Bio.Align.MultipleSeqAlignment(seqs))


def _msa_remove_gap_cols(msa):
   array,records = _msa_to_array(msa)
   array2 = array.T
   rows = 0
   gaps = ''.join(r.seq.alphabet.gap_char for r in records)
   for i in xrange(len(array2)):
      if array2[i].tostring() != gaps:
         rows += 1

   array3 = numpy.empty((rows,len(records)),dtype='|S1')

   rows = 0
   for i in xrange(len(array2)):
      if array2[i].tostring() != gaps:
         array3[rows,:] = array2[i,:]
         rows += 1
   assert(array3.shape[0] == rows)

   msa = _array_to_msa(array3.T,records)
   return msa

def _msa_to_array(msa):
   records = [copy.deepcopy(r) for r in msa]
   array = numpy.empty((len(records),len(records[0])),dtype='|S1')
   for i in xrange(len(records)):
      array[i,:] = list(records[i].seq.tostring())

   return array,records

def _array_to_msa(array,records):
   for i in xrange(len(records)):
      r = records[i]
      r.seq = Bio.Seq.Seq(array[i].tostring(),alphabet=r.seq.alphabet)

   return Bio.Align.MultipleSeqAlignment(records)

def _label_tree(tree):
   if not tree.is_terminal():
      leaves = []
      for child in tree:
         leaves += _label_tree(child)
      leaves.sort()
      tree.name = '('+','.join(leaves)+')'
   else:
      leaves = [tree.name]

   return leaves


def _evolve(seq,tree,events):
   info = dict()

   for child in tree:
      name = child.name
      mutated = Bio2.SeqEvolution.mutate_sequence(seq,child.branch_length,events,return_all=True)

      if not child.is_terminal():
         res = _evolve(mutated,child,events)
         info.update(res)

      info[name] = mutated

   return info

