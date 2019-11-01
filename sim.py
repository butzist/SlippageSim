from Bio2.SeqEvolution.Simulator import SimulationResult
import SlippageSim
from multiprocessing import Pool
import os
import sys
import Bio.AlignIO
import Bio.SeqIO
import Bio.Phylo
import Bio.Seq
import Bio.SeqRecord
import Bio2.SeqEvolution
import Bio2.Phylo.BDTree
import scipy.stats
import HMMERProfile
import pickle

def scale_tree(t,h):
	t.clade.branch_length = 0
	totlen = t.total_branch_length()
	adjustment = h/totlen
	for c in t.find_clades():
		c.branch_length *= adjustment

def do_sim(d,h,r):
   print(d)

   try:
      os.mkdir(d)
   except:
      print("%s already exists"%d)
      if os.path.isfile(d+"/repeats.full_info"):
         return
      else:
         print("recomputing "+d)

   tree = Bio2.Phylo.BDTree.simulate(5,1,6,1.0)
   scale_tree(tree,h)
   hmm = HMMERProfile.HMMERProfile("GALA_LRR.hmm",count=10000,binary="hmmemit")
   profiles = iter(hmm)
   length1 = scipy.stats.geom.rvs(1.0/100)
   length2 = scipy.stats.geom.rvs(1.0/100)
   rootseq1 = Bio2.SeqEvolution.random_sequence(length1,pi=Bio2.SeqEvolution.WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein)
   rootseq2 = Bio2.SeqEvolution.random_sequence(length2,pi=Bio2.SeqEvolution.WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein)
   repeats,homology = reduce(lambda x,y: (x[0]+y[0],x[1]+y[1]), [profiles.next() for j in xrange(6)], ("",[]))

   if False:
      nzix = 0
      while homology[nzix] == -1:
         nzix += 1

      repeats = repeats + repeats[nzix]
      homology = homology + [homology[nzix]]

   rootseq = rootseq1 + repeats + rootseq2
   rootseq = Bio2.SeqEvolution.SimulatorSequence(rootseq.tostring(),rootseq.alphabet)
   rootseq.tr_homology = length1*[-1] + homology + length2*[-1]
   rootseq.tr_parentmap = range(hmm.length)
   rootseq.tr_rootmap = range(hmm.length)

   events = [
               Bio2.SeqEvolution.SubstitutionFactory(rate=1.0,Q=Bio2.SeqEvolution.WAG.Q),
               #Bio2.SeqEvolution.RatesPerSiteSubst.GammaRatesSubstitutionFactory(Q=Bio2.SeqEvolution.WAG.Q,shape=1),
               SlippageSim.SlippageFactory(profiles,ins_rate=.5*r,del_rate=.5*r,length_distr=scipy.stats.geom(1.0/1.5,0)),
               Bio2.SeqEvolution.GeometricInsertionFactory(rate=0.005,alen=3.5),
               Bio2.SeqEvolution.GeometricDeletionFactory(rate=0.005,alen=3.5)
            ]

   result = Bio2.SeqEvolution.Simulator.SimulationResult(rootseq=rootseq,tree=tree,events=events)

   for s in result.leaf_seqs:
      s.description = ""

   #propagate tr_parentmap if present + merge
   SlippageSim.merge_tr_parentmap(tree.clade,result.info)

   trinfo = [Bio.SeqRecord.SeqRecord(id=n,seq=Bio.Seq.Seq(''.join('0' if s==-1 else '1' for s in i.tr_homology))) for n,i in result.info.iteritems()]

   Bio.SeqIO.write(result.leaf_seqs,open(d+"/input.fasta",'w'),"fasta")
   Bio.SeqIO.write(trinfo,open(d+"/input.trinfo",'w'),"fasta")
   Bio.Phylo.write(result.tree,open(d+"/input.tree",'w'),"newick",branch_length_only=True)
   Bio.AlignIO.write([result.msa],open(d+"/repeats.fasta",'w'),"fasta")
   Bio.AlignIO.write([result.ancestral_msa],open(d+"/repeats_anc.fasta",'w'),"fasta")
   SlippageSim.write_treks(open(d+"/repeats.treks",'w'),result)
   #SlippageSim.write_pairs(open(d+"/repeats.pairs",'w'),result)
   pickle.dump(result,open(d+"/repeats.full_info",'w'),pickle.HIGHEST_PROTOCOL)

if __name__=="__main__":
	r = float(sys.argv[3])
	h = float(sys.argv[2])
	d = sys.argv[1]
	do_sim(d,h,r)

