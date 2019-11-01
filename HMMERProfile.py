import subprocess
import Bio.AlignIO
import Bio.SeqIO
import Bio.Alphabet.IUPAC

def _annotate_homology(seq):
   s = seq.seq
   ann = []
   for i in xrange(len(s)):
      if s[i].isupper():
         ann += [i]
      elif s[i].islower():
         ann += [-1]
   return ann

class HMMERProfile:
   def __init__(self,profile_filename,count=2000,binary="hmmemit"):
      p = subprocess.Popen([binary,"-a","-n",str(count),profile_filename],stdout=subprocess.PIPE,stdin=None,stderr=None,close_fds=True)
      while p.stdout.readline().strip() != "":
         pass
      alphabet = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
      alignment = Bio.AlignIO.read(p.stdout,"stockholm",alphabet=alphabet)
      self.alignment = alignment
      self.seqs = [s.seq.upper().ungap() for s in alignment]
      self.anns = [_annotate_homology(s) for s in alignment]
      self.length = len(alignment[0])

   def __iter__(self):
      return iter(zip(self.seqs,self.anns))
