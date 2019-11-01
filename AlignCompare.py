import Bio.AlignIO
import Bio.Align
import Bio.Alphabet
import Bio.Alphabet.IUPAC
import collections
import itertools

#class MsaChar:
#   __slots__ = ("id","char","index")
#   def __init__(self,id,char,index):
#      self.id = id
#      self.char = char
#      self.index = index
#      self.__setattr__ = self.__getattr__ = self.noassign
#
#   def noassign(self,*ignored):
#      raise AttributeError("object is immutable")
#
#   def __repr__(self):
#      return "MsaChar"+repr((self.id,self.char,self.index))
#
#   def __str__(self):
#      return repr(self)

MsaChar = collections.namedtuple("MsaChar",["id","char","index"])

def as_tuples(msa,gaps=False):
   for s1,s2 in itertools.combinations(msa,2):
      i1 = 0
      i2 = 0
      for i in xrange(msa.get_alignment_length()):
         c1 = s1[i]
         if c1 == s1.seq.alphabet.gap_char:
            j1 = -1
         else:
            j1 = i1
            i1 += 1

         c2 = s2[i]
         if c2 == s2.seq.alphabet.gap_char:
            j2 = -1
         else:
            j2 = i2
            i2 += 1

         if gaps or j1 != -1 and j2 != -1:
            yield frozenset([MsaChar(s1.id,c1,j1),MsaChar(s2.id,c2,j2)])

def as_columns(msa):
      ii = len(msa)*[0]
      m = len(msa)*[None]
      for i in xrange(msa.get_alignment_length()):
         for j,s in enumerate(msa):
            c = s[i]
            if c == s.seq.alphabet.gap_char:
               jj = -1
            else:
               jj = ii[j]
               ii[j] += 1
            m[j]=MsaChar(s.id,c,jj)

         yield frozenset(m)

def to_msa(something,alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)):
   if isinstance(something,str):
      something = open(something,'r')
   if isinstance(something,file):
      something = Bio.AlignIO.read(something,"fasta",alphabet=alphabet)

   if not isinstance(something,Bio.Align.MultipleSeqAlignment):
      raise ValueError("could not convert to Bio.Align.MultipleSeqAlignment")

   return something

def CS(reference,test,alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)):
   if not isinstance(reference,set):
      reference = set(as_columns(to_msa(reference,alphabet)))

   if not isinstance(test,set):
      test = set(as_columns(to_msa(test)))

   return len(test & reference) / float(len(reference,alphabet))

def SP(reference,test,alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)):
   if not isinstance(reference,set):
      reference = set(as_tuples(to_msa(reference,alphabet)))

   if not isinstance(test,set):
      test = set(as_tuples(to_msa(test,alphabet)))

   return len(test & reference) / float(len(reference))

def SPtd1(reference_pairs,test,alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)):
   if not isinstance(test,set):
      test = set(as_tuples(to_msa(test,alphabet)))

   overlap = 0
   for p in reference_pairs:
      p0 = set(p)
      p0 &= test
      if len(p0) > 0:
         overlap += 1
      #else:
      #   print("no match:", p)

   return overlap / float(len(reference_pairs))

def SPtd2(reference_pairs,test,alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)):
   if not isinstance(test,set):
      test = set(as_tuples(to_msa(test,alphabet)))

   overlap = 0
   for p in reference_pairs:
      p0 = set(p)
      p0 &= test
      if len(p0) > 0:
         overlap += 1
         #if len(p0) > 1:
         #   print("reference:", p)
         #   print("intersection:", p0)

   return overlap / float(len(test))

if __name__ == "__main__":
   alphabet = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
   msa = Bio.AlignIO.read(open("sim2.fasta",'r'),"fasta",alphabet=alphabet)
   print set(as_tuples(msa))
   print set(as_columns(msa))

