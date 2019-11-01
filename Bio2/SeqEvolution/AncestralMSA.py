import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.Align
import numpy

__all__ = ["ancestral_msa"]

def _traverse_tree(t,info):
   if not t.is_bifurcating():
      raise ValueError("need binary tree")

   name = t.name
   seq = info[name].seq
   mapping = info[name].mapping
   labels = [name]

   if not t.is_terminal():
      (msa1,labels1,mapping1) = _traverse_tree(t[0],info)
      (msa2,labels2,mapping2) = _traverse_tree(t[1],info)

      msa_len = len(seq) + numpy.equal(mapping1,-1).sum() + numpy.equal(mapping2,-1).sum()
      labels = labels1 + labels + labels2
      msa = numpy.empty((len(labels),msa_len),dtype='|S1')
      msa.fill(Bio.Alphabet.Gapped(seq.alphabet).gap_char)
      new_mapping = msa_len*[-1]

      col = 0
      i1 = 0
      i2 = 0
      i = 0
      while col < msa_len:
         if i1 < len(mapping1) and mapping1[i1] < i:
            start = i1
            l = 0
            while i1 < len(mapping1) and mapping1[i1] < i:
               i1 += 1
            l = i1-start
            msa[:len(labels1),col:col+l] = msa1[:,start:i1]
            new_mapping[col:col+l] = l*[-1]
            col +=l
            continue

         if i2 < len(mapping2) and mapping2[i2] < i:
            start = i2
            l = 0
            while i2 < len(mapping2) and mapping2[i2] < i:
               i2 += 1
            l = i2-start
            msa[-len(labels2):,col:col+l] = msa2[:,start:i2]
            new_mapping[col:col+l] = l*[-1]
            col +=l
            continue

         if i >= len(mapping):
            assert(col == msa_len)
            break

         msa[len(labels1),col] = seq[i]
         new_mapping[col] = mapping[i]

         if i1 < len(mapping1) and mapping1[i1] == i:
            msa[:len(labels1),col] = msa1[:,i1]
            i1 += 1

         if i2 < len(mapping2) and mapping2[i2] == i:
            msa[-len(labels2):,col] = msa2[:,i2]
            i2 += 1

         i += 1
         col += 1

      mapping = new_mapping
   else:
      msa = numpy.empty((1,len(seq)),dtype='|S1')
      msa[0,:] = list(seq.tostring())

   assert(msa.shape[0] == len(labels))
   assert(msa.shape[1] == len(mapping))
   return(msa,labels,mapping)


def ancestral_msa(tree,info):
   (msa,labels,mapping) = _traverse_tree(tree.clade,info)
   seqs = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq.tostring(),Bio.Alphabet.Gapped(info[name].seq.alphabet)),
      id=name,name=name,description="simulated sequence") for name,seq in zip(labels,msa)]
   return Bio.Align.MultipleSeqAlignment(seqs)

