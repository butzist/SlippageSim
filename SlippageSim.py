import Bio.Alphabet.IUPAC
import Bio.Alphabet
import Bio.AlignIO
import Bio2.Phylo.BDTree
import Bio2.SeqEvolution
import Bio2.SeqEvolution.Simulator
import Bio2.SeqEvolution.RatesPerSiteSubst
import numpy
import itertools
import HMMERProfile
import random
import scipy.stats
import sys
import pickle
import bisect
import AlignCompare


#  ____  _ _                                ____             __ _ _
# / ___|| (_)_ __  _ __   __ _  __ _  ___  |  _ \ _ __ ___  / _(_) | ___
# \___ \| | | '_ \| '_ \ / _` |/ _` |/ _ \ | |_) | '__/ _ \| |_| | |/ _ \
#  ___) | | | |_) | |_) | (_| | (_| |  __/ |  __/| | | (_) |  _| | |  __/
# |____/|_|_| .__/| .__/ \__,_|\__, |\___| |_|   |_|  \___/|_| |_|_|\___|
#           |_|   |_|          |___/

class TRInsertionEvent(Bio2.SeqEvolution.InsertionEvent):
   def __init__(self,distance,start,length,sequence,homology):
      Bio2.SeqEvolution.InsertionEvent.__init__(self,distance,start,length,sequence)
      self._homology = homology

   def get_homology(self):
      return self._homology

   def apply(self,sequence):
      Bio2.SeqEvolution.InsertionEvent.apply(self,sequence)

      sequence.tr_homology = sequence.tr_homology[:self.get_start()] + self._homology + sequence.tr_homology[self.get_start():]
      assert(len(sequence.tr_homology) == len(sequence.seq))

class TRDeletionEvent(Bio2.SeqEvolution.DeletionEvent):
   pass

def _homiter(seq,start_index=0):
   i = start_index
   for i in xrange(start_index,len(seq.tr_homology)):
      h = seq.tr_homology[i]
      if h != -1:
         yield h,i

class RepeatUnit:
   def __init__(self,sequence,start_index,length):
      self.sequence = sequence
      self.start = start_index
      self.length = length

def _split_units(sequence):
   units = []
   it = _homiter(sequence)
   try:
      h,i = it.next()
      prev = h
      start = i
   except StopIteration:
      return []

   for h,i in it:
      if h <= prev:
         units += [RepeatUnit(sequence,start,i-start)]
         start = i
      prev = h

   units += [RepeatUnit(sequence,start,i+1-start)]
   return units

class SlippageFactory(Bio2.SeqEvolution.EventFactory):
   def __init__(self,profiles,ins_rate,del_rate,length_distr):
      self._profiles = profiles
      self._insrate = ins_rate
      self._delrate = del_rate
      self._len = length_distr
      self._units = []

   def get_rate(self):
      #no deletions in last unit
      return self._insrate * len(self._units) + self._delrate * (len(self._units)-1)

   def get_event(self,sequence,distance):
      if random.random()*self.get_rate() > self._insrate * len(self._units):
         #deletion
         start_unit_index = int(random.random()*(len(self._units)-1))
         start_unit = self._units[start_unit_index]
         start = start_unit.start + int(random.random()*start_unit.length)
         end_unit_index = min(start_unit_index+int(self._len.rvs()),len(self._units)-1)
         end_unit = self._units[end_unit_index]
         end = end_unit.start
         while end < end_unit.start+end_unit.length and (sequence.tr_homology[end] == -1 or sequence.tr_homology[end] < sequence.tr_homology[start]):
            end += 1
         return TRDeletionEvent(distance,start,end-start)
      else:
         #insertion
         new_seq,new_hom = self._profiles.next()
         for i in xrange(self._len.rvs()-1):
            more_seq,more_hom = self._profiles.next()
            new_seq += more_seq
            new_hom += more_hom

         imap = _invert_map(sequence.tr_rootmap)
         if -1 in imap:
            del imap[-1]
         imap = _Identdict(imap)
         new_hom = [imap[i] for i in new_hom]

         start_unit_index = int(random.random()*len(self._units))
         start_unit = self._units[start_unit_index]
         start = start_unit.start + int(random.random()*start_unit.length)

         #rotate inserted profile
         pos = start
         while pos < len(sequence.seq) and sequence.tr_homology[pos] == -1:
            pos += 1

         if pos < len(sequence.seq):
            phase = sequence.tr_homology[pos]
         else:
            phase = 1

         pos = len(new_seq)-1
         while pos >= 0 and (new_hom[pos] == -1 or new_hom[pos] >= phase):
            pos -= 1

         #insertions left or right?
         if random.random() > .5:
            while pos < len(new_seq) and new_hom[pos] == -1:
               pos += 1

         pos += 1
         new_seq = new_seq[pos:] + new_seq[:pos]
         new_hom = new_hom[pos:] + new_hom[:pos]

         return TRInsertionEvent(distance,start,len(new_seq),new_seq,new_hom)

   @staticmethod
   def _handle_gap(sequence,start):
      assert(start >= len(sequence.tr_homology) or sequence.tr_homology[start] == -1)
      if start >= 0:
         index = sequence.tr_homology[start-1]+1
      else:
         index = 0

      end = start
      while end < len(sequence.tr_homology) and sequence.tr_homology[end] == -1:
         end += 1
      length = end-start

      sequence.tr_parentmap = sequence.tr_parentmap[:index] + length*[-1] + sequence.tr_parentmap[index:]
      sequence.tr_rootmap = sequence.tr_rootmap[:index] + length*[-1] + sequence.tr_rootmap[index:]
      #tr_start = min(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
      #tr_end = max(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)

      for i in xrange(len(sequence.tr_homology)):
         if sequence.tr_homology[i] >= index:
            sequence.tr_homology[i] += length

      sequence.tr_homology[start:end] = range(index,index+length)

   def update(self,sequence,event):
      if isinstance(event,Bio2.SeqEvolution.IndelEvent) or event is None:
         if isinstance(event,TRInsertionEvent) or event is None:
            assert(len(sequence.tr_homology) == len(sequence.seq))

            if event is None:
               sequence.tr_homology = sequence.tr_homology[:]
               sequence.tr_parentmap = range(len(sequence.tr_parentmap))
               sequence.tr_rootmap = sequence.tr_rootmap[:]

            start = min(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
            end = max(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
            for i in xrange(start,end):
               if sequence.tr_homology[i] == -1:
                  self._handle_gap(sequence,i)

         else:
            #update tr homology
            if event.get_length() > 0:
               sequence.tr_homology = sequence.tr_homology[:event.get_start()] + \
                  event.get_length()*[-1] + sequence.tr_homology[event.get_start():]
               #new indices if inside tr
               if not event.get_start() == 0 and not sequence.tr_homology[event.get_start()-1] == -1 or not event.get_start()+event.get_length() == len(sequence.tr_homology) and not sequence.tr_homology[event.get_start()+event.get_length()] == -1:
                  self._handle_gap(sequence,event.get_start())
            else:
               sequence.tr_homology = sequence.tr_homology[:event.get_start()] + \
                  sequence.tr_homology[event.get_start()-event.get_length():]

         start = min(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
         end = max(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
         assert(min(sequence.tr_homology[start:end]) > -1)

         #count units
         self._units = _split_units(sequence)

#  ____  _ _
# / ___|| (_)_ __  _ __   __ _  __ _  ___
# \___ \| | | '_ \| '_ \ / _` |/ _` |/ _ \
#  ___) | | | |_) | |_) | (_| | (_| |  __/
# |____/|_|_| .__/| .__/ \__,_|\__, |\___|
#           |_|   |_|          |___/
#  ____              _ _           _   _
# |  _ \ _   _ _ __ | (_) ___ __ _| |_(_) ___  _ __
# | | | | | | | '_ \| | |/ __/ _` | __| |/ _ \| '_ \
# | |_| | |_| | |_) | | | (_| (_| | |_| | (_) | | | |
# |____/ \__,_| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
#             |_|

class TRDupInsertionEvent(TRInsertionEvent):
   def __init__(self,distance,start,length,sequence):
      seq = sequence.seq[start:start+length]
      hom = sequence.tr_homology[start:start+length]
      assert(not -1 in hom)
      TRInsertionEvent.__init__(self,distance,start+length,length,seq,hom)

   def apply(self,sequence):
      assert(len(sequence.tr_homology) == len(sequence.seq))
      assert(len(sequence.n21mapping) == len(sequence.seq))

      TRInsertionEvent.apply(self,sequence)
      sequence.n21mapping = sequence.n21mapping[:self.get_start()] + \
         sequence.n21mapping[self.get_start()-self.get_length():]

      if hasattr(sequence,"rates"):
         homrates = dict()
         for i in xrange(len(sequence.rates)):
            homrates[sequence.tr_homology[i]] = sequence.rates[i]
         del homrates[-1]

         old_rates = sequence.rates
         sequence.rates = numpy.empty((len(sequence.seq),))
         sequence.rates[:self.get_start()] = old_rates[:self.get_start()]
         sequence.rates[self.get_start()+self.get_length():] = old_rates[self.get_start():]

         sequence.rates[self.get_start():self.get_start()+self.get_length()] = map(lambda x:homrates[x], self._homology)
         assert(len(sequence.rates) == len(sequence.seq))

      assert(len(sequence.tr_homology) == len(sequence.seq))
      assert(len(sequence.n21mapping) == len(sequence.seq))

class SlippageDupFactory(Bio2.SeqEvolution.EventFactory):
   def __init__(self,ins_rate,del_rate,length_distr):
      self._insrate = ins_rate
      self._delrate = del_rate
      self._lendistr = length_distr
      self._pairs = []
      self._probs = []
      self._n_units = -1

   def _init_pairs(self,sequence):
      homs = dict()
      for i,h in enumerate(sequence.tr_homology):
         try:
            homs[h] += [i]
         except:
            homs[h] = [i]

      del homs[-1]

      self._pairs = []
      self._probs = []
      for p in homs.itervalues():
         for i1 in xrange(len(p)):
            for i2 in xrange(i1+1,len(p)):
               self._pairs += [(p[i1],p[i2])]
               self._probs += [self._lendistr.pmf(abs(p[i2]-p[i1]))]

      self._probs = numpy.array(self._probs)
      if len(self._probs > 0):
         for i in xrange(1,len(self._probs)):
            self._probs[i] += self._probs[i-1]
         self._probs /= float(self._probs[-1])
         assert(self._probs[-1] <= 1.0)

   def get_rate(self):
      if len(self._pairs) > 1:
         #no deletions in last unit
         #return self._insrate + self._delrate
         return self._insrate * self._n_units + self._delrate * (self._n_units-1)
      else:
         return 0.0

   def get_event(self,sequence,distance):
      assert(len(self._pairs) > 0)
      # select two homologous positions
      p = bisect.bisect_left(self._probs,random.random())
      start,end = self._pairs[p]
      length = end - start

      if random.random()*self.get_rate() < self._insrate or len(self._pairs) == 1:
         #insertion
         # copy seq,hom and rates
         return TRDupInsertionEvent(distance,start,length,sequence)
      else:
         # delete all in between
         return TRDeletionEvent(distance,start,end-start)

   def update(self,sequence,event):
      if isinstance(event,Bio2.SeqEvolution.IndelEvent) or event is None:
         if isinstance(event,TRInsertionEvent) or event is None:
            assert(len(sequence.tr_homology) == len(sequence.seq))
            if event is None:
               sequence.tr_homology = sequence.tr_homology[:]
               sequence.n21mapping = sequence.mapping[:]
               start = min(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
               end = max(i for i in xrange(len(sequence.tr_homology)) if sequence.tr_homology[i] != -1)
               for i in xrange(start,end):
                  if sequence.tr_homology[i] == -1:
                     SlippageFactory._handle_gap(sequence,i)
            assert(len(sequence.n21mapping) == len(sequence.seq))

         else:
            #update tr homology and n21mapping
            if event.get_length() > 0:
               sequence.tr_homology = sequence.tr_homology[:event.get_start()] + \
                  event.get_length()*[-1] + sequence.tr_homology[event.get_start():]
               sequence.n21mapping = sequence.n21mapping[:event.get_start()] + \
                  event.get_length()*[-1] + sequence.n21mapping[event.get_start():]
               #new indices if inside tr
               if not event.get_start() == 0 and not sequence.tr_homology[event.get_start()-1] == -1 or not event.get_start()+event.get_length() == len(sequence.tr_homology) and not sequence.tr_homology[event.get_start()+event.get_length()] == -1:
                  SlippageFactory._handle_gap(sequence,event.get_start())
            else:
               sequence.tr_homology = sequence.tr_homology[:event.get_start()] + \
                  sequence.tr_homology[event.get_start()-event.get_length():]
               sequence.n21mapping = sequence.n21mapping[:event.get_start()] + \
                  sequence.n21mapping[event.get_start()-event.get_length():]

               #find alternative homologs if "original" copy is deleted (greedy)
               homs = set(sequence.mapping)
               mmax = -1
               assert(len(sequence.mapping)==len(sequence.n21mapping))
               for i in xrange(len(sequence.n21mapping)):
                  c = sequence.n21mapping[i]
                  if c > mmax:
                     mmax = c
                     if not c in homs:
                        homs.add(c)
                        assert(sequence.mapping[i] == -1)
                        sequence.mapping[i] = c

         self._init_pairs(sequence)
         self._n_units = len(_split_units(sequence))

def write_pairs(out,result):
   pairs = get_pairs(result)
   pickle.dump(pairs,out,pickle.HIGHEST_PROTOCOL)

def get_pairs(result):
   pairs = list()
   for site in xrange(len(result.ancestral_msa[0])):
      pairs += _get_pairs_site(result,result.tree.clade,site,-1,None,None)[1]

   #move single pairs to separate set?
   return pairs

def _get_pairs_site(result,tree,site,ancpos,ancchar,ancname,onlyleaves=True):
   homologs = dict()
   pairs = list()
   msa = result.ancestral_msa
   seq = itertools.ifilter(lambda s:s.id == tree.name, msa).next().seq
   tr_hom = result.info[tree.name].tr_homology
   mapping = result.info[tree.name].n21mapping
   try:
      gap_char = seq.alphabet.gap_char
   except:
      gap_char = '-'

   new_homs = set()
   new_pairs = [] # pairs betw us and parent
   homologs = dict() # homologs between all children through current pos
   p = 0
   for i in xrange(len(seq)):
      if seq[i] != gap_char:
         if ancpos != -1 and mapping[p] == ancpos:
            #pairs between current seq and parent
            new_pairs += [frozenset([AlignCompare.MsaChar(ancname,ancchar,ancpos),AlignCompare.MsaChar(tree.name,seq[i],p)])]

         if ancpos != -1 and mapping[p] == ancpos or i == site:
            new_homs |= set([(p,seq[i])])

            if not tree.is_terminal():
               (homs1,pairs1) = _get_pairs_site(result,tree[0],site,p,seq[i],tree.name)
               (homs2,pairs2) = _get_pairs_site(result,tree[1],site,p,seq[i],tree.name)

               #collect all homs into homologs1 and homologs2
               for homs in [homs1,homs2]:
                  for s,h in homs.iteritems():
                     try:
                        homologs[s] |= h
                     except:
                        homologs[s] = set(h)

               if i == site:
                  pairs += pairs1 + pairs2
                  #compute all pairs between hom1 and hom2
                  for s1,hh1 in homs1.iteritems():
                     for s2,hh2 in homs2.iteritems():
                        pp = []
                        for h1 in hh1:
                           for h2 in hh2:
                               pp += [frozenset([AlignCompare.MsaChar(s1,h1[1],h1[0]),AlignCompare.MsaChar(s2,h2[1],h2[0])])]
                        assert(len(pp) > 0)
                        pairs += [pp]

         p += 1
      elif i == site and not tree.is_terminal():
         (homs1,pairs1) = _get_pairs_site(result,tree[0],site,-1,seq[i],tree.name)
         (homs2,pairs2) = _get_pairs_site(result,tree[1],site,-1,seq[i],tree.name)
         pairs += pairs1 + pairs2

   if not onlyleaves:
      if len(new_pairs) > 0:
         pairs += [new_pairs]
   if not onlyleaves or tree.is_terminal():
      if len(new_homs) > 0:
         homologs[tree.name] = new_homs

   return (homologs,pairs)

def write_treks(out,result):
   for name in result.leaves:
      info = result.info[name]
      seq = info.seq
      hom = info.tr_homology
      units = _split_units(info)

      all_sites = set(hom)
      all_sites.remove(-1)
      all_sites = list(all_sites)
      all_sites.sort()

      p = dict(zip(all_sites,range(len(all_sites))))
      assert(not -1 in p)

      for u in units:
         for i in xrange(u.start,u.start+u.length):
            hom[i] = p[hom[i]]

      all_sites = set(hom)
      all_sites.remove(-1)

      out.write(">%s\n"%name)
      unit_len = len(all_sites)
      start = units[0].start
      end = units[-1].start + units[-1].length - 1
      region_len = sum(u.length for u in units)
      out.write("Length: %d residues - nb: XXX  from  %d to %d - Psim:1.0 region Length:%d\n"%(unit_len,start+1,end+1,region_len))

      for u in units:
         s = bytearray(unit_len*Bio.Alphabet.Gapped(seq.alphabet).gap_char)
         for c,h in zip(seq[u.start:u.start+u.length],hom[u.start:u.start+u.length]):
            if h != -1:
               s[h] = c
         out.write(str(s+"\n"))

      out.write("**********************\n\n")

def _invert_map(mapping):
   out = dict()
   for i,v in enumerate(mapping):
      out[v] = i
   return out

class _Identdict:
   def __init__(self,dict):
      self.dict = dict

   def __getitem__(self,key):
      try:
         return self.dict[key]
      except:
         return key

def _counter(start=0,increment=1):
   i = start
   while True:
      ii = yield i
      if not ii is None:
         i = ii
      else:
         i += increment

def _alternateiter(*args):
   args = [iter(a) for a in args]
   while True:
      for a in args:
         yield a.next()

def merge_tr_parentmap(clade,info):
   ext_tr_mapping = _merge_tr_parentmap_up(clade,info,_counter(-100,-1))
   _merge_tr_parentmap_down(clade,info,ext_tr_mapping)

def _merge_tr_parentmap_up(clade,info,counter):
   if hasattr(info[clade.name],"tr_parentmap"):
      tr_parentmap = info[clade.name].tr_parentmap
   else:
      tr_parentmap = range()
   for i in xrange(len(tr_parentmap)):
      if tr_parentmap[i] == -1:
         tr_parentmap[i] = counter.next()
   ext_tr_parentmap = tr_parentmap[:]

   if clade.is_terminal():
      gaps = [list() for i in xrange(len(tr_parentmap+1))]
      for c in clade:
         chld_pmap = _merge_tr_parentmap_up(c,info,counter)
         index = -1
         for m in chld_pmap:
            if m >= 0:
               index = m
            else:
               gaps[index+1].append(m)

      ext_tr_parentmap = reduce(lambda x,y:x+y,_alternateiter(gaps,[[m] for m in ext_tr_parentmap]),[])

   return ext_tr_parentmap

def _merge_tr_parentmap_down(clade,info,ext_tr_mapping):
   tr_homology = info[clade.name].tr_homology
   imap = _Identdict(_invert_map(info[clade.name].tr_parentmap))

   ext_tr_mapping2 = [imap[i] for i in ext_tr_mapping]
   imap2 = _Identdict(_invert_map(ext_tr_mapping2))

   tr_homology[:] = [imap2[i] for i in tr_homology]

   if clade.is_terminal():
      for c in clade:
         _merge_tr_parentmap_down(c,info,ext_tr_mapping2)


if __name__ == "__main__":
   while True:
      tree = Bio2.Phylo.BDTree.simulate(5,1,10)
      hmm = HMMERProfile.HMMERProfile("LRR2.hmm",binary="/cluster/home/infk/sadam/WrapperBinaries/Other/hmmer-2.3.2/bin/hmmemit")
      profiles = iter(hmm)
      length1 = scipy.stats.geom.rvs(1.0/150)
      length2 = scipy.stats.geom.rvs(1.0/150)
      rootseq1 = Bio2.SeqEvolution.random_sequence(length1,pi=Bio2.SeqEvolution.WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein)
      rootseq2 = Bio2.SeqEvolution.random_sequence(length2,pi=Bio2.SeqEvolution.WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein)
      repeats,homology = reduce(lambda x,y: (x[0]+y[0],x[1]+y[1]), [profiles.next() for i in xrange(3)], ("",[]))

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
                  SlippageFactory(profiles,ins_rate=1,del_rate=1,length_distr=scipy.stats.geom(1.0/1.5,0)),
                  #SlippageDupFactory(ins_rate=1,del_rate=1,length_distr=scipy.stats.geom(1.0/50,0)),
                  Bio2.SeqEvolution.GeometricInsertionFactory(rate=0.005,alen=3.5),
                  Bio2.SeqEvolution.GeometricDeletionFactory(rate=0.005,alen=3.5)
               ]

      result = Bio2.SeqEvolution.Simulator.SimulationResult(rootseq=rootseq,tree=tree,events=events)

      print(result.msa)

      #propagate tr_parentmap if present + merge
      merge_tr_parentmap(tree.clade,result.info)

      Bio.AlignIO.write([result.msa],open("repeats.fasta",'w'),"fasta")
      Bio.AlignIO.write([result.ancestral_msa],open("repeats_anc.fasta",'w'),"fasta")
      write_treks(open("repeats.treks",'w'),result)
      #write_treks(open("/dev/null",'w'),result)

      if False:
         #write_pairs(open("repeats.pairs",'w'),result)
         write_pairs(open("/dev/null",'w'),result)
         pairs = get_pairs(result)

         sp1 = AlignCompare.SPtd1(pairs, result.msa)
         sp2 = AlignCompare.SPtd2(pairs, result.msa)

         print(sp1)
         print(sp2)
