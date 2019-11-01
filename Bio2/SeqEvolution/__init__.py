import Bio.Seq
import Bio.Alphabet.IUPAC
import Bio2.SeqEvolution
from Bio2.SeqEvolution import WAG
import scipy.stats
import random
import numpy
import copy

class SimulatorSequence:
   def __init__(self,seq=None,alphabet=None):
      if not seq is None:
         seq = str(seq)
      else:
         seq = ""

      self.seq = Bio.Seq.MutableSeq(seq,alphabet)
      self.mapping = range(len(self.seq))
      self.events = []

   def __deepcopy__(self,memo):
      new = copy.copy(self)
      new.seq = Bio.Seq.MutableSeq(new.seq.tostring(),new.seq.alphabet)
      return new


class Event:
   def __init__(self):
      raise NotImplemented()

   def apply(self,sequence):
      raise NotImplemented()

   def get_distance(self):
      return self._distance


class EventFactory:
   def __init__(self):
      raise NotImplemented()

   def get_rate(self):
      return 0.0;

   def get_event(self,sequence,distance):
      raise NotImplemented()

   def update(self,sequence,event):
      return None


class IndelEvent(Event):
   def get_start(self):
      return self._start

   def get_length(self):
      return self._length


class InsertionEvent(IndelEvent):
   def __init__(self,distance,start,length,sequence):
      self._distance = distance
      self._start = start
      self._length = length
      self._sequence = sequence

   def get_sequence(self):
      return self._sequence

   def apply(self,sequence):
      sequence.seq = sequence.seq[:self._start] + self._sequence + sequence.seq[self._start:]
      sequence.mapping = sequence.mapping[:self._start] + self._length*[-1] + sequence.mapping[self._start:]


class DeletionEvent(IndelEvent):
   def __init__(self,distance,start,length):
      self._distance = distance
      self._start = start
      self._length = -length

   def apply(self,sequence):
      sequence.seq = sequence.seq[:self._start] + sequence.seq[self._start-self._length:]
      sequence.mapping = sequence.mapping[:self._start] + sequence.mapping[self._start-self._length:]


class InsertionFactory(EventFactory):
   def __init__(self,rate,length_distribution,pi=WAG.pi):
      self._rate = rate
      self._lendist = length_distribution
      self._pi = pi
      self._seqlen = 0

   def get_rate(self):
      return self._seqlen * self._rate

   def get_event(self,sequence,distance):
      start = int(random.random()*(self._seqlen+1))
      length = int(self._lendist.rvs())
      sequence = random_sequence(length,self._pi)
      return InsertionEvent(distance,start,length,sequence)

   def update(self,sequence,event):
      self._seqlen = len(sequence.seq)

class GeometricInsertionFactory(InsertionFactory):
   def __init__(self,rate,alen,pi=WAG.pi):
      InsertionFactory.__init__(self,rate,scipy.stats.geom(1.0/alen,0),pi)


class DeletionFactory(EventFactory):
   def __init__(self,rate,length_distribution):
      self._rate = rate
      self._lendist = length_distribution
      self._seqlen = 0

   def get_rate(self):
      return self._seqlen * self._rate

   def get_event(self,sequence,distance):
      start = int(random.random()*self._seqlen)
      length = min(int(self._lendist.rvs()),self._seqlen-start)
      return DeletionEvent(distance,start,length)

   def update(self,sequence,event):
      self._seqlen = len(sequence.seq)

class GeometricDeletionFactory(DeletionFactory):
   def __init__(self,rate,alen):
      DeletionFactory.__init__(self,rate,scipy.stats.geom(1.0/alen,0))


class SubstitutionEvent(Event):
   def __init__(self,distance,position,character):
      self._distance = distance
      self._pos = position
      self._char = character

   def apply(self,sequence):
      sequence.seq[self._pos] = self._char


class SubstitutionFactory(EventFactory):
   def __init__(self,Q=WAG.Q,rate=1.0):
      self._Q = Q
      self._counts = numpy.empty((Q.dim,1),dtype=numpy.int32)
      self._rate = rate
      self._rates = None

   def get_rate(self):
      return self._rate * -(numpy.matrix(self._Q).diagonal() * self._counts).sum()

   def get_event(self,sequence,distance):
      pos = _sample_pos(self._rates)
      old = self._Q.indexof[sequence.seq[pos]]
      row = self._Q[old].A1.copy()
      row[old] = 0
      new =_sample_character(row,self._Q.alphabet)

      return SubstitutionEvent(distance,pos,new)

   def update(self,sequence,event):
      if isinstance(event,SubstitutionEvent):
         old = self._Q.indexof[sequence.seq[event._pos]]
         new = self._Q.indexof[event._char]
         self._counts[old] -= 1
         self._counts[new] += 1
         self._rates[event._pos] = -self._Q[new,new]
      elif isinstance(event,InsertionEvent):
         old_rates = self._rates
         startpos = event.get_start()
         endpos = event.get_start() + event.get_length()
         self._rates = numpy.empty((len(sequence.seq),))
         self._rates[:startpos] = old_rates[:startpos]
         self._rates[endpos:] = old_rates[startpos:]
         index = self._Q.indexof
         for i,c in enumerate(sequence.seq[startpos:endpos]):
            ic = index[c]
            self._counts[ic] += 1
            self._rates[startpos+i] = -self._Q[ic,ic]
      else:
         self._counts.fill(0)
         self._rates = numpy.empty((len(sequence.seq),))
         index = self._Q.indexof
         for i,c in enumerate(sequence.seq):
            ic = index[c]
            self._counts[ic] += 1
            self._rates[i] = -self._Q[ic,ic]


def _sample_pos(column):
   column = numpy.array(column)
   target = random.random()*column.sum()
   for i in xrange(len(column)):
      target -= column[i];
      if target <= 0:
         return i

def _sample_character(column,alphabet):
   i = _sample_pos(column)
   return alphabet.letters[i]

def random_sequence(length=None,pi=WAG.pi,alphabet=Bio.Alphabet.IUPAC.protein):
   if length is None:
      length = scipy.stats.geom.rvs(1.0/300.0)

   seq = numpy.zeros((length,),dtype='c')
   for i in xrange(length):
      seq[i] = _sample_character(pi,alphabet)

   return Bio.Seq.Seq(seq.tostring(),alphabet)

default_factories = [SubstitutionFactory(),GeometricInsertionFactory(rate=0.005,alen=3.5),GeometricDeletionFactory(rate=0.005,alen=3.5)]

def mutate_sequence(seq,dist,event_factories=default_factories,return_all=False):
   if not isinstance(seq,SimulatorSequence):
      seq = SimulatorSequence(seq.tostring(),seq.alphabet)
   else:
      seq = copy.deepcopy(seq)
      seq.mapping = range(len(seq.seq))
      seq.events = []

   d = 0
   last_event = None

   while True:
      map(lambda x:x.update(seq,last_event),event_factories)
      rates = [f.get_rate() for f in event_factories]
      d += scipy.stats.expon.rvs(0,1.0/sum(rates))

      if d >= dist:
         break

      fac = event_factories[_sample_pos(rates)]
      last_event = fac.get_event(seq,d)
      last_event.apply(seq)

      if return_all:
         seq.events += [last_event]

   seq.seq = seq.seq.toseq()

   if return_all:
      return seq
   else:
      return seq.seq

