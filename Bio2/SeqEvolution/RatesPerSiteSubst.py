import Bio2.SeqEvolution
from Bio2.SeqEvolution import WAG
import scipy.stats
import numpy

class RatesPerSiteSubstitutionFactory(Bio2.SeqEvolution.SubstitutionFactory):
   def __init__(self,rate_distr,Q=WAG.Q):
      self._Q = Q
      self._rate_distr = rate_distr
      self._rates = None

   def get_rate(self):
      return self._rates.sum()

   def update(self,sequence,event):
      if isinstance(event,Bio2.SeqEvolution.SubstitutionEvent):
         old = self._Q.indexof[sequence.seq[event._pos]]
         new = self._Q.indexof[event._char]
         rate = sequence.rates[event._pos]
         self._rates[event._pos] = -self._Q[new,new] * rate
      elif isinstance(event,Bio2.SeqEvolution.DeletionEvent):
         startpos = event.get_start()
         endpos = event.get_start() - event.get_length()
         old_rates = self._rates
         self._rates = numpy.empty((len(sequence.seq),))
         self._rates[:startpos] = old_rates[:startpos]
         self._rates[startpos:] = old_rates[endpos:]
         old_rates = sequence.rates
         sequence.rates = numpy.empty((len(sequence.seq),))
         sequence.rates[:startpos] = old_rates[:startpos]
         sequence.rates[startpos:] = old_rates[endpos:]
      elif isinstance(event,Bio2.SeqEvolution.InsertionEvent):
         startpos = event.get_start()
         endpos = event.get_start() + event.get_length()
         old_rates = self._rates
         self._rates = numpy.empty((len(sequence.seq),))
         self._rates[:startpos] = old_rates[:startpos]
         self._rates[endpos:] = old_rates[startpos:]
         old_rates = sequence.rates
         if len(sequence.rates) == len(sequence.seq):
            index = self._Q.indexof
            for i,c in enumerate(sequence.seq[startpos:endpos]):
               ic = index[c]
               rate = sequence.rates[startpos+i]
               self._rates[startpos+i] = -self._Q[ic,ic] * rate
         else:
            sequence.rates = numpy.empty((len(sequence.seq),))
            sequence.rates[:startpos] = old_rates[:startpos]
            sequence.rates[endpos:] = old_rates[startpos:]

            index = self._Q.indexof
            for i,c in enumerate(sequence.seq[startpos:endpos]):
               ic = index[c]
               rate = self._rate_distr.rvs()
               sequence.rates[startpos+i] = rate
               self._rates[startpos+i] = -self._Q[ic,ic] * rate
      elif event is None:
         self._rates = numpy.empty((len(sequence.seq),))
         if not hasattr(sequence,"rates"):
            sequence.rates = numpy.empty((len(sequence.seq),))
            for i in xrange(len(sequence.seq)):
               sequence.rates[i] = self._rate_distr.rvs()
         
         self._rates = numpy.empty((len(sequence.seq),))
         index = self._Q.indexof
         for i,c in enumerate(sequence.seq):
            ic = index[c]
            self._rates[i] = -self._Q[ic,ic]
         self._rates *= sequence.rates
      else:
         raise NotImplemented()

class GammaRatesSubstitutionFactory(RatesPerSiteSubstitutionFactory):
   def __init__(self,Q=WAG.Q,shape=1.0,rate=1.0):
      gamma_distr = scipy.stats.gamma(float(shape),scale=float(rate)/float(shape))
      RatesPerSiteSubstitutionFactory.__init__(self,gamma_distr,Q)

