import numpy

eps = 1e-6

__all__ = ['RateMatrix','S','Q']

class RateMatrix:
   pass

class S(RateMatrix):
   """
   This class (S lower triangular exchangeability matrix + aminoacid rates)
   is the common exchange format among the rate matrices defined in this
   module.
   """

   def r(self):
      return self.__r

   def pi(self):
      return self.__pi

   def S(self):
      return self

   def __init__(self,matrix,pi,alphabet):
      """
      Pass a symmetrical (or lower triangle) exchangeability matrix and aminoacid frequency vector
      as ndarrays.
      """

      self.alphabet = alphabet
      self.dim = len(alphabet.letters)
      self.__pi = numpy.array(pi).reshape((self.dim))/numpy.sum(pi)

      self.indexof = dict()
      for i in xrange(len(alphabet.letters)):
         self.indexof[alphabet.letters[i]] = i

      try:
	 self.__S = numpy.matrix(matrix,dtype=numpy.double).reshape((self.dim,self.dim))
      except:
	 raise TypeError("provided matrix has wrong shape")

      self.__normalize()

   def __normalize(self):
      """
      normalize exchangeability matrix by scaling the sum of off-diagonal elements to dim^2
      """

      self.__S = numpy.multiply(self.__S,numpy.tri(self.dim,k=-1))
      self.__r = numpy.sum(self.__S)*2.0 # because the upper triangular matrix is zero
      self.__S *= self.dim*self.dim/self.__r

   def __array__(self,dtype=None):
      if dtype:
	 return numpy.array(self.__S)
      else:
	 return numpy.array(self.__S,dtype=dtype)

   def __getitem__(self,item):
      return self.__S.__getitem__(item)


class Q(RateMatrix):
   """
   A (normalized) rate matix
   """

   def r(self):
      return self.__r

   def pi(self):
      if self.__pi == None:
	 (w,V) = numpy.linalg.eig(self.__Q.T)
	 (i,n) = max(enumerate(w),key=lambda x:x[1])
	 assert(abs(n) < eps)
	 self.__pi = numpy.array(V[:,i]).reshape((self.dim))
	 self.__pi /= numpy.sum(self.__pi)
      return self.__pi

   def S(self):
      pi = self.pi()
      s = self.__Q * numpy.diag(1/pi)
      return S(s,pi,self.alphabet)

   def __init__(self,matrix,alphabet=None):
      """
      Either pass an S exchangeability matrix object or a ndarray (the diagonal
      will be filled automatically)
      """
      self.__pi = None

      if isinstance(matrix,RateMatrix):
	 S = matrix.S()
	 pi = S.pi()
         self.dim = S.dim
	 Q = numpy.matrix(S)
	 Q = Q + Q.T
	 Q = Q * numpy.diag(pi)

         self.alphabet = S.alphabet
         self.indexof = S.indexof
	 self.__Q = Q
	 self.__pi = pi
      else:
         self.alphabet = alphabet
         self.dim = len(alphabet.letters)

         self.indexof = dict()
         for i in xrange(len(alphabet.letters)):
            self.indexof[alphabet.letters[i]] = i

	 try:
	    self.__Q = numpy.matrix(matrix,dtype=numpy.double).reshape((self.dim,self.dim))
	 except:
	    raise TypeError("provided matrix has wrong shape")

      self.__normalize()
      if not self.__check():
	 raise ValueError("provided matrix does not obey detailed balance (pi_i*Q_ij = pi_j*Q_ji)")

   def __array__(self,dtype=None):
      if dtype:
	 return numpy.array(self.__Q)
      else:
	 return numpy.array(self.__Q,dtype=dtype)

   def __getitem__(self,item):
      return self.__Q.__getitem__(item)

   def __pow__(self,t):
      (sigma,V) = numpy.linalg.eig(self.__Q)
      V = numpy.matrix(V)
      sigma = numpy.power(sigma,t)
      M = V*numpy.diag(sigma)*V.I
      return M

   def __rpow__(self,t):
      (sigma,V) = numpy.linalg.eig(self.__Q)
      V = numpy.matrix(V)
      sigma = numpy.power(t,sigma)
      M = V*numpy.diag(sigma)*V.I
      return M

   def exp(self,t=1.0):
      '''
      return matrix exponential of t*self
      '''

      (sigma,V) = numpy.linalg.eig(t*self.__Q)
      V = numpy.matrix(V)
      sigma = numpy.exp(sigma)
      M = V*numpy.diag(sigma)*V.I
      return M

#   def log(self,base=None):
#      '''
#      return matrix logrithm to given base (default=e)
#      '''
#
#      (sigma,V) = numpy.linalg.eig(self.__Q)
#      V = numpy.matrix(V)
#      sigma = numpy.log(sigma,base)
#      M = V*numpy.diag(sigma)*V.I
#      return M

   def __check(self):
      P = numpy.diag(self.pi())
      Q = self.__Q

      return numpy.all(numpy.abs(P*Q - Q.T*P) < eps)

   def __normalize(self):
      # compute diagonal elements
      for i in xrange(self.dim):
	 self.__Q[(i,i)] = 0
      r = numpy.sum(self.__Q,axis=1).A1
      self.__Q -= numpy.diag(r)

      # normalize sum of diagonal elements to minus one
      self.__r = numpy.sum(r*self.pi())
      self.__Q /= self.__r

#def R_to_M(R,pi,t):
#   P = diag(pi)
#   R = lower_triangular(R)
#   R = R + R.T
#   (sigma,V) = numpy.linalg.eigh(Q)
#   V = matrix(V)
#   M = P.I*V*diag(exp(sigma,t))*V.I*P
#   return M

def change_alphabet(initial,target,m):
   """
   change_alphabet(initial,target,matrix):

   initial: initial (source) alphabet as list of characters
   target:  target alphabet
   m:       input matrix with shape (len(initial),len(initial)) or (1,len(initial)) or (len(initial),1)
   """

   n = numpy.matrix(numpy.zeros((len(initial),len(target))))
   for x in xrange(len(target)):
      y = initial.index(target[x])

      n[y,x] = 1

   if len(m.shape) == 2 and m.shape[1] == len(initial):
      m = n.T * m
   if m.shape[0] == len(initial):
      m = m * n
   return m

