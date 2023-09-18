import os
import time
from fpylll import IntegerMatrix, LLL

"""
  Implementation of Coppersmith's method.
  Params:
  - polys, a set of polynomials
  - modulus, an integer
  - bounds, a set of integers.

  This function uses Coppersmith's method to try to compute all small common roots of the polynomials in polys
  modulo modulus, whose absolute values are upper bounded by the bounds in bounds.

  When given additional information about the small roots, this information may be encoded as polynomials into gbRelations.
  These polynomials are then added to the ideal, of which Coppersmith's method computes the Gröbner basis to extract the small roots.

  Set verbose = True to run in verbose mode.
"""
def coppersmithsMethod( polys, modulus, bounds, gbRelations = [], verbose=False ):
  
  R = polys[0].parent()
  
  for poly in polys:
    if poly.parent() != R:
      raise ValueError("Can't instantiate coppersmiths method with polynomials from different rings.")
  
  #Create basis matrix
  tt = cputime()
  
  monList = [] #Stores a list of all monomials in the basis.
  monDict = {} #Stores which monomial corresponds to which column in the basis.
  
  for poly in polys:
    for mon in poly.monomials():
      if mon not in monDict:
        index = len(monDict)
        monDict[mon] = index
        monList.append(mon)

  rows = len(polys)
  cols = len(monList)
  B = zero_matrix( ZZ, rows, cols )
  
  for i, poly in enumerate(polys):
    for mon in poly.monomials():
      B[i, monDict[mon]] = int( poly.monomial_coefficient(mon) * mon(*bounds) )
  
  if verbose:
    print("Finished basis generation. Polynomials: %d. Time: %fs." % (len(polys), cputime(tt)))
  
  #Reduce basis
  start = time.time() #Can't use cputime here, because it does not catch os.system(...).  

  with open("basis.tmp", "w+") as f:
    #Would prefer to simply use IntegerMatrix here, but IntegerMatrix.__str()__ results in crash.
    B_str = B.str()
    B_str = '\n'.join(' '.join(line.split()) for line in B_str.split('\n'))
    f.write( "[\n" + B_str + "\n]")
  
  success = os.system("flatter -v basis.tmp basis_out.tmp >/dev/null 2>&1")
  os.remove("basis.tmp")
  
  if success == 0:
    B_LLL = matrix(IntegerMatrix.from_file("basis_out.tmp"))
    os.remove("basis_out.tmp")
  else:
    if verbose:
      print("flatter not found. Resorting to FPLLL.")
    B_LLL = B.LLL()
  
  stop = time.time()
  
  #TODO: Catch non-full-rank matrices
  
  if verbose:
    print("Finished basis reduction. Time: %fs." % (stop-start))
 
  #Extract short polynomials
  tt = cputime()
  
  solutionPolynomials = gbRelations
  for v in B_LLL:
    sqNorm = sum( [v_i^2 for v_i in v] )
    norm = RR(sqrt(sqNorm))

    if norm < RR(modulus / sqrt( B_LLL.ncols() )):
      
      i = 0

      poly = R(0)
      while i < len(monList):
        mon = monList[i]
        poly += R( ZZ( v[i] / mon(*bounds) ) ) * mon
        i += 1
      
      solutionPolynomials.append(poly)
    
  if verbose:
    print("Found " + str(len(solutionPolynomials)) + " short polynomials. Time: %fs." % cputime(tt))
    
  #Extract solutions
  #We compute the Gröbner basis over finite fields, because Sage's default is too slow.
  tt = cputime()
  
  k = len(R.gens())
  
  if len(solutionPolynomials) < k:
    raise RuntimeError("LLL did not find enough short polynomials. Can't extract solution.")
  
  p = 0
  maxBound = max(bounds)
  gbModulus = 1
  gbFailCounter = 0
  
  crtResults = [ [] for _ in range(k) ]
  moduli = []
  
  while gbModulus < maxBound:
    p = next_prime(p+1)
    
    R = R.change_ring(GF(p))
    I = R*solutionPolynomials
    
    success = True
    
    try:
      solutions = I.variety()
    except ValueError:
      success = False
    
    if success and len(solutions) == 1:
      #TODO: Catch non-unique solutions.
      
      solutions = solutions[0]
      
      gbModulus *= p
      moduli.append(p)
      
      for i in range( k ):
        result_i = solutions[ R.gens()[i] ]
        crtResults[i].append( ZZ(result_i) )
      
    else:
      gbFailCounter +=1
      
      #We stop when the GB computation failed more than 100 times.
      #The value of 100 was chosen arbitrarily, and may be changed, if necessary.
      if gbFailCounter > 100:
        raise RuntimeError("Coppersmith heuristic failed. Could not extract solution from Gröbner basis.")

  solutions = [ crt( crtResults[i], moduli ) for i in range(k) ]
    
  if verbose:
    print("Finished extracting solutions. Time: %fs." % cputime(tt))
    
  return solutions