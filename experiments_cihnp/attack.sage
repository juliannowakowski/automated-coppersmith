import sys

load("../coppersmithsMethod.sage")
load("../optimalShiftPolys.sage")

n = int(sys.argv[1])
k = int(sys.argv[2])
m = int(sys.argv[3])


def attack(polys, bounds, modulus, i, solutions_verify=[]):
  print("")
  
  m = i*len(polys)
  M = ( prod(polys)^i ).monomials()
  
  tt = cputime()
  F = constructOptimalShiftPolys( polys, M, modulus, m )
  
  print("Finished construction of F. Time: %fs." % cputime(tt))

  solutions = coppersmithsMethod( F, modulus^m, bounds, verbose=True )
  if solutions_verify and solutions != solutions_verify:
    raise ValueError("Did not return correct result.")
  else:
    print("Found the following solutions:", solutions)
    
def split(x, bits):
  lsb = (x % 2^bits)
  msb = x - lsb
  
  return msb, lsb