import copy

#Tree-based approach to construct optimal shift-polynomial p, using the polynomials in polys, such that
# - p has leading monomial mon,
# - p is defined over M.
def getBestShiftPoly(mon, polys, M, poly = 1, label=0, best_label = 0, best_poly = 1, start=0):
    R = polys[0].parent()
    gens = R.gens()

    n = len(polys)
    
    if label==0:
        label = [0]*n
        best_label = [0]*n
        best_poly = R(1)
    
    shift_poly = poly*mon
    
    if set(shift_poly.monomials()).issubset(M):
        if sum(best_label) <= sum(label):
            best_label = label
            best_poly = shift_poly
    
        for i in range(start,n):
            lm = polys[i].lm()
            if mon % lm == 0:
                label_new = deepcopy(label)
                label_new[i] += 1

                poly_new = poly * polys[i]

                mon_new = R(mon/lm)

                best_label, best_poly = getBestShiftPoly(mon_new, polys, M, poly_new, label_new, best_label, best_poly, i)
            
    return best_label, best_poly

#Constructs an optimal (M,lex)-suitable set of shift-polynomials,
#using the polynomials in polys.
#All shift-polynomials have the desired small root modulo modulus^m.
def constructOptimalShiftPolys(polys,M,modulus,m):
  F = []
  
  for mon in M:
    label, poly = getBestShiftPoly( mon, polys, M )
    poly = poly * modulus^(m - sum(label))
    
    F.append(poly)
  
  
  return F