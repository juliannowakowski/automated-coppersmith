load("attack.sage")

p = next_prime( ZZ.random_element( 2^(n-1), 2^n ) )
while p%8 != 7:
  p = next_prime( p+1 )

C = ZZ( GF(p).random_element( ) )
C = (power_mod(C,2,p) - 2) % p #Ensure that C+2 is a square.
A = ( (C+6) * inverse_mod( 2 * power_mod(C+2,(p+1)/4,p), p  ) ) % p
B = ( (A+6) * inverse_mod( 2 * power_mod(A+2,(p+1)/4,p), p  ) ) % p

unknown_bits = n-k

A_MSB, A_LSB = split(A,unknown_bits)
B_MSB, B_LSB = split(B,unknown_bits)
C_MSB, C_LSB = split(C,unknown_bits)

R.<z,x,y> = PolynomialRing(QQ, order="lex")
f = (A_MSB+x)^2 + 12*(A_MSB+x) - 4*(A_MSB+x)*(B_MSB+y)^2 - 8*(B_MSB+y)^2 + 36
g = (C_MSB+z)^2 + 12*(C_MSB+z) - 4*(C_MSB+z)*(A_MSB+x)^2 - 8*(A_MSB+x)^2 + 36

X = 2^unknown_bits
Y = 2^unknown_bits
Z = 2^unknown_bits

polys = [f,g]
bounds = [X,Y,Z]
solutions_verify = [C_LSB,A_LSB,B_LSB]

print( "Starting attack on CSURF with n=%d and k=%d." % (n,k) )
attack(polys, bounds, p, m//2, solutions_verify)