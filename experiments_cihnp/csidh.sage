load("attack.sage")

p = next_prime( ZZ.random_element( 2^(n-1), 2^n ) )
while p%4 != 3:
  p = next_prime( p+1 )

A = ZZ( GF(p).random_element( ) )
B = (2*(A+6)*inverse_mod(-A+2,p)) % p
C = (2*(A-6)*inverse_mod(A+2,p)) % p

unknown_bits = n-k

A_MSB, A_LSB = split(A,unknown_bits)
B_MSB, B_LSB = split(B,unknown_bits)
C_MSB, C_LSB = split(C,unknown_bits)

R.<x,y,z> = PolynomialRing(QQ, order="lex")

f = (A_MSB+x)*(B_MSB+y) + 2*(A_MSB+x) - 2*(B_MSB+y) + 12
g = (C_MSB+z)*(B_MSB+y) + 2*(B_MSB+y) - 2*(C_MSB+z) + 12
h = (A_MSB+x)*(C_MSB+z) - 2*(A_MSB+x) + 2*(C_MSB+z) + 12


X = 2^unknown_bits
Y = 2^unknown_bits
Z = 2^unknown_bits

polys = [f,g,h]
bounds = [X,Y,Z]
solutions_verify = [A_LSB,B_LSB,C_LSB]

print( "Starting attack on CSIDH with n=%d and k=%d." % (n,k) )
attack(polys, bounds, p, m//3, solutions_verify)