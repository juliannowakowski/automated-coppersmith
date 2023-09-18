# Automated Coppersmith
An automated and efficient python implementation of Coppersmith's method, accompanying the paper

**Solving the Hidden Number Problem for CSIDH and CSURF via Automated Coppersmith.**


## Requirements

* [Sage](https://www.sagemath.org/)
* [flatter](https://github.com/keeganryan/flatter) (optional, for faster LLL-reduction)

If flatter is not installed, then our implementation resorts to Sage's native LLL implementation, which internally calls [fplll](https://github.com/fplll/).

Our code is written for Python 3.x / Sage 9.x.

## Tutorial

Suppose we are given the following polynomials $f,g,h$ and modulus $p$:
```py
f = x*y + 589027211857763577534285004407500599515877701606204166373378*x + 1210227469872595954120304353034891111272800242671483019591678*y + 712856912292730764373550860496219675510804490821800240884045436352648323240616556057658921760694390255205809419649548300
g = y*z + 180903529981528496775978209575327607639654632700732707438594*y + 589027211857763577534285004407500599515877701606204166373374*z + 106557101880247071038513607218438121606755644416752074440223180138737614542165434004098656165965878724132245889781596172
h = x*z + 180903529981528496775978209575327607639654632700732707438590*x + 1210227469872595954120304353034891111272800242671483019591682*z + 218934421380566537752529711732076738457083728873215369882112435160195713717062213345893111966722806427323100105889808396
p = 1255164847533299303722831198034117851714068657650841640180487
```
We want to use Coppersmith's method to compute all roots $(x_0,y_0,z_0)$ of $f$ modulo $p$ with $|x_0|,|y_0|,|z_0| \leq 2^{50}$.

We first include `coppersmithsMethod.sage` and `optimalShiftPolys.sage`:
```py
sage: load("coppersmithsMethod.sage")
sage: load("optimalShiftPolys.sage")
```
Additionally, we define $f,g,h$ and $p$:
```py
sage: R.<x,y,z> = PolynomialRing(QQ, order="lex")
sage: f = x*y + 589027211857763577534285004407500599515877701606204166373378*x + 1210227469872595954120304353034891111272800242671483019591678*y + 712856912292730764373550860496219675510804490821800240884045436352648323240616556057658921760694390255205809419649548300
sage: g = y*z + 180903529981528496775978209575327607639654632700732707438594*y + 589027211857763577534285004407500599515877701606204166373374*z + 106557101880247071038513607218438121606755644416752074440223180138737614542165434004098656165965878724132245889781596172
sage: h = x*z + 180903529981528496775978209575327607639654632700732707438590*x + 1210227469872595954120304353034891111272800242671483019591682*z + 218934421380566537752529711732076738457083728873215369882112435160195713717062213345893111966722806427323100105889808396
sage: p = 1255164847533299303722831198034117851714068657650841640180487
```

Following the strategy from our paper, we a pick a paramter $m$ and define the following set of monomials $\mathcal{M}$:
```py
sage: i = 1
sage: m = 3*i
sage: M = []
sage: for j_1 in range(i+1):
sage:   for j_2 in range(i+1):
sage:     for j_3 in range(i+1):
sage:       M += ( f^(j_1) * g^(j_2) * h^(j_3) ).monomials()
sage: M = list(set(M))   
```
Next, we construct an optimal set of shift-polynomials $\mathcal{F}$:
```py
sage: F = constructOptimalShiftPolys( [f,g,h], M, p, m )
```

Finally, we run Coppersmith's method and obtain the desired roots of $f,g,h$:
```py
sage: coppersmithsMethod( F, p^m, [2^50,2^50,2^50] )
[84100628905608, 269562626268817, 867486747849504]
```

The above example code can also be run using our `tutorial/tutorial.sage` file.

When running it on a laptop (without flatter) the output should look as follows:
```
Tutorial for Automated Coppersmith

Finished basis generation. Polynomials: 27. Time: 0.293104s.
flatter not found. Resorting to FPLLL.
Finished basis reduction. Time: 0.939303s.
Found 26 short polynomials. Time: 1.211586s.
Finished extracting solutions. Time: 0.565005s.

Found the following solutions: [84100628905608, 269562626268817, 867486747849504]
```

## Reproducing experiments from the paper

To reproduce our experiments for the Hidden Number Problem for CSIDH and CSURF please run `experiments_cihnp/csidh.sage` and `experiments_cihnp/csurf.sage`.

The parameters n, k and m have to be passed as command line arguments.

For instance, to run our experiments on CSIDH with n = 512, k = 318 and m = 3, please run
```console
sage experiments_cihnp/csidh.sage 512 318 3
```

Executing the above command should take less than 10 seconds on a laptop.