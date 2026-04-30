RR = RealField(1000)
F.<alpha> = GF(2^8)
Ring.<X> = PolynomialRing(F)
n = 255
k = 245
pex = 1 - (1 - 0.001)^8
p= RR(pex)
t = (n - k) // 2

def T(j, l, w):
    return sum ( binomial(w,z) * binomial(w-z, z+j-l) * (q - 2^(l+w-j-2*z)) * binomial(n-w, n-j-z) * (q-1)^(j+z-w) for z in range(0, w+1))

def Tbig(j, w, t):
    return sum(T(j, l, w) for l in range(0, t+1))

def A_w(n, k, w):
    d = n - k + 1
    if w < d:
        return 0
    return binomial(n, w) * sum( (-1)^j * binomial(w, j) * (q^(w - d - j + 1) - 1) for j in range(0, w - d + 1))

Perr = sum( (p/(q-1))**j * (1-p)**(n-j) * sum(A_w(n,k,w) * Tbig(j,w,t) for w in range(1, n+1)) for j in range(t+1, n+1))

print(Perr)