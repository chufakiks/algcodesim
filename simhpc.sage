RR = RealField(1000)
F.<alpha> = GF(2^8)
Ring.<X> = PolynomialRing(F)
n = 255
k = 245
p = RR(0.001)
x_pts = [alpha^i for i in range(n)]

def lagrange_basis(i,x_pts, rrange):
    return prod((X - x_pts[h]) / (x_pts[i] - x_pts[h]) for h in range(rrange) if h != i)

def buildG(k):
    def Aij(i,j,k,x_pts):
        Li = lagrange_basis(i,x_pts,k)
        return Li(x_pts[k + j])
    A = matrix(F, k, n - k, lambda i,j : Aij(i,j,k,x_pts))
    return block_matrix([[identity_matrix(F, k), A]])

def G(x_pts):
    return prod(X - x_pts[i] for i in range(n))

def R(x_pts, r):
    return sum(r[i] * lagrange_basis(i, x_pts, n) for i in range(n))

def buildP(r, x_pts):
    return matrix(Ring, [[G(x_pts), Ring(0)], [-R(x_pts, r), Ring(1)]])

def isitkreduced(P, k):
    return kminus1degree(P[0,0],P[0,1],k) < kminus1degree(P[1,0],P[1,1],k)

def kreduce(P):
    if isitkreduced(P,k):
        return P
    Q,R = P[0,0].quo_rem(P[1,0])
    P = matrix(Ring, [[Ring(0), Ring(1)],[Ring(1), -Q]]) * P
    return kreduce(P)

def kminus1degree(A,B,k):
    return max(A.degree(),B.degree() + k - 1)

def fastinterpolation(r, x_pts):
    P = buildP(r, x_pts)
    P = kreduce(P)
    if kminus1degree(P[0,0],P[0,1],k) <= kminus1degree(P[1,0],P[1,1],k):
        f = -P[0,0]/P[0,1]
    else:
        f = -P[1,0]/P[1,1]
    return vector([f(x_pts[i]) for i in range(k)])

def qarychannel(p, codeword):
    for i in range(n):
        if random() < p:
            e = F.random_element()
            while e == codeword[i]:
                e = F.random_element()
            codeword[i] = e
    return codeword

def sim(Gmatrix,k):
    m = vector([F.random_element() for _ in range(k)])
    codeword = list(m*Gmatrix)
    r = qarychannel(p, codeword)
    decoded = fastinterpolation(r,x_pts)
    if m == decoded:
        return True
    return False

num_trials = 1000

Gmatrix = buildG(k)

failures = 0
for i in range(num_trials):
    success = sim(Gmatrix, k)
    if not success:
        failures += 1 
    
print("Number of failures: ", failures)
print("Rate of failures: ", failures/num_trials)