#
# Finds the phi (3-body resonant angles) solutions for a system of n planets
#
# Enid M Cruz-Colón 
# 19 July 2017
#
# Department of Astronomy & Astrophysics -- UChicago
# PI Dan Fabrycky
#
# For a better understanding of the algorithm see
# Delisle 2017 arXiv:1706.09613 [astro-ph.EP]
#

from numpy import array, loadtxt, arange, empty, zeros, ones, dot
from pylab import plot, scatter, show, xlabel, ylabel, xlim, ylim, legend
from math import cos, sin, pi, sqrt
from numpy.linalg import solve
from time import time

t1 = time() # just to see how much time it takes to run the code

G = 2.95912203632e-4 # Newton's constant, AU^3*Msolar^(-1)*days^(-2)

# The next parameters are specific for a planetary system
n = 4
k = loadtxt("kvals.txt",int)
mass = loadtxt("masses.txt")
Cdir = loadtxt("Cdir.txt")
Cind = loadtxt("Cind.txt")

# Do not change anything from here down
c = ones(n-2,int)

def q(i,j): # order of resonance
    return k[j,i]-k[i,j]

def l(i,j):
    if j > i:
        return 0
    sigma = 0
    for r in range(j,i+1):
        p = 1
        for s in range(r+1,i+1):
            p *= k[s+1,s]/k[s-1,s]
        sigma += q(r,r+1)*p
    return sigma/k[i,i+1]

def phif(i):
    arr = zeros(n,int)
    if i == (n-2):
        arr[n-2] = 1
        arr[n-1] = -1
    elif i == (n-1):
        arr[n-2] = -k[n-2,n-1]/q(n-2,n-1)
        arr[n-1] = k[n-1,n-2]/q(n-2,n-1)
    else:
        arr[i] = k[i,i+1]/q(i,i+1)
        arr[i+2] = k[i+2,i+1]/q(i+1,i+2)
        arr[i+1] = -(k[i+1,i]/q(i,i+1)+k[i+1,i+2]/q(i+1,i+2))
    return arr

phi = []
for i in range(n):
    phi.append(phif(i))
    if i < n-2:
        m = min(abs(j) for j in phi[i] if abs(j) > 0)
        for j in range(m,0,-1):
            mod = phi[i] % j
            gcm = True
            for num in mod:
                if num != 0:
                    gcm = False
                    break
            if gcm:
                c[i] = j
                break
        phi[i] /= c[i]

def lambdaf(i):
    arr = zeros(n,float)
    for j in range(i,n-2):
        arr[j] = c[j]*l(j,i)
    p = 1
    for j in range(i,n-1):
        p *= k[j+1,j]/k[j,j+1]
    arr[n-2] = k[n-2,n-1]/q(n-2,n-1)*p
    arr[n-1] = 1

    ret = zeros(n,int)
    for i in range(len(ret)):
        ret[i] = round(arr[i])
    return ret

lam = []
for i in range(n):
    lam.append(lambdaf(i))

omegabar = []
for i in range(n):
    omegabar.append(zeros(n,int))
    omegabar[i][n-1] = 1

p = []
for i in range(n):
    p.append([])
    for j in range(n):
        p[i].append(zeros(2,int))

for i in range(n):
    for j in range(i+1,n):
        if q(i,j) == 1:
            curr = k[j,i]*lam[j] - k[i,j]*lam[i] - omegabar[i]
            for index in range(2):
                p[i][j][index]=p[j][i][index]=curr[index]

#Find C values
def mu(i):
    return G*(mass[0]+mass[i+1])

def alpha(i,j):
    return (k[i,j]/k[j,i])**(2/3) * (mu(i)/mu(j))**(1/3)

def beta(i):
    return mass[i+1]*mass[0]/(mass[0]+mass[i+1])

C = zeros([n,n],float)

for i in range(n):
    for j in range(i+1,n):
        if j == i:
            C[i,j] = 0
            continue
        
        c1 = -mass[i+1]*mass[j+1]*alpha(i,j)*2**(1/2)
        c2 = -G*mass[0]/sqrt(mu(i)*mu(j)*alpha(i,j))

        C[i,j] = c1*(1/beta(i)**(1/2))*(1/(mu(i)*alpha(i,j))**(1/4))*\
                 (Cdir[i,j] + c2*Cind[i,j])

        C[j,i] = c1*(1/beta(j)**(1/2))*(alpha(i,j)/mu(i))**(1/4)*\
                 (Cdir[j,i] + c2*Cind[j,i])

def f(phi):
    func = zeros(n-2,float)
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            for r in range(n):
                if r == i:
                    continue
                func += p[i][j]*C[i,j]*C[i,r]*sin(dot(p[i][j]-p[i][r],phi))
    return func

def phijac(phi):
    func = zeros([n-2,n-2],float)
    for l in range(n-2):
        for m in range(n-2):
            for i in range(n):
                for j in range(n):
                    if j == i:
                        continue
                    for r in range(n):
                        if r == i:
                            continue
                        func[l,m] += p[i][j][l]*C[i,j]*C[i,r]*\
                                     cos(dot(p[i][j]-p[i][r],phi))*\
                                     (p[i][j]-p[i][r])[m]
    return func

accuracy = 1e-8

a = 0.0
b = 2*pi
N = 15
h = (b-a)/N

phi_range = arange(a,b,h)
phi_solutions = []

dim = n-2
index = zeros(dim,int)

def update(i,pos,N):
    if i[pos]+1==N and pos-1>=0:
        i[pos]=0
        update(i,pos-1,N)
    else:
        i[pos]+=1

count = 0
"""
while index[0] < N:
    error = 1

    phi = zeros(dim,float)
    for i in range(len(index)):
        phi[i] = phi_range[index[i]]
    
    while error > accuracy:
        delta_phi = solve(phijac(phi),f(phi))
        phi = phi-delta_phi
        errorsq = 0
        for k in range(len(delta_phi)):
            errorsq += delta_phi[k]**2
        error = sqrt(errorsq)
    
    for k in range(len(phi)):
        while phi[k] < 0:
            phi[k] += 2*pi
        while phi[k] >= 2*pi:
            phi[k] -= 2*pi
    
    soln = True
    g = f(phi)
    for k in range(len(g)):
        if abs(g[k]) > accuracy: # not solution
            soln = False
            break
    if not soln:
        update(index,dim-1,N)
        continue   

    repeated = False
    for solution in phi_solutions:
        rep = abs(solution-phi)<=accuracy
        if False not in rep:
            repeated = True
            break
    if repeated:
        update(index,dim-1,N)
        continue

    phi_solutions.append(phi)
    update(index,dim-1,N)
    count += 1
"""    
for solution in phi_solutions:
    print(solution*180/pi)

t2 = time()

print(count,t2-t1)
# Do not change anything from here up


# The next part of the code is specific for a given planetary system                
##phi1=[]
##phi2=[]
##
##for i in range(len(phi_solutions)):
##    phi1.append(phi_solutions[i][0]*180/pi)
##    phi2.append(phi_solutions[i][1]*180/pi)
##
##plot(phi1,phi2,'+')
##xlabel(r'$\phi_1$ (deg)')
##ylabel(r'$\phi_2$ (deg)')
##xlim(-50,400)
##ylim(-50,400)
##show()

def delta_Lambdaf(i): #i >= 0
    arr = zeros(dim,int)
    eps_coef = 0
    if i == (n-1): #delLambda_n
        if (n-1)-1 >= 0:
            eps_coef = k[n-1,(n-1)-1]/q((n-1)-1,n-1)
        if (n-2)-1 >= 0:
            arr[(n-2)-1]=1/c[(n-2)-1]*k[n-1,(n-1)-1]/q((n-1)-1,n-1)
    elif i == (n-2): #deltaLambda_n-1
        if (n-1)-1 >= 0:
            eps_coef = -k[(n-1)-1,n-1]/q((n-1)-1,n-1)
        if (n-2)-1 >= 0:
            arr[(n-2)-1]=-1/c[(n-2)-1]*(k[(n-1)-1,(n-2)-1]/q((n-2)-1,(n-1)-1)+\
                                        k[(n-1)-1,n-1]/q((n-1)-1,n-1))
        if (n-3)-1 >= 0:
            arr[(n-3)-1]=1/c[(n-3)-1]*k[(n-1)-1,(n-2)-1]/q((n-2)-1,(n-1)-1)
    else:
        if i+1 < n:
            arr[i]=k[i,i+1]/(c[i]*q(i,i+1))
        if i-2 >= 0:
            arr[i-2]=1/c[i-2]*k[i,i-1]/q(i-1,i)
        if i-1 >= 0:
            arr[i-1]=-1/c[i-1]*(k[i,i-1]/q(i-1,i)+k[i,i+1]/q(i,i+1))
    return arr,int(eps_coef)
    
for i in range(n):
    print('delta Lambda_',(i+1),' = ',sep='',end='')
    deltaL,eps=delta_Lambdaf(i)
    for j in range(len(deltaL)):
        print(' + ',deltaL[j],'*delta L_',(j+1),sep='',end='')
    print(' + ',eps,'*epsilon',sep='')

delta_Lambda_l = []
epsilon_l = []
for i in range(n):
    deltaL,eps = delta_Lambdaf(i)
    delta_Lambda_l.append(deltaL)
    epsilon_l.append(eps)

delta_Lambda = array(delta_Lambda_l,int)
epsilon = array(epsilon_l,int)

def dH0_dDeltaL(j):#=0
    deltaL_coef = zeros(dim,float)
    epsilon_coef = 0
    for i in range(n):
        deltaL_coef += (k[0,i]/k[i,0])/beta(i)*sqrt(alpha(0,i)/mu(i))*\
                       delta_Lambda[i]*delta_Lambda[i,j]
        epsilon_coef += -(k[0,i]/k[i,0])/beta(i)*sqrt(alpha(0,i)/mu(i))*\
                       epsilon[i]*delta_Lambda[i,j]
    return deltaL_coef,epsilon_coef

deltaL_c = []
eps_c = []
for i in range(dim):
    Atemp,vtemp=dH0_dDeltaL(i)
    deltaL_c.append(Atemp)
    eps_c.append(vtemp)

deltaL_coef = array(deltaL_c,float)
eps_coef = array(eps_c,float)
deltaL_e = solve(deltaL_coef,eps_coef)
print(deltaL_e)

