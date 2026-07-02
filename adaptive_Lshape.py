import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

def theta(x):
    t=np.arctan2(x[1],x[0]); return t+2*np.pi if t<0 else t
def u_ex(p):
    r=np.hypot(p[0],p[1]); return 0.0 if r==0 else r**(2/3)*np.sin(2/3*theta(p))
def grad_u(p):
    r=np.hypot(p[0],p[1])
    if r==0: return np.array([0.,0.])
    th=theta(p); c=(2/3)*r**(-1/3); return np.array([-c*np.sin(th/3), c*np.cos(th/3)])
Ghat=np.array([[-1.,-1.],[1.,0.],[0.,1.]]).T

def initial_mesh():
    V=[np.array(p,float) for p in [(-1,1),(0,1),(1,1),(-1,0),(0,0),(1,0),(-1,-1),(0,-1)]]
    T=[]
    for ll,lr,ur,ul in [(3,4,1,0),(4,5,2,1),(6,7,4,3)]:
        T.append([ll,lr,ur]); T.append([ll,ur,ul])
    return V,T

def longest_edge(V,t):
    E=[(t[0],t[1]),(t[1],t[2]),(t[2],t[0])]
    return max(E,key=lambda e: (V[e[0]][0]-V[e[1]][0])**2+(V[e[0]][1]-V[e[1]][1])**2)

def refine(V,T,marked):
    mid={}
    def getmid(e):
        k=frozenset(e)
        if k not in mid:
            a,b=tuple(e); mid[k]=len(V); V.append(0.5*(V[a]+V[b]))
        return mid[k]
    em={}
    def add(ti):
        t=T[ti]
        for a,b in ((t[0],t[1]),(t[1],t[2]),(t[2],t[0])):
            em.setdefault(frozenset((a,b)),set()).add(ti)
    def rem(ti):
        t=T[ti]
        for a,b in ((t[0],t[1]),(t[1],t[2]),(t[2],t[0])):
            em[frozenset((a,b))].discard(ti)
    for ti in range(len(T)): add(ti)
    stack=list(marked)
    while stack:
        ti=stack.pop()
        e=longest_edge(V,T[ti]); efs=frozenset(e)
        nb=[tj for tj in em.get(efs,()) if tj!=ti]
        if nb and frozenset(longest_edge(V,T[nb[0]]))!=efs:
            stack.append(ti); stack.append(nb[0]); continue
        m=getmid(e); a,b=tuple(e)
        rem(ti); t=T[ti]; c=[v for v in t if v!=a and v!=b][0]
        T[ti]=[a,m,c]; add(ti); ni=len(T); T.append([m,b,c]); add(ni)
        if nb:
            tj=nb[0]
            rem(tj); t2=T[tj]; c2=[v for v in t2 if v!=a and v!=b][0]
            T[tj]=[a,m,c2]; add(tj); nj=len(T); T.append([m,b,c2]); add(nj)
    return V,T

def edge_tris(T):
    em={}
    for ti,t in enumerate(T):
        for a,b in ((t[0],t[1]),(t[1],t[2]),(t[2],t[0])):
            em.setdefault(frozenset((a,b)),[]).append(ti)
    return em

def solve(V,T):
    n=len(V); A=lil_matrix((n,n)); f=np.zeros(n)
    for t in T:
        P=[V[i] for i in t]; B=np.array([P[1]-P[0],P[2]-P[0]]).T
        area=abs(np.linalg.det(B))/2; Bi=np.linalg.inv(B); Aloc=area*(Ghat.T@(Bi@Bi.T)@Ghat)
        for i in range(3):
            for j in range(3): A[t[i],t[j]]+=Aloc[i,j]
    bnd=set(v for e,tr in edge_tris(T).items() if len(tr)==1 for v in e)
    for d in bnd:
        A.rows[d]=[d]; A.data[d]=[1.0]; f[d]=u_ex(V[d])
    return spsolve(A.tocsr(),f)

def estimator(V,T,uh):
    grads=[]; areas=[]
    for t in T:
        P=[V[i] for i in t]; B=np.array([P[1]-P[0],P[2]-P[0]]).T
        areas.append(abs(np.linalg.det(B))/2); grads.append(np.linalg.solve(B.T,Ghat@uh[t]))
    eta=np.zeros(len(T))
    for e,tr in edge_tris(T).items():
        if len(tr)!=2: continue
        a,b=tuple(e); ev=V[b]-V[a]; L=np.hypot(*ev); n=np.array([ev[1],-ev[0]])/L
        j=(grads[tr[0]]-grads[tr[1]])@n
        for ti in tr: eta[ti]+=np.sqrt(areas[ti])*L*j*j
    return np.sqrt(eta)

def energy_error(V,T,uh):
    e2=0.
    for t in T:
        P=[V[i] for i in t]; B=np.array([P[1]-P[0],P[2]-P[0]]).T
        area=abs(np.linalg.det(B))/2; gU=np.linalg.solve(B.T,Ghat@uh[t])
        for m in [(P[0]+P[1])/2,(P[1]+P[2])/2,(P[2]+P[0])/2]:
            d=gU-grad_u(m); e2+=(d@d)/3*area
    return np.sqrt(e2)

def dorfler(eta,th):
    idx=np.argsort(eta**2)[::-1]; tot=(eta**2).sum(); s=0.; M=[]
    for i in idx:
        s+=eta[i]**2; M.append(int(i))
        if s>=th*th*tot: break
    return M

V,T=initial_mesh(); rows=[]
for it in range(70):
    uh=solve(V,T); rows.append((len(V),energy_error(V,T,uh)))
    if len(V)>45000: break
    V,T=refine(V,T,dorfler(estimator(V,T,uh),0.5))
print(f"{'Np':>7} {'H1err':>12} {'rate_Np':>8}")
for k,(Np,err) in enumerate(rows):
    r = np.log(rows[k-1][1]/err)/np.log(Np/rows[k-1][0]) if k>0 else float('nan')
    print(f"{Np:>7} {err:>12.4e} {r:>8.3f}")
