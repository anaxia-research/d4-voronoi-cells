import numpy as np, json, sys
from math import comb
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack, identity
S=np.load(sys.argv[1] if len(sys.argv)>1 else 'samples.npz')
P,pk,T,tk=S['P'],S['pk']-1,S['T'],S['tk']-1
PT,TT,MM=S['PT'],S['TT'],S['MM']; nc=len(TT)
PK=[(a,b) for a in range(7) for b in range(a,7)]
TK=[(a,b,c) for a in range(7) for b in range(a,7) for c in range(b,7)]
def N2(n,a,b): return comb(n[a],2) if a==b else n[a]*n[b]
def N3(n,k):
    from collections import Counter
    r=1
    for b,m in Counter(k).items(): r*=comb(n[b],m)
    return r
def lin_cost(n): return TT/2 + (MM+PT)@np.array(n,float)
def solve(n, cols=None):
    cols=list(range(nc)) if cols is None else cols
    w2=np.array([N2(n,*k) for k in PK],float); w3=np.array([N3(n,k) for k in TK],float)
    a2=[i for i in range(len(PK)) if w2[i]>0]; a3=[i for i in range(len(TK)) if w3[i]>0]
    m2={k:j for j,k in enumerate(a2)}; m3={k:j for j,k in enumerate(a3)}
    rp=np.isin(pk,a2); rt=np.isin(tk,a3)
    nv=len(cols)+len(a2)+len(a3)
    c=np.concatenate([lin_cost(n)[cols], w2[a2], w3[a3]])
    Ap=np.hstack([P[rp][:,cols], -np.eye(len(a2))[[m2[k] for k in pk[rp]]], np.zeros((rp.sum(),len(a3)))])
    At=np.hstack([T[rt][:,cols], np.zeros((rt.sum(),len(a2))), -np.eye(len(a3))[[m3[k] for k in tk[rt]]]])
    A=np.vstack([Ap,At]); b=np.zeros(A.shape[0])
    Aeq=np.concatenate([np.ones(len(cols)),np.zeros(len(a2)+len(a3))])[None]
    bnd=[(0,None)]*len(cols)+[(None,None)]*(len(a2)+len(a3))
    r=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=[1],bounds=bnd,method='highs')
    lam=np.zeros(nc); lam[cols]=r.x[:len(cols)]
    return r.fun, lam
if __name__=='__main__':
    names=[l.strip() for l in open('certs.txt') if l.strip()]
    for v,nm in [((16,0,9,0,0,0,2),'combo27_k16f2v0902'),((16,8,2,0,0,1,0),'combo27_k16f2v8210')]:
        i=names.index(nm); val,_=solve(v,[i]); print(nm,'alone on its vector',val)
        val,lam=solve(v); print(' mixed', val, {names[j]:round(lam[j],3) for j in np.argsort(-lam)[:4]})
