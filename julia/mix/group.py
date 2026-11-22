import numpy as np, json, sys
from math import comb
from scipy.optimize import linprog
import lp
from lp import PK,TK,N2,N3,lin_cost
def solve_set(vecs, S=None):
    P,pk,T,tk,nc=lp.P,lp.pk,lp.T,lp.tk,lp.nc
    W2=np.array([[N2(n,*k) for k in PK] for n in vecs],float); W3=np.array([[N3(n,k) for k in TK] for n in vecs],float)
    a2=[i for i in range(len(PK)) if W2[:,i].max()>0]; a3=[i for i in range(len(TK)) if W3[:,i].max()>0]
    m2={k:j for j,k in enumerate(a2)}; m3={k:j for j,k in enumerate(a3)}
    rp=np.isin(pk,a2); rt=np.isin(tk,a3); n2,n3=len(a2),len(a3)
    nv=nc+n2+n3+1
    rows=[]
    Ap=np.hstack([P[rp], -np.eye(n2)[[m2[k] for k in pk[rp]]], np.zeros((rp.sum(),n3+1))])
    At=np.hstack([T[rt], np.zeros((rt.sum(),n2)), -np.eye(n3)[[m3[k] for k in tk[rt]]], np.zeros((rt.sum(),1))])
    Av=np.hstack([np.array([lin_cost(n) for n in vecs]), W2[:,a2], W3[:,a3], -np.ones((len(vecs),1))])
    A=np.vstack([Ap,At,Av]); b=np.zeros(A.shape[0])
    c=np.zeros(nv); c[-1]=1
    Aeq=np.concatenate([np.ones(nc),np.zeros(n2+n3+1)])[None]
    bnd=[(0,None)]*nc+[(None,None)]*(n2+n3+1)
    r=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=[1],bounds=bnd,method='highs')
    vals=Av[:,:-1]@r.x[:-1]
    return r.x[-1], r.x[:nc], vals
if __name__=='__main__':
    res=json.load(open('screen.json'))
    from collections import defaultdict
    fam=defaultdict(list)
    for val,v in res:
        v=tuple(v); fam[(v[0],v[6],min(v[5],3))].append(v)
    out={}
    for k in sorted(fam):
        s,lam,vals=solve_set(fam[k]); out[str(k)]=(s,len(fam[k]))
        print(k,len(fam[k]),round(s,5),flush=True)
    json.dump(out,open('families.json','w'))
