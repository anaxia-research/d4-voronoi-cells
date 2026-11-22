# Cutting-plane grouping: LP over mixtures of the 37 certificates with sample rows + cut rows + the exact check's margins.
import numpy as np, json, os, glob, time, sys
from scipy.optimize import linprog
from lp import PK,TK,N2,N3
MG2,MG3,MGM=2e-6,5e-6,1e-6
BASE=os.environ.get('SAMP','samples.npz')
Q='cutq'; os.makedirs(Q,exist_ok=True)
class Pool:
    def __init__(s): s.files=[]; S=np.load(BASE); s.P,s.pk,s.T,s.tk=S['P'],S['pk']-1,S['T'],S['tk']-1; s.PT,s.TT,s.MM=S['PT'],S['TT'],S['MM']; s.nc=len(s.TT); s.refresh()
    def refresh(s):
        fs=sorted(glob.glob(Q+'/*.npz'))+['cuts_h1.npz']
        for f in fs:
            if f in s.files or not os.path.exists(f): continue
            C=np.load(f); s.P=np.vstack([s.P,C['P']]); s.pk=np.concatenate([s.pk,C['pk']-1]); s.T=np.vstack([s.T,C['T']]); s.tk=np.concatenate([s.tk,C['tk']-1]); s.files.append(f)
pool=None
def W(vecs):
    W2=np.array([[N2(n,*k) for k in PK] for n in vecs],float); W3=np.array([[N3(n,k) for k in TK] for n in vecs],float); return W2,W3
def lin(n): return pool.TT/2+(pool.MM+pool.PT)@np.array(n,float)
def solve_set(vecs):
    P,pk,T,tk,nc=pool.P,pool.pk,pool.T,pool.tk,pool.nc
    W2,W3=W(vecs)
    a2=[i for i in range(len(PK)) if W2[:,i].max()>0]; a3=[i for i in range(len(TK)) if W3[:,i].max()>0]
    m2={k:j for j,k in enumerate(a2)}; m3={k:j for j,k in enumerate(a3)}
    rp=np.isin(pk,a2); rt=np.isin(tk,a3); n2,n3=len(a2),len(a3)
    Ap=np.hstack([P[rp], -np.eye(n2)[[m2[k] for k in pk[rp]]], np.zeros((rp.sum(),n3+1))])
    At=np.hstack([T[rt], np.zeros((rt.sum(),n2)), -np.eye(n3)[[m3[k] for k in tk[rt]]], np.zeros((rt.sum(),1))])
    Av=np.hstack([np.array([lin(n) for n in vecs]), W2[:,a2], W3[:,a3], -np.ones((len(vecs),1))])
    const=np.array([W2[i].sum()*MG2+W3[i].sum()*MG3+sum(n)*MGM for i,n in enumerate(vecs)])
    A=np.vstack([Ap,At,Av]); b=np.concatenate([np.zeros(Ap.shape[0]+At.shape[0]), -const])
    c=np.zeros(nc+n2+n3+1); c[-1]=1
    Aeq=np.concatenate([np.ones(nc),np.zeros(n2+n3+1)])[None]
    r=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=[1],bounds=[(0,None)]*nc+[(None,None)]*(n2+n3+1),method='highs')
    if r.status!=0: return 9.9,None,None
    vals=Av[:,:-1]@r.x[:-1]+const
    return r.x[-1], r.x[:nc], vals
def evaluate(lam, vecs):
    pv=pool.P@lam; tv=pool.T@lam
    c2=np.full(len(PK),-1e9); np.maximum.at(c2,pool.pk,pv); c3=np.full(len(TK),-1e9); np.maximum.at(c3,pool.tk,tv)
    W2,W3=W(vecs); out=[]
    for i,n in enumerate(vecs):
        m2=W2[i]>0; m3=W3[i]>0
        out.append(lin(n)@lam + W2[i][m2]@(c2[m2]+MG2) + W3[i][m3]@(c3[m3]+MG3) + sum(n)*MGM)
    return np.array(out)
def cuts(lam, tag):
    f=f'{Q}/{tag}.json'; json.dump(list(map(float,lam)),open(f,'w'))
    while not os.path.exists(f'{Q}/{tag}.npz'): time.sleep(5)
    pool.refresh()
