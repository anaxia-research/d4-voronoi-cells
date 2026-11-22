import numpy as np, json, sys, os
from lp import PK,TK,N2,N3
def load(base, cutfiles):
    S=dict(np.load(base)); P,pk,T,tk=S['P'],S['pk']-1,S['T'],S['tk']-1
    for f in cutfiles:
        Cf=np.load(f); P=np.vstack([P,Cf['P']]); pk=np.concatenate([pk,Cf['pk']-1]); T=np.vstack([T,Cf['T']]); tk=np.concatenate([tk,Cf['tk']-1])
    return P,pk,T,tk,S['PT'],S['TT'],S['MM']
MG2,MG3,MGM=2e-6,5e-6,1e-6
def ev(D, lam, vecs):
    P,pk,T,tk,PT,TT,MM=D; lam=np.array(lam); pv=P@lam; tv=T@lam
    c2=np.full(len(PK),-1e9); np.maximum.at(c2,pk,pv); c3=np.full(len(TK),-1e9); np.maximum.at(c3,tk,tv)
    out=[]
    for n in vecs:
        w2=np.array([N2(n,*k) for k in PK]); w3=np.array([N3(n,k) for k in TK])
        out.append(TT@lam/2+(MM+PT)@np.array(n,float)@lam + w2[w2>0]@(c2[w2>0]+MG2) + w3[w3>0]@(c3[w3>0]+MG3) + sum(n)*MGM)
    return max(out)
if __name__=='__main__':
    cf=[f for f in sys.argv[2:] if f.endswith('.npz')]; gs=[g for g in sys.argv[2:] if not g.endswith('.npz')]
    D=load(sys.argv[1], cf)
    for g in gs: print(g, ev(D, json.load(open(g+'_lam.json')), json.load(open(g+'_vecs.json'))))
