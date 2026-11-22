import json, numpy as np, sys, os
sys.argv=[sys.argv[0], 'samples.npz']
import lp
from group import solve_set
TARGET=float(os.environ.get('TARGET','3.072'))
res=json.load(open('screen.json'))
g1=set(tuple(v) for v in json.load(open('g1_vecs.json')))|{(15,0,7,0,5,0,0)}
ind={tuple(v):a for a,v in res}
left=sorted([v for v in ind if v not in g1 and ind[v]<TARGET], key=lambda v: ind[v])
out='groups2.json'
groups=json.load(open(out)) if os.path.exists(out) else []
done=set(tuple(v) for g in groups for v in g['vecs'])
left=[v for v in left if v not in done]
while left:
    G=[left[0]]; s,lam,_=solve_set(G); rest=[]
    for v in left[1:]:
        if v[0]!=G[0][0] or len(G)>=400: rest.append(v); continue
        s2,lam2,_=solve_set(G+[v])
        if s2<=TARGET: G.append(v); s,lam=s2,lam2
        else: rest.append(v)
    groups.append({'vecs':[list(v) for v in G],'lp':s,'lam':list(lam)}); json.dump(groups,open(out,'w'))
    print('group',len(groups),'size',len(G),'lp',round(s,5),'left',len(rest),flush=True)
    left=rest
