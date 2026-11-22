import json, os, sys, numpy as np, subprocess
import cg
cg.pool=cg.Pool()
TARGET=float(os.environ.get('TARGET','3.1005'))
res=json.load(open('screen.json'))
skip=set(tuple(v) for f in ['g1_vecs.json'] for v in json.load(open(f)))|{(15,0,7,0,5,0,0)}
cand=[tuple(v) for a,v in sorted(res) if tuple(v) not in skip]+[tuple(v) for f in ['g2_vecs.json','g4_vecs.json'] for v in json.load(open(f))]; cand=list(dict.fromkeys(cand))
# 1. individual screen with cuts and margins
sf='cscreen.json'
sc=dict((tuple(k),v) for k,v in json.load(open(sf))) if os.path.exists(sf) else {}
for i,v in enumerate(cand):
    if v in sc: continue
    s,lam,_=cg.solve_set([v]); sc[v]=s
    if i%100==0:
        json.dump([[list(k),x] for k,x in sc.items()],open(sf,'w')); print('screen',i,len(cand),round(s,5),flush=True)
json.dump([[list(k),x] for k,x in sc.items()],open(sf,'w'))
ok=sorted([v for v in cand if sc[v]<TARGET], key=lambda v: sc[v])
print('reachable individually', len(ok), 'of', len(cand), flush=True)
# 2. groups
gf='cgroups2.json'
groups=json.load(open(gf)) if os.path.exists(gf) else []
done=set(tuple(v) for g in groups for v in g['vecs'])
left=[v for v in ok if v not in done]
while left:
    G=[left[0]]; s,lam,_=cg.solve_set(G); rest=[]; print('building from',left[0],'with',len(left),'left',flush=True)
    for v in left[1:]:
        if len(G)>=int(os.environ.get('CAP','5000')): rest.append(v); continue
        if cg.evaluate(lam,[v]).max()<=TARGET-0.0015: G.append(v); continue
        s2,lam2,_=cg.solve_set(G+[v])
        if s2<=TARGET-0.0015: G.append(v); s,lam=s2,lam2
        else: rest.append(v)
    k=len(groups)+1
    for it in range(5):
        cg.cuts(lam, f'K{k}_{it}')
        s0=cg.evaluate(lam,G).max()
        s,lam,vals=cg.solve_set(G)
        while s>TARGET and len(G)>1:
            j=int(np.argmax(vals)); rest.append(G.pop(j)); s,lam,vals=cg.solve_set(G)
        print(f'  group {k} cut round {it}: old lam now {s0:.5f}, re-solved {s:.5f}, size {len(G)}',flush=True)
        if s0<=TARGET+2e-4: break
    if s>TARGET: print('  group',k,'dropped',flush=True); left=sorted(rest,key=lambda v: sc[v]); left=[v for v in left if v!=G[0]]; continue
    groups.append({'vecs':[list(v) for v in G],'lp':s,'lam':list(map(float,lam))}); json.dump(groups,open(gf,'w'))
    json.dump([list(v) for v in G],open(f'K{k}_vecs.json','w')); json.dump(list(map(float,lam)),open(f'K{k}_lam.json','w'))
    print('group',k,'size',len(G),'lp',round(s,5),'left',len(rest),flush=True)
    left=sorted(rest,key=lambda v: sc[v])
