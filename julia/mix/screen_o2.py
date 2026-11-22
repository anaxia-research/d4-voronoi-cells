import json
from lp import solve
O=[]
for n1 in range(17,25):
  for n7 in range(4):
    for n6 in range(28):
      for n5 in range(28):
        for n4 in range(28):
          for n3 in range(28):
            n2=27-n1-n3-n4-n5-n6-n7
            if n2<0: continue
            n=(n1,n2,n3,n4,n5,n6,n7)
            if n1+n2<14 or n1+n2+n3<20 or n1+n2+n3+n4<21 or n1+n2+n3+n4+n5<22 or 27-n7<24: continue
            O.append(n)
print('O2 vectors',len(O),flush=True)
json.dump(O,open('o2_vectors.json','w'))
from collections import defaultdict
best=defaultdict(lambda:(9,None)); res=[]
for i,v in enumerate(O):
    val,_=solve(v); res.append((val,v))
    if val<best[v[0]][0]: best[v[0]]=(val,v)
    if i%500==0: print(i,flush=True)
res.sort(); json.dump([[a,list(b)] for a,b in res],open('screen_o2.json','w'))
for k in sorted(best): print(k,best[k])
print('below 3.075:',sum(a<3.075 for a,_ in res),'of',len(res))
