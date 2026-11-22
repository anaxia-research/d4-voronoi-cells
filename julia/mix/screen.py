import json,numpy as np,sys
from lp import solve
O=[tuple(v) for v in json.load(open('../o1_vectors.json'))]
covered={(16,5,2,0,0,3,1),(16,4,2,0,0,4,1),(16,7,2,0,0,2,0),(16,4,2,0,0,3,2),(16,6,2,0,0,3,0),(16,4,2,0,0,5,0),(16,5,2,0,0,4,0),(16,6,2,0,0,2,1),(16,5,2,0,0,2,2),(16,8,2,0,0,1,0)}
O=[v for v in O if v not in covered]
res=[]
for i,v in enumerate(O):
    val,lam=solve(v); res.append((val,v))
res.sort()
json.dump([[a,list(b)] for a,b in res],open('screen.json','w'))
for a,b in res[:25]: print(round(a,5),b)
print('count <3.08:',sum(a<3.08 for a,_ in res),' <3.09:',sum(a<3.09 for a,_ in res), ' median', res[len(res)//2][0])
