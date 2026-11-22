# verify_group.py VECTORS.json LOG [LOG ...]
# Checks that the logs of the split exact check of one group prove it: every run used the same
# certificate thresholds (identical pair, triple and point-term lines and the same bound), the
# positivity, pair and bin inequalities passed in some run, and every triple kind that the
# count vectors use passed in some run. Prints PROVED or what is missing.
import sys, re, json, itertools
from math import comb
from collections import Counter
T=['A1','A2','B3','B4','B5','B6','F7']
def N3(n,k):
    r=1
    for b,m in Counter(k).items(): r*=comb(n[b],m)
    return r
V=json.load(open(sys.argv[1])); logs=sys.argv[2:]
need=set(''.join(T[i] for i in k) for k in itertools.combinations_with_replacement(range(7),3) if any(N3(v,k)>0 for v in V))
setup=None; passed=set(); pos=pairs=bins=False; bound=None; ok=True
for f in logs:
    s=open(f).read()
    if 'FAIL' in s: print('FAIL line in',f); ok=False
    st=[l for l in s.splitlines() if re.match(r'\s+(pair|triple) \w+: float largest',l) or l.startswith('  point terms')]
    if setup is None: setup=st
    elif st!=setup: print('thresholds differ in',f); ok=False
    b=re.findall(r'U\(Y\) <= ([0-9.]+)',s)
    if b:
        bound=bound or b[-1]
        if b[-1]!=bound: print('bound differs in',f); ok=False
    pos|=('PASS A_1..A_D' in s and 'PASS three-point blocks positive definite' in s)
    if 'PASS f <= m_b on every bin' in s:
        bins=True; pairs=pairs or all(f'PASS K + PAIR3_{a}{c} ' in s for a,c in re.findall(r'  pair (\w\d)(\w\d): float',s))
    passed|=set(re.findall(r'PASS TRIPLE3_(\w+) ',s))
miss=sorted(need-passed)
print(f'{len(V)} count vectors; {len(need)} triple kinds needed, {len(need&passed)} passed; positivity {pos}; pairs {pairs}; bins {bins}; bound {bound}')
if miss: print('missing kinds:',' '.join(miss)); ok=False
if not (pos and pairs and bins): ok=False
print('PROVED' if ok and bound and float(bound)<3.1033049 else 'NOT YET')
