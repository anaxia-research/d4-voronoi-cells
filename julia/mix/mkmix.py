import numpy as np, json, sys
import os; R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','combo')+os.sep
names=[l.strip() for l in open('certs.txt') if l.strip()]
cases=json.load(open('certcases.json')); mt=json.load(open('certmt.json'))
def tmap(nm):
    J=json.load(open(R+'cases/'+cases[nm])); T=J['types']
    return [T.index(b[0]) for b in J['bins']]
def make(lam, out, tol=1e-9):
    A=0; z=0; t=0; m=0; X=[]; G=[]
    for i,l in enumerate(lam):
        if l<tol: continue
        Z=np.load(R+names[i]+'_d10.npz')
        A=A+l*Z['A']; z=z+l*Z['z']; t=t+l*float(np.ravel(Z['t'])[0]); m=m+l*np.ravel(Z['m'])
        if mt[names[i]]:
            for g in range(Z['gtyp'].shape[0]): X.append(l*Z['x3g'][g]); G.append(Z['gtyp'][g])
        else:
            X.append(l*np.ravel(Z['x3'])); G.append(np.array(tmap(names[i])))
    np.savez(out, A=A, z=z, t=np.array([t]), m=m, x3g=np.array(X), gtyp=np.array(G,dtype=np.int64))
    return len(X)
if __name__=='__main__':
    lam=np.array(json.loads(open(sys.argv[1]).read())); print(make(lam, sys.argv[2]))
