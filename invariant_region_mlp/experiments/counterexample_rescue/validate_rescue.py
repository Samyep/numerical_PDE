"""Independent checks of Eq. (6), exact moments, and active-subspace reduction."""
import json, math
from pathlib import Path
import numpy as np
from counterexample_rescue import Solver

class Literal:
    def __init__(self,d,m,seed):
        self.d,self.m=d,m
        self.times,self.active,self.inactive=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(3)]
    def normal(self,b,k):
        z=np.empty((b,k,self.d)); z[:,:,0]=self.active.standard_normal((b,k))
        if self.d>1: z[:,:,1:]=self.inactive.standard_normal((b,k,self.d-1))
        return z
    def solve(self,n,h,x):
        b=len(h)
        if n==0: return np.zeros((b,self.d+1))
        k=self.m**n; w=np.sqrt(h)[:,None,None]*self.normal(b,k)
        diff=np.abs(x[:,None]+w[:,:,0])-np.abs(x[:,None])
        u=np.column_stack((np.abs(x),np.zeros((b,self.d))))
        weights=np.concatenate((np.ones((b,k,1)),w/h[:,None,None]),axis=2)
        u+=np.mean(diff[:,:,None]*weights,axis=1)
        for l in range(1,n):
            k=self.m**(n-l); r=self.times.random((b,k))**2
            dt=h[:,None]*r; w=np.sqrt(dt)[:,:,None]*self.normal(b,k)
            hp=h[:,None]*(1-r); xp=x[:,None]+w[:,:,0]
            a=self.solve(l,hp.ravel(),xp.ravel())
            df=np.linalg.norm(a[:,2:],axis=1)
            if l>1:
                a=self.solve(l-1,hp.ravel(),xp.ravel()); df-=np.linalg.norm(a[:,2:],axis=1)
            df=df.reshape(b,k); rho=1/(2*np.sqrt(h)[:,None]*np.sqrt(dt))
            weights=np.concatenate((np.ones((b,k,1)),w/dt[:,:,None]),axis=2)
            u+=np.mean((df/rho)[:,:,None]*weights,axis=1)
        return u

class Ranked(Solver):
    def __init__(self,d,m,seed,rank):
        super().__init__(d,m,seed,modes=('raw','recursive_subspace'))
        self.rank=rank
        self.extra=np.random.default_rng(np.random.SeedSequence([seed,984725]))
    def normals(self,b,k):
        z=np.empty((b,k,self.d)); z[:,:,0]=self.active.standard_normal((b,k))
        if self.rank>1: z[:,:,1:self.rank]=self.inactive.standard_normal((b,k,self.rank-1))
        if self.d>self.rank: z[:,:,self.rank:]=self.extra.standard_normal((b,k,self.d-self.rank))
        return z

def main():
    records=[]
    for d in [1,2,7]:
        for n in [1,2,3]:
            h=np.array([1.,.75,.01]); x=np.array([0.,.3,-1.])
            expected=Literal(d,2,91+n).solve(n,h,x)
            actual=Solver(d,2,91+n,modes=('raw',)).solve(n,h,x)[:,0]
            diff=float(np.max(np.abs(actual-expected)))
            assert np.allclose(actual,expected,rtol=2e-12,atol=2e-12)
            records.append({'d':d,'n':n,'max_abs_diff':diff,'check':'literal Eq.6'})
    moments=[]; rng=np.random.default_rng(7248532); R=65536
    for N in [1,4,27,256]:
        errors=[]
        for _ in range(0,R,256):
            z=rng.normal(size=(256,N)); u=np.abs(z).mean(1); v=(np.abs(z)*z).mean(1)
            errors.append((u-math.sqrt(2/math.pi))**2+v**2)
        e=np.concatenate(errors); emp=float(e.mean()); se=float(e.std(ddof=1)/math.sqrt(R)); theory=(4-2/math.pi)/N
        moments.append({'N':N,'replicates':R,'empirical_mse':emp,'exact_mse':theory,'mean_se':se,'ratio':emp/theory})
    print(json.dumps({'checks':records,'moments':moments},indent=2))
if __name__=='__main__': main()
