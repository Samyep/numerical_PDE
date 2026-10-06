import math, time
import numpy as np
from scipy.special import roots_laguerre


def seeded_rng(seed, path):
    vals=[int(seed)]+[int(v)&0xffffffff for v in path]
    return np.random.default_rng(np.random.SeedSequence(vals))


class HJBRosenbrock:
    def __init__(self,d=20,T=1.0,seed_A=0):
        self.d=d; self.T=T; self.sigma=math.sqrt(2.0)
        rng=np.random.default_rng(seed_A)
        self.c1=rng.uniform(0.5,1.5,d-1)
        self.c2=rng.uniform(0.5,1.5,d-1)
        A=np.zeros((d,d))
        for i in range(d-1):
            c=self.c1[i]
            A[i,i]+=c; A[i+1,i+1]+=c; A[i,i+1]-=c; A[i+1,i]-=c
            A[i+1,i+1]+=self.c2[i]
        self.A=A
        self.evals,self.evecs=np.linalg.eigh(A)
        self.lmax=float(self.evals.max()); self.trA=float(np.trace(A))
    def q(self,x): return float(x@self.A@x)
    def g(self,x): return math.log((1.0+self.q(x))/2.0)
    def f(self,t,x,v):
        z=np.asarray(v[1:]); return -0.5*float(z@z)
    def sample(self,t,s,x,rng):
        dt=s-t
        W=rng.normal(0.0,math.sqrt(dt),size=self.d)
        X=x+self.sigma*W
        dI=np.concatenate(([1.0],W/dt))
        return X,dI
    def z_radius(self,t,factor=1.0):
        return factor*math.sqrt(2.0*self.lmax)
    def u_bounds(self,t,x):
        tau=self.T-t
        lo=-math.log(2.0)
        hi=math.log((1.0+self.q(x)+2.0*tau*self.trA)/2.0)
        return lo,hi
    def exact(self,t,x,nquad=96):
        tau=self.T-t
        y=x@self.evecs
        s,w=roots_laguerre(nquad)
        den=1.0+4.0*tau*s[:,None]*self.evals[None,:]
        logdet=-0.5*np.sum(np.log(den),axis=1)
        expo=-s*np.sum((self.evals[None,:]*y[None,:]**2)/den,axis=1)
        integ=float(np.sum(w*np.exp(logdet+expo)))
        return -math.log(2.0*integ)


class FullHistoryMLP:
    def __init__(self,eq,M=3,alpha=0.5,mode='raw',radius_factor=1.0,seed=0,heuristic_clip=10.0):
        self.eq=eq; self.M=M; self.alpha=alpha; self.mode=mode; self.radius_factor=radius_factor; self.seed=seed; self.heuristic_clip=heuristic_clip
        self.g_evals=0; self.f_evals=0; self.proj_count=0; self.preproj_viol=0
    def project(self,v,t,x):
        v=np.asarray(v,dtype=float).copy(); self.proj_count+=1
        if self.mode=='raw': return v
        if self.mode=='heuristic': return np.clip(v,-self.heuristic_clip,self.heuristic_clip)
        y=float(v[0]); z=v[1:].copy()
        # certified u interval
        if self.mode in ('hard','hard_u','hard_box'):
            lo,hi=self.eq.u_bounds(t,x); y=min(max(y,lo),hi)
        # certified/relaxed gradient region
        if self.mode in ('hard','hard_z'):
            rad=self.eq.z_radius(t,self.radius_factor); nrm=float(np.linalg.norm(z))
            if nrm>rad:
                self.preproj_viol+=1; z*=rad/nrm
        elif self.mode=='hard_box':
            # coordinatewise necessary bound |z_i| <= radius; weaker than joint l2 ball
            rad=self.eq.z_radius(t,self.radius_factor)
            if np.any(np.abs(z)>rad): self.preproj_viol+=1
            z=np.clip(z,-rad,rad)
        return np.concatenate(([y],z))
    def solve(self,n,t,x,path=(0,)):
        d=self.eq.d
        if n==0: return np.zeros(d+1)
        M=self.M; T=self.eq.T; a=self.alpha
        gx=self.eq.g(x); self.g_evals+=1
        rhs_g=np.zeros(d+1)
        for i in range(M**n):
            rng=seeded_rng(self.seed,path+(11,n,i))
            XT,dI=self.eq.sample(t,T,x,rng)
            gT=self.eq.g(XT); self.g_evals+=1
            rhs_g+=(gT-gx)*dI
        rhs=np.concatenate(([gx],np.zeros(d)))+rhs_g/(M**n)
        for level in range(n):
            N=M**(n-level); rhs_f=np.zeros(d+1)
            for i in range(N):
                rng=seeded_rng(self.seed,path+(21,n,level,i))
                r=rng.power(a); R=t+(T-t)*r
                XR,dI=self.eq.sample(t,R,x,rng)
                vl=self.solve(level,R,XR,path+(100+level,i,1))
                fl=self.eq.f(R,XR,vl); self.f_evals+=1
                if level==0: fp=0.0
                else:
                    vp=self.solve(level-1,R,XR,path+(100+level,i,2))
                    fp=self.eq.f(R,XR,vp); self.f_evals+=1
                rhs_f+=(r**(1.0-a))*(fl-fp)*dI
            rhs+=(T-t)*rhs_f/(a*N)
        return self.project(rhs,t,x)


def sample_ball(rng,n,d,r=1.0):
    z=rng.normal(size=(n,d)); z/=np.linalg.norm(z,axis=1,keepdims=True)
    rr=rng.random(n)**(1.0/d)*r
    return z*rr[:,None]


def run_config(d=20,n=3,M=3,seeds=10,npts=16,modes=('raw','heuristic','hard','hard_box'),factors=(1.0,),seed_points=2026):
    eq=HJBRosenbrock(d=d)
    rng=np.random.default_rng(seed_points+d)
    xs=sample_ball(rng,npts,d,1.0); ts=rng.uniform(0.0,0.9,size=npts)
    truth=np.array([eq.exact(float(t),x) for t,x in zip(ts,xs)])
    rows=[]
    for mode in modes:
        use_factors=factors if mode=='hard' else (1.0,)
        for fac in use_factors:
            pred=[]; viol=[]; fevals=[]; gevals=[]; t0=time.time()
            for seed in range(seeds):
                vals=[]
                for j,(t,x) in enumerate(zip(ts,xs)):
                    sol=FullHistoryMLP(eq,M=M,mode=mode,radius_factor=fac,seed=seed)
                    out=sol.solve(n,float(t),x,path=(j,))
                    vals.append(out[0]); viol.append(sol.preproj_viol/max(sol.proj_count,1)); fevals.append(sol.f_evals); gevals.append(sol.g_evals)
                pred.append(vals)
            pred=np.asarray(pred)
            err=pred-truth[None,:]
            rels=np.linalg.norm(err,axis=1)/np.linalg.norm(truth)
            rows.append(dict(d=d,n=n,M=M,mode=mode,factor=fac,seeds=seeds,npts=npts,
                             rel_l2_mean=float(rels.mean()),rel_l2_std=float(rels.std()),
                             mae=float(np.mean(np.abs(err))),p95=float(np.quantile(np.abs(err),.95)),
                             bias=float(np.mean(err)),run_std=float(np.mean(np.std(pred,axis=0))),
                             proj_violation=float(np.mean(viol)),f_evals=int(np.mean(fevals)),g_evals=int(np.mean(gevals)),
                             seconds=time.time()-t0))
    return rows

if __name__=='__main__':
    import json
    allrows=[]
    for d in (20,40):
        for n,M in ((2,3),(2,6),(3,2),(3,3),(3,4)):
            allrows += run_config(d=d,n=n,M=M,seeds=10,npts=16,modes=('raw','heuristic','hard','hard_box'),factors=(1.0,2.0,4.0))
    print(json.dumps(allrows,indent=2))
