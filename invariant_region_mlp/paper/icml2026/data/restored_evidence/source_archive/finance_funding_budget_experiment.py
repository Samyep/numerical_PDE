import math, time, json
import numpy as np


def seeded_rng(seed,path):
    vals=[int(seed)]+[int(v)&0xffffffff for v in path]
    return np.random.default_rng(np.random.SeedSequence(vals))

class FundingMaxSpread:
    def __init__(self,d=100,T=0.5,sigma=0.2,mu=0.06,Rl=0.04,Rb=0.06):
        self.d=d; self.T=T; self.sigma=sigma; self.mu=mu; self.Rl=Rl; self.Rb=Rb
    def g(self,x):
        m=float(np.max(x)); return max(m-120.,0.)-2.*max(m-150.,0.)
    def f(self,t,x,v):
        y=float(v[0]); z=np.asarray(v[1:]); sz=float(np.sum(z))
        return -self.Rl*y - ((self.mu-self.Rl)/self.sigma)*sz + (self.Rb-self.Rl)*max(sz/self.sigma-y,0.)
    def sample(self,t,s,x,rng):
        dt=s-t; W=rng.normal(0.,math.sqrt(dt),size=self.d)
        X=x*np.exp((self.mu-.5*self.sigma**2)*dt+self.sigma*W)
        dI=np.concatenate(([1.],W/dt)); return X,dI
    def delta_radius(self,t,factor=1.0): return factor*math.exp(.5*self.sigma**2*(self.T-t))

class FullHistoryMLP:
    def __init__(self,eq,M=3,alpha=.5,mode='baseline',radius_factor=1.,seed=0,heuristic_delta_clip=1.0):
        self.eq=eq; self.M=M; self.alpha=alpha; self.mode=mode; self.radius_factor=radius_factor; self.seed=seed; self.heuristic_delta_clip=heuristic_delta_clip
        self.g_evals=0; self.f_evals=0; self.proj_count=0; self.preproj_viol=0
    def project(self,v,t,x):
        self.proj_count+=1; v=np.asarray(v,dtype=float).copy()
        if self.mode=='baseline': return v
        y=float(v[0]); z=v[1:].copy(); den=self.eq.sigma*np.maximum(x,1e-12); delta=z/den
        if self.mode=='hard':
            rad=self.eq.delta_radius(t,self.radius_factor); nr=float(np.linalg.norm(delta))
            if nr>rad: self.preproj_viol+=1; delta*=rad/nr
        elif self.mode=='box':
            rad=self.eq.delta_radius(t,self.radius_factor)
            if np.any(np.abs(delta)>rad): self.preproj_viol+=1
            delta=np.clip(delta,-rad,rad)
        elif self.mode=='heuristic':
            c=self.heuristic_delta_clip
            if np.any(np.abs(delta)>c): self.preproj_viol+=1
            delta=np.clip(delta,-c,c)
        return np.concatenate(([y],den*delta))
    def solve(self,n,t,x,path=(0,)):
        d=self.eq.d
        if n==0: return np.zeros(d+1)
        M=self.M; T=self.eq.T; a=self.alpha; gx=self.eq.g(x); self.g_evals+=1
        rg=np.zeros(d+1)
        for i in range(M**n):
            rng=seeded_rng(self.seed,path+(11,n,i)); XT,dI=self.eq.sample(t,T,x,rng); gt=self.eq.g(XT); self.g_evals+=1; rg+=(gt-gx)*dI
        rhs=np.concatenate(([gx],np.zeros(d)))+rg/(M**n)
        for level in range(n):
            N=M**(n-level); rf=np.zeros(d+1)
            for i in range(N):
                rng=seeded_rng(self.seed,path+(21,n,level,i)); r=rng.power(a); R=t+(T-t)*r; XR,dI=self.eq.sample(t,R,x,rng)
                vl=self.solve(level,R,XR,path+(100+level,i,1)); fl=self.eq.f(R,XR,vl); self.f_evals+=1
                if level==0: fp=0.
                else:
                    vp=self.solve(level-1,R,XR,path+(100+level,i,2)); fp=self.eq.f(R,XR,vp); self.f_evals+=1
                rf+=(r**(1-a))*(fl-fp)*dI
            rhs+=(T-t)*rf/(a*N)
        return self.project(rhs,t,x)

def run_config(n=3,M=3,seeds=100,d=100,reference=21.299,modes=('baseline','hard'),factors=(1.,2.,4.)):
    eq=FundingMaxSpread(d=d); x0=np.full(d,100.); rows=[]
    for mode in modes:
        usef=factors if mode=='hard' else (1.,)
        for fac in usef:
            vals=[]; norms=[]; viol=[]; fe=[]; ge=[]; t0=time.time()
            for seed in range(seeds):
                s=FullHistoryMLP(eq,M=M,mode=mode,radius_factor=fac,seed=seed)
                out=s.solve(n,0.,x0)
                vals.append(out[0]); norms.append(np.linalg.norm(out[1:]/(eq.sigma*x0))); viol.append(s.preproj_viol/max(s.proj_count,1)); fe.append(s.f_evals); ge.append(s.g_evals)
            vals=np.array(vals); ae=np.abs(vals-reference)
            rows.append(dict(n=n,M=M,mode=mode,factor=fac,seeds=seeds,mean=float(vals.mean()),mae=float(ae.mean()),median_ae=float(np.median(ae)),p90_ae=float(np.quantile(ae,.9)),std=float(vals.std()),bias=float(np.mean(vals-reference)),delta_med=float(np.median(norms)),delta_p90=float(np.quantile(norms,.9)),proj_violation=float(np.mean(viol)),f_evals=int(np.mean(fe)),g_evals=int(np.mean(ge)),seconds=time.time()-t0))
    return rows

if __name__=='__main__':
    rows=[]
    for n,M in ((2,2),(2,3),(2,4),(2,6),(2,10),(3,2),(3,3),(3,4),(3,5),(4,2),(4,3)):
        print('RUN',n,M,flush=True)
        rr=run_config(n=n,M=M,seeds=100,modes=('baseline','hard','box'),factors=(1.,2.,4.))
        rows+=rr
        for r in rr: print(r,flush=True)
    open('/mnt/data/finance_funding_budget_results.json','w').write(json.dumps(rows,indent=2))
