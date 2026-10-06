import numpy as np, math, json
from pathlib import Path

REF=21.299
class FundingEq:
    d=100; T=.5; sigma=.2; mu=.06; Rl=.04; Rb=.06
    def g(self,x):
        m=np.max(x,axis=-1)
        return np.maximum(m-120.,0.)-2*np.maximum(m-150.,0.)
    def f(self,y,z):
        sz=np.sum(z,axis=-1)
        return -self.Rl*y-((self.mu-self.Rl)/self.sigma)*sz+(self.Rb-self.Rl)*np.maximum(sz/self.sigma-y,0.)
    def radius(self,t):
        return np.exp(.5*self.sigma**2*(self.T-t))

class Solver:
    def __init__(self,M,seed,mode="raw",factor=1.0,c=1.0,a=.5):
        self.M=M; self.eq=FundingEq(); self.mode=mode; self.factor=factor; self.c=c; self.a=a
        self.rng=np.random.default_rng(seed)
    def terminal(self,t,x,n):
        b=len(t); h=self.eq.T-t; k=self.M**n
        normal=self.rng.normal(size=(b,k,self.eq.d))
        W=np.sqrt(h)[:,None,None]*normal
        XT=x[:,None,:]*np.exp((self.eq.mu-.5*self.eq.sigma**2)*h[:,None,None]+self.eq.sigma*W)
        gx=self.eq.g(x); diff=self.eq.g(XT)-gx[:,None]
        out=np.empty((b,self.eq.d+1))
        out[:,0]=gx+diff.mean(1)
        out[:,1:]=np.mean(diff[:,:,None]*W/np.maximum(h[:,None,None],1e-14),axis=1)
        return out
    def transform(self,state,t,x,B,K):
        out=state.copy()
        if self.mode=="raw": return out
        if self.mode=="scale":
            out[:,1:]*=self.c; return out
        z=out[:,1:]; delta=z/(self.eq.sigma*np.maximum(x,1e-12))
        rad=self.factor*self.eq.radius(t); nr=np.linalg.norm(delta,axis=1)
        if self.mode=="sample":
            s=np.minimum(1.,rad/np.maximum(nr,1e-30))
            out[:,1:]=self.eq.sigma*x*(delta*s[:,None]); return out
        ratio=(rad/np.maximum(nr,1e-30)).reshape(B,K)
        alpha=np.minimum(1.,np.min(ratio,axis=1))
        out[:,1:]=self.eq.sigma*x*(delta*alpha.repeat(K)[:,None]); return out
    def solve(self,n,t,x):
        t=np.asarray(t,float).reshape(-1); x=np.asarray(x,float).reshape(len(t),self.eq.d)
        b=len(t)
        if n<=0: return np.zeros((b,self.eq.d+1))
        out=self.terminal(t,x,n); h=self.eq.T-t
        for level in range(1,n):
            K=self.M**(n-level)
            r=self.rng.power(self.a,size=(b,K)); dt=h[:,None]*r
            normal=self.rng.normal(size=(b,K,self.eq.d)); W=np.sqrt(dt)[:,:,None]*normal
            XR=x[:,None,:]*np.exp((self.eq.mu-.5*self.eq.sigma**2)*dt[:,:,None]+self.eq.sigma*W)
            R=t[:,None]+dt
            high=self.solve(level,R.ravel(),XR.reshape(-1,self.eq.d)).reshape(b,K,self.eq.d+1)
            st=self.transform(high.reshape(-1,self.eq.d+1),R.ravel(),XR.reshape(-1,self.eq.d),b,K)
            fd=self.eq.f(st[:,0],st[:,1:]).reshape(b,K)
            if level>1:
                low=self.solve(level-1,R.ravel(),XR.reshape(-1,self.eq.d)).reshape(b,K,self.eq.d+1)
                st=self.transform(low.reshape(-1,self.eq.d+1),R.ravel(),XR.reshape(-1,self.eq.d),b,K)
                fd-=self.eq.f(st[:,0],st[:,1:]).reshape(b,K)
            w=h[:,None]*(r**(1-self.a))/self.a
            out[:,0]+=np.mean(w*fd,axis=1)
            out[:,1:]+=np.mean(w[:,:,None]*fd[:,:,None]*W/np.maximum(dt[:,:,None],1e-14),axis=1)
        return out

def run(n,M,mode="raw",factor=1.0,c=1.0,reps=100,block=10,seed=20261005):
    vals=[]
    for st in range(0,reps,block):
        b=min(block,reps-st)
        s=Solver(M,seed+100000*n+1000*M+st,mode=mode,factor=factor,c=c)
        vals.extend(s.solve(n,np.zeros(b),np.full((b,100),100.))[:,0])
    vals=np.asarray(vals)
    return float(np.mean(np.abs(vals-REF))), vals.tolist()

def fzero(n,M,reps=100,block=10,seed=20261005):
    vals=[]
    for st in range(0,reps,block):
        b=min(block,reps-st)
        rng=np.random.default_rng(seed+100000*n+1000*M+st)
        s=Solver(M,0)
        s.rng=rng
        vals.extend(s.terminal(np.zeros(b),np.full((b,100),100.),n)[:,0])
    vals=np.asarray(vals)
    return float(np.mean(np.abs(vals-REF))), vals.tolist()

if __name__=="__main__":
    settings=[(2,10),(3,8),(4,3)]
    out={}
    for n,M in settings:
        key=f"{n},{M}"
        out[key]={}
        for name,kw in [
            ("raw",dict(mode="raw")),
            ("sample",dict(mode="sample",factor=1)),
            ("batch",dict(mode="batch",factor=1)),
            ("z0",dict(mode="scale",c=0))]:
            out[key][name]=run(n,M,**kw)[0]
        out[key]["fzero"]=fzero(n,M)[0]
        out[key]["radius"]=[]
        for a in [0,.25,.5,.75,1]:
            out[key]["radius"].append([a,run(n,M,mode="sample",factor=a)[0],run(n,M,mode="batch",factor=a)[0]])
        out[key]["scale"]=[]
        for c in np.linspace(0,1,11):
            out[key]["scale"].append([float(c),run(n,M,mode="scale",c=float(c))[0]])
        print(key,out[key])
    Path("/mnt/data/funding_life_or_death_summary.json").write_text(json.dumps(out,indent=2))
