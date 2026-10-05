import math, json, time
from pathlib import Path
import numpy as np
from scipy.special import beta
from scipy.integrate import quad
from scipy.stats import norm

MODES=('raw','samplewise','batch')

def project_rows_ball(z,radii):
    n=np.linalg.norm(z,axis=1); s=np.minimum(1.,np.asarray(radii)/np.maximum(n,1e-30)); return z*s[:,None]
def uniform_shrink_ball(z,radii):
    n=np.linalg.norm(z,axis=1); a=min(1.,float(np.min(np.asarray(radii)/np.maximum(n,1e-30)))); return z*a,a

class NWBatchMLP:
    def __init__(self,d,M,N=12,T=.25,method='raw',seed=0):
        self.d=d; self.M=M; self.N=N; self.T=T; self.method=method; self.rng=np.random.default_rng(seed)
        self.mu0=.06; self.sigma=.2; self.L=10.; self.K0=25.; self.K1=95.; self.K2=120.; self.active=[]
    def g(self,x):
        m=np.max(x); return max(0.,m-self.K1)-2*max(0.,m-self.K2)
    def radius(self,t): return math.exp(self.mu0*(self.T-t))
    def f(self,z): return (self.L/self.d)*max(0.,np.max(np.abs(z))-self.K0)
    def transform(self,times,z):
        radii=np.array([self.radius(float(t)) for t in times]); nr=np.linalg.norm(z,axis=1)
        if self.method=='raw': return z
        if self.method=='samplewise': return project_rows_ball(z,radii)
        out,a=uniform_shrink_ball(z,radii); self.active.append(float(a<1)); return out
    def simulate(self,t,s,x):
        dt=max(self.T-t,1e-32)/self.N; S=self.N if abs(s-self.T)<1e-14 else int(np.floor((s-t)/dt)+1)
        W=self.rng.normal(size=(S,self.d),scale=math.sqrt(dt)); X=x.copy(); J=1.; V=np.zeros(self.d)
        for k in range(S):
            V += J*(1/self.sigma)*W[k]; J *= (1+self.mu0*dt); X=X+self.mu0*X*dt+self.sigma*W[k]
        return X,V/max(s-t,1e-32)
    def compute(self,t,x,n):
        if n==0: return np.zeros(self.d+1)
        gx=self.g(x); term=np.zeros(self.d+1)
        for _ in range(self.M**n):
            X,V=self.simulate(t,self.T,x); term+=(self.g(X)-gx)*np.r_[1.,V]
        u=term/(self.M**n)+np.r_[gx,np.zeros(self.d)]
        for l in range(n):
            Nout=self.M**(n-l); times=[]; vecs=[]; child=[]; prev=[]; rhos=[]
            for _ in range(Nout):
                R=t+(self.T-t)*self.rng.beta(.5,.5); Y,V=self.simulate(t,R,x); times.append(R); vecs.append(np.r_[1.,V])
                q=(R-t)/max(self.T-t,1e-32); rhos.append(math.sqrt(q*(1-q))*beta(.5,.5))
                child.append(self.compute(R,Y,l))
                if l>0: prev.append(self.compute(R,Y,l-1))
            times=np.asarray(times); child=np.asarray(child); fc=np.array([self.f(z) for z in self.transform(times,child[:,1:])])
            fp=np.zeros(Nout) if l==0 else np.array([self.f(z) for z in self.transform(times,np.asarray(prev)[:,1:])])
            u += (self.T-t)*np.sum((np.asarray(rhos)*(fc-fp))[:,None]*np.asarray(vecs),axis=0)/Nout
        return u

def nw_reference(d):
    mu=.06; sig=.2; T=.25; x0=100.; mean=x0*math.exp(mu*T); sd=sig*math.sqrt((math.exp(2*mu*T)-1)/(2*mu))
    def call(K): return quad(lambda y:1.0-norm.cdf((y-mean)/sd)**d,K,mean+12*sd+10,epsabs=1e-10,limit=200)[0]
    return call(95)-2*call(120)

def run(d,M=2,n=2,reps=30):
    ref=nw_reference(d); vals={m:[] for m in MODES}; act=[]; t0=time.time()
    for rep in range(reps):
        seed=20262000+d*10+rep
        for mode in MODES:
            s=NWBatchMLP(d,M,method=mode,seed=seed); v=s.compute(0.,np.full(d,100.),n)[0]; vals[mode].append(v)
            if mode=='batch': act.extend(s.active)
    rows={}
    for m in MODES:
        a=np.array(vals[m]); err=np.abs(a-ref); rows[m]=dict(mae=float(err.mean()),rmse=float(np.sqrt(np.mean((a-ref)**2))),mean=float(a.mean()),std=float(a.std(ddof=1)),values=a.tolist())
    return dict(d=d,n=n,M=M,reps=reps,reference=ref,results=rows,batch_activation_rate=float(np.mean(act)) if act else 0.,
                batch_win_fraction=float(np.mean(np.abs(np.array(vals['batch'])-ref)<np.abs(np.array(vals['samplewise'])-ref))),elapsed_seconds=time.time()-t0)

if __name__=='__main__':
    import sys
    out=[]
    for d in (100,200,300):
        r=run(d);out.append(r);print('NW',d,{m:round(r['results'][m]['mae'],6) for m in MODES},'win',r['batch_win_fraction'],flush=True)
    Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/neufeld_batch_m2_results.json').write_text(json.dumps(dict(settings=out),indent=2))
