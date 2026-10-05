import math, json, time, argparse
from pathlib import Path
import numpy as np
from scipy.special import roots_laguerre, beta
from scipy.integrate import quad
from scipy.stats import norm

METHODS = ['raw','samplewise_ir','mpbc','uniform_shrink']

def project_rows_ball(z, radii):
    z=np.asarray(z,float); radii=np.asarray(radii,float).reshape(-1)
    n=np.linalg.norm(z,axis=1)
    s=np.minimum(1.0, radii/np.maximum(n,1e-30))
    return z*s[:,None]

def uniform_shrink_ball(z, radii):
    z=np.asarray(z,float); radii=np.asarray(radii,float).reshape(-1)
    n=np.linalg.norm(z,axis=1)
    a=min(1.0,float(np.min(radii/np.maximum(n,1e-30))))
    return z*a, a

def mpbc_ball(z, radii):
    z=np.asarray(z,float); radii=np.asarray(radii,float).reshape(-1)
    m=np.mean(z,axis=0)
    if np.any(np.linalg.norm(m) > radii + 1e-12):
        return z.copy(), 1.0, False
    alpha=1.0
    for zi,R in zip(z,radii):
        d=zi-m; A=float(np.dot(d,d))
        if A <= 1e-30: continue
        B=2*float(np.dot(m,d)); C=float(np.dot(m,m)-R*R)
        disc=max(B*B-4*A*C,0.0)
        alpha=min(alpha,(-B+math.sqrt(disc))/(2*A))
    alpha=min(1.0,max(0.0,alpha))
    return m[None,:]+alpha*(z-m[None,:]),alpha,True

def transform_ball_batch(z,radii,method):
    if method=='raw':
        return z.copy(),{'feasible':False,'alpha':1.0,'mean_shift':0.0}
    if method=='samplewise_ir':
        out=project_rows_ball(z,radii)
        return out,{'feasible':True,'alpha':np.nan,'mean_shift':float(np.linalg.norm(out.mean(0)-z.mean(0)))}
    if method=='uniform_shrink':
        out,a=uniform_shrink_ball(z,radii)
        return out,{'feasible':True,'alpha':a,'mean_shift':float(np.linalg.norm(out.mean(0)-z.mean(0)))}
    if method=='mpbc':
        out,a,feas=mpbc_ball(z,radii)
        return out,{'feasible':feas,'alpha':a,'mean_shift':float(np.linalg.norm(out.mean(0)-z.mean(0)))}
    raise ValueError(method)

class HJBLog:
    def __init__(self,d,T=1.0): self.d=d; self.T=T; self.sigma=math.sqrt(2.0)
    def g(self,x): return np.log((1+np.sum(x*x,axis=-1))/2.0)
    def fz(self,z): return -0.5*np.sum(z*z,axis=-1)
    def exact(self,t,x,nquad=160):
        r2=np.sum(x*x,axis=1); tau=self.T-t
        nodes,weights=roots_laguerre(nquad); s=nodes[None,:]; den=1+4*tau[:,None]*s
        integ=np.sum(weights[None,:]*np.exp((-self.d/2)*np.log(den)-s*r2[:,None]/den),axis=1)
        return -np.log(2*integ)

def sample_ball(rng,n,d,radius=1.0):
    z=rng.normal(size=(n,d)); z/=np.linalg.norm(z,axis=1,keepdims=True)
    return z*(rng.random(n)**(1/d)*radius)[:,None]

def hjb_one_rep(eq,t,x,M,seed):
    B,d=x.shape; T=eq.T; dt=T-t; rng=np.random.default_rng(seed)
    ZT=rng.normal(size=(B,M*M,d)); XT=x[:,None,:]+eq.sigma*np.sqrt(dt)[:,None,None]*ZT
    uterm=np.mean(eq.g(XT),axis=1)
    tau=rng.random((B,M)); ds=tau*dt[:,None]; sr=np.sqrt(np.maximum(ds,1e-12))
    ZR=rng.normal(size=(B,M,d)); XR=x[:,None,:]+eq.sigma*sr[:,:,None]*ZR; R=t[:,None]+ds
    dt1=T-R; sq1=np.sqrt(np.maximum(dt1,1e-12))
    Z1=rng.normal(size=(B,M,M,d)); XT1=XR[:,:,None,:]+eq.sigma*sq1[:,:,None,None]*Z1
    g1=eq.g(XT1)
    z1=np.mean(g1[:,:,:,None]*Z1/np.maximum(sq1[:,:,None,None],1e-12),axis=2)
    outs={}; meta={}; rad=math.sqrt(2.0)
    for method in METHODS:
        fvals=np.zeros((B,M)); feas=[]; alphas=[]; shifts=[]
        for b in range(B):
            zz,st=transform_ball_batch(z1[b],np.full(M,rad),method)
            fvals[b]=eq.fz(zz); feas.append(st['feasible']); alphas.append(st['alpha']); shifts.append(st['mean_shift'])
        outs[method]=uterm+dt*np.mean(fvals,axis=1)
        meta[method]={'mpbc_feasible_rate':float(np.mean(feas)) if method=='mpbc' else None,
                      'alpha_mean':float(np.nanmean(alphas)) if method!='samplewise_ir' else None,
                      'mean_shift':float(np.mean(shifts))}
    return outs,meta

class Funding:
    def __init__(self,d=100,T=.5,sigma=.2,mu=.06,Rl=.04,Rb=.06):
        self.d=d; self.T=T; self.sigma=sigma; self.mu=mu; self.Rl=Rl; self.Rb=Rb
    def g(self,x):
        m=np.max(x,axis=-1); return np.maximum(m-120,0)-2*np.maximum(m-150,0)
    def f(self,y,z):
        sz=np.sum(z,axis=-1)
        return -self.Rl*y-((self.mu-self.Rl)/self.sigma)*sz+(self.Rb-self.Rl)*np.maximum(sz/self.sigma-y,0)
    def radius(self,t): return np.exp(.5*self.sigma**2*(self.T-t))

def funding_one_seed(eq,M,seed):
    d=eq.d; T=eq.T; x0=np.full(d,100.0); rng=np.random.default_rng(seed)
    W=rng.normal(size=(M*M,d))*math.sqrt(T)
    XT=x0[None,:]*np.exp((eq.mu-.5*eq.sigma**2)*T+eq.sigma*W)
    gg=eq.g(XT); gx=float(eq.g(x0[None,:])[0]); uterm=gx+np.mean(gg-gx)
    a=.5; r=rng.power(a,size=M); R=T*r
    W0=rng.normal(size=(M,d))*np.sqrt(R)[:,None]
    XR=x0[None,:]*np.exp((eq.mu-.5*eq.sigma**2)*R[:,None]+eq.sigma*W0)
    tau=T-R; W1=rng.normal(size=(M,M,d))*np.sqrt(tau)[:,None,None]
    XT1=XR[:,None,:]*np.exp((eq.mu-.5*eq.sigma**2)*tau[:,None,None]+eq.sigma*W1)
    gT=eq.g(XT1); gXR=eq.g(XR)
    y1=gXR+np.mean(gT-gXR[:,None],axis=1)
    z1=np.mean((gT-gXR[:,None])[:,:,None]*W1/np.maximum(tau[:,None,None],1e-12),axis=1)
    delta=z1/(eq.sigma*np.maximum(XR,1e-12)); radii=eq.radius(R)
    outs={}; meta={}
    for method in METHODS:
        dp,st=transform_ball_batch(delta,radii,method); zp=eq.sigma*XR*dp
        fval=eq.f(y1,zp); outs[method]=uterm+T*np.mean((r**(1-a))*fval)/a
        meta[method]={'mpbc_feasible':st['feasible'] if method=='mpbc' else None,
                      'alpha':st['alpha'],'mean_shift_delta':st['mean_shift'],
                      'raw_mean_delta_norm':float(np.linalg.norm(delta.mean(0)))}
    return outs,meta

class NWBatchMLP:
    def __init__(self,d,M,N=12,T=.25,method='raw',seed=0):
        self.d=d; self.M=M; self.N=N; self.T=T; self.method=method; self.rng=np.random.default_rng(seed)
        self.mu0=.06; self.sigma=.2; self.L=10.; self.K0=25.; self.K1=95.; self.K2=120.; self.feas=[]
    def g(self,x):
        m=np.max(x); return max(0.,m-self.K1)-2*max(0.,m-self.K2)
    def radius(self,t): return math.exp(self.mu0*(self.T-t))
    def f(self,z): return (self.L/self.d)*max(0.,np.max(np.abs(z))-self.K0)
    def transform(self,times,z):
        out,st=transform_ball_batch(z,np.array([self.radius(float(t)) for t in times]),self.method)
        if self.method=='mpbc': self.feas.append(float(st['feasible']))
        return out
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
    def call(K):
        return quad(lambda y:1.0-norm.cdf((y-mean)/sd)**d,K,mean+12*sd+10,epsabs=1e-10,limit=200)[0]
    return call(95)-2*call(120)

if __name__=='__main__':
    print('Use this module from batch_shrink_paired.py or a notebook; see docs/batch_contraction_report.md.')
