import numpy as np, math, json, time
from pathlib import Path

MODES=('raw','samplewise','batch')
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

class FundingSolver:
    def __init__(self,M,seed,a=.5):
        self.M=M; self.eq=FundingEq(); self.a=a; self.rng=np.random.default_rng(seed)
        self.batch_groups=0; self.batch_active=0; self.alpha_sum=0.; self.sample_viol=0; self.sample_total=0
    def simulate(self,t,x,dt,normal):
        return x[:,None,:]*np.exp((self.eq.mu-.5*self.eq.sigma**2)*dt[:,:,None]+self.eq.sigma*np.sqrt(dt)[:,:,None]*normal)
    def terminal_state(self,t,x,n):
        b=len(t); d=self.eq.d; h=self.eq.T-t; k=self.M**n
        normal=self.rng.normal(size=(b,k,d)); W=np.sqrt(h)[:,None,None]*normal
        XT=x[:,None,:]*np.exp((self.eq.mu-.5*self.eq.sigma**2)*h[:,None,None]+self.eq.sigma*W)
        gx=self.eq.g(x); gT=self.eq.g(XT); diff=gT-gx[:,None]
        term=np.empty((b,d+1),float); term[:,0]=gx+diff.mean(1)
        term[:,1:]=np.mean(diff[:,:,None]*W/np.maximum(h[:,None,None],1e-14),axis=1)
        return term
    def transform(self,state,t,x,mode,group_shape):
        # state flattened B*K x (d+1), t B*K, x B*K x d. group_shape=(B,K)
        if mode=='raw': return state
        out=state.copy(); z=out[:,1:]; delta=z/(self.eq.sigma*np.maximum(x,1e-12)); rad=self.eq.radius(t); nr=np.linalg.norm(delta,axis=1)
        self.sample_total += len(nr) if mode=='samplewise' else 0
        self.sample_viol += int(np.sum(nr>rad)) if mode=='samplewise' else 0
        if mode=='samplewise':
            s=np.minimum(1.,rad/np.maximum(nr,1e-30)); out[:,1:]=self.eq.sigma*x*(delta*s[:,None]); return out
        B,K=group_shape; ratio=(rad/np.maximum(nr,1e-30)).reshape(B,K)
        alpha=np.minimum(1.,np.min(ratio,axis=1)); self.batch_groups+=B; self.batch_active+=int(np.sum(alpha<1)); self.alpha_sum+=float(alpha.sum())
        out[:,1:]=self.eq.sigma*x*(delta*alpha.repeat(K)[:,None]); return out
    def solve(self,n,t,x):
        t=np.asarray(t,float).reshape(-1); x=np.asarray(x,float).reshape(len(t),self.eq.d); b=len(t); d=self.eq.d
        if n<=0: return np.zeros((b,len(MODES),d+1))
        term=self.terminal_state(t,x,n); out=np.repeat(term[:,None,:],len(MODES),axis=1); h=self.eq.T-t
        for level in range(1,n):
            K=self.M**(n-level)
            r=self.rng.power(self.a,size=(b,K)); dt=h[:,None]*r; normal=self.rng.normal(size=(b,K,d)); W=np.sqrt(dt)[:,:,None]*normal
            XR=x[:,None,:]*np.exp((self.eq.mu-.5*self.eq.sigma**2)*dt[:,:,None]+self.eq.sigma*W); R=t[:,None]+dt
            high=self.solve(level,R.ravel(),XR.reshape(-1,d)).reshape(b,K,len(MODES),d+1)
            fd=np.empty((b,K,len(MODES)))
            for j,mode in enumerate(MODES):
                st=self.transform(high[:,:,j,:].reshape(-1,d+1),R.ravel(),XR.reshape(-1,d),mode,(b,K))
                fd[:,:,j]=self.eq.f(st[:,0],st[:,1:]).reshape(b,K)
            if level>1:
                low=self.solve(level-1,R.ravel(),XR.reshape(-1,d)).reshape(b,K,len(MODES),d+1)
                for j,mode in enumerate(MODES):
                    st=self.transform(low[:,:,j,:].reshape(-1,d+1),R.ravel(),XR.reshape(-1,d),mode,(b,K))
                    fd[:,:,j]-=self.eq.f(st[:,0],st[:,1:]).reshape(b,K)
            # importance weight inverse Beta(a,1) density
            w=h[:,None]*(r**(1-self.a))/self.a
            out[:,:,0]+=np.mean(w[:,:,None]*fd,axis=1)
            V=W/np.maximum(dt[:,:,None],1e-14)
            out[:,:,1:]+=np.mean(w[:,:,None,None]*fd[:,:,:,None]*V[:,:,None,:],axis=1)
        return out

def run_setting(n,M,reps=100,block=10,seed=20261005):
    vals=[]; metas=[]
    start=time.time()
    for st in range(0,reps,block):
        b=min(block,reps-st); s=FundingSolver(M,seed+100000*n+1000*M+st)
        t=np.zeros(b); x=np.full((b,100),100.)
        y=s.solve(n,t,x)[:,:,0]
        vals.append(y)
        metas.append(dict(batch_groups=s.batch_groups,batch_active=s.batch_active,alpha_sum=s.alpha_sum,
                          sample_viol=s.sample_viol,sample_total=s.sample_total))
    vals=np.concatenate(vals,axis=0)
    rows={}
    for j,m in enumerate(MODES):
        err=np.abs(vals[:,j]-REF); rows[m]=dict(mean=float(vals[:,j].mean()),std=float(vals[:,j].std(ddof=1)),
                                                mae=float(err.mean()),rmse=float(np.sqrt(np.mean((vals[:,j]-REF)**2))),values=vals[:,j].tolist())
    bg=sum(m['batch_groups'] for m in metas); ba=sum(m['batch_active'] for m in metas); asum=sum(m['alpha_sum'] for m in metas)
    sv=sum(m['sample_viol'] for m in metas); stot=sum(m['sample_total'] for m in metas)
    return dict(n=n,M=M,reps=reps,results=rows,
                batch_activation_rate=ba/bg if bg else 0.,batch_alpha_mean=asum/bg if bg else 1.,samplewise_violation_rate=sv/stot if stot else 0.,
                batch_win_fraction=float(np.mean(np.abs(vals[:,2]-REF)<np.abs(vals[:,1]-REF))),elapsed_seconds=time.time()-start)

if __name__=='__main__':
    import sys
    settings=[(2,10),(3,8),(4,3)]
    out=[]
    for n,M in settings:
        r=run_setting(n,M,reps=100,block=10);out.append(r)
        print('funding',n,M,{k:round(r['results'][k]['mae'],6) for k in MODES},'batchwin',r['batch_win_fraction'],'sec',round(r['elapsed_seconds'],2),flush=True)
    Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/funding_batchwise_deep_results.json').write_text(json.dumps(dict(reference=REF,settings=out),indent=2))
