import numpy as np, math, json, time
from pathlib import Path

MODES=('raw','recursive_subspace','batch_structural')

class Solver:
    def __init__(self,d,m,seed):
        self.d=d; self.m=m
        ss=np.random.SeedSequence(seed).spawn(3)
        self.times,self.active,self.inactive=[np.random.default_rng(s) for s in ss]
        self.batch_groups=0; self.batch_active=0; self.alpha_sum=0.
    def normals(self,b,k):
        g=np.empty((b,k,self.d)); g[:,:,0]=self.active.standard_normal((b,k))
        if self.d>1:g[:,:,1:]=self.inactive.standard_normal((b,k,self.d-1))
        return g
    def finish(self,out):
        # only the samplewise structural method returns a projected state.
        out[:,1,2:]=0.
        return out
    def driver_grouped(self,child,b,k):
        # child: (b*k, modes, d+1); reshape siblings to apply one common factor.
        c=child.reshape(b,k,len(MODES),self.d+1)
        f=np.linalg.norm(c[:,:,:,2:],axis=-1)
        j=MODES.index('batch_structural')
        z=c[:,:,j,1:].copy(); z[:,:,1:]=0.
        maxact=np.max(np.abs(z[:,:,0]),axis=1)
        alpha=np.minimum(1.,1./np.maximum(maxact,1e-30))
        z*=alpha[:,None,None]
        # Q kills all inactive coordinates, hence the generator is exactly zero.
        f[:,:,j]=np.linalg.norm(z[:,:,1:],axis=-1)
        self.batch_groups+=b; self.batch_active+=int(np.sum(alpha<1));self.alpha_sum+=float(alpha.sum())
        return f
    def solve(self,n,h,x1):
        b=len(h)
        if n<=0:return np.zeros((b,len(MODES),self.d+1))
        k=self.m**n; normal=self.normals(b,k); sh=np.sqrt(h)
        diff=np.abs(x1[:,None]+sh[:,None]*normal[:,:,0])-np.abs(x1[:,None])
        term=np.empty((b,self.d+1));term[:,0]=np.abs(x1)+diff.mean(1);term[:,1:]=np.einsum('bi,bij->bj',diff,normal)/k/sh[:,None]
        out=np.repeat(term[:,None,:],len(MODES),axis=1)
        for level in range(1,n):
            k=self.m**(n-level);r=self.times.random((b,k));
            while np.any(r==0):r[r==0]=self.times.random(np.count_nonzero(r==0))
            r=r*r;normal=self.normals(b,k);xp=x1[:,None]+sh[:,None]*np.sqrt(r)*normal[:,:,0];hp=h[:,None]*(1-r)
            high=self.solve(level,hp.ravel(),xp.ravel());fd=self.driver_grouped(high,b,k)
            if level>1:
                low=self.solve(level-1,hp.ravel(),xp.ravel());fd-=self.driver_grouped(low,b,k)
            out[:,:,0]+=np.mean(2*h[:,None,None]*np.sqrt(r)[:,:,None]*fd,axis=1)
            out[:,:,1:]+=2*sh[:,None,None]*np.einsum('bik,bid->bkd',fd,normal)/k
        return self.finish(out)

def run(d,n=3,m=3,reps=512,batch=16,seed=20261005):
    vals=[]; active=groups=0; asum=0.; start=time.time()
    for block in range(reps//batch):
        block_seed=int(np.random.SeedSequence([seed,n,m,block]).generate_state(1)[0])
        s=Solver(d,m,block_seed);y=s.solve(n,np.ones(batch),np.zeros(batch));vals.append(y[:,:,0]);active+=s.batch_active;groups+=s.batch_groups;asum+=s.alpha_sum
    vals=np.concatenate(vals);truth=math.sqrt(2/math.pi);rows={}
    for j,mode in enumerate(MODES):
        e=vals[:,j]-truth;rows[mode]=dict(value_rmse=float(np.sqrt(np.mean(e*e))),value_bias=float(e.mean()),value_mean=float(vals[:,j].mean()))
    return dict(d=d,n=n,m=m,reps=reps,results=rows,
                max_abs_batch_vs_subspace=float(np.max(np.abs(vals[:,2]-vals[:,1]))),
                batch_activation_rate=active/groups if groups else 0.,batch_alpha_mean=asum/groups if groups else 1.,elapsed_seconds=time.time()-start)

if __name__=='__main__':
    import sys
    out=[]
    for d in (10,100,1000):
        r=run(d);out.append(r);print('counter',d,{k:round(v['value_rmse'],6) for k,v in r['results'].items()},'diff',r['max_abs_batch_vs_subspace'],'act',r['batch_activation_rate'],flush=True)
    Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/counterexample_batchwise_sanity_results.json').write_text(json.dumps(dict(settings=out),indent=2))
