import json, math, time
from pathlib import Path
import numpy as np
from scipy.special import roots_laguerre, logsumexp
import jax.random as jr

DIMS=(100,120,140,160)
M=10
N=2
REPS=10
RADIUS=math.sqrt(15.0)

class RosenbrockHJB:
    T=1.0
    def __init__(self,d):
        self.d=d
        self.c1=np.asarray(jr.uniform(jr.PRNGKey(0),(d-1,),minval=.5,maxval=1.5),dtype=np.float64)
        self.c2=np.asarray(jr.uniform(jr.PRNGKey(1),(d-1,),minval=.5,maxval=1.5),dtype=np.float64)
        A=np.zeros((d,d),dtype=np.float64); idx=np.arange(d-1)
        A[idx,idx]+=self.c1; A[idx+1,idx+1]+=self.c1+self.c2
        A[idx,idx+1]-=self.c1; A[idx+1,idx]-=self.c1
        self.A=A
        self.evals,self.evecs=np.linalg.eigh(A)
        self.trA=float(np.trace(A))
    def q(self,x):
        x=np.asarray(x)
        dx=x[..., :-1]-x[...,1:]
        return np.sum(self.c1*dx*dx+self.c2*x[...,1:]*x[...,1:],axis=-1)
    def g(self,x):
        return np.log((1.0+self.q(x))/2.0)
    def f(self,z):
        return -0.5*np.sum(z*z,axis=-1)

def stable_u_z(eq,t,x,nquad=32,chunk=256,need_z=True):
    t=np.asarray(t,dtype=np.float64); x=np.asarray(x,dtype=np.float64)
    nodes,weights=roots_laguerre(nquad); logw=np.log(weights); lam=eq.evals; E=eq.evecs
    u=np.empty(len(t)); zout=np.empty((len(t),eq.d)) if need_z else None
    for st in range(0,len(t),chunk):
        en=min(st+chunk,len(t)); tc=t[st:en]; xc=x[st:en]; y=xc@E
        tau=eq.T-tc; q=eq.q(xc); c=1.0+q+2.0*tau*eq.trA
        s=nodes[None,:,None]/c[:,None,None]
        den=1.0+4.0*tau[:,None,None]*s*lam[None,None,:]
        logj=-0.5*np.sum(np.log(den),axis=2)-s[:,:,0]*np.sum(lam[None,None,:]*(y[:,None,:]**2)/den,axis=2)
        logterms=logw[None,:]-np.log(c)[:,None]+nodes[None,:]*(1.0-1.0/c[:,None])+logj
        logI=logsumexp(logterms,axis=1)
        u[st:en]=-math.log(2.0)-logI
        if need_z:
            alpha=np.exp(logterms-logI[:,None])
            numer=np.sum(alpha[:,:,None]*(s/den),axis=1)
            grad_y=2.0*(lam[None,:]*y)*numer
            zout[st:en]=math.sqrt(2.0)*(grad_y@E.T)
    return (u,zout) if need_z else u

def test_points(d,seed=20261005):
    # SCaSML-style scale/geometry: 1000 interior unit-ball points + 200 unit-sphere points,
    # each with a time drawn uniformly on [0,1]. Fixed across paired repetitions.
    rng=np.random.default_rng(seed+d)
    ni,nb=1000,200
    xi=rng.normal(size=(ni,d)); xi/=np.linalg.norm(xi,axis=1,keepdims=True)
    xi*=rng.random(ni)[:,None]**(1.0/d)
    xb=rng.normal(size=(nb,d)); xb/=np.linalg.norm(xb,axis=1,keepdims=True)
    x=np.vstack([xi,xb]); t=rng.random(ni+nb)
    return t,x

def project_sample(z,r=RADIUS):
    n=np.linalg.norm(z,axis=-1)
    return z*np.minimum(1.0,r/np.maximum(n,1e-30))[...,None]

def project_batch(z,r=RADIUS):
    # z shape B x M x d; one common scale per parent Monte Carlo sibling batch.
    n=np.linalg.norm(z,axis=-1)
    a=np.minimum(1.0,r/np.maximum(np.max(n,axis=1),1e-30))
    return z*a[:,None,None],a

def one_rep(eq,t,x,seed,diag_parents=64,chunk=96):
    rng=np.random.default_rng(seed); B=len(t); d=eq.d
    outs={k:np.empty(B) for k in ('raw','samplewise','batch')}
    corrections={k:np.empty(B) for k in outs}
    diag=[]
    batch_alpha=[]
    batch_active=[]
    sample_viols=[]
    parents_seen=0
    for st in range(0,B,chunk):
        en=min(B,st+chunk); tc=t[st:en]; xc=x[st:en]; b=en-st; dt=1.0-tc
        # Root terminal block (M^2 samples).
        ZT=rng.normal(size=(b,M*M,d)); XT=xc[:,None,:]+math.sqrt(2.0)*np.sqrt(dt)[:,None,None]*ZT
        uterm=eq.g(XT).mean(axis=1)
        # l=1 nonlinear block (M intermediate states). The centered terminal difference
        # is the variance-reduced Elworthy--Bismut--Li estimator used in the corrected protocol.
        ds=rng.random((b,M))*dt[:,None]; sr=np.sqrt(np.maximum(ds,1e-14))
        ZR=rng.normal(size=(b,M,d)); XR=xc[:,None,:]+math.sqrt(2.0)*sr[:,:,None]*ZR; R=tc[:,None]+ds
        rem=1.0-R; sq=np.sqrt(np.maximum(rem,1e-14))
        Z1=rng.normal(size=(b,M,M,d)); XT1=XR[:,:,None,:]+math.sqrt(2.0)*sq[:,:,None,None]*Z1
        gT=eq.g(XT1); gR=eq.g(XR)
        zraw=np.mean((gT-gR[:,:,None])[:,:,:,None]*Z1/sq[:,:,None,None],axis=2)
        zsample=project_sample(zraw)
        zbatch,alpha=project_batch(zraw)
        nr=np.linalg.norm(zraw,axis=-1)
        sample_viols.extend((nr>RADIUS).ravel().tolist())
        batch_alpha.extend(alpha.tolist()); batch_active.extend((alpha<1.0).tolist())
        for mode,z in [('raw',zraw),('samplewise',zsample),('batch',zbatch)]:
            corr=dt*eq.f(z).mean(axis=1)
            corrections[mode][st:en]=corr
            outs[mode][st:en]=uterm+corr
        # Keep complete sibling groups for generator-MSE diagnostics.
        take=min(diag_parents-parents_seen,b)
        if take>0:
            diag.append(dict(t=R[:take].copy(),x=XR[:take].copy(),zraw=zraw[:take].copy(),
                             zsample=zsample[:take].copy(),zbatch=zbatch[:take].copy()))
            parents_seen+=take
    return outs,corrections,diag,dict(
        raw_violation_rate=float(np.mean(sample_viols)),
        batch_activation_rate=float(np.mean(batch_active)),
        batch_alpha_mean=float(np.mean(batch_alpha)),
        batch_alpha_p10=float(np.quantile(batch_alpha,.1)),
        batch_alpha_median=float(np.median(batch_alpha)))

def flatten_diag(diags,key):
    return np.concatenate([x[key].reshape((-1, x[key].shape[-1])) for x in diags],axis=0)

def flatten_scalar(diags,key):
    return np.concatenate([x[key].ravel() for x in diags],axis=0)

def run_dim(d):
    eq=RosenbrockHJB(d); t,x=test_points(d); truth=stable_u_z(eq,t,x,need_z=False)
    denom=np.linalg.norm(truth)
    rel={k:[] for k in ('raw','samplewise','batch')}; corr={k:[] for k in rel}
    all_diag=[]; meta=[]
    for rep in range(REPS):
        outs,cs,diag,mm=one_rep(eq,t,x,seed=910000+d*100+rep)
        for k in rel:
            rel[k].append(float(np.linalg.norm(outs[k]-truth)/denom)); corr[k].append(cs[k])
        all_diag.extend(diag); meta.append(mm)
    # Stable reference gradient on the diagnostic intermediate states.
    td=flatten_scalar(all_diag,'t'); xd=flatten_diag(all_diag,'x')
    _,ztrue=stable_u_z(eq,td,xd,need_z=True,chunk=128)
    ftrue=eq.f(ztrue)
    gmse={}; incvar={}
    for mode,key in [('raw','zraw'),('samplewise','zsample'),('batch','zbatch')]:
        z=flatten_diag(all_diag,key); ff=eq.f(z)
        gmse[mode]=float(np.mean((ff-ftrue)**2)); incvar[mode]=float(np.var(ff,ddof=1))
    # Variance across paired repetitions of the actual level correction, averaged over fixed test points.
    corrvar={k:float(np.mean(np.var(np.stack(corr[k],axis=0),axis=0,ddof=1))) for k in corr}
    return dict(
        d=d,n=N,M=M,reps=REPS,n_test=len(t),radius=RADIUS,
        lambda_max=float(eq.evals[-1]),trace_A=eq.trA,
        relative_l2={k:dict(mean=float(np.mean(v)),std=float(np.std(v,ddof=1)),values=v) for k,v in rel.items()},
        paired_batch_vs_sample=dict(
            mean_difference=float(np.mean(np.array(rel['batch'])-np.array(rel['samplewise']))),
            win_fraction=float(np.mean(np.array(rel['batch'])<np.array(rel['samplewise'])))),
        mechanism=dict(generator_mse=gmse,individual_increment_variance=incvar,correction_variance=corrvar,
                       generator_mse_reduction_batch_vs_raw=1-gmse['batch']/gmse['raw'],
                       correction_variance_reduction_batch_vs_raw=1-corrvar['batch']/corrvar['raw']),
        activation={k:float(np.mean([m[k] for m in meta])) for k in meta[0]}
    )

def main(out):
    start=time.time(); rows=[]
    for d in DIMS:
        r=run_dim(d); rows.append(r)
        print('HJB',d,'raw/sample/batch',*[round(r['relative_l2'][k]['mean'],6) for k in ('raw','samplewise','batch')],
              'gmse batch red',round(r['mechanism']['generator_mse_reduction_batch_vs_raw'],5),flush=True)
    payload=dict(protocol='SCaSML-scale corrected Rosenbrock HJB; fixed A per dimension; z-only corrections; final u never clipped',
                 rows=rows,elapsed_seconds=time.time()-start)
    Path(out).write_text(json.dumps(payload,indent=2))

if __name__=='__main__':
    import sys
    main(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/hjb_batchwise_headline_results.json')
