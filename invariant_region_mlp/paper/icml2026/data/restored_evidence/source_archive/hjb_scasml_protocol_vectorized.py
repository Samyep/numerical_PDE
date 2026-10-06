import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.special import roots_laguerre


class RosenbrockHJB:
    def __init__(self, d, T=1.0, seed_A=0):
        self.d = int(d); self.T=float(T); self.sigma=math.sqrt(2.0)
        rng=np.random.default_rng(seed_A)
        self.c1=rng.uniform(0.5,1.5,self.d-1)
        self.c2=rng.uniform(0.5,1.5,self.d-1)
        A=np.zeros((self.d,self.d),dtype=np.float64)
        for i in range(self.d-1):
            c=self.c1[i]
            A[i,i]+=c; A[i+1,i+1]+=c; A[i,i+1]-=c; A[i+1,i]-=c
            A[i+1,i+1]+=self.c2[i]
        self.A=A
        self.evals,self.evecs=np.linalg.eigh(A)
        self.lmax=float(self.evals.max()); self.trA=float(np.trace(A))
        self.zrad=math.sqrt(2.0*self.lmax)

    def q_batch(self, x):
        # exploit tridiagonal Rosenbrock form, avoids Bxdxd temporaries
        dx=x[...,:-1]-x[...,1:]
        return np.sum(self.c1*dx*dx + self.c2*x[...,1:]*x[...,1:],axis=-1)

    def g_batch(self,x):
        return np.log((1.0+self.q_batch(x))/2.0)

    def u_bounds(self,t,x):
        q=self.q_batch(x)
        lo=-math.log(2.0)*np.ones_like(q)
        hi=np.log((1.0+q+2.0*(self.T-t)*self.trA)/2.0)
        return lo,hi

    def exact_batch(self,t,x,nquad=64,chunk=256):
        s,w=roots_laguerre(nquad)
        out=np.empty(len(t),dtype=np.float64)
        # eigenbasis transform can be chunked
        for st in range(0,len(t),chunk):
            en=min(st+chunk,len(t)); tc=t[st:en]; xc=x[st:en]
            y=xc@self.evecs
            tau=(self.T-tc)[:,None,None]                  # B,1,1
            ss=s[None,:,None]                            # 1,Q,1
            lam=self.evals[None,None,:]                  # 1,1,d
            den=1.0+4.0*tau*ss*lam                       # B,Q,d
            logdet=-0.5*np.sum(np.log(den),axis=2)       # B,Q
            expo=-s[None,:]*np.sum((self.evals[None,None,:]*(y[:,None,:]**2))/den,axis=2)
            integ=np.sum(w[None,:]*np.exp(logdet+expo),axis=1)
            out[st:en]=-np.log(2.0*integ)
        return out


def sample_test_points(d, num_domain=1000, num_boundary=200, seed=2026):
    rng=np.random.default_rng(seed+d)
    # mimic LQG.test_geometry: unit hypersphere x [0,1]
    z=rng.normal(size=(num_domain,d)); z/=np.linalg.norm(z,axis=1,keepdims=True)
    rr=rng.random(num_domain)**(1.0/d)
    xin=z*rr[:,None]
    zb=rng.normal(size=(num_boundary,d)); zb/=np.linalg.norm(zb,axis=1,keepdims=True)
    x=np.concatenate([xin,zb],axis=0)
    t=rng.uniform(0.0,0.95,size=len(x))  # avoid t=T singular EBL weight; repo tests hit continuous time almost surely
    return t.astype(np.float64),x.astype(np.float64)


def project_z(z, mode, zrad):
    if mode=='heuristic':
        return np.clip(z,-10.0,10.0), np.zeros(z.shape[:-1],dtype=bool)
    norms=np.linalg.norm(z,axis=-1)
    if mode=='hard':
        scale=np.minimum(1.0,zrad/np.maximum(norms,1e-30))
        return z*scale[...,None], norms>zrad
    if mode=='hard_box':
        return np.clip(z,-zrad,zrad), np.any(np.abs(z)>zrad,axis=-1)
    if mode=='raw':
        return z, np.zeros(z.shape[:-1],dtype=bool)
    raise ValueError(mode)


def solve_n2_chunk(eq,t,x,M,rng,mode='heuristic',terminal_norm='corrected',radius_factor=1.0,dtype=np.float32):
    # Batched full-history n=2 MLP with uniform random time, matching SCaSML's MLP_full_history sampling.
    B,d=x.shape; T=eq.T; sig=eq.sigma; zrad=radius_factor*eq.zrad
    tf=t.astype(dtype); xf=x.astype(dtype)
    dt=(T-tf).astype(dtype)
    sqdt=np.sqrt(dt)
    gx=eq.g_batch(xf).astype(dtype)

    # Top-level terminal term: M^2 samples.
    ZT=rng.normal(size=(B,M*M,d)).astype(dtype)
    XT=xf[:,None,:]+sig*sqdt[:,None,None]*ZT
    gT=eq.g_batch(XT).astype(dtype)
    u0=np.mean(gT,axis=1)
    if terminal_norm=='corrected':
        z0=np.mean((gT-gx[:,None])[:,:,None]*ZT/np.maximum(sqdt[:,None,None],1e-6),axis=1)
    elif terminal_norm=='repo':
        # SCaSML MLP_full_history.py as written: mean(g*std_normal)/(T-t)
        z0=np.mean(gT[:,:,None]*ZT,axis=1)/np.maximum(dt[:,None],1e-6)
    else: raise ValueError(terminal_norm)

    # l=1 contribution: M uniformly sampled intermediate times.
    tau=rng.uniform(size=(B,M)).astype(dtype)
    dtr=np.maximum(tau*dt[:,None],1e-6).astype(dtype)
    sr=np.sqrt(dtr)
    ZR=rng.normal(size=(B,M,d)).astype(dtype)
    XR=xf[:,None,:]+sig*sr[:,:,None]*ZR
    R=tf[:,None]+dtr
    gxR=eq.g_batch(XR).astype(dtype)

    # U_1(R, X_R): terminal term with M samples. Shape B,M,M,d.
    dt1=np.maximum(T-R,1e-6).astype(dtype)
    s1=np.sqrt(dt1)
    Z1=rng.normal(size=(B,M,M,d)).astype(dtype)
    XT1=XR[:,:,None,:]+sig*s1[:,:,None,None]*Z1
    gT1=eq.g_batch(XT1).astype(dtype)
    u1=np.mean(gT1,axis=2)
    if terminal_norm=='corrected':
        z1=np.mean((gT1-gxR[:,:,None])[:,:,:,None]*Z1/np.maximum(s1[:,:,None,None],1e-6),axis=2)
    else:
        z1=np.mean(gT1[:,:,:,None]*Z1,axis=2)/np.maximum(dt1[:,:,None],1e-6)

    z1p,viol1=project_z(z1,mode,zrad)
    if mode=='heuristic':
        u1p=np.clip(u1,-10.0,10.0)
    elif mode in ('hard','hard_box'):
        lo1,hi1=eq.u_bounds(R,XR)
        u1p=np.minimum(np.maximum(u1,lo1),hi1)
    else: u1p=u1
    # HJB f only depends on z.
    f1=-0.5*np.sum(z1p*z1p,axis=-1)  # B,M

    uint=dt*np.mean(f1,axis=1)
    zint=dt[:,None]*np.mean(f1[:,:,None]*ZR/np.maximum(sr[:,:,None],1e-6),axis=1)
    u2=u0+uint; z2=z0+zint

    z2p,viol2=project_z(z2,mode,zrad)
    if mode=='heuristic':
        u2p=np.clip(u2,-10.0,10.0)
    elif mode in ('hard','hard_box'):
        lo2,hi2=eq.u_bounds(tf,xf)
        u2p=np.minimum(np.maximum(u2,lo2),hi2)
    else: u2p=u2

    stats={
        'u1_grad_violation': float(np.mean(viol1)),
        'u2_grad_violation': float(np.mean(viol2)),
        'u1_preproj_grad_norm_median': float(np.median(np.linalg.norm(z1,axis=-1))),
        'u2_preproj_grad_norm_median': float(np.median(np.linalg.norm(z2,axis=-1))),
    }
    return u2p.astype(np.float64),stats


def run_dim(d, M=10, reps=10, num_domain=1000, num_boundary=200, modes=('heuristic','hard','hard_box'), terminal_norm='corrected', chunk=32, seed_points=2026):
    eq=RosenbrockHJB(d=d)
    t,x=sample_test_points(d,num_domain,num_boundary,seed_points)
    truth=eq.exact_batch(t,x,nquad=64)
    rows=[]
    predictions={m:[] for m in modes}
    aggstats={m:[] for m in modes}
    tstart=time.time()
    # Pair methods by using identical random streams per rep/chunk: reconstruct rng with same seed per mode.
    for rep in range(reps):
        for mode in modes:
            vals=np.empty(len(t),dtype=np.float64); stlist=[]
            for st in range(0,len(t),chunk):
                en=min(st+chunk,len(t))
                seed=np.random.SeedSequence([98731,d,rep,st])
                rng=np.random.default_rng(seed)
                v,ss=solve_n2_chunk(eq,t[st:en],x[st:en],M,rng,mode=mode,terminal_norm=terminal_norm)
                vals[st:en]=v; stlist.append(ss)
            predictions[mode].append(vals)
            aggstats[mode].append(stlist)
    for mode in modes:
        pred=np.asarray(predictions[mode]); err=pred-truth[None,:]
        rel=np.linalg.norm(err,axis=1)/np.linalg.norm(truth)
        ae=np.abs(err)
        flatstats=[s for repstats in aggstats[mode] for s in repstats]
        rows.append({
            'd':d,'n':2,'M':M,'reps':reps,'npts':len(t),'mode':mode,'terminal_norm':terminal_norm,
            'rel_l2_mean':float(rel.mean()),'rel_l2_std':float(rel.std(ddof=1) if reps>1 else 0),
            'mae':float(ae.mean()),'p95_abs_error':float(np.quantile(ae,0.95)),
            'bias':float(err.mean()),'pointwise_run_std':float(np.mean(np.std(pred,axis=0,ddof=1))) if reps>1 else 0.0,
            'u1_grad_violation':float(np.mean([s['u1_grad_violation'] for s in flatstats])),
            'u2_grad_violation':float(np.mean([s['u2_grad_violation'] for s in flatstats])),
            'u1_preproj_grad_norm_median':float(np.median([s['u1_preproj_grad_norm_median'] for s in flatstats])),
            'u2_preproj_grad_norm_median':float(np.median([s['u2_preproj_grad_norm_median'] for s in flatstats])),
            'z_radius':float(eq.zrad),'elapsed_total':float(time.time()-tstart),
            'f_evals_nominal_per_point':220,'g_evals_nominal_per_point':211,
        })
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--dims',nargs='+',type=int,default=[100,120,140,160])
    ap.add_argument('--M',type=int,default=10); ap.add_argument('--reps',type=int,default=10)
    ap.add_argument('--domain',type=int,default=1000); ap.add_argument('--boundary',type=int,default=200)
    ap.add_argument('--terminal-norm',choices=['corrected','repo'],default='corrected')
    ap.add_argument('--chunk',type=int,default=32); ap.add_argument('--out',default='/mnt/data/hjb_scasml_protocol_results.json')
    args=ap.parse_args(); allrows=[]
    for d in args.dims:
        print('RUN',d,args.terminal_norm,flush=True)
        rows=run_dim(d,M=args.M,reps=args.reps,num_domain=args.domain,num_boundary=args.boundary,terminal_norm=args.terminal_norm,chunk=args.chunk)
        allrows+=rows
        for r in rows: print(json.dumps(r),flush=True)
        Path(args.out).write_text(json.dumps(allrows,indent=2))
    Path(args.out).write_text(json.dumps(allrows,indent=2))

if __name__=='__main__': main()
