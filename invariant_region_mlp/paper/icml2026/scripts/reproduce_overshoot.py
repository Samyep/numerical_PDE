#!/usr/bin/env python3
"""Fresh 36-setting n=2 value-channel MLP diagnostic, not the lost original samples.
Only child gradients are projected; final u is not clipped. Corrected scaled
Hopf--Cole quadrature. Uniform time is used for the n=2 value channel, not as
a claim about L2 convergence of arbitrary gradient-weighted recursions.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.special import roots_laguerre,logsumexp
from scipy.stats import pearsonr,spearmanr


def reference(t,x,lam,V,q,Q):
 r,w=roots_laguerre(Q);out=np.empty(len(t))
 for st in range(0,len(t),64):
  tc=t[st:st+64];xc=x[st:st+64];h=1-tc;y=xc@V;c=1+q(xc)+2*h*lam.sum()
  s=r[None,:,None]/c[:,None,None];den=1+4*h[:,None,None]*s*lam[None,None,:]
  logs=np.log(w)[None,:]-np.log(c)[:,None]+r[None,:]*(1-1/c[:,None])-.5*np.log(den).sum(2)-s[:,:,0]*(lam[None,None,:]*y[:,None,:]**2/den).sum(2)
  out[st:st+64]=-math.log(2)-logsumexp(logs,axis=1)
 return out


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 d=100;reps=10;Ms=[2,3,4,6,8,10];fac=np.array([1,1.25,1.5,2,3,4.]);rng=np.random.default_rng(0)
 c1=rng.uniform(.5,1.5,d-1);c2=rng.uniform(.5,1.5,d-1);A=np.zeros((d,d))
 for i,c in enumerate(c1):A[i,i]+=c;A[i+1,i+1]+=c+c2[i];A[i,i+1]-=c;A[i+1,i]-=c
 lam,V=np.linalg.eigh(A);radii=fac*math.sqrt(2*lam[-1])
 def q(x):return np.sum(c1*(x[...,:-1]-x[...,1:])**2+c2*x[...,1:]**2,axis=-1)
 def g(x):return np.log((1+q(x))/2)
 rng=np.random.default_rng(2026+d);xi=rng.normal(size=(1000,d));xi/=np.linalg.norm(xi,axis=1,keepdims=True);xi*=rng.random(1000)[:,None]**(1/d)
 xb=rng.normal(size=(200,d));xb/=np.linalg.norm(xb,axis=1,keepdims=True);x=np.concatenate([xi,xb]);t=rng.uniform(0,.95,len(x))
 truth=reference(t,x,lam,V,q,64);truth128=reference(t,x,lam,V,q,128)
 raw=np.empty((6,reps,len(t)));proj=np.empty((6,6,reps,len(t)));omega=np.zeros((6,6,reps));viol=np.zeros_like(omega);rows=[];start=time.time()
 for mi,M in enumerate(Ms):
  for rep in range(reps):
   for st in range(0,len(t),32):
    en=min(st+32,len(t));xc=x[st:en];h=1-t[st:en];B=en-st;rng=np.random.default_rng(np.random.SeedSequence([98731,d,rep,st]))
    ZT=rng.normal(size=(B,M*M,d));uterm=g(xc[:,None,:]+math.sqrt(2)*np.sqrt(h)[:,None,None]*ZT).mean(1)
    tau=rng.random((B,M));zr=rng.normal(size=(B,M,d));xr=xc[:,None,:]+math.sqrt(2)*np.sqrt(h[:,None]*tau)[:,:,None]*zr
    hc=h[:,None]*(1-tau);zz=rng.normal(size=(B,M,M,d));xt=xr[:,:,None,:]+math.sqrt(2)*np.sqrt(hc)[:,:,None,None]*zz
    dg=g(xt)-g(xr)[:,:,None];z1=(dg[:,:,:,None]*zz/np.sqrt(hc)[:,:,None,None]).mean(2);norm=np.linalg.norm(z1,axis=-1)
    raw[mi,rep,st:en]=uterm-.5*h*(norm**2).mean(1)
    for fi,R in enumerate(radii):
     proj[mi,fi,rep,st:en]=uterm-.5*h*(np.minimum(norm,R)**2).mean(1)
     omega[mi,fi,rep]+=(np.maximum(norm-R,0)**2).sum();viol[mi,fi,rep]+=(norm>R).sum()
   omega[mi,:,rep]/=len(t)*M;viol[mi,:,rep]/=len(t)*M
  re=np.linalg.norm(raw[mi]-truth,axis=1)/np.linalg.norm(truth)
  for fi,f in enumerate(fac):
   ie=np.linalg.norm(proj[mi,fi]-truth,axis=1)/np.linalg.norm(truth)
   rows.append(dict(d=d,n=2,M=M,radius_factor=float(f),radius=float(radii[fi]),raw_rel_l2_mean=float(re.mean()),ir_rel_l2_mean=float(ie.mean()),gain_fraction=float(1-ie.mean()/re.mean()),overshoot_energy=float(omega[mi,fi].mean()),violation_rate=float(viol[mi,fi].mean()),raw_rel_l2_sd=float(re.std(ddof=1)),ir_rel_l2_sd=float(ie.std(ddof=1))))
  print('completed M',M,'seconds',round(time.time()-start,2),flush=True)
  (a.out/'overshoot_36_results.json').write_text(json.dumps(rows,indent=2))
 O=[r['overshoot_energy'] for r in rows];G=[r['gain_fraction'] for r in rows];P=[r['violation_rate'] for r in rows]
 np.savez_compressed(a.out/'overshoot_36_raw.npz',t=t,x=x,A=A,truth=truth,truth128=truth128,M_values=Ms,radius_factors=fac,raw_value=raw,projected_value=proj,overshoot_by_rep=omega,violation_by_rep=viol)
 report=dict(source='fresh reproduction, not original Oct 4 samples or digitized points',d=d,n=2,reps=reps,n_domain=1000,n_boundary=200,chunk=32,M_values=Ms,radius_factors=fac.tolist(),radius_base=float(radii[0]),dtype='float64',reference='scaled Gauss-Laguerre 64 nodes',max_reference_64_128_abs_diff=float(np.max(np.abs(truth-truth128))),projection='child z only; final u not clipped',random_time='uniform; only n=2 value channel studied',time_denominator_floor=False,spearman_overshoot_gain=float(spearmanr(O,G).statistic),pearson_overshoot_gain=float(pearsonr(O,G).statistic),spearman_violation_gain=float(spearmanr(P,G).statistic),finite=bool(np.isfinite(raw).all() and np.isfinite(proj).all()),correlations_descriptive=True,elapsed_seconds=time.time()-start)
 (a.out/'overshoot_36_protocol.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
