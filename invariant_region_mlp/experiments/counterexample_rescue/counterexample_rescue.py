"""Hutzenthaler--Nguyen (2025), arXiv:2506.23969v1, Eq. (6).

Independent full-history MLP trees; float64; r=Uniform(0,1)^2.
PDE: u_t + 0.5 Delta u + ||grad u[1:]||_2 = 0; g(x)=|x[0]|.
All methods share every terminal/time/Brownian draw within a replicate.
The zero l=0 generator term is skipped identically for every method.
No clipping, antithetic sampling, or time floor in the raw baseline.

The deterministic x[1:] coordinates need not be stored: coefficients and
payoff do not depend on them. Their independent Brownian weights ARE sampled
and their noisy gradient estimates ARE propagated by the raw recursion.

Counts are scalar f/g evaluations PER REPLICA and per method (not batch calls).
The experiment traverses the full tree even for IR, to match nominal budgets.
A separately justified optimized IR implementation may prune zero corrections.
"""
from __future__ import annotations
import argparse, hashlib, json, math, platform, time
from dataclasses import dataclass
from pathlib import Path
import numpy as np

MODES = ('raw', 'recursive_subspace', 'recursive_ball', 'recursive_box', 'generator_only')

@dataclass
class Counts:
    terminal_samples: int = 0
    terminal_centers: int = 0
    driver_evals: int = 0
    random_times: int = 0
    normal_scalars: int = 0

class Solver:
    def __init__(self, d: int, m: int, seed: int, modes=MODES, beta: float=0.0):
        if d < 1 or m < 1: raise ValueError('d and m must be positive')
        self.d, self.m, self.modes, self.beta = d, m, tuple(modes), beta
        ss = np.random.SeedSequence(seed).spawn(3)
        self.times, self.active, self.inactive = [np.random.default_rng(s) for s in ss]
        self.counts = Counts()
        self.root_terminal = None
        self.max_ir_driver = 0.0
        self.min_remaining_time = 1.0

    def normals(self, b: int, k: int) -> np.ndarray:
        g = np.empty((b, k, self.d), dtype=np.float64)
        g[:, :, 0] = self.active.standard_normal((b,k))
        if self.d > 1: g[:, :, 1:] = self.inactive.standard_normal((b,k,self.d-1))
        self.counts.normal_scalars += b*k*self.d
        return g

    def finish(self, out):
        for j, mode in enumerate(self.modes):
            z = out[:,j,1:]
            if mode == 'recursive_subspace': z[:,1:] = 0.0
            elif mode == 'recursive_ball':
                norms = np.linalg.norm(z, axis=1, keepdims=True)
                z *= 1.0/np.maximum(1.0, norms)
            elif mode == 'recursive_box': np.clip(z, -1.0, 1.0, out=z)
            elif mode not in ('raw','generator_only'): raise ValueError(mode)
        return out

    def driver(self, child):
        z = child[:,:,1:]
        f = np.linalg.norm(z[:,:,1:], axis=-1)
        for j, mode in enumerate(self.modes):
            if mode == 'generator_only': f[:,j] = 0.0
        if self.beta: f = f + self.beta*np.sin(z[:,:,0])
        self.counts.driver_evals += child.shape[0]
        if self.beta == 0 and 'recursive_subspace' in self.modes:
            j=self.modes.index('recursive_subspace')
            self.max_ir_driver = max(self.max_ir_driver, float(np.max(np.abs(f[:,j]),initial=0)))
        return f

    def solve(self, n: int, h: np.ndarray, x1: np.ndarray, root=True):
        """h is remaining time 1-t, kept directly to avoid cancellation near T."""
        b = len(h); k_modes = len(self.modes)
        if n <= 0: return np.zeros((b,k_modes,self.d+1))
        if np.any(h <= 0): raise ValueError('h must be strictly positive')
        self.min_remaining_time=min(self.min_remaining_time,float(h.min()))
        k=self.m**n
        normal=self.normals(b,k)
        sh=np.sqrt(h)
        diff=np.abs(x1[:,None]+sh[:,None]*normal[:,:,0])-np.abs(x1[:,None])
        term=np.empty((b,self.d+1))
        term[:,0]=np.abs(x1)+diff.mean(axis=1)
        term[:,1:]=np.einsum('bi,bij->bj', diff, normal)/k/sh[:,None]
        self.counts.terminal_samples += b*k
        self.counts.terminal_centers += b
        if root: self.root_terminal=term.copy()
        out=np.repeat(term[:,None,:],k_modes,axis=1)
        del normal, diff
        for level in range(1,n):
            k=self.m**(n-level)
            r=self.times.random((b,k))
            # PCG64 can return exactly zero. Resample that event; no time floor.
            while np.any(r==0): r[r==0]=self.times.random(np.count_nonzero(r==0))
            r=r*r
            self.counts.random_times += b*k
            normal=self.normals(b,k)
            xp=x1[:,None]+sh[:,None]*np.sqrt(r)*normal[:,:,0]
            hp=h[:,None]*(1-r)
            high=self.solve(level,hp.ravel(),xp.ravel(),root=False)
            fd=self.driver(high).reshape(b,k,k_modes)
            del high
            if level > 1:
                low=self.solve(level-1,hp.ravel(),xp.ravel(),root=False)
                fd-=self.driver(low).reshape(b,k,k_modes)
                del low
            else:
                # F(V_0)=0 evaluated conceptually; don't alter any random stream.
                self.counts.driver_evals += b*k
            out[:,:,0] += np.mean(2*h[:,None,None]*np.sqrt(r)[:,:,None]*fd,axis=1)
            # (1/rho)*(W_R-W_t)/(R-t)=2 sqrt(h) G, exactly.
            out[:,:,1:] += 2*sh[:,None,None]*np.einsum('bik,bid->bkd',fd,normal)/k
        return self.finish(out)


def bootstrap_rmse(sq, rng, nboot=1000):
    ix=rng.integers(0,len(sq),(nboot,len(sq)))
    return np.sqrt(np.quantile(np.mean(sq[ix],axis=1),[.025,.975])).tolist()


def run_config(d,n,m,reps=512,batch=16,seed=20261004,beta=0.0):
    if reps % batch: raise ValueError('Use reps divisible by batch for this protocol')
    states=[]; terms=[]; counts=None; err=0.0; ir_driver=0.0; min_h=1.
    start=time.perf_counter()
    for block in range(reps//batch):
        # These seeds deliberately do NOT depend on d, for paired dimension sweeps.
        block_seed=int(np.random.SeedSequence([seed,n,m,block]).generate_state(1)[0])
        s=Solver(d,m,block_seed,beta=beta)
        out=s.solve(n,np.ones(batch),np.zeros(batch))
        root=s.root_terminal
        if beta==0:
            p=root.copy();p[:,2:]=0.
            err=max(err,float(np.max(np.abs(out[:,1,:]-p))))
        ir_driver=max(ir_driver,s.max_ir_driver);min_h=min(min_h,s.min_remaining_time)
        states.append(out);terms.append(root)
        c={k:v/batch for k,v in vars(s.counts).items()}
        if counts is not None and c!=counts: raise AssertionError('counts vary by block')
        counts=c
    states=np.concatenate(states);term=np.concatenate(terms)
    elapsed=time.perf_counter()-start
    post=states[:,0,:].copy();post[:,2:]=0
    reduced=term.copy();reduced[:,2:]=0
    all_states=np.concatenate((states,post[:,None,:],reduced[:,None,:]),axis=1)
    modes=MODES+('final_only','reduced_mc')
    truth=np.zeros(d+1);truth[0]=math.sqrt(2/math.pi)
    sq=np.sum((all_states-truth[None,None,:])**2,axis=2)
    usq=(all_states[:,:,0]-truth[0])**2
    rows=[];rng=np.random.default_rng(829421)
    for j,mode in enumerate(modes):
        rows.append(dict(d=d,n=n,m=m,reps=reps,mode=mode,
            state_mse=float(sq[:,j].mean()),state_rmse=float(np.sqrt(sq[:,j].mean())),
            state_rmse_boot95=bootstrap_rmse(sq[:,j],rng),
            value_rmse=float(np.sqrt(usq[:,j].mean())),
            value_bias=float(all_states[:,j,0].mean()-truth[0]),
            value_mean=float(all_states[:,j,0].mean()),
            gradient_rmse=float(np.sqrt(np.sum(all_states[:,j,1:]**2,axis=1).mean())),
            finite_fraction=float(np.isfinite(all_states[:,j,:]).all(axis=1).mean()),
            full_tree_counts_per_replica=counts,
            ir_exact_mse=(4-2/math.pi)/(m**n),
            generator_only_exact_mse=(d+3-2/math.pi)/(m**n)))
    compact=dict(d=d,n=n,m=m,elapsed_seconds_all_methods=elapsed,
        ir_terminal_max_abs_diff=err,max_ir_generator_abs=ir_driver,
        min_remaining_time=min_h,nominal_counts_per_method=counts)
    # Full per-replica states are saved; compression benefits exact-zero IR coords.
    return rows, compact, all_states, sq, usq


def tests():
    checks=[]
    for d in [1,2,7]:
        for n in [1,2,3]:
            s=Solver(d,2,77+n)
            y=s.solve(n,np.array([1.,.75,.02]),np.array([0.,.3,-1.]))
            term=s.root_terminal.copy();term[:,2:]=0
            assert np.array_equal(y[:,1,:],term)
            assert np.array_equal(y[:,4,:],s.root_terminal)
            assert s.max_ir_driver==0
            checks.append(f'collapse d={d}, n={n}')
    # Nontrivial ACTIVE generator: beta sin(z_1). No collapse; dimension reduction
    # must still commute with the recursion when projected Brownian streams match.
    for d in [2,7,50]:
        full=Solver(d,2,12345,beta=.3).solve(3,np.ones(3),np.array([.2,.4,-.3]))[:,1,:2]
        low=Solver(1,2,12345,beta=.3).solve(3,np.ones(3),np.array([.2,.4,-.3]))[:,0,:2]
        assert np.allclose(full,low,rtol=0,atol=5e-14)
        checks.append(f'active-generator commuting d={d}')
    # Norm conventions audit: a ball projector is Euclidean nonexpansive,
    # not necessarily max-norm nonexpansive. No automatic dimension-free theorem.
    return {'checks':checks,'n_passed':len(checks),'status':'passed'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--dims',default='2,5,10,20,50,100,200,500,1000')
    ap.add_argument('--levels',default='1,2,3,4')
    ap.add_argument('--reps',type=int,default=512)
    ap.add_argument('--batch',type=int,default=16)
    ap.add_argument('--seed',type=int,default=20261004)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    test=tests();(a.out/'tests.json').write_text(json.dumps(test,indent=2))
    print('Tests passed:',test['n_passed'],flush=True)
    metadata=dict(seed=a.seed,reps=a.reps,batch=a.batch,dims=a.dims,levels=a.levels,
        numpy=np.__version__,python=platform.python_version(),platform=platform.platform(),
        source='Hutzenthaler and Nguyen (2025), arXiv:2506.23969v1, Eq. (6)',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        coupling='same independent time/active/inactive streams within replicas; active stream shared across d',
        raw_notes='no projection, no time flooring; l=0 zero term skipped in every mode',
        timing_notes='elapsed is combined vectorized time for all methods, NOT separate speed benchmarking')
    (a.out/'metadata.json').write_text(json.dumps(metadata,indent=2))
    rows=[];audit=[]
    for n in map(int,a.levels.split(',')):
        for d in map(int,a.dims.split(',')):
            label=f'd{d}_n{n}_m{n}'
            rowpath=a.out/(label+'.json');rawpath=a.out/(label+'.npz')
            if rowpath.exists() and rawpath.exists():
                record=json.loads(rowpath.read_text());r,c=record['rows'],record['audit']
            else:
                r,c,st,sq,usq=run_config(d,n,n,a.reps,a.batch,a.seed)
                np.savez_compressed(rawpath,states=st,state_squared_errors=sq,value_squared_errors=usq,modes=np.array(MODES+('final_only','reduced_mc')))
                rowpath.write_text(json.dumps({'rows':r,'audit':c},indent=2))
            rows.extend(r);audit.append(c)
            print(label, 'raw',round(r[0]['state_rmse'],6),'IR',round(r[1]['state_rmse'],6),'secs',round(c['elapsed_seconds_all_methods'],2),flush=True)
            (a.out/'summary.json').write_text(json.dumps({'metadata':metadata,'rows':rows,'audit':audit},indent=2))

if __name__=='__main__':main()
