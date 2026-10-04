"""Recreate counterexample figures from results/summary_all.csv."""
from pathlib import Path
import argparse, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--summary',type=Path,default=Path('../../results/counterexample_rescue/summary_all.csv'))
    ap.add_argument('--validation',type=Path,default=Path('../../results/counterexample_rescue/independent_validation.json'))
    ap.add_argument('--out',type=Path,default=Path('../../results/counterexample_rescue/figures_png'))
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.summary)
    def series(n,mode): return df[(df.n==n)&(df['mode']==mode)].sort_values('d')
    for metric,ylabel,name,modes in [
      ('state_rmse','RMSE of the full (value, gradient) state','counterexample_state_rmse.png',
       [('raw','Original full-history MLP'),('final_only','Final-only subspace'),('generator_only','Generator-only subspace'),('recursive_ball','Certified unit ball'),('recursive_subspace','Recursive subspace IR')]),
      ('value_rmse','RMSE of u(0,0)','counterexample_value_rmse.png',
       [('raw','Original MLP'),('recursive_box','Coordinate box'),('recursive_ball','Certified unit ball'),('recursive_subspace','Recursive subspace IR')])]:
        fig,ax=plt.subplots(figsize=(6.6,4.3))
        for mode,label in modes:
            r=series(3,mode); ax.loglog(r.d,r[metric],marker='o',label=label)
        ax.set_xlabel('Ambient dimension d'); ax.set_ylabel(ylabel); ax.grid(True,which='both',alpha=.2); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(a.out/name,dpi=200); plt.close(fig)
    val=json.loads(a.validation.read_text())['moments']; N=np.array([x['N'] for x in val]); emp=np.array([x['empirical_mse'] for x in val]); se=np.array([x['mean_se'] for x in val])
    fig,ax=plt.subplots(figsize=(6.2,3.8)); ax.errorbar(N,emp,yerr=1.96*se,fmt='o',capsize=4,label='empirical'); ax.plot(N,(4-2/np.pi)/N,'--',label='exact')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('N'); ax.set_ylabel('Full-state MSE'); ax.grid(True,which='both',alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(a.out/'counterexample_exact_mse.png',dpi=200); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.2,3.8))
    for n in [1,2,3,4]:
        r=series(n,'raw'); ax.loglog(r.d,r.state_rmse,marker='o',label=f'n=m={n}')
    ax.set_xlabel('Ambient dimension d'); ax.set_ylabel('Full-state RMSE'); ax.grid(True,which='both',alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(a.out/'counterexample_depth_sweep.png',dpi=200)
if __name__=='__main__': main()
