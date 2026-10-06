#!/usr/bin/env python3
"""Render recovered and reproduced evidence; no implicit PDE reruns.

Every plot is a distinct Matplotlib figure. Existing figure assets are untouched.
Data provenance and limitations are in data/restored_evidence/PROVENANCE.json.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'restored_evidence'
OUT=ROOT/'figures'/'restored'
OUT.mkdir(parents=True,exist_ok=True)

def load(name):
    return json.loads((DATA/name).read_text())

def canvas():
    fig=plt.figure(figsize=(4.55,3.25))
    ax=fig.add_subplot(111)
    ax.tick_params(labelsize=10)
    return fig,ax

def save(fig,ax,name,legend=True,loc='best',ncol=1):
    ax.set_axisbelow(True)
    ax.grid(True,axis='y',alpha=.22)
    if legend:
        ax.legend(loc=loc,fontsize=8.2,frameon=False,ncol=ncol)
    fig.tight_layout(pad=.6)
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=200,bbox_inches='tight')
    plt.close(fig)

O=load('overshoot_36_results.json')
for metric,name,xlabel in [('overshoot_energy','overshoot_energy_reproduction','Invariant-region overshoot energy'),('violation_rate','violation_gain_reproduction','Child gradient violation rate (%)')]:
    fig,ax=canvas()
    for M,mk in zip([2,3,4,6,8,10],['o','s','^','D','v','P']):
        rows=sorted([r for r in O if r['M']==M],key=lambda r:r[metric])
        mult=100 if metric=='violation_rate' else 1
        ax.plot([mult*r[metric] for r in rows],[100*r['gain_fraction'] for r in rows],marker=mk,markersize=4.5,linewidth=1.3,label=f'$M={M}$')
    ax.set_xlabel(xlabel,fontsize=10.5)
    ax.set_ylabel('Final relative-$L_2$ error reduction (%)',fontsize=10)
    ax.set_ylim(0,100)
    save(fig,ax,name,loc='lower right',ncol=2)

fund=json.loads((DATA/'source_archive'/'finance_higher_levels.json').read_text())
fig,ax=canvas()
settings=[(4,2),(4,3),(5,2)]
for a,label,mk in [('baseline','Raw','o'),('hard','Joint ellipsoid','s'),('box','Coordinate box','^')]:
    vals=[next(r['mae'] for r in fund if r['n']==n and r['M']==M and r['mode']==a) for n,M in settings]
    ax.plot(range(3),vals,marker=mk,markersize=6,linewidth=1.4,label=label)
ax.set_xticks(range(3),[r'$(4,2)$',r'$(4,3)$',r'$(5,2)$'])
ax.set_xlabel('Depth / sampling-base setting $(n,M)$',fontsize=10)
ax.set_ylabel('Value MAE',fontsize=11)
ax.set_ylim(bottom=0)
save(fig,ax,'funding_deep_geometry')

# Keep all budget points, not only the cumulative-best envelope in the old PNG.
fundall=json.loads((DATA/'source_archive'/'finance_selected_results.json').read_text())
fig,ax=canvas()
for a,label,mk in [('baseline','Raw','o'),('hard','Joint ellipsoid','s')]:
    rs=sorted([r for r in fundall if r['mode']==a],key=lambda r:r['f_evals'])
    ax.plot([r['f_evals'] for r in rs],[r['mae'] for r in rs],marker=mk,markersize=5,linestyle='none',label=label)
ax.set(xscale='log',yscale='log')
ax.set_xlabel('Nonlinear-generator evaluations',fontsize=10.5)
ax.set_ylabel('Recorded MAE at each setting',fontsize=10.5)
save(fig,ax,'funding_all_budget_points')

sur=json.loads((DATA/'source_archive'/'surrogate_quality_sweep_multibudget.json').read_text())
for M in [1,2,4]:
    fig,ax=canvas()
    for a,label,mk in [('raw','Raw defect MLP','o'),('hard_total','Total-state IR','s'),('defect_clip','Fixed defect clip','^')]:
        rs=sorted([r for r in sur if r['M']==M and r['mode']==a],key=lambda r:r['surrogate_rel_l2'])
        ax.errorbar([100*r['surrogate_rel_l2'] for r in rs],[100*r['rel_l2_mean'] for r in rs],yerr=[100*r['rel_l2_std'] for r in rs],marker=mk,markersize=4.5,linewidth=1.2,capsize=2,label=label)
    ax.set_yscale('symlog',linthresh=.1)
    ax.set_yticks([0,.1,1,10,100],['0','0.1','1','10','100'])
    ax.set_ylim(-.02,150)
    ax.set_xticks([0,25,50,75,100])
    ax.set_xlabel('Controlled surrogate error (%)',fontsize=10.5)
    ax.set_ylabel('Corrected relative-$L_2$ error (%)',fontsize=10)
    save(fig,ax,f'surrogate_quality_M{M}',loc='lower right')

B=load('batch_additional_report_values.json')
fig,ax=canvas()
for a,label,mk,ls,ms in [('raw','Raw','o','-',7),('mean_preserving','Mean-preserving','s','--',5),('samplewise','Samplewise IR','^','-',7),('batch','Batch-IR','x',':',5)]:
    ax.plot([r['d'] for r in B['nw_m3']],[r[a] for r in B['nw_m3']],marker=mk,fillstyle='none',markersize=ms,linestyle=ls,linewidth=1.4,label=label)
ax.set_xlabel('Dimension $d$',fontsize=11)
ax.set_ylabel('Value MAE',fontsize=11)
ax.set_xticks([100,200,300]);ax.set_ylim(bottom=0)
save(fig,ax,'batch_neufeld_M3')
fig,ax=canvas()
ax.plot([r['M'] for r in B['funding_mp_feasible']],[100*r['feasible_fraction'] for r in B['funding_mp_feasible']],marker='o',linewidth=1.4)
ax.set_xlabel('Funding sampling base $M$',fontsize=11)
ax.set_ylabel('Runs with feasible batch mean (%)',fontsize=10)
ax.set_xticks([10,20,40]);ax.set_ylim(-3,103)
save(fig,ax,'mean_preserving_feasibility',legend=False)

neg=load('negative_controls_extracted.json')
for b,name,ykey,factor in [('Allen-Cahn 100D','negative_allen_batch','rel_mae',100),('Counterparty credit risk 100D','negative_credit_batch','mae',1),('Linear convection-diffusion 100D','negative_linear_batch','mae',1)]:
    fig,ax=canvas()
    for a,label,mk,ls,ms in [('raw','Raw','o','-',8),('samplewise','Samplewise IR','s','--',6),('batch','Batch-IR','x',':',5)]:
        rs=[r for r in neg if r['benchmark']==b and r['method']==a]
        xx=range(len(rs)) if 'credit' in b else [r['M'] for r in rs]
        ax.plot(xx,[factor*r[ykey] for r in rs],marker=mk,markersize=ms,fillstyle='none',linestyle=ls,linewidth=1.3,label=label)
    if 'credit' in b:
        ax.set_xticks(range(len(rs)),[f"$({r['n']},{r['M']})$" for r in rs]);ax.set_xlabel('Setting $(n,M)$',fontsize=11)
    else:
        ax.set_xticks([r['M'] for r in rs]);ax.set_xlabel('Sampling base $M$',fontsize=11)
    ax.set_ylabel('Relative value MAE (%)' if factor==100 else 'Value MAE',fontsize=11)
    ax.set_ylim(bottom=0)
    save(fig,ax,name)
fig,ax=canvas();rs=[r for r in neg if 'Linear ' in r['benchmark'] and r['method']=='batch']
for key,label,mk,ls,ms in [('violation_rate','Gradient violation rate','o','-',7),('batch_activation_rate','Batch activation rate','x','--',5)]:
    ax.plot([r['M'] for r in rs],[100*r[key] for r in rs],marker=mk,fillstyle='none',markersize=ms,linestyle=ls,linewidth=1.4,label=label)
ax.set_xlabel('Sampling base $M$',fontsize=11);ax.set_ylabel('Recorded rate (%)',fontsize=11)
ax.set_xticks([4,8,16]);ax.set_ylim(0,105)
save(fig,ax,'negative_linear_activity')

E=load('elliptic_report_values.json')
for key,name,ylabel in [('driver_abs_error','elliptic_driver_diagnostic','Mean absolute quadratic-driver error'),('terminal_loss','elliptic_terminal_diagnostic','Mean squared terminal BSDE residual')]:
    fig,ax=canvas();rs=E['seed0'];xx=np.arange(4)
    ax.bar(xx,[r[key] for r in rs])
    ax.set_xticks(xx,['No scale','Joint ball','Box','Upstream\n$/d$'])
    ax.set_yscale('log');ax.set_ylabel(ylabel,fontsize=10)
    for i,r in enumerate(rs):
        ax.annotate(f"{r[key]:.4g}",(i,r[key]),xytext=(0,4),textcoords='offset points',ha='center',fontsize=8)
    ax.set_ylim(top=max(r[key] for r in rs)*2.1)
    save(fig,ax,name,legend=False)

print(f'Wrote {len(list(OUT.glob("*.pdf")))} new plot PDFs and PNGs to {OUT}; no old figure was changed.')
