"""Typeset archived, verified results. Does not run new PDE experiments.

Original PDF assets remain byte-for-byte in figures/originals/. The three
counterexample SVG plots retain the exact published SVG vertex coordinates;
only text/legend layout is changed. Those vertices are graphics, not raw data.
"""
from pathlib import Path
import json
import math
import html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter
import cairosvg

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / 'figures'
D = ROOT / 'data'
F.mkdir(exist_ok=True)
D.mkdir(exist_ok=True)

# Archived working-draft summaries; these arrays do not imply unreported precision.
DATA = {
    'hjb': {'d':[100,120,140,160],
        'heuristic':[2.080,2.277,2.414,2.492], 'exact_ball':[.561,.555,.547,.565],
        'uniform_ball':[.634,.641,.635,.643], 'box':[1.110,1.117,1.121,1.121],
        'generator_raw':[3589.1,4709.2,7687.4,10009.1],
        'generator_ir':[39.75,41.03,42.44,43.61],
        'variance_raw':[19.261,30.352,46.107,59.050],
        'variance_ir':[.223,.243,.223,.195]},
    'nw': {'d':[100,200,300], 'raw_m2':[.05687,.04989,.03813],
        'ir_m2':[.01888,.01505,.01887], 'raw_m3':[.01637,.01103,.01397],
        'ir_m3':[.00469,.00580,.00668]},
    'batch_hjb': {'d':[20,40], 'raw':[314.136,1034.262],
        'sample':[.2412,.1993], 'batch':[.0966,.0833]},
    'batch_funding': {'M':[10,20,40], 'raw':[1.4042,.8848,.5221],
        'sample':[.5102,.2965,.1993], 'mp':[1.4042,.3350,.1773],
        'batch':[.4807,.2539,.1486]},
    'allen_cahn': {'d':100, 'n':3, 'paired_seeds':100, 'M':[2,3,4],
        'raw_relative_mae_percent':[3.30,1.96,1.42],
        'ir_relative_mae_percent':[3.30,1.96,1.42], 'violation_percent':[0,0,0]},
    'primary_rmse': {'reps':512, 'depth':[1,2,3,4], 'N':[1,4,27,256],
        'exact':[1.833952,.916976,.352944,.114622],
        'empirical':[1.797422,.945043,.356391,.117464]},
    'independent_mse': {'reps':65536, 'N':[1,4,27,256],
        'empirical':[3.380034450498047,.8389663936455724,.12525171685975256,.013202442726108567],
        'mean_se':[.04056729388335935,.006264968739407201,.0006716141939526456,.00006603391824740008],
        'ci95':[[3.3005225544866628,3.4595463465094314],[.8266870549163342,.8512457323748105],
                [.12393535303960537,.12656808067989975],[.013073016246343664,.01333186920587347]]}
}
(D/'figure_data.json').write_text(json.dumps(DATA, indent=2)+'\n')

# Font sizes are specified in physical points at the intended column width.
# No custom colour palette or matplotlib style is used.
def newfig(height=2.55):
    fig, ax = plt.subplots(figsize=(3.25, height))
    ax.tick_params(labelsize=8)
    ax.xaxis.label.set_size(8.5)
    ax.yaxis.label.set_size(8.5)
    return fig, ax

def finish(fig, ax, name, legend=True, loc='best', ncol=1):
    if legend:
        ax.legend(fontsize=7.5, loc=loc, framealpha=.95, ncol=ncol,
                  handlelength=1.6, borderpad=.4, labelspacing=.25)
    fig.tight_layout(pad=.7)
    fig.savefig(F/f'{name}.pdf', metadata={'Creator':'Archived-data figure restoration'})
    plt.close(fig)

h=DATA['hjb']
fig,ax=newfig()
for k,lab,m in [('heuristic','Heuristic','o'),('exact_ball','Exact-matrix ball','s'),
                ('uniform_ball',r'Uniform $\sqrt{15}$ ball','^'),('box','Coordinate box','D')]:
    ax.plot(h['d'], h[k], marker=m, markersize=3.8, linewidth=1.25, label=lab)
ax.set(xlabel='Dimension $d$',ylabel='Mean relative $L_2$ error',xticks=h['d'])
finish(fig,ax,'hjb_corrected',loc='center right')
for a,b,yl,nm in [('generator_raw','generator_ir','Generator MSE','hjb_generator_mse'),
                  ('variance_raw','variance_ir','Variance of level correction','hjb_correction_variance')]:
    fig,ax=newfig()
    ax.plot(h['d'],h[a],'o-',markersize=4,linewidth=1.25,label='Raw')
    ax.plot(h['d'],h[b],'s--',markersize=4,linewidth=1.25,label='Samplewise IR')
    ax.set(xlabel='Dimension $d$',ylabel=yl,yscale='log',xticks=h['d'])
    finish(fig,ax,nm,loc='center right')

nw=DATA['nw']; fig,ax=newfig(); x=np.arange(3);w=.2
for i,(k,lab) in enumerate([('raw_m2','Raw, $M=2$'),('ir_m2','IR, $M=2$'),
                            ('raw_m3','Raw, $M=3$'),('ir_m3','IR, $M=3$')]):
    ax.bar(x+(i-1.5)*w,nw[k],width=w,label=lab)
ax.set(xlabel='Dimension $d$',ylabel='MAE',xticks=x,xticklabels=nw['d'],ylim=(0,.072))
finish(fig,ax,'neufeld_wu',loc='upper right',ncol=2)

bh=DATA['batch_hjb']; fig,ax=newfig()
for k,lab,m in [('raw','Raw / mean-preserving','o'),('sample','Samplewise IR','s'),('batch','Batch-IR','^')]:
    ax.plot(bh['d'],bh[k],marker=m,markersize=4,linewidth=1.25,label=lab)
ax.set(xlabel='Dimension $d$',ylabel='Relative $L_2$ error',yscale='log',xticks=bh['d'])
finish(fig,ax,'batch_hjb',loc='center right')
bf=DATA['batch_funding'];fig,ax=newfig()
for k,lab,m in [('raw','Raw','o'),('sample','Samplewise IR','s'),('mp','Mean-preserving','D'),('batch','Batch-IR','^')]:
    ax.plot(bf['M'],bf[k],marker=m,markersize=4,linewidth=1.25,label=lab)
ax.set(xlabel='Samples parameter $M$',ylabel='MAE',xticks=bf['M'],ylim=(0,1.52))
finish(fig,ax,'batch_funding',loc='upper right')

pr=DATA['primary_rmse'];fig,ax=newfig()
ax.plot(pr['depth'],pr['exact'],'o-',markersize=4,label='Exact RMSE')
ax.plot(pr['depth'],pr['empirical'],'s--',markersize=4,label='512-replica RMSE')
ax.set(xlabel='Depth $n=m$',ylabel='Full-state RMSE',yscale='log',xticks=pr['depth'])
finish(fig,ax,'counterexample_exact_rmse',loc='upper right')

im=DATA['independent_mse'];fig,ax=newfig()
ns=np.array(im['N']); exact=(4-2/math.pi)/ns; yy=np.array(im['empirical']); ci=np.array(im['ci95'])
ax.plot(ns,exact,'--',linewidth=1.25,label=r'Exact $(4-2/\pi)/N$')
ax.errorbar(ns, yy, yerr=np.stack([yy-ci[:,0],ci[:,1]-yy]),fmt='o',
            capsize=2.5,markersize=4,label='65,536-replica MSE')
ax.set(xscale='log',yscale='log',xlabel='Terminal samples $N$',ylabel='Full-state MSE')
ax.set_xticks(ns);ax.xaxis.set_major_formatter(ScalarFormatter())
finish(fig,ax,'counterexample_independent_mse',loc='lower left')

ac=DATA['allen_cahn'];fig,ax=newfig()
ax.plot(ac['M'],ac['raw_relative_mae_percent'],'o-',markersize=7,linewidth=1.25,label='Raw MLP')
ax.plot(ac['M'],ac['ir_relative_mae_percent'],'x--',markersize=6,linewidth=1.25,label='Samplewise IR')
ax.set(xlabel='Samples parameter $M$',ylabel='Relative MAE (%)',xticks=ac['M'],ylim=(1.1,3.65))
ax.text(.03,.07,'0% constraint violations\nat every tested budget', transform=ax.transAxes, fontsize=8)
finish(fig,ax,'allen_cahn_control',loc='upper right')

# Exact graphic vertices from repository SVGs. Never treat these as raw numeric observations.
common_x=[92,186.8,258.5,330.2,425,496.8,568.5,663.3,735]
SVGDATA={
 'counterexample_state_full': {
 'source':'counterexample_state_rmse.svg', 'ylabel':'Full-state RMSE',
 'ticks':[[364.47886031402726,'1e0'],[289.6671515849008,'1e1'],[214.8554428557744,'1e2'],[140.04373412664802,'1e3'],[65.23202539752162,'1e4']],
 'curves':[
 ['Original MLP',[358.7,316.1,281.7,244.7,197.6,162.3,127.6,82.2,48.0]],
 ['Final-only subspace',[366.5,335.2,309.6,284.4,252.2,228.3,205.0,174.7,151.7]],
 ['Generator-only subspace',[393.3,385.9,377.3,367.5,353.4,342.6,331.7,316.9,305.7]],
 ['Certified unit ball',[372.7,357.8,354.5,353.7,353.1,352.9,352.9,352.8,352.8]],
 ['Recursive subspace IR',[398]*9]]},
 'counterexample_value_rmse': {
 'source':'counterexample_value_rmse.svg', 'ylabel':'RMSE of u(0,0)',
 'ticks':[[306.2111709262824,'1e0'],[206.85203543098876,'1e1'],[107.4928999356951,'1e2']],
 'curves':[
 ['Original / final-only',[333.3,283.6,251.,218.9,178.1,147.8,117.6,78.1,48.]],
 ['Coordinate box',[338.8,297.1,272.8,252.4,227.9,210.8,194.3,173.3,157.7]],
 ['Certified unit ball',[347.4,317.1,309.8,307.4,306.,305.5,305.3,305.1,305.]],
 ['Subspace IR / reduced MC',[398]*9]]},
 'counterexample_depth_sweep': {
 'source':'counterexample_depth_sweep.svg', 'ylabel':'Full-state RMSE',
 'ticks':[[388.1344582438935,'1e0'],[324.902276557692,'1e1'],[261.6700948714905,'1e2'],[198.437913185289,'1e3'],[135.2057314990875,'1e4'],[71.97354981288595,'1e5']],
 'curves':[
 ['n = m = 1',[368.4,360.4,353.5,345.7,333.9,324.4,315.4,302.7,293.3]],
 ['n = m = 2',[371.8,346.6,328.1,309.1,283.8,264.8,245.6,220.5,201.3]],
 ['n = m = 3',[383.2,347.2,318.2,286.9,247.1,217.2,187.9,149.5,120.6]],
 ['n = m = 4',[398.,359.7,324.1,279.6,221.8,180.,139.6,87.,48.]]]}
}
(D/'source_svg_vertices.json').write_text(json.dumps(SVGDATA,indent=2)+'\n')

# Retain paths; move legends below axes so none occludes a curve.
for name,dd in SVGDATA.items():
    # 26px at 3.25in / 760px is approximately 8pt in the manuscript.
    ss=['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="590" viewBox="0 0 760 590">',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.axis{stroke:#111;stroke-width:2}.grid{stroke:#bbb;stroke-width:.7;stroke-dasharray:3 4}.line{fill:none;stroke:#111;stroke-width:3}.m{fill:#fff;stroke:#111;stroke-width:2}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<line class="axis" x1="92" y1="410" x2="746" y2="410"/><line class="axis" x1="92" y1="36" x2="92" y2="410"/>']
    for x,l in [(92,'2'),(258.52193691149506,'10'),(496.76096845574756,'100'),(735,'1000')]:
        ss.append(f'<line class="grid" x1="{x}" y1="36" x2="{x}" y2="410"/><text x="{x}" y="443" text-anchor="middle" font-size="26">{l}</text>')
    for y,l in dd['ticks']:
        ss.append(f'<line class="grid" x1="92" y1="{y}" x2="746" y2="{y}"/><text x="82" y="{y+8}" text-anchor="end" font-size="26">{l}</text>')
    ss.append(f'<text x="414" y="479" text-anchor="middle" font-size="28">Ambient dimension d</text><text transform="translate(22 224) rotate(-90)" text-anchor="middle" font-size="28">{html.escape(dd["ylabel"])}</text>')
    dashes=['','10 5','3 5','12 4 3 4','16 5']
    # Main path geometry is exactly that of the archived SVG file.
    for k,(lab,ys) in enumerate(dd['curves']):
        pts=' '.join(f'{x:.1f},{y:.1f}' for x,y in zip(common_x,ys))
        dash=f' stroke-dasharray="{dashes[k]}"' if dashes[k] else ''
        ss.append(f'<polyline class="line" points="{pts}"{dash}/>')
        for x,y in zip(common_x,ys):ss.append(f'<circle class="m" cx="{x:.1f}" cy="{y:.1f}" r="{4+k*.2}"/>')
        col=k%2;row=k//2;lx=50+col*355;ly=516+row*30
        ss.append(f'<line class="line" x1="{lx}" y1="{ly}" x2="{lx+38}" y2="{ly}"{dash}/><text x="{lx+48}" y="{ly+7}" font-size="22">{html.escape(lab)}</text>')
    ss.append('</svg>')
    svg='\n'.join(ss); (F/f'{name}.svg').write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(F/f'{name}.pdf'))

print('Restored/re-typeset 12 vector graphics. Original PDFs remain in figures/originals.')

# Re-typeset the archived vector geometry at column size. This is a graphics
# transformation: x/y below are SVG vertices, NOT recovered simulation values.
# Default matplotlib line colours reproduce the familiar coloured-curve layout.
for name, dd in SVGDATA.items():
    fig=plt.figure(figsize=(3.25,2.65))
    ax=fig.add_axes([.17,.32,.80,.65])
    marker_seq=['o','s','^','D','v']
    line_seq=['-','--',':','-.','--']
    for k,(lab,ys) in enumerate(dd['curves']):
        ax.plot(common_x,ys,marker=marker_seq[k],markersize=3.0,
                linestyle=line_seq[k],linewidth=1.1,label=lab)
    ax.set_xlim(75,756); ax.set_ylim(412,30)
    ax.set_xticks([92,258.52193691149506,496.76096845574756,735],['2','10','100','1000'])
    ax.set_yticks([v for v,_ in dd['ticks']],
                 [r'$10^{'+lab.split('e')[1]+'}$' for _,lab in dd['ticks']])
    ax.tick_params(labelsize=7.5,pad=2)
    ax.set_xlabel('Ambient dimension $d$',fontsize=8.3,labelpad=2)
    ax.set_ylabel(dd['ylabel'].replace('u(0,0)','$u(0,0)$'),fontsize=8.3,labelpad=2)
    ax.grid(True,alpha=.22,linewidth=.5)
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.015),ncol=2,
               fontsize=6.5,handlelength=1.8,columnspacing=.6,labelspacing=.25,
               handletextpad=.45,frameon=False)
    fig.savefig(F/f'{name}.pdf',metadata={'Creator':'Re-typeset archived SVG geometry; no new PDE runs'})
    plt.close(fig)
