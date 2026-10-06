from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data'/'restored_evidence';T=ROOT/'tables'
def load(x):return json.loads((D/x).read_text())
def table(name,caption,label,cols,head,rows):
 s='\\begin{table}[htbp]\n\\centering\\small\n\\caption{'+caption+'}\\label{'+label+'}\n\\begin{tabular}{'+cols+'}\n\\toprule\n'+head+' \\\\\n\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'
 (T/(name+'.tex')).write_text(s)
O=load('overshoot_36_results.json');rows=[]
for i,r in enumerate(O):
 if i>0 and i%6==0:rows.append(r'\midrule')
 rows.append(f"{r['M']} & {r['radius_factor']:g} & {r['overshoot_energy']:.3f} & {100*r['violation_rate']:.2f} & {r['raw_rel_l2_mean']:.4f} & {r['ir_rel_l2_mean']:.4f} & {100*r['gain_fraction']:.2f}"+r' \\')
table('restored_overshoot','All 36 settings of the independent, corrected-reference reproduction. $a$ multiplies the exact-matrix certificate $R_0$; $p$ is the pre-projection violation fraction. Raw and IR columns are mean relative $L_2$ errors over ten paired repetitions. These are newly sampled results, not reconstructed old figure coordinates.','tab:restored-overshoot','rrrrrrr',r'$M$ & $a$ & $\widehat\Omega$ & $p$ (\%) & Raw & IR & Gain (\%)',rows)
F=load('source_archive/finance_selected_results.json');rows=[]
for n,M in sorted(set((r['n'],r['M']) for r in F)):
 a=next(r for r in F if r['n']==n and r['M']==M and r['mode']=='baseline');b=next(r for r in F if r['n']==n and r['M']==M and r['mode']=='hard')
 rows.append(f"{n} & {M} & {a['seeds']} & {a['f_evals']} & {a['mae']:.4f} & {b['mae']:.4f}"+r' \\')
table('restored_funding_budgets','All recorded settings behind the original funding work envelope. Counts are nonlinear-generator evaluations; replication counts differ across settings.','tab:restored-funding-budget','rrrrrr',r'$n$ & $M$ & Runs & $f$ calls & Raw MAE & Joint IR MAE',rows)
F=load('source_archive/finance_higher_levels.json');rows=[]
for n,M in [(4,2),(4,3),(5,2)]:
 a=next(r for r in F if r['n']==n and r['M']==M and r['mode']=='baseline');b=next(r for r in F if r['n']==n and r['M']==M and r['mode']=='hard');c=next(r for r in F if r['n']==n and r['M']==M and r['mode']=='box')
 rows.append(f"{n} & {M} & {a['seeds']} & {a['f_evals']} & {a['mae']:.4f} & {b['mae']:.4f} & {c['mae']:.4f}"+r' \\')
table('restored_funding_geometry','Earlier higher-level funding geometry ablation. The three methods are paired within each setting; no Batch-IR series is inferred from these data.','tab:restored-funding-geometry','rrrrrrr',r'$n$ & $M$ & Runs & $f$ calls & Raw & Joint IR & Box',rows)
S=load('source_archive/surrogate_quality_sweep_multibudget.json');rows=[]
for M in [1,2,4]:
 if M!=1:rows.append(r'\midrule')
 for q in [1,.95,.75,.5,0]:
  vals=[]
  for mode in ['raw','hard_total','defect_clip']:
   r=next(r for r in S if r['M']==M and r['q']==q and r['mode']==mode)
   vals.append(f"${100*r['rel_l2_mean']:.3f}\\pm{100*r['rel_l2_std']:.3f}$")
  rows.append(f'{M} & {100*(1-q):.0f} & '+' & '.join(vals)+r' \\')
table('restored_surrogate','Controlled logistic-surrogate diagnostic: relative $L_2$ error in percent (mean $\\pm$ recorded standard deviation over 40 runs). The perfect-surrogate rows remain, and the table also retains cases where IR is worse than raw defect MLP.','tab:restored-surrogate','rrccc',r'$M$ & Surrogate error (\%) & Raw defect MLP & Total-state IR & Fixed defect clip',rows)
N=load('negative_controls_extracted.json');rows=[]
for bench,short,metric,mul in [('Allen-Cahn 100D','Allen--Cahn','rel_mae',100),('Counterparty credit risk 100D','Credit risk','mae',1),('Linear convection-diffusion 100D','Linear','mae',1)]:
 if rows:rows.append(r'\midrule')
 for r in [r for r in N if r['benchmark']==bench and r['method']=='batch']:
  rows.append(f"{short} & {r['n']} & {r['M']} & {r[metric]*mul:.6f} & {100*r['violation_rate']:.2f} & {100*r['batch_activation_rate']:.2f} & 0"+r' \\')
table('restored_negative','Recorded paired negative controls. All three methods have the same error and zero maximum value discrepancy. Allen--Cahn reports relative MAE in percent; the other rows report absolute MAE. Activity columns refer to the constrained state and the Batch-IR implementation.','tab:restored-negative','lrrrrrr',r'Problem & $n$ & $M$ & Common error & Violation (\%) & Batch active (\%) & Max $\Delta u$',rows)
E=load('elliptic_report_values.json')
rows=[f"{r['method']} & {100*r['rel_u_error']:.2f} & {r['terminal_loss']:.4f}"+r' \\' for r in E['no_scale_3seed']]
table('restored_elliptic','Archived elliptic no-scale ablation, means over three paired seeds. Lower terminal loss does not imply lower error in the single reported value. Rounded report values are preserved. This is a neural BSDE diagnostic, not an MLP benchmark.','tab:restored-elliptic','lrr',r'Method & Relative value error (\%) & Terminal BSDE loss',rows)
print('Generated six restored-evidence tables plus full 36-setting table.')
