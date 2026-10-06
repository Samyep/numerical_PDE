from pathlib import Path
import json
P=Path(__file__).resolve().parents[1];D=P/'data/batchwise_extended';T=P/'tables'
s=json.loads((D/'batchwise_extended_suite_summary.json').read_text());h=json.loads((D/'hjb_headline_extracted.json').read_text())['rows']
def table(name,caption,label,align,header,rows):
 text='\\begin{table}[htbp]\n\\centering\\small\n\\caption{'+caption+'}\\label{'+label+'}\n\\begin{tabular}{@{}'+align+'@{}}\n\\toprule\n'+header+'\\\\\n\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'
 (T/(name+'.tex')).write_text(text)
def sci(x):
 a,e=f'{x:.2e}'.split('e'); return f'${a}\\times10^{{{int(e)}}}$'
rows=[]
for r,su in zip(h,s['headline_hjb']):
 vals=[]
 for key in ['raw','samplewise','batch']:
  v=r['relative_l2'][key];vals.append(f"${v['mean']:.5f}\\pm{v['std']:.5f}$")
 rows.append(f"{r['d']} & "+' & '.join(vals)+f" & {su['batch_vs_sample_reduction']*100:.1f}\\% & 10/10 & {sci(su['paired_t_p'])}\\\\")
table('extended_hjb','Headline HJB: mean relative $L_2$ value error $\\pm$ sample SD over ten paired repetitions. Reduction is relative to Samplewise IR. The raw method has no correction.','tab:extended-hjb','rccc rcc','$d$ & Raw & Samplewise IR & Batch-IR & Reduction & Wins & Paired $p$',rows)
rows=[]
for r in h:
 vals=[r['mechanism'][metric][key] for metric in ['generator_mse','correction_variance'] for key in ['raw','samplewise','batch']]
 rows.append(str(r['d'])+' & '+' & '.join(f'{x:.4f}' if j>=3 else f'{x:.2f}' for j,x in enumerate(vals))+r'\\')
table('extended_mechanism','Matched HJB mechanism diagnostics. Generator MSE and averaged-correction variance use the definitions and parent sets given in the text.','tab:extended-mechanism','rrrrrrr',r'& \multicolumn{3}{c}{Generator MSE} & \multicolumn{3}{c}{Averaged-correction variance}\\ \cmidrule(lr){2-4}\cmidrule(l){5-7} $d$ & Raw & Samplewise IR & Batch-IR & Raw & Samplewise IR & Batch-IR',rows)
rows=[]
for r in h:
 v=r['mechanism']['individual_increment_variance']; rows.append(str(r['d'])+' & '+' & '.join(f'{v[k]:.4f}' for k in ['raw','samplewise','batch'])+' & '+f"{100*(v['batch']/v['samplewise']-1):+.1f}\\%"+r'\\')
table('extended_increments','Pooled individual-increment variance on the diagnostic child states. Positive change means Batch-IR is worse on this particular statistic; no uniform dominance is claimed.','tab:extended-increments','rrrrr','$d$ & Raw & Samplewise IR & Batch-IR & Change vs. sample',rows)
rows=[f"{r['d']} & {100*r['batch_activation_rate']:.3f}\\% & {r['batch_alpha_mean']:.4f}"+r'\\' for r in s['headline_hjb']]
table('extended_activity','Common-scale activity in the headline HJB study. These are sibling-batch statistics, not scalar-coordinate violation rates.','tab:extended-activity','rrr','$d$ & Batch activation & Mean $\\alpha_B$',rows)
rows=[f"{r['n']} & {r['M']} & {r['raw_mae']:.5f} & {r['samplewise_mae']:.5f} & {r['batch_mae']:.5f} & {100*r['batch_vs_sample_reduction']:.1f}\\% & {100*r['batch_win_fraction']:.0f}/100 & {sci(r['paired_t_p'])}"+r'\\' for r in s['funding_deep']]
table('extended_funding','Extended 100D funding: 100 paired root estimates per setting; all MAEs use reference 21.299. Paired tests are reported by the archived summary; no replica-level confidence intervals are inferred.','tab:extended-funding','rrrrrrrr','$n$ & $M$ & Raw MAE & Samplewise IR & Batch-IR & Reduction & Wins & Reported $p$',rows)
rows=[f"{r['d']} & {r['reference']:.8f} & {r['raw_mae']:.5f} & {r['samplewise_mae']:.5f} & {r['batch_mae']:.5f}"+r'\\' for r in s['neufeld_m2']]
table('extended_neufeld','Neufeld--Wu extension at $n=M=2$, 30 repetitions per dimension. Samplewise IR and Batch-IR tie; this table must not be pooled with the older runs in Table~\\ref{tab:nw-extra}.','tab:extended-nw','rrrrr','$d$ & Reference & Raw MAE & Samplewise IR MAE & Batch-IR MAE',rows)
rows=[f"{r['d']} & {r['raw_value_rmse']:.6f} & {r['samplewise_subspace_value_rmse']:.6f} & {r['batch_structural_value_rmse']:.6f} & 0"+r'\\' for r in s['counterexample_sanity']]
table('extended_sanity','Separate batchwise structural sanity check: value RMSE over 512 paired replicas at $n=m=3$. Maximum discrepancy compares paired values, not gradients or full states.','tab:extended-sanity','rrrrr','$d$ & Raw & Subspace IR & Batch structural & Max. paired discrepancy',rows)
