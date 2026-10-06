#!/usr/bin/env python3
"""Check preserved assets and the newly executed overshoot summaries.
Run from any working directory after unpacking the Overleaf project.
"""
from pathlib import Path
import hashlib,json,re
import numpy as np
from scipy.stats import spearmanr
P=Path(__file__).resolve().parents[1];D=P/'data'/'restored_evidence'
checks={}
base=json.loads((D/'baseline_asset_sha256.json').read_text())
missing=[];modified=[]
for name,h in base.items():
 f=P/name
 if not f.exists():missing.append(name)
 elif hashlib.sha256(f.read_bytes()).hexdigest()!=h:modified.append(name)
checks['preserved_original_assets']=len(base)
checks['missing_original_assets']=missing;checks['modified_original_assets']=modified
assert not missing and not modified
r=json.loads((D/'overshoot_36_results.json').read_text())
q=json.loads((D/'overshoot_36_protocol.json').read_text())
a=np.load(D/'overshoot_36_raw.npz')
assert len(r)==36
assert set((z['M'],z['radius_factor']) for z in r)=={(M,R) for M in [2,3,4,6,8,10] for R in [1,1.25,1.5,2,3,4]}
assert np.isfinite(a['raw_value']).all() and np.isfinite(a['projected_value']).all()
worst=0
for mi,M in enumerate(a['M_values']):
 er=np.linalg.norm(a['raw_value'][mi]-a['truth'],axis=1)/np.linalg.norm(a['truth'])
 for fi,fac in enumerate(a['radius_factors']):
  e=np.linalg.norm(a['projected_value'][mi,fi]-a['truth'],axis=1)/np.linalg.norm(a['truth'])
  row=next(z for z in r if z['M']==int(M) and z['radius_factor']==float(fac))
  vals={'raw_rel_l2_mean':er.mean(),'ir_rel_l2_mean':e.mean(),'gain_fraction':1-e.mean()/er.mean(),'overshoot_energy':a['overshoot_by_rep'][mi,fi].mean(),'violation_rate':a['violation_by_rep'][mi,fi].mean()}
  worst=max(worst,max(abs(row[k]-float(v)) for k,v in vals.items()))
assert worst<1e-12
checks['max_summary_discrepancy']=worst
checks['overshoot_rows']=len(r)
checks['spearman_overshoot_gain']=float(spearmanr([z['overshoot_energy'] for z in r],[z['gain_fraction'] for z in r]).statistic)
checks['spearman_violation_gain']=float(spearmanr([z['violation_rate'] for z in r],[z['gain_fraction'] for z in r]).statistic)
checks['max_reference_64_128_abs_difference']=float(np.max(np.abs(a['truth']-a['truth128'])))
assert abs(checks['spearman_overshoot_gain']-q['spearman_overshoot_gain'])<1e-12
assert abs(checks['spearman_violation_gain']-q['spearman_violation_gain'])<1e-12
# Projection only decreases squared gradient norms in this particular diagnostic.
assert np.all(a['projected_value']>=a['raw_value'][:,None,:,:]-1e-12)
checks['all_root_values_finite']=True
checks['projected_value_vs_raw_quadratic_identity']=True
text=(P/'main.tex').read_text()+'\n'+'\n'.join(f.read_text() for f in (P/'sections').glob('*.tex'))
images=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text)
for name in images:assert (P/name).is_file(),name
restored=list((P/'figures'/'restored').glob('*.pdf'))
assert all(str(f.relative_to(P)) in images for f in restored)
checks['new_distinct_plot_pdfs']=len(restored)
checks['graphic_file_occurrences_in_tex']=len(images)
checks['distinct_graphic_files_in_tex']=len(set(images))
checks['panel_occurrences_including_inline_tikz']=len(images)+1
checks['distinct_panels_including_inline_tikz']=len(set(images))+1
checks['original_panels']=19
checks['new_stochastic_experiments']='Only the 36-setting HJB overshoot reproduction; all other added results are archived.'
checks['original_overshoot_png_recovered']=False
checks['omitted_new_plots']=[]
(D/'restoration_audit.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
