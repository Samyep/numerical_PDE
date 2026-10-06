"""Check imported summary arithmetic and available HJB repetition-level records.

Requires NumPy and SciPy. Does not run any PDE experiment or use the network.
"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import ttest_rel

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'batchwise_extended'
h = json.loads((DATA / 'hjb_headline_extracted.json').read_text())['rows']
s = json.loads((DATA / 'batchwise_extended_suite_summary.json').read_text())
for r,summary in zip(h,s['headline_hjb']):
    assert r['d'] == summary['d']
    for key, record in r['relative_l2'].items():
        a = np.asarray(record['values'],dtype=float)
        assert len(a) == r['reps'] == 10
        np.testing.assert_allclose(a.mean(),record['mean'],rtol=1e-13)
        np.testing.assert_allclose(a.std(ddof=1),record['std'],rtol=1e-11)
        np.testing.assert_allclose(record['mean'],summary[f'{key}_rel_l2'],rtol=1e-13)
    a = np.asarray(r['relative_l2']['samplewise']['values'])
    b = np.asarray(r['relative_l2']['batch']['values'])
    assert np.all(b<a)
    np.testing.assert_allclose(ttest_rel(a,b).pvalue,summary['paired_t_p'],rtol=1e-10,atol=0)
    np.testing.assert_allclose(1-b.mean()/a.mean(),summary['batch_vs_sample_reduction'])
for r in s['funding_deep']:
    np.testing.assert_allclose(1-r['batch_mae']/r['samplewise_mae'],r['batch_vs_sample_reduction'])
for r in s['neufeld_m2']:
    assert r['samplewise_mae']==r['batch_mae']
for r in s['counterexample_sanity']:
    assert r['samplewise_subspace_value_rmse']==r['batch_structural_value_rmse']
    assert r['max_abs_batch_vs_subspace']==0
print('All imported numerical consistency checks passed.')
print('Funding paired tests were not recomputed: per-replica arrays are absent from the imported summary.')
