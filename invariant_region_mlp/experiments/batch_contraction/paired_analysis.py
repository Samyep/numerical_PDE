import json, numpy as np
from scipy.stats import ttest_rel

def paired_summary(raw, method):
    raw=np.asarray(raw,float); method=np.asarray(method,float)
    diff=raw-method
    return {
        "raw_mean_error": float(raw.mean()),
        "method_mean_error": float(method.mean()),
        "paired_improvement_mean": float(diff.mean()),
        "paired_improvement_se": float(diff.std(ddof=1)/np.sqrt(len(diff))),
        "win_rate": float(np.mean(method<raw)),
        "p_ttest": float(ttest_rel(raw,method).pvalue),
    }

# The committed result file was produced with paired random paths.
# Re-run the benchmark module with the seed schedules documented in the report
# to regenerate per-run arrays, then call paired_summary on absolute/relative
# errors. This file is intentionally small so the statistical comparison is
# transparent.
