# Counterexample rescue results

- `summary_all.csv`: all 252 method/configuration aggregate rows from 36 (d,n) settings.
- `raw_n3_d{2,10,100,1000}_replicas.csv`: 512-replica observations for every compared method at the four headline dimensions.
- `independent_validation.json`: literal Eq. (6) comparisons and 65,536-replica Gaussian-moment checks.
- `figures/*.svg`: repository-native versions of the paper figures.

The full temporary NPZ arrays are not committed because they exceed 80 MB and contain redundant full gradient vectors. The committed raw replica tables contain the quantities used for all headline RMSE/error plots and are sufficient to audit those reported numbers.
