# Three-seed compute scaling

This directory implements the frozen study in `PREREGISTRATION.md` with a
small sequence: timing-only calibration, independent formal runs, fixed-bank
evaluation, and mean±standard-deviation aggregation.  C0 reuses checkpoints
from one trajectory per seed; C1/C2 train each parameter value independently.
