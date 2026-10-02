# SCaSML in-place integration

The complete patch files are in the source archive. The minimal change adds an optional `state_projector` callback to `MLP_full_history.py` and `ScaSML_full_history.py`. Recursive calls use the same solver instance, so the projected state is fed back inside the Picard recursion.

For the public Rosenbrock coefficient ranges c1,c2 in [0.5,1.5], a conservative bound is lambda_max(A) <= 7.5, hence with sigma=sqrt(2): `||z||_2 <= sqrt(15)`.

For ScaSML, the recursive state is a defect. The projector therefore forms the total gradient `z_hat + z_breve`, projects the total state, then subtracts `z_hat` again.

Main experiments use the corrected EBL terminal-gradient normalization; the public-repo normalization is kept only as a diagnostic.
