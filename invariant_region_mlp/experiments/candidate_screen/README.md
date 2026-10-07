# Candidate PDE screen (rough: 240 points, 2 paired repetitions, d in {20, 100})

Reuses the mechanism-suite `MechanismMLP` unchanged; only `_project_segment` is generalised (in `screen.py`)
to the new families. Results: `results/candidate_screen/` (`screen.jsonl`, `summary.md`).

| id | PDE | mechanism | reference | certificate |
|---|---|---|---|---|
| C1 | l1-control HJB `u_t + Lap u - lam ||grad u||_1 = 0`, `lam = 1/||w||_1` | convex, noise aggregated in l1 | same ridge solution as P4 (cached 1-D reference) | segment `|psi_s| <= 1`, box, ball, span |
| C2 | zero-sum LQ game `u_t + Lap u - a|grad_A u|^2 + b|grad_B u|^2 = 0`, `a - b = 2`, w has mass 1/2 per block | nonconvex (Isaacs); one-step orthogonal Jensen bias `(c/2M)[-a(d_A-1/2) + b(d_B-1/2)]` | same closed form as P1 | P1 segment |
| C3 | cubic viscous HJ `u_t + Lap u - 0.5|grad u|^3 = 0`, `g = 2 log cosh(2s)/2`, T=0.25 | superlinear growth | 1-D Godunov + Richardson (screening grade, residual ~1e-5) | segment `|psi_s| <= 2` |

C2 configurations: convex (a=2, b=0, d_A=d/2), cancel (a=3, b=1, d_A=d/4, predicted factor 0.5),
flip (a=3, b=1, d_A=d/8, predicted factor > 0).

Run: `python screen.py ../../results/candidate_screen/screen.jsonl 2` (about 10 min on 2 CPUs), then `python analyze.py`.
Not a pre-registered study; use it to choose candidates for the next pre-registered round.
