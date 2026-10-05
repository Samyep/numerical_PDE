import math
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent


def summary(vals, ref):
    vals = np.asarray(vals, dtype=float)
    err = np.abs(vals - ref)
    return {
        "mean": float(vals.mean()),
        "std": float(vals.std(ddof=1)),
        "mae": float(err.mean()),
        "rel_mae": float(err.mean() / max(abs(ref), 1e-15)),
    }


AC_REF = 0.052802
AC_T = 0.3
AC_D = 100


def ac_g(x):
    return 1.0 / (2.0 + 0.4 * np.dot(x, x))


def ac_f(u):
    return u - u**3


class AllenMLP:
    def __init__(self, M, mode, seed):
        self.M = M
        self.mode = mode
        self.rng = np.random.default_rng(seed)
        self.n_checked = 0
        self.n_viol = 0
        self.n_batch = 0
        self.n_batch_active = 0

    def transform(self, vals):
        vals = np.asarray(vals, dtype=float)
        self.n_checked += vals.size
        self.n_viol += int(np.sum((vals < 0.0) | (vals > 1.0)))
        if self.mode == "raw":
            return vals
        if self.mode == "samplewise":
            return np.clip(vals, 0.0, 1.0)
        self.n_batch += 1
        if np.any(vals < 0):
            alpha = 0.0
        else:
            vmax = float(np.max(vals)) if vals.size else 0.0
            alpha = min(1.0, 1.0 / vmax) if vmax > 0 else 1.0
        if alpha < 1.0:
            self.n_batch_active += 1
        return alpha * vals

    def comp(self, t, x, n):
        if n == 0:
            return 0.0
        dt = AC_T - t
        Mn = self.M**n
        Z = self.rng.normal(size=(Mn, AC_D))
        XT = x[None, :] + math.sqrt(2.0 * dt) * Z
        u = float(np.mean([ac_g(xx) for xx in XT]))
        for ell in range(n):
            m = self.M ** (n - ell)
            up, um = [], []
            for _ in range(m):
                R = t + dt * self.rng.uniform()
                Y = x + math.sqrt(2.0 * (R - t)) * self.rng.normal(size=AC_D)
                up.append(self.comp(R, Y, ell))
                if ell > 0:
                    um.append(self.comp(R, Y, ell - 1))
            up = self.transform(up)
            corr = np.array([ac_f(v) for v in up])
            if ell > 0:
                um = self.transform(um)
                corr -= np.array([ac_f(v) for v in um])
            u += dt * float(corr.mean())
        return u


CR_REF = 2.626
CR_T = 2.0
CR_D = 100
CR_SIGMA = 0.2
CR_BETA = 0.03
CR_K1, CR_K2, CR_L = 30.0, 60.0, 15.0


def cr_g(x):
    m = float(np.min(x))
    return max(m - CR_K1, 0.0) - max(m - CR_K2, 0.0) - CR_L


def cr_f(u):
    return CR_BETA * (max(u, 0.0) - u)


class CreditMLP:
    def __init__(self, M, mode, seed):
        self.M = M
        self.mode = mode
        self.rng = np.random.default_rng(seed)
        self.n_checked = 0
        self.n_viol = 0
        self.n_batch = 0
        self.n_batch_active = 0

    def gbm_step(self, t, s, x):
        dt = s - t
        z = self.rng.normal(size=CR_D)
        return x * np.exp(-0.5 * CR_SIGMA**2 * dt + CR_SIGMA * math.sqrt(dt) * z)

    def transform(self, vals):
        vals = np.asarray(vals, dtype=float)
        self.n_checked += vals.size
        self.n_viol += int(np.sum((vals < -15.0) | (vals > 15.0)))
        if self.mode == "raw":
            return vals
        if self.mode == "samplewise":
            return np.clip(vals, -15.0, 15.0)
        self.n_batch += 1
        vmax = float(np.max(np.abs(vals))) if vals.size else 0.0
        alpha = min(1.0, 15.0 / vmax) if vmax > 0 else 1.0
        if alpha < 1.0:
            self.n_batch_active += 1
        return alpha * vals

    def comp(self, t, x, n):
        if n == 0:
            return 0.0
        dt = CR_T - t
        Mn = self.M**n
        Z = self.rng.normal(size=(Mn, CR_D))
        XT = x[None, :] * np.exp(-0.5 * CR_SIGMA**2 * dt + CR_SIGMA * math.sqrt(dt) * Z)
        u = float(np.mean([cr_g(xx) for xx in XT]))
        for ell in range(n):
            m = self.M ** (n - ell)
            up, um = [], []
            for _ in range(m):
                R = t + dt * self.rng.uniform()
                Y = self.gbm_step(t, R, x)
                up.append(self.comp(R, Y, ell))
                if ell > 0:
                    um.append(self.comp(R, Y, ell - 1))
            up = self.transform(up)
            corr = np.array([cr_f(v) for v in up])
            if ell > 0:
                um = self.transform(um)
                corr -= np.array([cr_f(v) for v in um])
            u += dt * float(corr.mean())
        return u


LC_D = 100
LC_T = 0.5
A_DIR = np.ones(LC_D) / np.sqrt(LC_D)
B_DRIFT = 0.20 * np.ones(LC_D) / np.sqrt(LC_D)
LC_REF = math.exp(-0.5 * LC_T * np.dot(A_DIR, A_DIR)) * math.sin(
    LC_T * np.dot(A_DIR, B_DRIFT)
)


def run():
    rows = []

    for M in [2, 3, 4]:
        vals_by_mode, meta_by_mode = {}, {}
        for mode in ["raw", "samplewise", "batch"]:
            vals, vrs, bars = [], [], []
            for rep in range(100):
                alg = AllenMLP(M, mode, 10000 + 100 * M + rep)
                vals.append(alg.comp(0.0, np.zeros(AC_D), 3))
                vrs.append(alg.n_viol / alg.n_checked if alg.n_checked else 0.0)
                bars.append(alg.n_batch_active / alg.n_batch if alg.n_batch else 0.0)
            vals_by_mode[mode] = vals
            meta_by_mode[mode] = {
                "violation_rate": float(np.mean(vrs)),
                "batch_activation_rate": float(np.mean(bars)),
            }
        for mode in ["raw", "samplewise", "batch"]:
            row = {"benchmark": "Allen-Cahn 100D", "n": 3, "M": M, "method": mode}
            row.update(summary(vals_by_mode[mode], AC_REF))
            row.update(meta_by_mode[mode])
            row["max_abs_difference_from_raw"] = float(
                np.max(np.abs(np.asarray(vals_by_mode[mode]) - np.asarray(vals_by_mode["raw"])))
            )
            rows.append(row)

    for n, M in [(2, 10), (3, 6), (4, 3)]:
        vals_by_mode, meta_by_mode = {}, {}
        for mode in ["raw", "samplewise", "batch"]:
            vals, vrs, bars = [], [], []
            for rep in range(100):
                alg = CreditMLP(M, mode, 20000 + 1000 * n + 100 * M + rep)
                vals.append(alg.comp(0.0, np.ones(CR_D) * 100.0, n))
                vrs.append(alg.n_viol / alg.n_checked if alg.n_checked else 0.0)
                bars.append(alg.n_batch_active / alg.n_batch if alg.n_batch else 0.0)
            vals_by_mode[mode] = vals
            meta_by_mode[mode] = {
                "violation_rate": float(np.mean(vrs)),
                "batch_activation_rate": float(np.mean(bars)),
            }
        for mode in ["raw", "samplewise", "batch"]:
            row = {"benchmark": "Counterparty credit risk 100D", "n": n, "M": M, "method": mode}
            row.update(summary(vals_by_mode[mode], CR_REF))
            row.update(meta_by_mode[mode])
            row["max_abs_difference_from_raw"] = float(
                np.max(np.abs(np.asarray(vals_by_mode[mode]) - np.asarray(vals_by_mode["raw"])))
            )
            rows.append(row)

    for M in [4, 8, 16]:
        N = M**2
        raw_vals, sample_vals, batch_vals = [], [], []
        z_viol_rates, batch_active = [], []
        for rep in range(200):
            rng = np.random.default_rng(30000 + 100 * M + rep)
            W = math.sqrt(LC_T) * rng.normal(size=(N, LC_D))
            XT = B_DRIFT[None, :] * LC_T + W
            gv = np.sin(XT @ A_DIR)
            uhat = float(gv.mean())
            zhat = np.mean(gv[:, None] * W / LC_T, axis=0)
            z_norm = float(np.linalg.norm(zhat))
            z_viol_rates.append(float(z_norm > 1.0))
            alpha = min(1.0, 1.0 / z_norm) if z_norm > 0 else 1.0
            batch_active.append(float(alpha < 1.0))
            raw_vals.append(uhat)
            sample_vals.append(uhat)
            batch_vals.append(uhat)

        for mode, vals in [("raw", raw_vals), ("samplewise", sample_vals), ("batch", batch_vals)]:
            row = {"benchmark": "Linear convection-diffusion 100D", "n": 2, "M": M, "method": mode}
            row.update(summary(vals, LC_REF))
            row["violation_rate"] = float(np.mean(z_viol_rates))
            row["batch_activation_rate"] = float(np.mean(batch_active)) if mode == "batch" else 0.0
            row["max_abs_difference_from_raw"] = float(
                np.max(np.abs(np.asarray(vals) - np.asarray(raw_vals)))
            )
            rows.append(row)

    (OUT / "negative_control_batchir_results.json").write_text(json.dumps(rows, indent=2))
    return rows


if __name__ == "__main__":
    run()
