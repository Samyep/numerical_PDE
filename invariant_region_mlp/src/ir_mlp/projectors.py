"""Certified projection operators for IR-MLP."""
import math
import jax.numpy as jnp

def project_l2_ball(z, radius: float):
    norms = jnp.linalg.norm(z, axis=-1, keepdims=True)
    scale = jnp.minimum(1.0, float(radius) / jnp.maximum(norms, 1e-12))
    return z * scale

class HJBGradientBallProjector:
    def __init__(self, radius: float = math.sqrt(15.0)):
        self.radius = float(radius)
    def __call__(self, output_uz, x_t, solver):
        del x_t, solver
        output_uz = jnp.asarray(output_uz)
        return jnp.concatenate((output_uz[:, :1], project_l2_ball(output_uz[:, 1:], self.radius)), axis=-1)

class ScaSMLHJBDefectProjector:
    def __init__(self, radius: float = math.sqrt(15.0)):
        self.radius = float(radius)
    def __call__(self, defect_uz, x_t, solver):
        defect_uz = jnp.asarray(defect_uz)
        grad_u_hat = jnp.asarray(solver.model.predict(x_t, operator=solver.equation.grad))
        z_hat = jnp.asarray(solver.equation.sigma(x_t)) * grad_u_hat
        z_total_proj = project_l2_ball(z_hat + defect_uz[:, 1:], self.radius)
        return jnp.concatenate((defect_uz[:, :1], z_total_proj - z_hat), axis=-1)
