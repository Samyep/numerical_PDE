"""Read-only re-export of the validated joint nonseparable benchmark."""

from invariant_region_mlp.experiments.single_point_joint_extension.benchmark import (
    JointNonseparableBenchmark,
    frozen_basis,
    make_benchmark,
    reduced_to_ambient_z,
)

__all__ = [
    "JointNonseparableBenchmark", "frozen_basis", "make_benchmark",
    "reduced_to_ambient_z",
]
