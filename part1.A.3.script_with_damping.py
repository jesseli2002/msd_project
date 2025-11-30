"""
Modal analysis
"""

import numpy as np
import control as ct
import sympy as sp
from numpy import linalg as LA


# Define constants
m1 = 2
m2 = 0.2
m3 = 0.03
k1 = 1e4
k2 = 3e4
k3 = 4e4
c1 = 0.1
c2 = 0.1
c3 = 0.1

mass_matrix = np.array([[m1, 0, 0], [0, m2, 0], [0, 0, m3]])
damping_matrix = np.array(
    [
        [c1 + c2, -c2, 0],
        [-c2, c2 + c3, -c3],
        [0, -c3, c3],
    ]
)
stiffness_matrix = np.array(
    [
        [k1 + k2, -k2, 0],
        [-k2, k2 + k3, -k3],
        [0, -k3, k3],
    ]
)

sp.print_latex(sp.Matrix(mass_matrix))

eigvals, eigvecs = LA.eig(LA.inv(mass_matrix) @ stiffness_matrix)
# print(f"{eigvals=}")
print(f"{eigvecs=}")

phi_to_x = eigvecs # conversion from modal coordinates to [x1, x2, x3]
modal_mass = phi_to_x.T @ mass_matrix @ phi_to_x
modal_damping = phi_to_x.T @ damping_matrix @ phi_to_x
modal_stiffness = phi_to_x.T @ stiffness_matrix @ phi_to_x

# Sanity check - make sure these matrices are diagonal
assert np.isclose(modal_mass, np.diag(np.diag(modal_mass))).all()
assert np.isclose(modal_stiffness, np.diag(np.diag(modal_stiffness))).all()

modal_mass = np.diag(modal_mass)
modal_stiffness = np.diag(modal_stiffness)
modal_damping = np.diag(modal_damping) # ignore off-diagonal terms

modal_tfs = np.array([ct.tf([1], [m, c, k]) for m, c, k in zip(modal_mass, modal_damping, modal_stiffness)])

forces = np.array([[1, 0, 0], [-1, 1, 0]]).T
# warning - associativity is no joke; at least one of these parentheses is load bearing
sys = phi_to_x @ (modal_tfs[:, None] * (phi_to_x.T @ forces))

N_OUTPUTS = sys.shape[0]
N_INPUTS = sys.shape[1]

common_den = np.squeeze(np.array(sys[0, 0].den))
print(f"Gij denominator: ")
print(common_den / common_den[-1])
for out_i in range(N_OUTPUTS):
    for in_i in range(N_INPUTS):
        print(f"\n\n====================\n")
        print(f"G{out_i}{in_i} numerator:")
        
        # Normalize so constant term in denominator is 1
        num = np.squeeze(np.array(sys[out_i, in_i].num))
        den = np.squeeze(np.array(sys[out_i, in_i].den))
        if not np.all(common_den == den):
            raise RuntimeError(f"Unequal denominators {out_i=} {in_i=}")

        den_const = den[-1]
        
        print(num/den_const)



