# %%
import sympy as sp
import control as ct
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt


m1, m2, m3, k1, k2, k3, c1, c2, c3, F1, F2, F3 = sp.symbols(
    "m1 m2 m3 k1 k2 k3 c1 c2 c3 F1 F2 F3"
)
s = sp.symbols("s")


# F_eqn = sp.Matrix([[1, -1], [0, 1], [0, 0]]) # original
F_eqn = sp.Matrix(
    [
        [1, -1, 0],
        [0, 1, -1],
        [0, 0, 1],
    ]
)  # modified

mass_matrix = sp.Matrix([[m1, 0, 0], [0, m2, 0], [0, 0, m3]])
damping_matrix = sp.Matrix(
    [
        [c1 + c2, -c2, 0],
        [-c2, c2 + c3, -c3],
        [0, -c3, c3],
    ]
)
stiffness_matrix = sp.Matrix(
    [
        [k1 + k2, -k2, 0],
        [-k2, k2 + k3, -k3],
        [0, -k3, k3],
    ]
)
print("Mass matrix: ")
sp.print_latex(mass_matrix)
print("Damping matrix: ")
sp.print_latex(damping_matrix)
print("stiffness matrix: ")
sp.print_latex(stiffness_matrix)
print("Force matrix: ")
sp.print_latex(F_eqn)
print("\n\n")

x_eqn = mass_matrix * s**2 + damping_matrix * s + stiffness_matrix

N_INPUTS = F_eqn.cols
N_OUTPUTS = F_eqn.rows

full_tf = x_eqn.inv() @ F_eqn


# %% [markdown]
# Demonstrate that all elements are rational functions and the denominator is the same:

# %%
reference_denom = sp.collect(full_tf[0, 0], s).as_numer_denom()[1]
error = False

for out_i in range(N_OUTPUTS):
    for in_i in range(N_INPUTS):
        if not full_tf[out_i, in_i].is_rational_function:
            print(f"{out_i=}, {in_i=} is not rational")

        elem_denom = sp.collect(full_tf[out_i, in_i], s).as_numer_denom()[1]
        if elem_denom != reference_denom:
            print(f"{out_i=}, {in_i=} denominator is wrong")
            error = True

if not error:
    print("All denominators are the same")


# %% [markdown]
# Print out all the transfer functions as Latex, evaluated

# %%
subs_dict = {
    m1: 2,
    m2: 0.2,
    m3: 0.03,
    k1: 1e4,
    k2: 3e4,
    k3: 4e4,
    c1: 0.1,
    c2: 0.1,
    c3: 0.1,
}
full_tf_eval = full_tf.subs(subs_dict)

# Divide out constant term in denominator to normalize
ref_constant = reference_denom.subs(subs_dict).as_poly(s).all_coeffs()[-1]

print(f"Denominator:")
print(sp.latex(reference_denom.subs(subs_dict) / ref_constant))

for out_i in range(N_OUTPUTS):
    for in_i in range(N_INPUTS):
        numer = sp.collect(full_tf_eval[out_i, in_i], s).as_numer_denom()[0]
        print(f"Numerator G{out_i+1}{in_i+1}:")
        sp.print_latex(numer / ref_constant)

# %% [markdown]
# Convert into control TransferFunction - use Sympy to get everything in to polynomials for numerator + denominator:

# %%
all_tfs = []
for out_i in range(N_OUTPUTS):
    tfs_row = []
    for in_i in range(N_INPUTS):
        numer, denom = sp.collect(full_tf_eval[out_i, in_i], s).as_numer_denom()
        numer_coeffs = [float(obj) for obj in numer.as_poly(s).all_coeffs()]
        denom_coeffs = [float(obj) for obj in denom.as_poly(s).all_coeffs()]

        tf = ct.TransferFunction(numer_coeffs, denom_coeffs)
        tfs_row.append(tf)

        print(f"--------------------")
        print(f"G{out_i}{in_i}:")

        # Normalize so constant term in denominator is 1
        num = np.squeeze(np.array(tf.num))
        den = np.squeeze(np.array(tf.den))
        den_const = den[-1]
        print(ct.tf(num / den_const, den / den_const))

    all_tfs.append(tfs_row)


# %%
# Merge all TFs into a MIMO system
system_tf = ct.combine_tf(all_tfs)


omega_max = 0
omega_min = np.inf

for out_i in range(N_OUTPUTS):
    for in_i in range(N_INPUTS):
        _, _, omega = ct.frequency_response(system_tf[out_i, in_i])
        omega_min = min(omega_min, omega[0])
        omega_max = max(omega_max, omega[-1])

N_POINTS_PER_DECADE = 100
N_POINTS = int(np.log10(omega_max / omega_min) * N_POINTS_PER_DECADE)
omega = np.geomspace(omega_min, omega_max, N_POINTS)

magnitude, phase, omega = ct.frequency_response(system_tf, omega=omega)

# Make bode plots of original system =======
fig, ax = plt.subplots(N_OUTPUTS * 2, 2, sharex=True, height_ratios=[3,2,3,2,3,2])
SCALE = 1.3# A4 size with margins
fig.set_size_inches((210 - 50) / 25.4 * SCALE, (297 - 50) / 25.4 * SCALE * 5/6)  

for out_i in range(N_OUTPUTS):
    for in_i in range(2):
        ax_mag = ax[out_i * 2, in_i]
        ax_phase = ax[out_i * 2 + 1, in_i]

        # mag_dB = 20 * np.log10(magnitude[out_i, in_i])
        freq = omega / (2 * pi)  # convert to Hz
        phase_deg = np.unwrap(phase[out_i, in_i]) * (180 / np.pi)
        ax_mag.loglog(freq, magnitude[out_i, in_i])
        ax_phase.semilogx(freq, phase_deg)

        ax_mag.grid(which="both")
        ax_phase.grid(which="both")
        ax_mag.set_ylabel(f"|$G_{{{out_i+1}{in_i+1}}}$(jω)|")
        ax_phase.set_ylabel(f"∠$G_{{{out_i+1}{in_i+1}}}$(jω) [deg]")

        # Set phase Y-ticks to be multiples of 90 degrees
        ylim = ax_phase.get_ylim()  # setting ticks changes limits
        ytick_vals = np.arange(
            np.floor(ylim[0] / 90) * 90, np.ceil(ylim[1] / 90) * 90 + 1, 90
        )
        ax_phase.set_yticks(ytick_vals)
        ax_phase.set_ylim(ylim)
ax[-1, 0].set_xlabel("Frequency (Hz)")
ax[-1, 1].set_xlabel("Frequency (Hz)")

fig.tight_layout()
fig.savefig("img/part1.A.2.bode.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# Make bode plot of 3rd input to system ========
fig, axs = plt.subplots(4, 2, sharex=True, height_ratios=[3,2,3,2])
ax_mags = [axs[0, 0], axs[0, 1], axs[2, 0]]
ax_phases = [axs[1, 0], axs[1, 1], axs[3, 0]]
for out_i in range(N_OUTPUTS):
    in_i = 2
    ax_mag = ax_mags[out_i]
    ax_phase = ax_phases[out_i]

    freq = omega / (2 * pi)  # convert to Hz
    phase_deg = np.unwrap(phase[out_i, in_i]) * (180 / np.pi)
    ax_mag.loglog(freq, magnitude[out_i, in_i])
    ax_phase.semilogx(freq, phase_deg)

    ax_mag.grid(which="both")
    ax_phase.grid(which="both")
    ax_mag.set_ylabel(f"|$G_{{{out_i+1}{in_i+1}}}$(jω)|")
    ax_phase.set_ylabel(f"∠$G_{{{out_i+1}{in_i+1}}}$(jω) [deg]")

    # Set phase Y-ticks to be multiples of 90 degrees
    ylim = ax_phase.get_ylim()  # setting ticks changes limits
    ytick_vals = np.arange(
        np.floor(ylim[0] / 90) * 90, np.ceil(ylim[1] / 90) * 90 + 1, 90
    )
    ax_phase.set_yticks(ytick_vals)
    ax_phase.set_ylim(ylim)

axs[-1, 0].set_xlabel("Frequency (Hz)")
axs[1, 1].set_xlabel("Frequency (Hz)")
axs[2, 1].axis('off')
axs[3, 1].axis('off')

# A4 size with margins, and scaling
fig.set_size_inches((210 - 50) / 25.4 * SCALE, (297 - 50) / 25.4* 2/3 * SCALE * 5/6)
fig.tight_layout()
axs[1, 1].tick_params(labelbottom=True)
fig.savefig("img/part1.A.5.bode.png", bbox_inches="tight", dpi=300)

# plt.show()
