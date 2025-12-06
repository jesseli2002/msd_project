"""
Modal analysis
"""

import numpy as np
import control as ct
import sympy as sp
from numpy import pi, linalg as LA
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

phi_to_x = eigvecs  # conversion from modal coordinates to [x1, x2, x3]
modal_mass = phi_to_x.T @ mass_matrix @ phi_to_x
modal_damping = phi_to_x.T @ damping_matrix @ phi_to_x
modal_stiffness = phi_to_x.T @ stiffness_matrix @ phi_to_x

# Sanity check - make sure these matrices are diagonal
assert np.isclose(modal_mass, np.diag(np.diag(modal_mass))).all()
assert np.isclose(modal_stiffness, np.diag(np.diag(modal_stiffness))).all()

modal_mass = np.diag(modal_mass)
modal_stiffness = np.diag(modal_stiffness)
modal_damping = np.diag(modal_damping)  # ignore off-diagonal terms

modal_tfs = np.array(
    [ct.tf([1], [m, 0, k]) for m, k in zip(modal_mass, modal_stiffness)]
)

# forces = np.array([[1, 0, 0], [-1, 1, 0]]).T # 
forces = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 1]]).T

# warning - associativity is no joke; at least one of these parentheses is load bearing
sys = phi_to_x @ (modal_tfs[:, None] * (phi_to_x.T @ forces))

N_OUTPUTS = sys.shape[0]
N_INPUTS = sys.shape[1]

common_den = np.squeeze(np.array(sys[0, 0].den))
print(f"Gij denominator: ")
print(common_den / common_den[-1])
for out_i in range(N_OUTPUTS):
    for in_i in range(N_INPUTS):
        print(f"\n\n====================")
        print(f"G{out_i+1}{in_i+1} numerator:")

        # Normalize so constant term in denominator is 1
        num = np.squeeze(np.array(sys[out_i, in_i].num))
        den = np.squeeze(np.array(sys[out_i, in_i].den))
        if not np.all(common_den == den):
            raise RuntimeError(f"Unequal denominators {out_i=} {in_i=}")

        den_const = den[-1]

        print(num / den_const)

# Section A.4 =============================
# Section A.4 =============================
# Section A.4 =============================
# Add damping
modal_tfs_damp = np.array(
    [ct.tf([1], [m, c, k]) for m, c, k in zip(modal_mass, modal_damping, modal_stiffness)]
)
sys_damp = phi_to_x @ (modal_tfs_damp[:, None] * (phi_to_x.T @ forces))

def plot_mode_shape(mode_vector, title="Mode Shape"):
    """
    Visualize a mode shape for a 3-mass system.
    
    Parameters:
    mode_vector : array-like with 3 elements
        The mode shape displacement for each of the 3 masses
    title : str
        Title for the plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Normalize mode vector for better visualization
    mode = np.array(mode_vector)
    mode = mode / np.max(np.abs(mode))  # Normalize to [-1, 1]
    
    # Base positions for the three masses
    base_positions = np.array([1, 2, 3])
    
    # Scale factor for displacement visibility
    scale = 0.4
    
    offsets = [-scale * mode, np.zeros(3), scale * mode]
    
    # Mass size
    mass_width = 0.3
    mass_height = 0.3
    
    # Wall size
    wall_height = 0.35

    # Draw the wall on the left
    wall_x = 0
    
    # Spring color
    spring_color = 'gray'

    # Three rows: -mode, 0 (equilibrium), +mode
    y_positions = np.arange(3) * (wall_height + 0.1)
    
    y_margin = 0.05
    
    # Draw each row
    for row_idx, (y_pos, offset) in enumerate(zip(y_positions, offsets)):
        ax.plot([wall_x, wall_x], [y_pos - wall_height/2, y_pos + wall_height / 2], 'k-', linewidth=3)
        
        # Draw masses and springs
        for i in range(3):
            x_pos = base_positions[i] + offset[i]
            
            # Draw spring from wall to first mass
            if i == 0:
                spring_start = wall_x
                spring_end = x_pos - mass_width/2
                spring_y = y_pos
                n_coils = 8
                spring_x = np.linspace(spring_start, spring_end, n_coils*2)
                spring_y_coords = np.zeros_like(spring_x)
                spring_y_coords[1:-1:2] = spring_y + 0.1
                spring_y_coords[2:-1:2] = spring_y - 0.1
                spring_y_coords[0] = spring_y
                spring_y_coords[-1] = spring_y
                ax.plot(spring_x, spring_y_coords, '-', color=spring_color, linewidth=1.5)
            
            # Draw spring to next mass
            if i < 2:
                spring_start = x_pos + mass_width/2
                spring_end = base_positions[i+1] + offset[i+1] - mass_width/2
                spring_y = y_pos
                n_coils = 8
                spring_x = np.linspace(spring_start, spring_end, n_coils*2)
                spring_y_coords = np.zeros_like(spring_x)
                spring_y_coords[1:-1:2] = spring_y + 0.1
                spring_y_coords[2:-1:2] = spring_y - 0.1
                spring_y_coords[0] = spring_y
                spring_y_coords[-1] = spring_y
                ax.plot(spring_x, spring_y_coords,'-', color=spring_color, linewidth=1.5)
            
            # Draw mass
            rect = patches.Rectangle((x_pos - mass_width/2, y_pos - mass_height/2),
                                     mass_width, mass_height,
                                     linewidth=2, edgecolor='black', facecolor='lightblue')
            ax.add_patch(rect)
            
            # Label the mass
            ax.text(x_pos, y_pos, f'$m_{i+1}$', ha='center', va='center',
                   fontsize=22, fontweight='bold')
    
    # Set axis properties
    ax.set_xlim(wall_x, base_positions[-1] + 0.8)
    ax.set_ylim(y_positions[0] - mass_height / 2 - y_margin, y_positions[-1] + mass_height / 2 + y_margin)
    ax.set_aspect('equal')
    ax.axis('off')  # Remove all axes, ticks, and labels
    
    fig.tight_layout()
    return fig, ax

# Create visualization of each mode
for mode_i in range(3):
    mode = eigvecs[:, mode_i]
    freq = np.sqrt(eigvals[mode_i]) / (2 * pi)
    fig, ax = plot_mode_shape(mode)
    ax.set_title(f"Mode {mode_i+1}, {freq:.1f} Hz", fontsize=22)
    
    fig.savefig(f'img/part.A.4.mode{mode_i+1}.png', bbox_inches='tight')


# Plot contributions from each mode, to G11 and G21
plt.close('all')
freq_min = 1
freq_max = 1e6
N_POINTS_PER_DECADE = 100
N_POINTS = int(np.log10(freq_max / freq_min) * N_POINTS_PER_DECADE)
freq = np.geomspace(freq_min, freq_max, N_POINTS + 1)
omega = freq * 2 * pi

# G11
for out_i in [0, 1]:  # output index (x1 and x2)
    in_i = 0  # input index (f1)

    fig, (ax_mag, ax_phase) = plt.subplots(nrows=2, sharex=True,layout='tight', height_ratios=[2, 1])

    def add_plot(sys_tf, **plot_kwargs):
        mag, phase, _ = ct.frequency_response(sys_tf, omega)
        ax_mag.loglog(freq, mag, **plot_kwargs)
        ax_phase.semilogx(freq, np.unwrap(phase) * 180 / pi, **plot_kwargs)

    # shape (N_MODES, N_INPUTS)
    mode_contributions = modal_tfs_damp[:, None] * (phi_to_x.T @ forces)
    mode_colors = ["red", "green", "blue"]

    for mode_i in range(3):
        mode_contribution = phi_to_x[out_i, mode_i] * mode_contributions[mode_i, in_i]
        add_plot(mode_contribution, label=f"Mode {mode_i + 1}", color=mode_colors[mode_i])

    add_plot(sys_damp[out_i, in_i], label="Total", color="black", alpha=0.6)
    ax_mag.grid()
    ax_phase.grid()
    ax_mag.legend()
    ax_phase.set_xlim(freq_min, freq_max)

    ylim = ax_phase.get_ylim()  # setting ticks changes limits
    ytick_spacing = 180
    ytick_vals = np.arange(
        np.floor(ylim[0] / ytick_spacing) * ytick_spacing, np.ceil(ylim[1] / ytick_spacing) * ytick_spacing + 1, ytick_spacing
    )
    ax_phase.set_yticks(ytick_vals)
    ax_phase.set_ylim(ylim)
    
    ax_mag.set_ylabel("Magnitude")
    ax_phase.set_ylabel("Phase [deg]")
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_mag.set_title(f"$G_{{{out_i + 1}{in_i + 1}}}$ Bode Plot")

    fig.set_size_inches(6 * 0.8, 5 * 0.8)
    fig.savefig(f"img/part.A.4.G{out_i}{in_i}.png", bbox_inches='tight')

    # plt.show()
