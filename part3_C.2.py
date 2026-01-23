"""
Feedforward controller design

Basically reused from Part 2, A.8
"""
# %%
import numpy as np
import control as ct
from scipy.io import loadmat, savemat
from scipy.signal import butter
import matplotlib.pyplot as plt

from colorama import Fore, Style

from numpy import pi, sin, cos, arctan2, arctan, sqrt

from control_utils import (
    impulse_response_freq_domain,
    make_double_pole,
    make_lpf_butter,
)

USE_DB_IN_PLOTS = False

mat = loadmat("dat/MSD2025_P2_Plant_numden.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])

# Normalize to reduce the range of coefficients
num /= den[-1]
den /= den[-1]

plant = ct.tf(num, den)

# %%
poles = ct.poles(plant)
poles = poles[np.argsort(np.abs(poles))]
zeros = ct.zeros(plant)
zeros = zeros[np.argsort(np.abs(zeros))]

poles_gen = (pole for pole in poles)
plant_num_tfs = []
plant_den_tfs = []


def decompose_polynomial(roots):
    """
    Pass poles or zeros into roots; get decomposed into product of first and second order transfer functions, in order of frequency
    """
    root_tfs = []
    roots_gen = (root for root in roots)
    while True:
        try:
            root = next(roots_gen)
            if np.isreal(root):
                omega = -np.real(root)
                root_tfs.append(ct.tf([1, omega], omega))
            else:  # double pole
                omega_n = np.abs(root)
                zeta = -np.real(root) / omega_n
                root_tfs.append(1 / make_double_pole(omega_n, zeta))

                next_pole = next(roots_gen)
                assert next_pole == np.conj(root)

        except StopIteration:
            break
    return root_tfs


plant_den_tfs = decompose_polynomial(poles)
plant_num_tfs = decompose_polynomial(zeros)
plant_dc = ct.dcgain(plant)

for i, tf in enumerate(plant_num_tfs):
    print(f"Numerator TF {i}:")
    print(tf)

for i, tf in enumerate(plant_den_tfs):
    print(f"Denominator TF {i}:")
    print(tf)

# For feedforward
omega_lpf = 700 * 2 * pi
num, den = butter(4, omega_lpf, analog=True, output="ba")
lpf_worst = ct.tf(num, den)

omega_zero = float(sqrt(plant_num_tfs[0].den[0][0][0]))
# feedforward = 1 / plant_dc * plant_den_tfs[0] * plant_den_tfs[1] * plant_den_tfs[2] * make_double_pole(omega_zero, 0.025) * lpf
feedforward = (
    1
    / plant_dc
    * plant_den_tfs[0]
    * plant_den_tfs[1]
    * plant_den_tfs[2]
    * make_double_pole(omega_zero, 0.023)
    * lpf_worst
)

# For max sin feedforward - we can keep most of original system
omega_lpf = 1000 * 2 * pi
num, den = butter(4, omega_lpf, analog=True, output="ba")
lpf_maxsin = ct.tf(num, den)
feedforward_sin = (
    1
    / plant_dc
    * plant_den_tfs[0]
    * plant_den_tfs[1]
    * plant_den_tfs[2]
    * make_double_pole(omega_zero, 0.016)
    * lpf_maxsin
)


# # Previous feedforward for comparison
# def get_old_feedforward():
#     omega_lpf = 1000 * 2 * pi
#     feedforward_lpf = make_lpf_butter(omega_lpf, order=4)
#     feedforward = feedforward_lpf / plant
#     return feedforward


#  Feedforward control bode plot -------------
def add_bode_plot(tf, label, fig, freq):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    if USE_DB_IN_PLOTS :
        axm.semilogx(freq, mag_db, label=label)
    else:
        axm.loglog(omega, np.squeeze(mag), label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


# Plot feedforward control
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
freq = np.logspace(1, 5, 3600)

add_bode_plot(feedforward_sin, "Sinusoidal gain limited", fig, freq)
add_bode_plot(feedforward, "Worst-case limited", fig, freq)

axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
ylim = axp.get_ylim()  # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 180))
axp.set_ylim(ylim)

if USE_DB_IN_PLOTS:
    axm.set_ylabel("Magnitude [dB]")
else:
    axm.set_ylabel("Magnitude")
axm.set_xlim([freq[0], freq[-1]])

fig.set_size_inches(8, 6)
fig.savefig("img3/C.2.ff.png", bbox_inches="tight")
plt.close(fig)

# Show time response of impulse response
fig, axs = plt.subplots(3, height_ratios=[1, 1, 2])
time, yout = ct.impulse_response(feedforward, np.linspace(0, 0.02, 10000))
axs[0].plot(time, yout / np.max(np.abs(yout)), label="Impulse response scaled")

worst_case = np.zeros(len(time)) # worst case input
worst_case[yout > 0] = 1
worst_case[yout <= 0] = -1
worst_case = worst_case[::-1] * 10 # input scale
axs[0].set_ylabel("Impulse response\n(scaled)")
axs[0].grid()
axs[0].xaxis.set_inverted(True)
axs[0].set_xlabel("Time [s] (plotted backwards)")

ax = axs[1]
ax.plot(time, worst_case)
ax.set_ylabel("Worst case input [V]")

ax = axs[2]
time, yout = ct.forced_response(feedforward, time, worst_case)
ax.plot(time, yout, label="Max response")
ax.set_ylabel("Response [V]")
ax.set_xlabel("Time [s]")


# # find maximum magnitude and corresponding frequency
# mag, phase, _ = ct.frequency_response(feedforward, freq * 2 * pi)
# worst_freq = freq[np.argmax(mag)]
# max_sine = sin(worst_freq * 2 * pi * time)
# time, yout = ct.forced_response(feedforward, time, max_sine)
# ax.plot(time, yout, label="Sin response")
# ax.legend()

# axs[0].set_title("Worst case behaviour")
axs[1].grid()
axs[2].grid()
fig.set_size_inches(8, 8)
fig.tight_layout()
fig.savefig("img3/C.2.ff_worst_case.png", bbox_inches="tight")

# plt.show()

# %% Plot prefilter for completeness
fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
freq = np.logspace(1, 5, 3600)
omega = freq * 2 * pi
DELAY_TIME = 0.18e-3 # 0.18ms

for lpf, label in zip([lpf_maxsin, lpf_worst], [ "Sinusoidal gain limited",  "Worst-case limited"]):
    mag, phase, _ = ct.frequency_response(lpf, omega)
    mag_db = 20 * np.log10(mag)
    if USE_DB_IN_PLOTS :
        axm.semilogx(freq, mag_db, label=label)
    else:
        axm.loglog(freq, np.squeeze(mag), label=label)

    phase = np.unwrap(np.squeeze(phase))
    delay_phase = omega * DELAY_TIME
    phase -= delay_phase
    
    axp.semilogx(freq, phase * 180 / pi)

axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
axp.set_yticks(np.arange(-720, 1, 180))
axp.set_ylim(-720, 90)
axp.set_xlim(freq[0], 1e4)
axm.set_ylim(1e-5, 10**0.5)
fig.set_size_inches([8, 6])
fig.savefig('img3/C.3.prefilter.png', bbox_inches='tight')

plt.show()
