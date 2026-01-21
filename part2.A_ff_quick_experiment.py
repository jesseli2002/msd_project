"""
Hypothesis: Overshoot from feedforward is due to overactive integral control
Test: Reduce low freuqency integrator
Result: No, since low frequency integrartor is in different frequency band; this is high frequency behaviour
"""

# %%
import numpy as np
import control as ct
from scipy.io import loadmat, savemat
from scipy.signal import butter
import matplotlib.pyplot as plt

from colorama import Fore, Style

from numpy import pi, sin, cos, arctan2, arctan, sqrt

from control_utils import impulse_response_freq_domain

mat = loadmat("dat/MSD2025_P2_Plant_numden.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])

# Normalize to reduce the range of coefficients
num /= den[-1]
den /= den[-1]

plant = ct.tf(num, den)


# Plot plant ----------------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
freq = np.logspace(2, 4, 1200)


def add_bode_plot(tf, label, fig):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    # axm.loglog(omega, np.squeeze(mag), label=label)
    axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


add_bode_plot(plant, "Plant", fig)

axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
axm.set_ylabel("Magnitude [dB]")
axm.set_xlim([freq[0], freq[-1]])
axp.set_yticks(np.arange(-270, 1, 90))

fig.set_size_inches(8, 4)
fig.savefig("img2/A.1.plant.png", bbox_inches="tight")
plt.close(fig)


def make_integrator(omega_corner):
    """
    Create integrator action
    :param omega_corner: Corner frequency
    """
    return ct.tf([1, omega_corner], [1, 0])


def make_double_pole(omega_n, zeta):
    """
    Create double pole (second order system)
    """
    return ct.tf([omega_n**2], [1, 2 * zeta * omega_n, omega_n**2])


def make_ll(omega_ll, phi_m):
    """
    Create lead controller

    :param omega_ll: Center frequency [rad/s]
    :param phi_m: Phase at center frequency
    """
    alpha = (1 - sin(phi_m)) / (1 + sin(phi_m))
    T = 1 / (omega_ll * sqrt(alpha))
    C_ll = ct.tf([T, 1], [alpha * T, 1])
    return C_ll


def make_notch(omega_n, gain=10, Q2=10):
    """
    Docstring for make_notch

    :param omega_n: Description
    :param Q: Q factor; gain of notch.
    :param Q2: Measure of width. Higher number = sharper peak
    """
    Q1 = gain * Q2
    C_notch = ct.tf([1, omega_n / Q1, omega_n**2], [1, omega_n / Q2, omega_n**2])

    return C_notch


def make_lpf_butter(omega_lpf, order=2):
    num, den = butter(order, omega_lpf, analog=True, output="ba")
    return ct.tf(num, den)


# Define controller
C1 = 2750 / ct.tf("s")
C2 = 1 / make_double_pole(4610, 0.01) * make_double_pole(8000, 1)
C3 = make_notch(6300, gain=5)
C4 = make_lpf_butter(28000)
C5 = make_ll(300 * 2 * pi, 20 * pi / 180)

integral_adjust_w = 100 * 2 * pi
integral_adjust = ct.tf(
    [1, integral_adjust_w / 10],
    [1, integral_adjust_w],
)
C = (
    ct.tf2ss(C1)
    * ct.tf2ss(C2)
    * ct.tf2ss(C3)
    * ct.tf2ss(C4)
    * ct.tf2ss(C5)
    * ct.tf2ss(integral_adjust)
)

print(f"Skewed notch:\n{C2}")
print(f"Notch2:\n{C3}")
print(f"LPF:\n{C4}")
print(f"Lead lag:\n{C5}")

loop = C * ct.tf2ss(plant)


# Plot open loop results ----------------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
freq = np.logspace(2, 4.5, 1200)
omega = freq * 2 * pi

add_bode_plot(plant, "Plant", fig)
add_bode_plot(C, "Controller", fig)
add_bode_plot(loop, "Loop", fig)


axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
axm.set_ylabel("Magnitude [dB]")
axm.set_ylim(bottom=-100)
axm.set_xlim([freq[0], freq[-1]])
ylim = axp.get_ylim()  # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 90))
axp.set_ylim(ylim)

fig.set_size_inches(8, 6)
fig.savefig("img2/A.1.open_loop.png", bbox_inches="tight")
# plt.show()
plt.close(fig)

# Plot sensitivity ----------------------------
# Evaluate controller modulus margin at same time
# S = 1 / (1 + loop)
S = ct.feedback(1, loop)
T = ct.feedback(loop)

freq_test = np.logspace(1, 5, 300 * 5 + 1)
omega_test = 2 * pi * freq_test
mag, _, _ = ct.frequency_response(S, omega_test)
modulus_margin = np.max(mag)
modulus_margin_db = 20 * np.log10(modulus_margin)

print(f"max(S(jw)) = {modulus_margin_db=} dB @ w = {freq_test[np.argmax(mag)]:.1f} Hz")

plt.close("all")
ct.bode_plot(S, label="Sensitivity")
ct.bode_plot(T, label="Complementary sensitivity")
fig = plt.gcf()
fig.set_size_inches(8, 6)
fig.axes[0].legend(loc="lower left")
fig.axes[0].set_xlim([10**2, 10**5])
fig.axes[0].set_title(
    f"Sensitivity functions; max(S(j$\\omega$))={float(modulus_margin_db):.2f} dB @ $\\omega$ = {freq_test[np.argmax(mag)]:.1f} Hz"
)
fig.savefig("img2/A.2.sensitivity.png", bbox_inches="tight", dpi=200)
plt.close(fig)


# Plot margins ----------------------------
gm, pm, wg, wp = ct.margin(loop)
ct.bode_plot(loop, display_margins=True, dB=True, Hz=True, title="")
fig = plt.gcf()
fig.axes[0].set_title(
    f"Gain margin: {20*np.log10(gm):.2f} dB at {wg/(2*pi):.1f} Hz; Phase margin: {pm:.2f}° at {wp/(2*pi):.1f} Hz"
)
fig.set_size_inches(8, 6)
fig.savefig("img2/A.2.margins.png", bbox_inches="tight")
# plt.show()
plt.close(fig)


# Step response --------------------------
# Things are better numerically in state space
C_ss = ct.tf2ss(C)
plant_ss = ct.tf2ss(plant)
loop_ss = C_ss * plant_ss
feedback_ss = ct.feedback(loop_ss)

time, yout = ct.step_response(feedback_ss)
fig, ax = plt.subplots()
ax.plot(time, yout, label="Closed loop")
ax.set_xlim([0, 5e-3])
ax.grid()
ax.set_xlabel("Time [s]")
ax.set_ylabel("Response")

fig.set_size_inches(6, 4)
fig.savefig("img2/A.3.step_response.png")
plt.close(fig)

# Get results
step_info = ct.step_info(
    feedback_ss, RiseTimeLimits=(0.1, 0.9), SettlingTimeThreshold=0.02
)
print(f"{Fore.CYAN}Feedback control design (A.3){Style.RESET_ALL}")
print(f"{Fore.CYAN}Rise time (10-90%): {step_info['RiseTime']}{Style.RESET_ALL}")
print(f"{Fore.CYAN}Settling time (2%): {step_info['SettlingTime']}{Style.RESET_ALL}")
print(f"{Fore.CYAN}Overshoot: {step_info['Overshoot']:.2f}%{Style.RESET_ALL}")

#  Feedforward control ----------------
omega_lpf = 1000 * 2 * pi
feedforward_lpf = make_lpf_butter(omega_lpf, order=4)
feedforward = feedforward_lpf / plant
print(f"Feedforward low pass filter: \n{feedforward_lpf}")
# feedforward = 1/ct.dcgain(plant) * make_notch(4650, gain=40, Q2=1) * make_notch(6300, gain=10, Q2=2.5) * make_lpf_butter(omega_lpf, order=2)

# Plot feedforward bode plot
# fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
# omega = np.logspace(1, 5, 1200)
# add_bode_plot(feedforward, 'feedforward', fig)
# add_bode_plot(plant, 'plant', fig)
# add_bode_plot(feedforward * plant, 'combined', fig)
# axm.legend()
# plt.show()
# plt.close(fig)

feedforward_ss = ct.tf2ss(feedforward)

# [f]eed[f]orward and [f]eed[b]ack [s]tate [s]pace model
fffb_ss = ct.feedback(plant_ss, C_ss) * feedforward_ss + feedback_ss

# feedforward only control
ff_control_ss = plant_ss * feedforward_ss

fig, ax = plt.subplots()
time, yout = ct.step_response(feedback_ss)
ax.plot(time, yout, label="Feedback only")
time, yout = ct.step_response(fffb_ss, time)
ax.plot(time, yout, label="Feedforward and feedback")
time, yout = ct.step_response(ff_control_ss, time)
ax.plot(time, yout, label="Feedforward only")
ax.legend()
ax.set_xlim([0, 5e-3])
ax.grid()
fig.set_size_inches(8, 6)
fig.savefig("img2/A.5.step_response_ff.png")
plt.show()
plt.close(fig)

# Get results --------------
step_info = ct.step_info(fffb_ss, RiseTimeLimits=(0.1, 0.9), SettlingTimeThreshold=0.02)
print(f"{Fore.LIGHTGREEN_EX}Feedforward control design (A.6){Style.RESET_ALL}")
print(f"{Fore.LIGHTGREEN_EX}Rise time: {step_info['RiseTime']}{Style.RESET_ALL}")
print(
    f"{Fore.LIGHTGREEN_EX}Settling time: {step_info['SettlingTime']}{Style.RESET_ALL}"
)
print(f"{Fore.LIGHTGREEN_EX}Overshoot: {step_info['Overshoot']:.2f}%{Style.RESET_ALL}")

# %% A.8 new performance limits
C2 = C * 0.9

sensitivity = 1 / (1 + C2 * plant)
omega_test = np.logspace(1, 5, 300 * 4 + 1)
freq_test = omega_test / (2 * pi)
mag, _, _ = ct.frequency_response(sensitivity, omega_test)
sensitivity_hinf = np.max(mag)

fig, axs = plt.subplots(2)
add_bode_plot(sensitivity, label="Sensitivity", fig=fig)
axs[0].set_title(f"A.8, S max={sensitivity_hinf:.3f}")

# Design limited-magnitude feedforward
plant_dc_gain = ct.dcgain(plant)
lpf_num, lpf_den = butter(2, 6000, analog=True, output="ba")
feedforward2 = (
    ct.tf(lpf_num, lpf_den) / plant_dc_gain * make_notch(4650, gain=3 * plant_dc_gain)
)

fig, ax = plt.subplots(2)
ct.bode_plot(feedforward2, title="A.8 new feedforward", Hz=True)

# Show time response of impulse response
fig, ax = plt.subplots()
time, yout = ct.impulse_response(feedforward2, np.linspace(0, 0.1, 10000))
ax.plot(time, yout, label="Impulse response")

worst_case = np.zeros(len(time))
worst_case[yout > 0] = 1
worst_case[yout <= 0] = -1
time, yout = ct.forced_response(feedforward2, time, worst_case)
ax.plot(time, yout, label="Max response")

max_sine = sin(4650 * time)
time, yout = ct.forced_response(feedforward2, time, max_sine)
ax.plot(time, yout, label="Sin response")
ax.legend()
ax.set_title("Problem Maximizer")

# plt.show()

# %%
