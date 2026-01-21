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


C_ll = make_ll(6200, 50 * pi / 180)

# Define controller
# C = 1/make_notch(6100, gain=10, Q2=3)
# C = make_notch(4650, gain=20, Q2=2.5) * make_integrator(3000) * make_ll(4800, 50 * pi / 180) / 3

# OK controller but not great
# C = 1000 / 0.3 / ct.tf('s') * make_notch(4650, gain=20, Q2=2.5)

# Gets a bandwidth around 1000 rad/s
C = (
    1000
    * 2
    / ct.tf("s")
    * make_notch(4650, gain=40, Q2=1)
    * make_notch(6300, gain=10, Q2=2.5)
)

# Try to push bandwidth higher; this gets to 200 Hz
C = (
    2000 
    / ct.tf("s")
    * make_notch(4610, gain=50, Q2=1)
    * make_notch(6300, gain=10, Q2=2.5)
    * make_ll(500 * 2 * pi, 30 * pi / 180)
)



loop = C * plant

# Evaluate controller modulus margin
sensitivity = 1 / (1 + loop)
omega_test = np.logspace(1, 5, 300 * 4 + 1)
freq_test = omega_test / (2 * pi)
mag, _, _ = ct.frequency_response(sensitivity, omega_test)
modulus_margin = np.max(mag)
modulus_margin_db = 20 * np.log10(modulus_margin)
print(f"max(S(jw)) = {modulus_margin_db=} dB @ w = {freq_test[np.argmax(mag)]:.1f} Hz")

# Plot open loop results ----------------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
omega = np.logspace(1, 5, 1200)
freq = omega / (2 * pi)


def add_bode_plot(tf, label, fig):
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    # axm.loglog(omega, np.squeeze(mag), label=label)
    axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


add_bode_plot(plant, "Plant", fig)
add_bode_plot(C, "Controller", fig)
add_bode_plot(loop, "Loop", fig)


axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
axm.set_ylabel("Magnitude [dB]")

fig.set_size_inches(6, 4)
fig.savefig("img2/A.1.open_loop.png", bbox_inches="tight")

plt.show()
plt.close(fig)
print("Exiting early!!")
exit()

# Plot sensitivity ----------------------------
S = 1 / (1 + loop)
T = ct.feedback(loop)

plt.close("all")
ct.bode_plot(S, title="Sensitivity functions", label="Sensitivity")
ct.bode_plot(T, label="Complementary sensitivity")
fig = plt.gcf()
fig.set_size_inches(6, 4)
# fig.set_size_inches(15, 10)
fig.savefig("img2/A.2.sensitivity.png", bbox_inches="tight", dpi=200)
# plt.close(fig)


plt.show() # ======================================== TO REMOVE
exit()

# Plot margins ----------------------------
# gm, pm, wg, wp = ct.margin(loop)
ct.bode_plot(loop, display_margins=True, dB=True, Hz=True)
fig = plt.gcf()

fig.set_size_inches(6, 4)
fig.savefig("img2/A.2.margins.png", bbox_inches="tight")
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

fig.set_size_inches(6, 4)
fig.savefig("img2/A.3.step_response.png")
plt.close(fig)

# Get results
step_info = ct.step_info(
    feedback_ss, RiseTimeLimits=(0.1, 0.9), SettlingTimeThreshold=0.02
)
print(f"{Fore.CYAN}Feedback control design (A.3){Style.RESET_ALL}")
print(f"{Fore.CYAN}Rise time: {step_info['RiseTime']}{Style.RESET_ALL}")
print(f"{Fore.CYAN}Settling time: {step_info['SettlingTime']}{Style.RESET_ALL}")
print(f"{Fore.CYAN}Overshoot: {step_info['Overshoot']:.2f}%{Style.RESET_ALL}")

#  Feedforward control ----------------
omega_lpf = 6000
num, den = butter(4, omega_lpf, analog=True, output="ba")
feedforward = ct.tf(num, den) / plant
# num = [omega_lpf ** 2] # critically damped system
# den = [1, omega_lpf * 1, omega_lpf ** 2]
# feedforward = ct.tf(num, den) * ct.tf(num, den) / plant

# fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
# omega = np.logspace(1, 5, 1200)
# add_bode_plot(feedforward, 'feedforward')
# add_bode_plot(plant, 'plant')
# plt.show()

feedforward_ss = ct.tf2ss(feedforward)

# [f]eed[f]orward and [f]eed[b]ack state space model
fffb_ss = ct.feedback(plant_ss, C_ss) * feedforward_ss + feedback_ss

# feedforward only control
ff_control_ss = plant_ss * feedforward_ss

fig, ax = plt.subplots()
time, yout = ct.step_response(feedback_ss)
ax.plot(time, yout, label="Feedback only")
time, yout = ct.step_response(fffb_ss, time)
ax.plot(time, yout, label="Feedforward and feedback")
# time, yout = ct.step_response(ff_control_ss, time)
# ax.plot(time, yout, label='Feedforward only')
ax.legend()
ax.set_xlim([0, 5e-3])
ax.grid()

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

plt.show()

# %%
