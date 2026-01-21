"""
Just trying to mess around with feedforward
"""
# %%
import numpy as np
import control as ct
from scipy.io import loadmat, savemat
from scipy.signal import butter
import matplotlib.pyplot as plt

from colorama import Fore, Style

from numpy import pi, sin, cos, arctan2, arctan, sqrt

from control_utils import impulse_response_freq_domain, make_lpf_butter, make_integrator, make_double_pole, make_ll, make_notch

plt.rcParams['savefig.dpi'] = 600

mat = loadmat("dat/MSD2025_P2_Plant_numden.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])

# Normalize to reduce the range of coefficients
num /= den[-1]
den /= den[-1]

plant = ct.tf(num, den)


# Plot plant ----------------------
freq = np.logspace(2, 4, 1200) 

def add_bode_plot(tf, label, fig):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    # axm.loglog(omega, np.squeeze(mag), label=label)
    axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


# Define controller
C1 = 2750 / ct.tf('s')
C2 = 1/make_double_pole(4610, 0.01) * make_double_pole(8000, 1)
C3 = make_notch(6300, gain=5)
C4 = make_lpf_butter(28000)
C5 = make_ll(300 * 2 * pi, 20 * pi / 180)
C = ( # for better numerical stability
    ct.tf2ss(C1)
    * ct.tf2ss(C2)
    * ct.tf2ss(C3)
    * ct.tf2ss(C4)
    * ct.tf2ss(C5)
)

print(f"Skewed notch:\n{C2}")
print(f"Notch2:\n{C3}")
print(f"LPF:\n{C4}")
print(f"Lead lag:\n{C5}")

loop = C * ct.tf2ss(plant)


freq = np.logspace(2, 4.5, 1200) 
omega = freq * 2 * pi

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

# Step response --------------------------
# Things are better numerically in state space
C_ss = ct.tf2ss(C)
plant_ss = ct.tf2ss(plant)
loop_ss = C_ss * plant_ss
feedback_ss = ct.feedback(loop_ss)

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

# try modifying with prefilter
fffb_ss_prefilter = ct.feedback(plant_ss, C_ss) * feedforward_ss + ct.series(feedforward_lpf, feedback_ss)

# feedforward only control
ff_control_ss = plant_ss * feedforward_ss

fig, ax = plt.subplots()
time = np.linspace(0, 0.02, 2000 + 1)
time, yout = ct.step_response(feedback_ss, time)
ax.plot(time, yout, label="Feedback (FB) only")
time, yout = ct.step_response(fffb_ss, time)
ax.plot(time, yout, label="Feedforward (FF) and FB")
time, yout = ct.step_response(fffb_ss_prefilter, time)
ax.plot(time, yout, label='FF + FB w/ prefilter')
# time, yout = ct.impulse_response(feedback_ss, time)
# ax.plot(time, yout / np.max(np.abs(yout)), label="Feedback impulse, scaled")
ax.legend()
ax.set_xlim([0, 5e-3])
ax.grid()
ax.set_xlabel("Time [s]")
ax.set_ylabel("Response")
fig.set_size_inches(6, 4)
fig.savefig('img2/A.5.step_response_ff.png')
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

# %% A.7 Transfer functions
# Recall all transfer functions and rename for alignment with LaTeX
C_ss = C_ss
P_ss = plant_ss
F_ss = feedforward_ss

# [f]eed[b]ack transfer functions
ref_track_fb = ct.feedback(P_ss * C_ss)
dist_rej_fb = ct.feedback(P_ss, C_ss)
noise_rej_fb = ct.feedback(1, P_ss * C_ss) # n -> y
noise_rej_fb2 = ct.feedback(P_ss * C_ss) # n -> x

# [f]eed[f]orward and [f]eed[b]ack reference tracking
ref_track_fffb = ref_track_fb + F_ss * ct.feedback(P_ss, C_ss)

# Plot reference tracking 
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
freq = np.logspace(1, 5, 1200)

add_bode_plot(ref_track_fb, "Feedback only", fig)
add_bode_plot(ref_track_fffb, "Feedforward and feedback", fig)

axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
axm.set_ylabel("Magnitude [dB]")
axm.set_xlim([freq[0], freq[-1]])
ylim = axp.get_ylim() # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 180))
axp.set_ylim(ylim)
axm.set_ylim(bottom=-100)

fig.set_size_inches(8, 6)
fig.savefig("img2/A.7.reference_tracking.png", bbox_inches="tight")
plt.close(fig)


# Plot other transfer functions
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
freq = np.logspace(1, 5, 1200)

add_bode_plot(dist_rej_fb, r"d$\to$y", fig)
add_bode_plot(noise_rej_fb, r"n$\to$y", fig)
add_bode_plot(noise_rej_fb2, r"n$\to$-x", fig)

axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
ylim = axp.get_ylim() # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 180))
axp.set_ylim(ylim)
axm.set_ylabel("Magnitude [dB]")
axm.set_xlim([freq[0], freq[-1]])
axm.set_ylim(bottom=-100)

fig.set_size_inches(8, 6)
fig.savefig("img2/A.7.disturb_noise.png", bbox_inches="tight")
plt.close(fig)



# %% A.8 new performance limits
# TODO: feedforward design is done separately in part2.A.8_ff.py
C2 = C * 0.65

# Compare controller to previous one in frequency domain
fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
freq = np.logspace(1, 5, 1200)

add_bode_plot(C_ss, "Original controller", fig)
add_bode_plot(C2, "New controller", fig)

axm.legend()
axm.grid(which="both")
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
ylim = axp.get_ylim() # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 90))
axp.set_ylim(ylim)
axm.set_ylabel("Magnitude [dB]")
axm.set_xlim([freq[0], freq[-1]])
fig.set_size_inches(8, 6)
fig.savefig("img2/A.8.fb_freq.png", bbox_inches="tight")
plt.close(fig)

C2_ss = ct.tf2ss(C2)
S = ct.feedback(1, C2_ss * plant_ss)
T = ct.feedback(C2_ss * plant_ss)

omega_test = np.logspace(1, 5, 300 * 4 + 1)
freq_test = omega_test / (2 * pi)
mag, _, _ = ct.frequency_response(S, omega_test)
sensitivity_hinf = np.max(mag)

fig, (axm, axp) = plt.subplots(2, sharex=True, height_ratios=[2,1])
add_bode_plot(S, label="Sensitivity", fig=fig)
add_bode_plot(T, label="Complementary sensitivity", fig=fig)
axm.set_title(f"max(S(jω))={sensitivity_hinf:.3f}")
axm.legend()
axm.grid(which="both")
axm.set_ylim(bottom=-150)
axp.grid(which="both")
axp.set_xlabel("Frequency [Hz]")
axp.set_ylabel("Phase [deg]")
ylim = axp.get_ylim() # store so set yticks doesn't affect it
axp.set_yticks(np.arange(-720, 1, 180))
axp.set_ylim(ylim)
axm.set_ylabel("Magnitude [dB]")
axm.set_xlim([freq[0], freq[-1]])
fig.set_size_inches(8, 6)
fig.savefig("img2/A.8.S_T_freq.png", bbox_inches="tight")
plt.close(fig)
