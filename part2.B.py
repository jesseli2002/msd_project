import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.signal import welch, csd
import scipy.signal as signal
from numpy.fft import rfft, rfftfreq
from numpy import pi, sin, cos
import control as ct

plt.rcParams['savefig.dpi'] = 600
# B.1: Commercial device; controller not accessible?

# %%
# B.2 Signals in time domain
data = loadmat("dat/MSD2025_P2_signals.mat")

r = np.squeeze(data["r"])
e = np.squeeze(data["e"])
u = np.squeeze(data["u"])
y = np.squeeze(data["y"])
ts = 30e-6

time = np.arange(len(r)) * ts
fig, axs = plt.subplots(4, sharex=True)
axs[0].plot(time, r, label="r")
axs[1].plot(time, e, label="e")
axs[2].plot(time, u, label="u")
axs[3].plot(time, y, label="y")
axs[0].set_ylabel("r")
axs[1].set_ylabel("e")
axs[2].set_ylabel("u")
axs[3].set_ylabel("y")
# axs[0].set_title("Time domain signals")
axs[-1].set_xlabel("Time [s]")
fig.set_size_inches(6, 6)
plt.savefig("img2/B.2.time_domain.png", bbox_inches="tight")

# Create zoomed version
axs[0].set_xlim(2.5, 3.1)
plt.savefig("img2/B.2.zoomed.png", bbox_inches="tight")
plt.close(fig)


# %%
# B.3 Closed loop frequency responses


def tf_cohere(input_t, output_t, ts):
    """
    Estimate transfer function and report coherence

    :param input_t: Input in time domain
    :param output_t: Output in time domain
    :param ts: Sample time
    :return: (freq, tf, coherence)
    """
    # u_fft = rfft(input_t)
    # y_fft = rfft(output_t)
    # freq = rfftfreq(len(u), ts)  # [Hz]

    window = signal.get_window("hann", 2048)
    # window = signal.get_window("hann", len(input_t))
    freq, Suu = welch(input_t, fs=1 / ts, window=window)
    _, Syy = welch(output_t, fs=1 / ts, window=window)
    _, Suy = csd(output_t, input_t, fs=1 / ts, window=window)

    tf = Syy / Suy

    # Syu = Suy.conj() # direct math implementation
    # coherence = Syu * Suy / (Syy * Suu)
    coherence = np.abs(Suy) ** 2 / (Syy * Suu)  # more practical numerically

    return freq, tf, coherence


freq, r2e_tf, r2e_cohere = tf_cohere(r, e, ts)
_, r2u_tf, r2u_cohere = tf_cohere(r, u, ts)
_, r2y_tf, r2y_cohere = tf_cohere(r, y, ts)

# Skip DC term
freq = freq[1:]
r2e_tf = r2e_tf[1:]
r2u_tf = r2u_tf[1:]
r2y_tf = r2y_tf[1:]
r2e_cohere = r2e_cohere[1:]
r2u_cohere = r2u_cohere[1:]
r2y_cohere = r2y_cohere[1:]


fig, (axm, axp, axc) = plt.subplots(3, sharex=True, height_ratios=[1, 1, 0.4])

tfs = [r2e_tf, r2u_tf, r2y_tf]
coheres = [r2e_cohere, r2u_cohere, r2y_cohere]
labels = ["r to e", "r to u", "r to y"]
for label, tf, coherance in zip(labels, tfs, coheres):
    mag = np.abs(tf)
    phase = np.angle(tf)
    axm.loglog(freq, mag, label=label)
    axp.semilogx(freq, phase * 180 / pi)
    axc.semilogx(freq, coherance)

axm.legend()
axm.set_xlim(10, 8e3)
axm.grid(which="both")
axp.grid()
axc.grid()
ylim = axp.get_ylim()
axp.set_yticks(np.arange(-180, 181, 90))
axp.set_ylim(ylim)
axm.set_ylim([7e-3, 3e0])
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axc.set_ylabel("Coherence")
axc.set_ylim(0, 1.05)
fig.axes[-1].set_xlabel("Frequency [Hz]")
fig.set_size_inches(8, 6)
fig.savefig("img2/B.3.freq_resp.png", bbox_inches="tight")
plt.close(fig)

# plt.show()

# %%
# B.4
_, e2u_tf_direct, e2u_cohere = tf_cohere(e, u, ts)
_, u2y_tf_direct, u2y_cohere = tf_cohere(u, y, ts)

# Skip DC term
e2u_tf_direct = e2u_tf_direct[1:]
u2y_tf_direct = u2y_tf_direct[1:]
e2u_cohere = e2u_cohere[1:]
u2y_cohere = u2y_cohere[1:]

C_tf = r2u_tf / r2e_tf
P_tf = r2y_tf / r2u_tf

fig, (axm, axp) = plt.subplots(2, sharex=True)

tfs = [C_tf, e2u_tf_direct]
labels = ["C = CS / S", "C = u / e"]
for label, tf in zip(labels, tfs):
    mag = np.abs(tf)
    phase = np.angle(tf)
    axm.loglog(freq, mag, label=label)
    axp.semilogx(freq, phase * 180 / pi)

axm.legend()
axm.set_xlim(10, 8e3)
axm.grid(which="both")
axp.grid()
ylim = axp.get_ylim()
axp.set_yticks(np.arange(-180, 181, 90))
axp.set_ylim(ylim)
axm.set_ylim([7e-3, 3e1])
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axp.set_xlabel("Frequency [Hz]")
fig.set_size_inches(6, 4.5)
fig.savefig("img2/B.4.C.png", bbox_inches="tight")
fig.axes[-1].set_xlabel("Frequency [Hz]")


fig, (axm, axp) = plt.subplots(2, sharex=True)
tfs = [P_tf, u2y_tf_direct]
labels = ["P = PCS / CS", "P = y / u"]
for label, tf in zip(labels, tfs):
    mag = np.abs(tf)
    phase = np.angle(tf)
    axm.loglog(freq, mag, label=label)
    axp.semilogx(freq, phase * 180 / pi)

axm.legend()
axm.set_xlim(10, 8e3)
axm.grid(which="both")
axp.grid()
ylim = axp.get_ylim()
axp.set_yticks(np.arange(-180, 181, 90))
axp.set_ylim(ylim)
axm.set_ylim([7e-3, 3e1])
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axp.set_xlabel("Frequency [Hz]")
fig.set_size_inches(6, 4.5)
fig.savefig("img2/B.4.P.png", bbox_inches="tight")
fig.axes[-1].set_xlabel("Frequency [Hz]")
plt.close(fig)
# plt.show()

# %%
# B.5 Compare plant model with measurements
plt.close('all')
mat = loadmat("dat/MSD2025_P2_Plant_numden.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])

# Normalize to reduce the range of coefficients
num /= den[-1]
den /= den[-1]
plant_model = ct.tf(num, den)

# Plot plant ----------------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])

def add_bode_plot(tf, label, fig):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    axm.loglog(freq, np.squeeze(mag), label=label)
    # axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)
add_bode_plot(plant_model, "Plant model", fig)

mag = np.abs(u2y_tf_direct)
phase = np.angle(u2y_tf_direct)
axm.loglog(freq, mag, label='Plant as measured')
axp.semilogx(freq, np.unwrap(phase) * 180 / pi)
axm.legend()
axm.set_xlim(10, 8e3)
axm.grid(which="both")
axp.grid()
ylim = axp.get_ylim()
axp.set_yticks(np.arange(-720, 181, 180))
axp.set_ylim(-720, 45)
axm.set_ylim([7e-3, 3e1])
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axp.set_xlabel("Frequency [Hz]")
fig.set_size_inches(8, 4.5)
fig.savefig("img2/B.5.P.png", bbox_inches="tight")

# %%
# B.7 recreate controller

# done in different file
