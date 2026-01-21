import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from scipy.signal import welch, csd
import scipy.signal as signal
from numpy.fft import rfft, rfftfreq
from numpy import pi, sin, cos, linalg as LA
import control as ct

from control_utils import (
    impulse_response_freq_domain,
    make_lpf_butter,
    make_integrator,
    make_double_pole,
    make_ll,
    make_notch,
)

plt.rcParams['savefig.dpi'] = 600

ts = 30e-6

def tf_cohere(input_t, output_t, ts, nperseg=2048):
    window = signal.get_window("hann", nperseg)
    fs = 1.0 / ts
    freq, Suu = signal.welch(input_t, fs=fs, window=window, nperseg=nperseg)
    _, Syy = signal.welch(output_t, fs=fs, window=window, nperseg=nperseg)
    _, Suy = signal.csd(output_t, input_t, fs=fs, window=window, nperseg=nperseg)
    tf = Syy / Suy
    coherence = np.abs(Suy)**2 / (Syy * Suu)
    return freq, tf, coherence

def compute_controller_tf(mat_path="dat/MSD2025_P2_signals.mat"):
    data = loadmat(mat_path)
    r = np.squeeze(data["r"])
    e = np.squeeze(data["e"])
    u = np.squeeze(data["u"])

    freq, r2e_tf, _ = tf_cohere(r, e, ts)
    _, r2u_tf, _ = tf_cohere(r, u, ts)
    _, e2u_tf_direct, e2u_cohere = tf_cohere(e, u, ts)

    # remove DC
    freq = freq[1:]
    C_tf = r2u_tf[1:] / r2e_tf[1:]
    C_tf_direct = e2u_tf_direct[1:]

    # return freq, C_tf
    return freq, C_tf_direct

freq, C_tf = compute_controller_tf()

# Outsource fitting to MATLAB
mat = loadmat("dat/part2_B7_C_tf_fit.done.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])
C_tf_model = ct.tf(num, den)

print(C_tf_model)

# Plot results --------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True)

def add_bode_plot(tf, label, fig, freq):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    axm.loglog(freq, np.squeeze(mag), label=label)
    # axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)

mag = np.abs(C_tf)
phase = np.angle(C_tf)
axm.loglog(freq, mag, label="Measured controller", marker='.', linestyle='None')
axp.semilogx(freq, np.unwrap(phase) * 180 / pi, marker='.', linestyle='None')

add_bode_plot(C_tf_model, "Recreated controller", fig, np.logspace(1, 4, 401))

axm.legend()
axm.set_xlim([10, 1e4])
axm.set_ylim(top=1e2)
axm.grid(which="both")
axp.grid()
axp.set_yticks(np.arange(-180, 181, 45))
axp.set_ylim(-180, 90)
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axp.set_xlabel("Frequency [Hz]")
fig.set_size_inches(8, 6)
fig.savefig('img2/B.7.fit.png', bbox_inches='tight')



# Compare to model
mat = loadmat("dat/MSD2025_P2_Plant_numden.mat")
num = np.squeeze(mat["num"])
den = np.squeeze(mat["den"])

plant = ct.tf(num, den)

fig, (axm, axp) = plt.subplots(2, 1, sharex=True)
ct.bode_plot(C_tf_model * plant, title='actual loop', Hz=True, dB=True, ax=(axm, axp))
ct.bode_plot(C_tf_model, title='cots controller', Hz=True, dB=True, ax=(axm, axp))
ct.bode_plot(plant, title='plant', Hz=True, dB=True, ax=(axm, axp))


plt.show()
