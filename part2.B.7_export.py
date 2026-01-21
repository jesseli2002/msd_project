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

# Limit frequency to 9kHz for fitting
idx = np.searchsorted(freq, 9000)
freq_fit = freq[:idx]
omega_fit = 2 * pi * freq_fit
C_tf_fit = C_tf[:idx]


# Outsource fitting to MATLAB
savemat("dat/part2_B7_C_tf_fit.mat", {"omega_fit": omega_fit, "C_tf_fit": C_tf_fit})


exit()
