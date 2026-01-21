import numpy as np
import control as ct

from numpy import pi, sin, cos, arctan2, arctan, sqrt
from numpy.fft import rfft, rfftfreq, irfft
from scipy.signal import butter


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

def impulse_response_freq_domain(sys: ct.TransferFunction, end_time, n_points = 10000):
    """
    Impulse response with hopefully better numerical behaviour, which works in frequency domain, for transfer functions
    Works better on some transfer functions but fails on others (most notably, 1/s)
    
    Implemented by taking two impulses in series, with opposite sign, delayed by `end_time`; this is to ensure result is periodic (and CTFT approximation with DFT is valid). But this causes other problems...
    
    :param sys: TransferFunction
    :param end_time: Time to simulate until
    :param n_points: Number of points in time to consider
    """
    # window_len = n_points * 2
    sampling_time = end_time / n_points
    
    # Some different behaviour happens when n_opints is odd
    assert n_points % 2 == 0
    n_freqs = n_points + 1
    max_freq = 0.5 / sampling_time

    freq_test = np.linspace(0, max_freq, n_freqs)
    omega_test = freq_test * 2 * pi

    with np.errstate(divide='ignore', invalid='ignore'):
        mag, phase, _ = sys.frequency_response(omega=freq_test * 2 * pi)
    freq_resp_up = mag * np.exp(1j * phase)

    # frequency response of down impulse (not sign flipped yet)
    freq_resp_down = mag * np.exp(1j * (phase - omega_test * end_time))
    
    # total frequency response
    freq_resp = freq_resp_up - freq_resp_down

    # DC term is zero; manually set it due to poles at 0 making nan
    freq_resp[0] = 0 

    # Equation:
    time_irfft = np.linspace(0, sampling_time * n_points, n_points, endpoint=False)
    time_response = irfft(freq_resp, norm='backward') / sampling_time
    
    # roughly estimate DC offset
    dc_offset = time_response[n_points]
    return time_irfft, time_response[:n_points] + dc_offset
    