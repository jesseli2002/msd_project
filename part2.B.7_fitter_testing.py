import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
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

# test_tf = make_double_pole(1000, 1) * ct.tf([1, 100], [1])
test_tf = ct.tf([1], [1, 0])

omega_test = np.logspace(1, 4, 201)
mag, phase, _ = ct.frequency_response(test_tf, omega_test)

C_tf = mag * np.exp(1j * phase)
C_tf = np.squeeze(C_tf)

omega_fit = omega_test

rng = np.random.default_rng()
# gain_noise = rng.normal(1, 0.05, size=C_tf.shape)
gain_noise = rng.uniform(0.9, 1.1, size=C_tf.shape)
C_tf_fit = C_tf * gain_noise


def levy_fit(omega, tf_data, n_poles, n_zeros, max_iter=30, rtol=1e-6):
    """
    Use Levy method to approximate transfer function.
    See http://dx.doi.org/10.5772/intechopen.71358, Rational Fitting Techniques for the Modeling of Electric Power Components and Systems Using MATLAB Environment, summarized in https://dsp.stackexchange.com/a/73817

    Keyword:  Sanathanan-Koerner (SK) iteration (section 3.3.2 in paper)

    Problem: Right now, this doesn't work for poles at zero, since in the formulation, denominator has constant term 1. 

    :param omega: Frequencies that tf_data corresponds to
    :param tf_data: Transfer function data (complex values)
    :param n_poles: Number of poles in the model
    :param n_zeros: Number of zeros in the model
    """
    # coefficients stored as all numerators, then all denominators, lowest order first, denominator leading coeff = 1
    coeffs = np.zeros(n_poles + n_zeros + 1)

    def eval_residual(coeffs):
        numer_coeffs = coeffs[: n_zeros + 1]
        denom_coeffs = np.r_[1, coeffs[n_zeros + 1 :]]  # leading 1
        numer = np.polynomial.polynomial.polyval(1j * omega, numer_coeffs)
        denom = np.polynomial.polynomial.polyval(1j * omega, denom_coeffs)
        residual = tf_data * denom - numer
        return residual

    def eval_denom(coeffs):
        # Used in SK iteration
        denom_coeffs = np.r_[1, coeffs[n_zeros + 1 :]]  # leading 1
        denom = np.polynomial.polynomial.polyval(1j * omega, denom_coeffs)
        return denom

    # Construct A matrix for least squares, or at least baseline
    s = 1j * omega
    s_poly_num = s[:, None] ** np.arange(n_zeros + 1)[None]
    s_poly_den = s[:, None] ** np.arange(1, n_poles + 1)[None]
    A_base = np.concatenate([s_poly_num, -tf_data[:, None] * s_poly_den], axis=1)

    b_base = tf_data

    # Want only real coefficients => split into real and imaginary parts

    converged = False
    for i in range(max_iter):  # TODO - until convergence
        denom = eval_denom(coeffs)
        A_complex = A_base / denom[:, None]
        b_complex = b_base / denom

        # Want only real coefficients => split into real and imaginary parts
        A = np.concatenate([A_complex.real, A_complex.imag], axis=0)
        b = np.concatenate([b_complex.real, b_complex.imag], axis=0)

        new_coeffs, _, _, _ = LA.lstsq(A, b)

        if LA.norm(new_coeffs - coeffs) < LA.norm(coeffs) * rtol:
            converged = True  # may as well use new coeffs since we calculated them
        coeffs = new_coeffs

        if converged:
            break
    else:
        print("Warning: levy_fit did not converge within max_iter")

    numer_coeffs = coeffs[: n_zeros + 1]
    denom_coeffs = np.r_[1, coeffs[n_zeros + 1 :]]
    # recall control library has highest order first
    return ct.tf(numer_coeffs[::-1], denom_coeffs[::-1])



fit_tf = levy_fit(omega_test, C_tf_fit, n_poles=1, n_zeros=0)


# Plot results --------------
fig, (axm, axp) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
freq = omega_test / (2 * pi)


def add_bode_plot(tf, label, fig):
    omega = freq * 2 * pi
    axm, axp = fig.axes
    mag, phase, _ = ct.frequency_response(tf, omega)
    mag_db = 20 * np.log10(mag)
    axm.loglog(freq, np.squeeze(mag), label=label)
    # axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


mag = np.abs(C_tf_fit)
phase = np.angle(C_tf_fit)
axm.loglog(freq, mag, label="Measured controller", linestyle="none", marker='.')
axp.semilogx(freq, np.unwrap(phase) * 180 / pi, linestyle="none", marker='.')

add_bode_plot(fit_tf, "Recreated controller", fig)

axm.legend()
axm.grid(which="both")
axp.grid()
ylim = axp.get_ylim()
axm.set_ylabel("Magnitude")
axp.set_ylabel("Phase [deg]")
axp.set_xlabel("Frequency [Hz]")
fig.set_size_inches(8, 4.5)
plt.show()
