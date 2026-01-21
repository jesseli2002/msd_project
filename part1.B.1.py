import scipy
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import rfft, rfftfreq
from numpy import pi, sin, cos
import control as ct

# Load data
data = scipy.io.loadmat("dat/part1.B.chirp.mat")

u = np.squeeze(data["u"])
t = np.squeeze(data["t"])
ts = np.squeeze(data["ts"])
y = np.squeeze(data["y"])
fmin = np.squeeze(data["fmin"])
fmax = np.squeeze(data["fmax"])

# fig, ax = plt.subplots()
# ax.plot(t, y)

# plt.show()

u_fft = rfft(u)
y_fft = rfft(y)
freq = rfftfreq(len(u), ts)  # [Hz]

omega = freq * 2 * pi

Syy = y_fft * y_fft.conj()
Suu = u_fft * u_fft.conj()
Suy = u_fft * y_fft.conj()
Syu = y_fft * u_fft.conj()
tf = Syy / Suy 

coherence = Syu * Suy / (Syy * Suu)


# # Plot results
SUBSAMPLE_RATE = 50
# fig, (axmag, axphase) = plt.subplots(2, sharex=True)
# axmag.loglog(freq[::SUBSAMPLE_RATE], np.abs(tf)[::SUBSAMPLE_RATE], marker='o')
# axphase.semilogx(freq[::SUBSAMPLE_RATE], np.unwrap(np.angle(tf[::SUBSAMPLE_RATE])) * 180 / np.pi)
# axphase.set_xlim(fmin, fmax)
# axmag.grid()
# axphase.grid()
# axphase.set_xlabel("Frequency [Hz]")
# axphase.set_ylabel("Phase [deg]")
# axmag.set_ylabel("Magnitude")
# plt.show()


# Estimate transfer function
# From inspection, looks like double pole, double zero, double pole, plus some time delay
def tf_model(params, omega):
    """ """
    kp, pomega0, pzeta0, pomega1, pzeta1, pomega2, zomega0, zzeta0, td = params

    s = 1j * omega
    # fmt: off
    pred_response = (
        kp # gain
        * pomega0 ** 2 / (s ** 2 + 2 * pzeta0 * pomega0 * s + pomega0 ** 2) # double pole 0
        * pomega1 ** 2 / (s ** 2 + 2 * pzeta1 * pomega1 * s + pomega1 ** 2) # double pole 1
        * pomega2 / (s + pomega2) # pole 2
        * (s ** 2 + 2 * zzeta0 * zomega0 * s + zomega0 ** 2) / zomega0 ** 2 # double zero 0
    ) * np.exp(-s * td)
    return pred_response

# Estimate transfer function
# From inspection, looks like double pole, double zero, double pole, plus some time delay
def tf_model_ct(params):
    """ """
    kp, pomega0, pzeta0, pomega1, pzeta1, pomega2, zomega0, zzeta0, td = params

    s = ct.tf('s')
    # fmt: off
    pred_response = (
        kp # gain
        * pomega0 ** 2 / (s ** 2 + 2 * pzeta0 * pomega0 * s + pomega0 ** 2) # double pole 0
        * pomega1 ** 2 / (s ** 2 + 2 * pzeta1 * pomega1 * s + pomega1 ** 2) # double pole 1
        * pomega2 / (s + pomega2) # pole 2
        * (s ** 2 + 2 * zzeta0 * zomega0 * s + zomega0 ** 2) / zomega0 ** 2 # double zero 0
    ) 
    return pred_response


# Reference set of parameters
x0_ref = np.array(
    [
        4.53853710e-01,
        4.59784851e03,
        1.03797057e-02,
        6.24658235e03,
        9.51312451e-03,
        6e3,
        6.12734133e03,
        3.23902434e-03,
        # 2.10633875e-04,
        1.95e-04,
    ]
)


def tf_eval(x, omega, response):
    """
    Evaluation function for scipy.optimize.minimize()
    :param x: All parameters to model
    :param omega: Frequencies in data (rad/s)
    :param response: Measured response
    """
    # Meaning of parameters:
    # :param kp: Overall gain
    # :param pomega0: Natural frequency of pole 0
    # :param pzet0: Damping ratio of pole 0
    # :param pomega1: Natural frequency of pole 1
    # :param pzet1: Damping ratio of pole 1
    # :param zomega1: Natural frequency of zero 0
    # :param zzeta1: Damping ratio of zero 0
    # :param td: Time delay
    # All frequencies in rad/s

    pred_response = tf_model(x * x0_ref, omega)
    error_tf = pred_response / response
    gain_err = np.abs(np.log(np.abs(error_tf)))
    # phase_err = np.abs(np.unwrap(np.angle(error_tf)))

    total_err = gain_err
    # total_err = gain_err + phase_err

    # Crude integration over log-omega space
    log_omega = np.log(omega)
    log_omega_interval = log_omega[1:] - log_omega[:-1]
    error = np.sum(total_err[:-1] ** 2 * log_omega_interval)

    # error = np.sum(np.log(np.abs(pred_response / response)) ** 2)
    # error = np.sum(np.abs(pred_response - response) ** 2)
    return error


# Skip past low frequencies and high freq with bad data
START_IDX = np.searchsorted(freq, 10)
END_IDX = np.searchsorted(freq, 5000)
res = scipy.optimize.minimize(
    tf_eval,
    x0=np.ones_like(x0_ref),
    args=(
        omega[START_IDX:END_IDX:SUBSAMPLE_RATE],
        tf[START_IDX:END_IDX:SUBSAMPLE_RATE],
    ),
)
if not res.success:
    print(f"Failed to optimize!")
    print(res.message)
    exit()

print(f"final params {res.x * x0_ref}")
tf_predicted = tf_model(res.x * x0_ref, omega)

tf_model_expr = tf_model_ct(res.x * x0_ref)
print(f"TF model expr: ")
tf_model_expr.num
print(tf_model_expr)

# tf_predicted = tf_model(x0_ref, omega)
# tf_predicted2 = tf_model(x0_ref, omega)

fig, (axmag, axphase, axcohere) = plt.subplots(3, sharex=True, height_ratios=[1, 1, 0.3])
axmag.loglog(freq, np.abs(tf), label="data")
axphase.semilogx(freq, np.unwrap(np.angle(tf)) * 180 / np.pi)

axmag.loglog(freq, np.abs(tf_predicted), label="model")
axphase.semilogx(freq, np.unwrap(np.angle(tf_predicted)) * 180 / np.pi, label="model")
axcohere.semilogx(freq, coherence)

axphase.set_xlim(fmin, fmax)
axcohere.set_ylim(0, 1.1)
axmag.grid()
axphase.grid()
axcohere.set_xlabel("Frequency [Hz]")
axcohere.set_ylabel("Coherence")
axphase.set_ylabel("Phase [deg]")
axmag.set_ylabel("Magnitude")
axmag.legend()
plt.show()
