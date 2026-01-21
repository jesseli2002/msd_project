import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.signal import welch, csd
import scipy.signal as signal
from numpy.fft import rfft, rfftfreq
from numpy import pi, sin, cos, squeeze
import control as ct
from control_utils import (
    impulse_response_freq_domain,
    make_lpf_butter,
    make_integrator,
    make_double_pole,
    make_ll,
    make_notch,
)


class HackedFRD:
    """
    Hacky frequency response data class
    frequency range is hard coded
    """

    # defined frequency data
    omega_data = np.geomspace(1e1, 1e5, 4 * 200 + 1)

    freq_data = omega_data / (2 * pi)

    def __init__(self, resp):
        """
        Docstring for __init__

        :param self: Description
        :param resp: Complex array giving transfer function
        """
        self._resp = resp

    @classmethod
    def from_complex(cls, mag, phase):
        return HackedFRD(mag * np.exp(1j * phase))

    @classmethod
    def from_lti(cls, sys: ct.LTI):
        mag, phase, _ = sys.frequency_response(cls.omega_data, True)
        return cls.from_complex(mag, phase)

    def get_mag_phase(self):
        return np.abs(self._resp), np.angle(self._resp)

    def get_response(self):
        return self._resp

    def to_control_frd(self):
        return ct.frd(self._resp, self.omega_data, smooth=True)

    # Arithmetic operations
    def __add__(self, other):
        """Addition: parallel connection for transfer functions"""
        if isinstance(other, HackedFRD):
            return HackedFRD(self._resp + other._resp)
        else:
            return HackedFRD(self._resp + other)

    def __radd__(self, other):
        """Right addition"""
        return self.__add__(other)

    def __sub__(self, other):
        """Subtraction"""
        if isinstance(other, HackedFRD):
            return HackedFRD(self._resp - other._resp)
        else:
            return HackedFRD(self._resp - other)

    def __rsub__(self, other):
        """Right subtraction"""
        return HackedFRD(other - self._resp)

    def __mul__(self, other):
        """Multiplication: series connection for transfer functions"""
        if isinstance(other, HackedFRD):
            return HackedFRD(self._resp * other._resp)
        else:
            return HackedFRD(self._resp * other)

    def __rmul__(self, other):
        """Right multiplication"""
        return self.__mul__(other)

    def __truediv__(self, other):
        """Division"""
        if isinstance(other, HackedFRD):
            return HackedFRD(self._resp / other._resp)
        else:
            return HackedFRD(self._resp / other)

    def __rtruediv__(self, other):
        """Right division"""
        return HackedFRD(other / self._resp)

    def __neg__(self):
        """Negation"""
        return HackedFRD(-self._resp)

    def __pow__(self, power):
        """Power operation"""
        return HackedFRD(self._resp**power)


def make_delay(delay_time):
    """
    Returns delay FRD
    :param delay_time: Delay time
    """
    return HackedFRD(np.exp(-delay_time * HackedFRD.omega_data * 1j))


def add_bode_plot(fig, tf: HackedFRD, label):
    axm, axp = fig.axes
    freq = HackedFRD.omega_data / (2 * pi)

    mag, phase = tf.get_mag_phase()
    mag_db = 20 * np.log10(mag)
    axm.loglog(freq, np.squeeze(mag), label=label)
    # axm.semilogx(freq, mag_db, label=label)
    axp.semilogx(freq, np.unwrap(np.squeeze(phase)) * 180 / pi)


def make_bode_fig_ax(height_ratios=None):
    if height_ratios is None:
        height_ratios = [1, 1]
    fig, ax = plt.subplots(2, sharex=True, height_ratios=height_ratios)

    axm, axp = ax
    axm.grid()
    axp.grid()
    axp.set_xlabel("Frequency [Hz]")
    axp.set_ylabel("Phase [deg]")
    # axm.set_ylabel("Magnitude [dB]")
    axm.set_ylabel("Magnitude")
    return fig, ax


mat = loadmat("dat/MSD2025_P3_Plant_numdendelay.mat")
num = np.squeeze(np.squeeze(mat["num"]).item())
den = np.squeeze(np.squeeze(mat["den"]).item())
delay = np.squeeze(mat["delay"])

plant = HackedFRD.from_lti(ct.tf(num, den)) * make_delay(delay)

C1 = 1770 / ct.tf("s")
C2 = make_double_pole(8000, 1) / make_double_pole(4610, 0.01)
controller = HackedFRD.from_lti(C1 * C2)

loop = controller * plant


fig, (axm, axp) = make_bode_fig_ax()
add_bode_plot(fig, controller, "controller")
add_bode_plot(fig, plant, "plant")
add_bode_plot(fig, loop, "loop")
axm.legend()
axm.set_title("Open loop transfer functions")


# Plot sensitivity ----------------------------
# Evaluate controller modulus margin at same time
S = 1 / (1 + loop)
T = loop * S

mag, _ = S.get_mag_phase()
modulus_margin = np.max(mag)
modulus_margin_db = 20 * np.log10(modulus_margin)

print(f"max(S(jw)) = {modulus_margin_db=} dB @ w = {HackedFRD.freq_data[np.argmax(mag)]:.1f} Hz")

fig, (axm, axp) = make_bode_fig_ax(height_ratios=[2, 1])
add_bode_plot(fig, S, 'sensitivity')
add_bode_plot(fig, T, 'complementary sensitivity')
fig.set_size_inches(8, 6)
fig.axes[0].legend(loc='lower left')
fig.axes[0].set_xlim([20, 2e4])
fig.axes[0].set_ylim(bottom=1e-8)
fig.axes[0].set_title(f"Sensitivity functions; max(S(j$\\omega$))={float(modulus_margin_db):.2f} dB @ $\\omega$ = {HackedFRD.freq_data[np.argmax(mag)]:.1f} Hz")


# Plot margins ----------------------------
fig, (axm, axp) = make_bode_fig_ax(height_ratios=[2, 1])
loop_frd = loop.to_control_frd()
gm, pm, wg, wp = ct.margin(loop_frd)
ct.bode_plot(loop_frd, display_margins=True, dB=False, Hz=True, title='')
fig = plt.gcf()
fig.axes[0].set_title(f'Gain margin: {20*np.log10(gm):.2f} dB at {wg/(2*pi):.1f} Hz; Phase margin: {pm:.2f}° at {wp/(2*pi):.1f} Hz')
fig.set_size_inches(8, 6)


#  C.2 =================================
#  C.2 =================================
#  C.2 =================================


plt.show()
plt.close('all')
