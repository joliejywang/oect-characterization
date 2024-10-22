# import packages
import numpy as np
import matplotlib.pyplot as plt

# # Plots transconductance
# def plotTransconductance(delta_Id, Vg_amp, freqs):
#     gm = np.divide(delta_Id, Vg_amp)
#     plt.loglog(freqs, gm)
#     plt.xlabel("Frequency (Hz)")
#     plt.ylabel("gm (S)")
#     return gm

# # Plots impedance graph
# def plotImpedance(Ig_fit, Vg_amp, freqs, points_omitted):
#     Z = (np.divide(Vg_amp, Ig_fit[:, 0]))
#     fig, ax1 = plt.subplots()
#     color = 'tab:blue'
#     ax1.set_xlabel('Frequency (Hz)')
#     ax1.set_ylabel('|Z| (Ohms)', color=color)
#     ax1.loglog(freqs, Z, 'x-')
#     m, b = np.polyfit(np.log10(freqs[:len(freqs) - points_omitted]), np.log10(Z[:len(Z) - points_omitted]), 1)
#     Z_fit = np.power(10, m * np.log10(freqs) + b)
#     ax1.loglog(freqs, Z_fit, 'r--')
#     ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#     color = 'tab:red'
#     ax2.set_ylabel('Phase (deg)', color=color)  # we already handled the x-label with ax1
#     phase = Ig_fit[:, 1]
#     ax2.semilogx(freqs, phase, color=color)
#     ax2.tick_params(axis='y', labelcolor=color)
#     fig.tight_layout()  # otherwise the right y-label is slightly clipped
#     capacitance = np.divide(np.power(10, -1 * b), 2 * np.pi)
#     print("Capacitance = " + str(capacitance) + " F")
#     plt.show()
#     return capacitance

# Plots curve for mobility fitting
def plotMobility(Ig_fit, Id_fit, freqs, points_omitted, length, Vd):
    freq_D_Id = np.multiply(Id_fit[:, 0], freqs)[:len(freqs) - points_omitted]
    plt.plot(freq_D_Id, Ig_fit[:len(freqs) - points_omitted, 0])
    m, b = np.polyfit(freq_D_Id, Ig_fit[:len(freqs) - points_omitted, 0], 1)
    plt.plot(freq_D_Id, m * freq_D_Id + b, 'rx--')
    plt.ylabel("∆Ig (A)")
    plt.xlabel("∆Id*freq (A s^-1)")
    tau = np.divide(m, 2 * np.pi)
    L = length * 1e-4
    mu = np.divide(np.square(L), np.multiply(abs(Vd), tau))
    print("Mobility = " + str(mu) + " cm^2 V^-1 s^-1")
    return mu

# Plots pulse transient
def plotPulse(pulse):
    fig, ax1 = plt.subplots()
    color = 'tab:blue'
    ax1.set_ylabel('I$_{DS}$ (A)')
    ax1.set_xlabel('time (s)', color=color)
    ax1.plot(pulse[:, 1], pulse[:, 2])
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('I$_{GS}$ (A)', color=color)  # we already handled the x-label with ax1
    ax2.plot(pulse[:, 1], pulse[:, 0], color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()