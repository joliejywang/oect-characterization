# import packages
import tkinter
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit

# get freqs from ac measurement voltage sweep
from ac_measurement_voltage_sweep_gui import get_freqs, get_points_per_cycle, get_sweeps, get_cycles, get_vg_amp
freqs = get_freqs()
points_per_cycle = get_points_per_cycle()
sweeps = get_sweeps()
cycles = get_cycles()
Vg_amplitude = get_vg_amp()

# get Vd from repeat transfer curve measurement
from transfer_curve_repeat_gui import get_vd_for_plot
Vd = get_vd_for_plot()

# get measurement
from math_functions import fitSineV3

# get plot function
from plot import plotMobility

# gui window
window = tkinter.Tk()

# get fit cycles entry
def get_fit_cycles():
    fit_cycles = int(fit_cycles_entry.get())
    return fit_cycles

# get plot data entry values
def get_points_omitted():
    points_omitted = int(points_omitted_entry.get())
    return points_omitted

def get_length():
    length = int(length_entry.get())
    return length

# extract freq data
def extract_freq_dep_data():
    fit_cycles = get_fit_cycles()

    # global variable
    global Id_fit, Ig_fit

    # Initialize variables for the fit results (format amp, phi, I0)
    Ig_fit = np.empty((len(freqs), 4))
    # Ig_cov = np.empty((len(freqs),3,3))
    Id_fit = np.empty((len(freqs), 4))
    # Id_cov = np.empty((len(freqs),3,3))

    # Setup figure for checking fits
    fig, axs = plt.subplots(len(freqs), 2)
    fig = plt.gcf()
    fig.set_size_inches(12, 10)
    freq = np.empty((len(freqs)))

    ## Add fit bounds

    for i in range(0, len(freqs)):
        fit_start_ind = round((cycles - fit_cycles) * points_per_cycle[i])

        # Get current data to fit from data
        Ig = sweeps[i][:, 0]
        Id = sweeps[i][:, 2]
        t = sweeps[i][:, 1]

        # Set global variable for frequency used
        f_freq = freqs[i]

        # Initialize the fitting for the gate amp, phi, I0
        p0 = [np.mean(Ig), freqs[i], 0.25 * (np.amax(Ig) - np.amin(Ig)), 0]
        # bounds = [[0,0,np.amin(Ig)],[(np.amax(Ig)-np.amin(Ig)),2*np.pi,np.amax(Ig)]]

        # Run curve optimisation for the sine fitting
        # Ig_fit[i],_ = scipy.optimize.curve_fit(fitSineV3, t[fit_start_ind:], Ig[fit_start_ind:], p0=p0)
        Ig_fit[i], _ = scipy.optimize.curve_fit(fitSineV3, t[fit_start_ind:], Ig[fit_start_ind:], p0=p0)

        # Initialize the fitting for the drain
        p0 = [np.mean(Id), freqs[i], 0.5 * (np.amax(Id) - np.amin(Id)), 0]
        # bounds = [[0,0,np.amin(Id)],[1.1*(np.amax(Id)-np.amin(Id)),2*np.pi,np.amax(Id)]]

        Id_fit[i], _ = scipy.optimize.curve_fit(fitSineV3, t[fit_start_ind:], Id[fit_start_ind:], p0=p0)

        # Plot fits to verify optimization is working
        axs[i, 0].plot(t, Ig, 'x')
        axs[i, 0].plot(t, fitSineV3(t, Ig_fit[i, 0], Ig_fit[i, 1], Ig_fit[i, 2], Ig_fit[i, 3]))
        axs[i, 1].plot(t, Id, 'x')
        axs[i, 1].plot(t, fitSineV3(t, Id_fit[i, 0], Id_fit[i, 1], Id_fit[i, 2], Id_fit[i, 3]))

def plot_extract_values_data():
    points_omitted = get_points_omitted()
    length = get_length()

    # global variable
    global Id_fit, Ig_fit

    plt.figure(0)
    plt.loglog(freqs[1:], abs(Id_fit[1:, 2]) / Vg_amplitude, 'o')
    plt.ylabel("gm")

    plt.figure(1)
    plt.loglog(freqs[1:], Vg_amplitude / abs(Ig_fit[1:, 2]), 'o')
    plt.ylabel("Z")

    plt.figure(2)
    mobility = plotMobility(Ig_fit, Id_fit, freqs, points_omitted, length, Vd)

# how many cycles to fit entry
Label(window, text="Number of cycles to fit: ").grid(row=0, column=0)
fit_cycles_entry = Entry(window)
fit_cycles_entry.grid(row=0, column=1)

# plotting variables entry
Label(window, text="Points omitted: ").grid(row=1, column=0)
points_omitted_entry = Entry(window)
points_omitted_entry.grid(row=1, column=1)

Label(window, text="Channel length in microns: ").grid(row=2, column=0)
length_entry = Entry(window)
length_entry.grid(row=2, column=1)

# extract data button
Button(window, text="Extract Frequency Dependent Data: ", command=extract_freq_dep_data).grid(row=3, column=0)

# plot button
Button(window, text="Plot Data and Extracted Values: ", command=plot_extract_values_data).grid(row=4, column=0)

# run gui
window.mainloop()