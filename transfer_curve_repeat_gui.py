# import packages
import tkinter
from tkinter import *
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt

# get inst from connect gui
from connect_gui import use_inst
inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
path = get_folder_path()

# import math functions
from measurement_functions import generateLinSweep, runListSweep

# gui window
window = tkinter.Tk()

# get individual variables
def get_vg_init():
    Vg_init = float(vg_init_entry.get())
    return Vg_init

def get_vg_final():
    Vg_final = float(vg_final_entry.get())
    return Vg_final

def get_vd():
    Vd = float(vd_entry.get())
    return Vd

def get_time_per_point():
    time_per_point = float(time_per_point_entry.get())
    return time_per_point

def get_points_per_sweep():
    points_per_sweep = int(points_per_sweep_entry.get())
    return points_per_sweep

def get_reverse_sweep():
    reverse_sweep = reverse_sweep_is_checked.get()
    return reverse_sweep

def get_ig_range():
    Ig_range = float(ig_range_entry.get())
    return Ig_range

def get_id_range():
    Id_range = float(id_range_entry.get())
    return Id_range

def get_t_hold_before_meas():
    t_hold_before_meas = float(time_hold_entry.get())
    return t_hold_before_meas

def get_n_repeat():
    n_repeat = int(n_repeat_entry.get())
    return n_repeat

# check settings
def set_settings():
    # get entered labels
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    Vd = get_vd()
    time_per_point = get_time_per_point()
    points_per_sweep = get_points_per_sweep()
    reverse_sweep = get_reverse_sweep()
    Ig_range = get_ig_range()
    Id_range = get_id_range()
    t_hold_before_meas = get_t_hold_before_meas()
    n_repeat = get_n_repeat()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nVg init: {Vg_init}"
                                  f"\nVg final: {Vg_final}"
                                  f"\nVd: {Vd}"
                                  f"\nTime per point: {time_per_point}"
                                  f"\nPoints per sweep: {points_per_sweep}"
                                  f"\nReverse sweep: {reverse_sweep}"
                                  f"\nIg range: {Ig_range}"
                                  f"\nId range: {Id_range}"
                                  f"\nTime Hold: {t_hold_before_meas}"
                                  f"\nNumber of repeat: {n_repeat}")

# do measurement
def run_repeat_transfer_curve():
    # get entered values
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    global Vd
    Vd = get_vd()
    time_per_point = get_time_per_point()
    points_per_sweep = get_points_per_sweep()
    reverse_sweep = get_reverse_sweep()
    Ig_range = get_ig_range()
    Id_range = get_id_range()
    t_hold_before_meas = get_t_hold_before_meas()
    n_repeat = get_n_repeat()

    # Generate sweep for gate voltage based on user input
    vg_list = generateLinSweep(points_per_sweep, Vg_init, Vg_final, reverse=reverse_sweep)
    vd_list = "{:.1E}".format(Vd)

    # parameters format (list of strings, units of seconds, amps, volts):
    # [(0) measurment points,(1) time per point,(2) measurement delay,(3) aquisition time,
    #  (4) gate current range,(5) drain current range,(6) gate points,(7) drain points,
    #  (8) pre-measurement voltage hold time]

    parameters = []
    parameters.append("{:.0f}".format(len(vg_list.split(','))))
    parameters.append("{:.1E}".format(time_per_point))
    parameters.append("{:.1E}".format(time_per_point / 2))
    parameters.append("{:.1E}".format((time_per_point / 2) * 0.8))
    parameters.append("{:.0E}".format(Ig_range))
    parameters.append("{:.0E}".format(Id_range))
    parameters.append("{:.0f}".format(len(vg_list.split(','))))
    parameters.append("{:.0f}".format(len(vd_list.split(','))))
    parameters.append(t_hold_before_meas)

    # Create empty array to store output data
    transfer = np.empty((n_repeat, int(parameters[0]), 4))
    Vg = np.asarray([float(i) for i in vg_list.split(',')])
    # Run transfer curve for each drain voltage and plot results
    # Run output curve for each gate voltage and plot results

    for i in range(0, n_repeat):
        transfer[i] = runListSweep(parameters, vg_list, vd_list, wait=False)
        plt.plot(Vg, transfer[0, :, 2])
        plt.show()

    plt.gca().invert_yaxis()
    plt.xlabel("$V_{GS}$ (V)")
    plt.ylabel("$I_{DS}$ (A)")
    plt.title("Transfer curves")
    # plt.legend()

    # Turn SMU output off after all measurements
    inst.write(":outp1 off")
    inst.write(":outp2 off")

def get_vd_for_plot():
    global Vd
    return Vd

# create widgets
# create input fields for setting
Label(window, text="Initial Vg: ").grid(row=0, column=0)
vg_init_entry = Entry(window)
vg_init_entry.grid(row=0, column=1)

Label(window, text="Final Vg: ").grid(row=1, column=0)
vg_final_entry = Entry(window)
vg_final_entry.grid(row=1, column=1)

Label(window, text="Vd: ").grid(row=2, column=0)
vd_entry = Entry(window)
vd_entry.grid(row=2, column=1)

Label(window, text="Time Per Point: ").grid(row=3, column=0)
time_per_point_entry = Entry(window)
time_per_point_entry.grid(row=3, column=1)

Label(window, text="Points Per Sweep: ").grid(row=4, column=0)
points_per_sweep_entry = Entry(window)
points_per_sweep_entry.grid(row=4, column=1)

reverse_sweep_is_checked = BooleanVar()
Checkbutton(window, text="Reverse Sweep: ", variable=reverse_sweep_is_checked).grid(row=5, columnspan=2)

Label(window, text="Ig Range (in 1e-X format with X between 0 and 7 inclusive): ").grid(row=6, column=0)
ig_range_entry = Entry(window)
ig_range_entry.grid(row=6, column=1)

Label(window, text="Id Range (in 1e-X format with X between 0 and 7 inclusive): ").grid(row=7, column=0)
id_range_entry = Entry(window)
id_range_entry.grid(row=7, column=1)

Label(window, text="Time Hold Before Measurement: ").grid(row=8, column=0)
time_hold_entry = Entry(window)
time_hold_entry.grid(row=8, column=1)

Label(window, text="Number of Repeats: ").grid(row=9, column=0)
n_repeat_entry = Entry(window)
n_repeat_entry.grid(row=9, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=10, column=0)

# run transfer curve button
Button(window, text="Run Transfer Curves: ", command=run_repeat_transfer_curve).grid(row=11, column=0)

# run gui
window.mainloop()